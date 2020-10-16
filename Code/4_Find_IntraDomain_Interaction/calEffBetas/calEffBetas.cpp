#ifndef FINDINTRA_CPP
#define FINDINTRA_CPP

#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <vector>

#include <gsl/gsl_cdf.h>

#include <scythestat/rng/mersenne.h>
#include <scythestat/distributions.h>
#include <scythestat/ide.h>
#include <scythestat/la.h>
#include <scythestat/matrix.h>
#include <scythestat/rng.h>
#include <scythestat/smath.h>
#include <scythestat/stat.h>
#include <scythestat/optimize.h>

#include <IRLS_glm/IRLS.h>

//#include <RInside.h>

#include <mlpack/methods/lars/lars.hpp>
#include <mlpack/methods/linear_regression/linear_regression.hpp>
#include <boost/test/unit_test.hpp>

#define INTTAG 0
#define BOOLTAG 1
#define DIETAG 2

#define TESTFREQTHRES 1
#define NEIGHBDIS 5
//#define MORANI 0.001
#define PCUTOFF 0.004
#define PRECISION 10e48
#define READBLOCK 409600
#define OPTDIST 2000
#define MAXDIST 100000

//using namespace std;
using std::abs;
using std::basic_string;
using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::exception;
using std::flush;
using std::ifstream;
using std::istringstream;
using std::make_pair;
using std::map;
using std::ofstream;
using std::ostringstream;
using std::pair;
using std::sort;
using std::string;
using std::vector;


int reader_dif(const string&, vector< vector< vector<int> > >&);
int reader_dsm(const string&, vector< vector<int> >&);
void reader_ictm(const string&, map<int, int>&);
int reader_coef(const string&, map<int, double>&);
int findIntraDomainInteraction(map<int, int>&, string, vector< vector<int> >&, vector< vector< vector<int> > >&, map<int, int>&, map<int, double>&, map< int, map< pair<int, int>, double > >&, double, double, double);
int rmDependentCol(vector< vector<double> >&, vector<double>&, vector< pair<int, int> >&, double, int&, vector< vector<double> >&);
int calStanScores(vector< vector<double> >&);
int calPearsonCorr(vector< vector<double> >&, vector< pair<int, int> >&, const double&, int, map<int, int>&);
double calBgFrags(map<int, double>&, int);
int writer_corrFile(map< pair<int, int>, double >&, string&);

void usage();

inline void usage()
{
    cout <<"Usage: findIntraDomainInteraction empericalDistributionFile interChrFreqFile domainSitesFile domainInteractionFreqFile chrNum MORANI CorrelationThresholdForScreen"<<endl;
}

//convtIndex is now 0 based
inline int convtIndex(int row_fun, int col_fun, int siteNum_fun)
{
    if (row_fun==col_fun || row_fun>=siteNum_fun || col_fun>=siteNum_fun || row_fun<0 || col_fun<0)
    {//cout<<"row or column exceed bound."<<endl;
        return -1;
    }
    else
    {
        int max_fun = 0, min_fun = 0;
        if (row_fun<col_fun)
        {
            //switch to 1 based
            max_fun = col_fun+1;
            min_fun = row_fun+1;
        }
        else
        {
            max_fun = row_fun+1;
            min_fun = col_fun+1;
        }
        return((2*siteNum_fun-min_fun)*(min_fun-1)/2+max_fun-min_fun) - 1;//-1 put the output to be 0 based
    }
}


//l=sigma(yi(theta*xi+offseti)-exp(theta*xi+offseti))
class PoissonModel {
    public:
    double operator() (const scythe::Matrix<double> beta){
        const int n = y_.rows();
        const int p = X_.cols();

        scythe::Matrix<double> eta = X_ * beta + offset_;
        scythe::Matrix<double> m = exp(eta);
        double loglike = 0.0;
        for (int i=0; i<n; ++i)
        loglike += y_(i) * log(m(i)) - m(i);
        return -1.0 * loglike;
    }
    scythe::Matrix<double> y_;
    scythe::Matrix<double> X_;
    scythe::Matrix<double> offset_;
};


int main(int argc, char* argv[])
{
    if (argc < 6)
    {
        usage();
        exit(1);
    }

    //string empDisFile = argv[1];
    string ictmFile = argv[1];
    string dsmFile = argv[2];
    string difFile = argv[3];
    string chrNum = argv[4];
    //double MORANI = atof(argv[6]);
    double corrThres = 0.0;
    string disCoefFile = argv[5];
    double disIntcp = 0.0;
    double disCoeff = 0.0;

    map<int, int> empDis_map;
    map<int, int> ictm_map;
    vector< vector<int> > dsm_map;
    vector< vector< vector<int> > > dif_map;
    map<int, double> disPolyCoef;

    //reader_ictm(empDisFile, empDis_map);
    //cout<<"\nReading empDis file done."<<endl;
    reader_ictm(ictmFile, ictm_map);
    cout<<"\nRead ictm file done."<<endl;
    reader_dsm(dsmFile, dsm_map);
    cout<<"\nRead dsm file done."<<endl;
    reader_dif(difFile, dif_map);
    cout<<"\nRead dif file done."<<endl;
    reader_coef(disCoefFile, disPolyCoef);
    cout<<"\nRead distance polynomial fit coefficients done."<<endl;

    map< int, map< pair<int, int>, double > > pVal_map;

    findIntraDomainInteraction(empDis_map, chrNum, dsm_map, dif_map, ictm_map, disPolyCoef, pVal_map, corrThres, disIntcp, disCoeff);

    return 0;
}


int reader_dif(const string& fileToread, vector< vector< vector<int> > >& dif_local)
{
    ifstream inputFile(fileToread.c_str());
    if (!inputFile)
    {
        //char a;
        //cin >>a;
        cout <<"\n"<< "Error opening " << fileToread << "." << endl;
        exit(1);
    }

    char lineStr[READBLOCK];
    char freq_str[15];
    int freq_int = 0;
    int domainNum = 0;

    int column = 0;

    char* it_lineStr;
    char* it_token;

    int peakNum = 0;

    //copy data to a map, which chrNo is the key, the value is a vector of ints(read's starting point)
    //ofstream testBadFile("test");
    vector< vector<int> > tempDomn;
    while(inputFile.getline(lineStr,READBLOCK))
    {
        it_lineStr = lineStr;

        column = 0;

        if (lineStr != NULL)
        {
            vector<int> tempFreq;
            //testBadFile<<lineStr<<endl;

            //pass head white space
            if(*it_lineStr == 'd')
            {
                if (tempDomn.begin() != tempDomn.end())
                {
                    ++domainNum;
                    dif_local.push_back(tempDomn);
                    tempDomn.clear();
                }
                continue;
            }

            while(*it_lineStr != '\0' && *it_lineStr != '\n')
            {
                if (*it_lineStr == ' ' || *it_lineStr == '\t')
                {
                    //remove white space
                    while(*it_lineStr == ' ' || *it_lineStr == '\t')
                    {
                        ++it_lineStr;
                    }

                    ++column;
                }
                else
                {
                    it_token = freq_str;
                    while(*it_lineStr != ' ' && *it_lineStr != '\t' && *it_lineStr != '\0' && *it_lineStr != '\n')
                    {
                        *it_token = *it_lineStr;
                        ++it_token;
                        ++it_lineStr;
                    }
                    *it_token = '\0';
                    freq_int = atoi(freq_str);
                    tempFreq.push_back(freq_int);
                }
            }

            tempDomn.push_back(tempFreq);
            ++peakNum;
        }
    }

    if(tempDomn.begin() != tempDomn.end())
    {
        dif_local.push_back(tempDomn);
    }

    inputFile.close();

    return peakNum;
}




int reader_dsm(const string& fileToread, vector< vector<int> >& dsm_local)
{
    ifstream inputFile(fileToread.c_str());
    if (!inputFile)
    {
        //char a;
        //cin >>a;
        cout <<"\n"<< "Error opening " << fileToread << "." << endl;
        exit(1);
    }

    char lineStr[READBLOCK];
    char pos_str[15];
    int pos_int = 0;

    int column = 0;

    char* it_lineStr;
    char* it_token;

    int peakNum = 0;

    //copy data to a map, which chrNo is the key, the value is a vector of ints(read's starting point)
    //ofstream testBadFile("test");
    while(inputFile.getline(lineStr,READBLOCK))
    {
        it_lineStr = lineStr;

        column = 0;

        if (lineStr != NULL)
        {
            vector<int> tempPos;
            //testBadFile<<lineStr<<endl;

            //pass head white space
            while(*it_lineStr == ' '|| *it_lineStr == '\t')
            {
                ++it_lineStr;
            }

            while(*it_lineStr != '\0' && *it_lineStr != '\n')
            {
                if (*it_lineStr == ' ' || *it_lineStr == '\t')
                {
                    //remove white space
                    while(*it_lineStr == ' ' || *it_lineStr == '\t')
                    {
                        ++it_lineStr;
                    }

                    ++column;
                }
                else
                {
                    it_token = pos_str;
                    while(*it_lineStr != ' ' && *it_lineStr != '\t' && *it_lineStr != '\0' && *it_lineStr != '\n')
                    {
                        *it_token = *it_lineStr;
                        ++it_token;
                        ++it_lineStr;
                    }
                    *it_token = '\0';
                    pos_int = atoi(pos_str);
                    tempPos.push_back(pos_int);
                }
            }

            dsm_local.push_back(tempPos);
            ++peakNum;
        }
    }

    inputFile.close();

    return peakNum;
}



int reader_coef(const string& fileToread, map<int, double>& disPolyCoef_local)
{
    ifstream inputFile(fileToread.c_str());
    if (!inputFile)
    {
        cout <<"\n"<< "Error opening " << fileToread << "." << endl;
        exit(1);
    }

    char degree[15];
    int degree_int = 0;
    char coef[50];
    double coef_dbl = 0.0;

    char* it_lineStr;
    char* it_token;

    char lineStr[512];
    while(inputFile.getline(lineStr,512))
    {
        if (lineStr != NULL)
        {
            it_lineStr = lineStr;

            it_token = degree;
            while(*it_lineStr != '\t' && *it_lineStr != ' ' && *it_lineStr != '\n' && *it_lineStr != '\0')
            {
                *it_token = *it_lineStr;
                ++it_token;
                ++it_lineStr;
            }
            *it_token = '\0';
            degree_int = atoi(degree);
	    ++it_lineStr;
	    while(*it_lineStr == ' ' || *it_lineStr == '\t')
	    {
		++it_lineStr;
	    }

            it_token = coef;
            while(*it_lineStr != '\t' && *it_lineStr != ' ' && *it_lineStr != '\n' && *it_lineStr != '\0')
            {
                *it_token = *it_lineStr;
                ++it_token;
                ++it_lineStr;
            }
            *it_token = '\0';
            coef_dbl = atof(coef);

            disPolyCoef_local[degree_int] = coef_dbl;
        }
    }
    return 0;
}
void reader_ictm(const string& fileToread, map<int, int>& ictm_local)
{
    ifstream inputFile(fileToread.c_str());
    if (!inputFile)
    {
        cout <<"\n"<< "Error opening " << fileToread << "." << endl;
        exit(1);
    }

    char pos[15];
    int pos_int = 0;
    char freq[15];
    int freq_int = 0;

    char* it_lineStr;
    char* it_token;

    char lineStr[512];
    while(inputFile.getline(lineStr,512))
    {
        if (lineStr != NULL)
        {
            it_lineStr = lineStr;

            it_token = pos;
            while(*it_lineStr != '\t')
            {
                *it_token = *it_lineStr;
                ++it_token;
                ++it_lineStr;
            }
            *it_token = '\0';
            ++it_lineStr;
            pos_int = atoi(pos);

            it_token = freq;
            while(*it_lineStr != '\t' && *it_lineStr != '\n' && *it_lineStr != '\0')
            {
                *it_token = *it_lineStr;
                ++it_token;
                ++it_lineStr;
            }
            *it_token = '\0';
            freq_int = atoi(freq);

            ictm_local[pos_int] = freq_int;
        }
    }
}



int findIntraDomainInteraction(map<int, int>& empDisMap_fun, string jobChr_fun, vector< vector<int> >& domainSitesMap_fun, vector< vector< vector<int> > >& domainCSinterFreq_fun, map<int, int>& csInterChromTotalMap_fun, map<int, double>& disPolyCoef_local, map< int, map< pair<int, int>, double > >& pVal_fun, double corrThres_fun, double disIntcp_fun, double disCoeff_fun)
{
    //string file_effBetas =  jobChr_fun + "_domainEffBetas";
    //ofstream outEffBetas(file_effBetas.c_str());

    //srand((unsigned)time(0));
    //double rndNum = rand()/(RAND_MAX + 1.0);

    int domainNum = domainSitesMap_fun.size();

    vector<double> freq_convert;
    vector<double> offset;//distance effect
    vector<double> effVec;//efficiency of cutting site and mappability, measured by the total inter chromosome hybrid frags
    //#pragma omp parallel for
    for (int domainIt = 0; domainIt < domainNum; ++domainIt)
    {
        map< pair<int, int>, double>& domainPval = pVal_fun[domainIt];
        vector<int>& sitesMap = domainSitesMap_fun[domainIt];
        vector< vector<int> >& csInterFreqMap = domainCSinterFreq_fun[domainIt];
        int siteNum = sitesMap.size();

        int size_convert = convtIndex(siteNum-2,siteNum-1,siteNum) + 1;//the convtIndex function's first two input is 0 based and the output is 0 based too
        for (int rowIt=0; rowIt<siteNum; ++rowIt)
        {
            for (int colIt=rowIt+1; colIt<siteNum; ++colIt)
            {
                int index_convert=convtIndex(rowIt,colIt,siteNum);

		if(abs(sitesMap[rowIt]-sitesMap[colIt])>OPTDIST && abs(sitesMap[rowIt]-sitesMap[colIt])<MAXDIST)//the index of freq_convert need to be changed because of this
		{
		    //log(y/ydist)=mu+betaeff*effvec
		    freq_convert.push_back(csInterFreqMap[rowIt][colIt]);
		    offset.push_back(calBgFrags(disPolyCoef_local, abs(sitesMap[rowIt]-sitesMap[colIt])));
		    //the digestion and ligation efficient of a cutting site should take log because, its log should proportional to log(freq)
		    effVec.push_back(log(sqrt(csInterChromTotalMap_fun[sitesMap[rowIt]])*sqrt(csInterChromTotalMap_fun[sitesMap[colIt]])+1));
		}
            }
        }
    }

    int size_convert_s = freq_convert.size();
    string file_freq =  jobChr_fun + "_freq";
    string file_offset =  jobChr_fun + "_offset";
    string file_effvec =  jobChr_fun + "_effvec";
    //ofstream outEffBetas(file_effBetas.c_str());
    ofstream outFreq(file_freq);
    ofstream outOffset(file_offset);
    ofstream outEffvec(file_effvec);
    for(int i=0;i<size_convert_s;++i)
    {
	outFreq<<freq_convert[i]<<endl;
	outOffset<<offset[i]<<endl;
	outEffvec<<effVec[i]<<endl;
    }
    cout<<"Out matrix done."<<endl;
    outFreq.close();
    outOffset.close();
    outEffvec.close();

    return 0;
}




int rmDependentCol(vector< vector<double> >& distMatrix_local, vector<double>& freqVec_local, vector< pair<int, int> >& oriCoor_local, double thres_local, int& testIndex_local, vector< vector<double> >& distMatRmDep_local)
{
    int rowNum_fun = distMatrix_local.size();
    if (rowNum_fun != freqVec_local.size())
    {
        cerr<<"Error: row number of the dist matrix is different from the size of the frequency vector."<<endl;
        exit(1);
    }
    int colNum_fun = (distMatrix_local.begin())->size();
    vector< vector<double> > distMatrix_std;
    for (int i = 0; i < colNum_fun; ++i)
    {
        vector<double> tempVec;
        for (int j = 0; j < rowNum_fun; ++j)
        {
            tempVec.push_back(distMatrix_local[j][i]);
        }
        distMatrix_std.push_back(tempVec);
    }
    distMatrix_std.push_back(freqVec_local);
    calStanScores(distMatrix_std);

    map< pair<int, int>, double > corrMap;
    map<int, int> rmIndex;
    int testFlag_local = calPearsonCorr(distMatrix_std, oriCoor_local, thres_local, testIndex_local, rmIndex);
    cout<<"testFlag: "<<testFlag_local<<endl;

    int oriTestIndex = testIndex_local;
    for (int m = 0; m < oriTestIndex; ++m)
    {
        if (rmIndex[m]==1)
        {
            testIndex_local--;
        }
    }

    for (int i = 0; i < rowNum_fun; ++i)
    {
        vector<double> tempVec;
        for (int j = 0; j < colNum_fun; ++j)
        {
            if (rmIndex[j]!=1)
            {
                tempVec.push_back(distMatrix_local[i][j]);
            }
        }
        distMatRmDep_local.push_back(tempVec);
    }

    return testFlag_local;
}




int calStanScores(vector< vector<double> >& distMatrix_std_fun)
{
    for (vector< vector<double> >::iterator colIt = distMatrix_std_fun.begin(), colIt_end = distMatrix_std_fun.end(); colIt != colIt_end; ++colIt) //colomn as in original matrix before transpose
    {
        int rowNum_fun = colIt->size();
        double sum = 0.0;
        for (vector<double>::const_iterator rowIt = colIt->begin(), rowIt_end = colIt->end(); rowIt != rowIt_end; ++rowIt)
        {
            sum += *rowIt;
        }
        double mean = sum / rowNum_fun;

        double var = 0.0;
        for (vector<double>::const_iterator rowIt = colIt->begin(), rowIt_end = colIt->end(); rowIt != rowIt_end; ++rowIt)
        {
            double diff = *rowIt - mean;
            var += diff * diff;
        }
        var = var / (rowNum_fun - 1);
        double std = sqrt(var);

        for (vector<double>::iterator rowIt = colIt->begin(), rowIt_end = colIt->end(); rowIt != rowIt_end; ++rowIt)
        {
            *rowIt = (*rowIt - mean) / std;
        }
    }

    return 0;
}




int calPearsonCorr(vector< vector<double> >& distMatrix_std_fun, vector< pair<int, int> >& oriCoor_fun, const double& thres_fun, int testIndex_fun, map<int, int>& rmIndex_fun)
{
    int sigNum = 0;

    map< pair<int, int>, double > corrMap_fun;
    vector<double> xyCorr;
    int colNum_fun = distMatrix_std_fun.size(); // colomn as in original matrix
    int rowNum_fun = (distMatrix_std_fun.begin())->size();
    int freedom = rowNum_fun - 1;
    for (int i = 0; i < colNum_fun; ++i)
    {
        for (int j = i+1; j < colNum_fun; ++j)
        {
            double pearsonCorr = 0.0;
            for (int m = 0; m < rowNum_fun; ++m)
            {
                pearsonCorr += distMatrix_std_fun[i][m] * distMatrix_std_fun[j][m];
            }
            pearsonCorr /= freedom;

            //if (pearsonCorr >= pos_thres || pearsonCorr <= neg_thres)

            if (j==colNum_fun-1)
            {
                xyCorr.push_back(pearsonCorr);
            }
            else if (pearsonCorr >= thres_fun || pearsonCorr <= -thres_fun)
            {
                //cout<<pearsonCorr<<endl;
                corrMap_fun[make_pair(i, j)] = pearsonCorr;
            }
        }
    }

    //string corrFile = "corrFile";
    //writer_corrFile(corrMap_fun, corrFile);

    for (map< pair<int, int>, double >::const_iterator sigIt = corrMap_fun.begin(), sigIt_end = corrMap_fun.end(); sigIt != sigIt_end; ++sigIt)
    {
        int col1 = sigIt->first.first, col2 = sigIt->first.second;
        pair<int, int>& oriCol1 = oriCoor_fun[col1], oriCol2 = oriCoor_fun[col2];
        cout<<testIndex_fun<<"\t"<<col1<<"\t"<<col2<<"\t"<<oriCol1.first<<"\t"<<oriCol1.second<<"\t"<<xyCorr[col1]<<"\t"<<oriCol2.first<<"\t"<<oriCol2.second<<"\t"<<xyCorr[col2]<<endl;

        if(col1==testIndex_fun)
        {
            //if within the 3*3 square centered on the testing spot, there is a spot with higher freq than the testing spot, and this high freq spot has high corr with the testing spot, then skip the testing spot, since distMatrix_std_fun[colNum_fun-1][col1] is not the frequency of col1 (need the origin col1), use corr instead
            if(abs(oriCol1.first-oriCol2.first)>1 || abs(oriCol1.second-oriCol2.second)>1 || xyCorr[col1]>=xyCorr[col2])
            {
                rmIndex_fun[col2] = 1;
            }
            else
            {
                return 0;
            }
        }
        else if(col2==testIndex_fun)
        {
            if(abs(oriCol1.first-oriCol2.first)>1 || abs(oriCol1.second-oriCol2.second)>1 || xyCorr[col2]>=xyCorr[col1])
            {
                rmIndex_fun[col1] = 1;
            }
            else
            {
                return 0;
            }
        }
        else
        {
            //always keep the one with higher positive correlation
            //cout<<"break 1"<<endl;
            rmIndex_fun[xyCorr[col1]>xyCorr[col2]?col2:col1] = 1;
            //cout<<"break 2"<<endl;
        }
    }

    return 1;
}




int writer_corrFile(map< pair<int, int>, double >& corrMap_local, string& fileName)
{
    ofstream outputFile(fileName.c_str());

    for (map< pair<int, int>, double>::const_iterator corrIt = corrMap_local.begin(), corrIt_end = corrMap_local.end(); corrIt != corrIt_end; ++corrIt)
    {
        outputFile<<corrIt->first.first<<"\t"<<corrIt->first.second<<"\t"<<corrIt->second<<endl;
    }
    outputFile.close();

    return 0;
}



double calBgFrags(map<int, double>& coef_fun, int dist_fun)
{
    double logDist=log(dist_fun);
    double bgFrags=0.0;
    map<int, double>::const_iterator coefIt=coef_fun.begin(), coefIt_end=coef_fun.end();
    for( ; coefIt != coefIt_end; ++coefIt)
    {
	bgFrags+=pow(logDist, coefIt->first)*coefIt->second;
    }
    return bgFrags;
}


#endif
