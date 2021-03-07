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

#include <../IRLS_glm/IRLS.h>

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
#define PCUTOFF 0.05
#define PRECISION 10e48
#define READBLOCK 409600
#define OPTDIST 1000
#define DISTHRES 20000

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
using std::pow;
using std::sort;
using std::string;
using std::vector;


int reader_dif(const string&, vector< vector< vector<int> > >&);
int reader_dsm(const string&, vector< vector<int> >&);
int reader_coef(const string&, map<int, double>&);
void reader_ictm(const string&, map<int, int>&);
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
    if (argc < 8)
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
    double mu_eff = atof(argv[6]);
    double beta_eff = atof(argv[7]);

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

    findIntraDomainInteraction(empDis_map, chrNum, dsm_map, dif_map, ictm_map, disPolyCoef, pVal_map, corrThres, mu_eff, beta_eff);

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
    vector< vector<int> > tempDomn;
    while(inputFile.getline(lineStr,READBLOCK))
    {
        it_lineStr = lineStr;

        column = 0;

        if (lineStr != NULL)
        {
            vector<int> tempFreq;

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




int findIntraDomainInteraction(map<int, int>& empDisMap_fun, string jobChr_fun, vector< vector<int> >& domainSitesMap_fun, vector< vector< vector<int> > >& domainCSinterFreq_fun, map<int, int>& csInterChromTotalMap_fun, map<int, double>& disPolyCoef_local, map< int, map< pair<int, int>, double > >& pVal_fun, double corrThres_fun, double mu_eff_fun, double beta_eff_fun)
{
    //string file_debug =  jobChr_fun + "_debug";
    //ofstream outDebug(file_debug.c_str());

    int domainNum = domainSitesMap_fun.size();

    //#pragma omp parallel for
    for (int domainIt = 0; domainIt < domainNum; ++domainIt)
    {
	//if(domainIt==1)
	//{exit(0);}
        map< pair<int, int>, double>& domainPval = pVal_fun[domainIt];
        //outDebug<<"domain "<<domainIt<<endl;
        vector<int>& sitesMap = domainSitesMap_fun[domainIt];
        vector< vector<int> >& csInterFreqMap = domainCSinterFreq_fun[domainIt];
        int siteNum = sitesMap.size();

        int size_convert = convtIndex(siteNum-2,siteNum-1,siteNum) + 1;//the convtIndex function's first two input is 0 based and the output is 0 based too
	vector<double> freq_convert;
        vector<double> offset;//distance effect
        vector<double> effVec;//efficiency of cutting site and mappability, measured by the total inter chromosome hybrid frags

	ostringstream domainIt_str;
	domainIt_str<<domainIt;
	//string domainFile = "domainData_" + domainIt_str.str();
	//ofstream outDomainData(domainFile);

        for (int rowIt=0; rowIt<siteNum; ++rowIt)
        {
            for (int colIt=rowIt+1; colIt<siteNum; ++colIt)
            {
                int index_convert=convtIndex(rowIt,colIt,siteNum);

		//if(abs(sitesMap[rowIt]-sitesMap[colIt])>OPTDIST)//the index of freq_convert need to be changed because of this
		{
		    //log(y/ydist)=mu+betaeff*effvec
		    freq_convert.push_back(csInterFreqMap[rowIt][colIt]);
		    offset.push_back(calBgFrags(disPolyCoef_local, abs(sitesMap[rowIt]-sitesMap[colIt])));
		    //the digestion and ligation efficient of a cutting site should take log because, its log should proportional to log(freq)
		    effVec.push_back(log(sqrt(csInterChromTotalMap_fun[sitesMap[rowIt]])*sqrt(csInterChromTotalMap_fun[sitesMap[colIt]])+1));
		    //outDomainData<<rowIt<<"\t"<<colIt<<"\t"<<csInterFreqMap[rowIt][colIt]<<"\t"<<csInterChromTotalMap_fun[sitesMap[rowIt]]<<"\t"<<csInterChromTotalMap_fun[sitesMap[colIt]]<<"\t"<<log(sqrt(csInterChromTotalMap_fun[sitesMap[rowIt]])*sqrt(csInterChromTotalMap_fun[sitesMap[colIt]])+1)<<"\t"<<sitesMap[rowIt]<<"\t"<<sitesMap[colIt]<<"\t"<<calBgFrags(disPolyCoef_local, abs(sitesMap[rowIt]-sitesMap[colIt])+OPTDIST)<<endl;
		}
            }
        }
	//outDomainData.close();

	vector<double> constTerms(size_convert, 0.0);
        vector<double> effConsts(size_convert, 0.0);//efficiency constand mu+beta_eff_fun*effVec[i]
	for (int i = 0; i < size_convert; ++i)
	{
	    //constTerms[i] = mu+dist_ij[i]*alpha+effVec[i]*beta;
	    effConsts[i] = mu_eff_fun+effVec[i]*beta_eff_fun;
	    constTerms[i] = exp(effConsts[i]+offset[i]);
	    //constTerms[i] = mu+effVec[i]*beta+offset[i];
	}

        int counter_fit=0;
        int counter_skip=0;

	int regionIt=0;
	int testIt=0;

	//string regionFile = "regionData_" + domainIt_str.str();
	//ofstream outRegionData(regionFile);
	bool flag_test = false;

	//string distMatrixFile = "distMatrix_" + domainIt_str.str();
	//ofstream outDistMatrix(distMatrixFile);

	string oneColTestFile = "oneCol_" + domainIt_str.str();
	ofstream oneColTest(oneColTestFile);
	string testPosPfile= "testPosP_" + domainIt_str.str();
	ofstream testPosP(testPosPfile);

        //#pragma omp parallel for
        for (int rowIt=0; rowIt<siteNum; ++rowIt)
        {
            for (int colIt=rowIt+1; colIt<siteNum; ++colIt)
            {
                if(csInterFreqMap[rowIt][colIt]>TESTFREQTHRES && abs(sitesMap[rowIt]-sitesMap[colIt])>=DISTHRES)
                {
                    cout<< "colIt\t"<<colIt<<"\trowIt\t"<<rowIt<<endl;
                    vector<double> neighbFreq;
                    vector<double> neighbConstTerms;
                    vector<double> fullExp;

                    int sideLen = 2*NEIGHBDIS + 1;
                    int colNum_mat = 0;
                    int rowNum_mat = 0;

		    //ostringstream regionIt_str;
		    //regionIt_str<<regionIt;
		    //string regionFile = "regionData_" + domainIt_str.str() + "_" + regionIt_str.str();
		    //ofstream outRegionData(regionFile);
                    vector< pair<int, int> > oriCoor;
                    for (int neighbRowIt=(rowIt-NEIGHBDIS), squareRowIt = 0; squareRowIt < sideLen; ++neighbRowIt, ++squareRowIt)
                    {
												cout <<"neighbRowIt\t"<<neighbRowIt<<endl;
                        for (int neighbColIt=colIt-NEIGHBDIS, squareColIt = 0; squareColIt < sideLen; ++neighbColIt, ++squareColIt)
                        {
														cout <<"neighbColIt\t"<<neighbColIt<<endl;
                            int distMatrix_colIt=convtIndex(neighbRowIt, neighbColIt, siteNum);
                            //if (distMatrix_colIt==-1 || neighbColIt<=neighbRowIt || colIt<=neighbRowIt || neighbColIt<=rowIt || (neighbRowIt==rowIt && neighbColIt==colIt) )
                            if (distMatrix_colIt==-1 || neighbColIt<=neighbRowIt || colIt<=neighbRowIt || neighbColIt<=rowIt)
                            {continue;}

                            //for future get back the original coordinates
                            if(freq_convert[distMatrix_colIt] > TESTFREQTHRES)
                            {
                                oriCoor.push_back(make_pair(neighbRowIt, neighbColIt));
                                ++colNum_mat;
                            }

                            neighbFreq.push_back(freq_convert[distMatrix_colIt]-constTerms[distMatrix_colIt]);
			    //outRegionData<<neighbRowIt<<"\t"<<neighbColIt<<"\t"<<distMatrix_colIt<<"\t"<<freq_convert[distMatrix_colIt]<<"\t"<<constTerms[distMatrix_colIt]<<"\t"<<mu_eff_fun<<"\t"<<beta_eff_fun<<"\t"<<effVec[distMatrix_colIt]<<"\t"<<sitesMap[rowIt]<<"\t"<<sitesMap[neighbRowIt]<<"\t"<<sitesMap[colIt]<<"\t"<<sitesMap[neighbColIt]<<"\t"<<mu_eff_fun<<"\t"<<beta_eff_fun<<"\t"<<calBgFrags(disPolyCoef_local, abs(sitesMap[rowIt]-sitesMap[neighbRowIt])+abs(sitesMap[colIt]-sitesMap[neighbColIt])+OPTDIST)<<endl;
                            //neighbConstTerms.push_back(constTerms[distMatrix_colIt]);
			    //
														cout<<"rowNum_mat\t"<<rowNum_mat<<endl;
                            fullExp.push_back(exp(effConsts[distMatrix_colIt]+calBgFrags(disPolyCoef_local, abs(sitesMap[rowIt]-sitesMap[neighbRowIt])+abs(sitesMap[colIt]-sitesMap[neighbColIt])+OPTDIST)));
                            ++rowNum_mat;
                        }
                    }
		    //exit(0);

                    //int neighbFreq_size = neighbFreq.size();
		    arma::vec neighbFreq_lm(neighbFreq);
				cout<<"neighbFreq.size\t"<<neighbFreq_lm.size()<<endl;
 
		    int fullExpMat_size = fullExp.size();

		    arma::mat fullExpMat_lm(1,fullExpMat_size);
		    //arma::mat fullExpMat_lm(fullExpMat_size,1);

		    for(int i=0; i<fullExpMat_size; ++i)
		    {
			//fullExpMat_lm(i,0)=1;
			//fullExpMat_lm(i,0)=fullExp[i];
			fullExpMat_lm(0,i)=fullExp[i];
						cout<<"i\t"<<i<<"\t"<<fullExpMat_lm(0,i)<<"\t"<<neighbFreq_lm[i]<<endl;
						//cout<<"i\t"<<i<<"\t"<<fullExpMat_lm(i,0)<<"\t"<<neighbFreq_lm[i]<<endl;
		    }
				cout<<"fullExpMat_lm.size\t"<<fullExpMat_lm.size()<<endl;

		    //mlpack::regression::LinearRegression lm(fullExpMat_lm, neighbFreq_lm,0.0,false);
		    mlpack::regression::LinearRegression lm(fullExpMat_lm, neighbFreq_lm,0.0,true);
		    cout<<"construct region model done."<<endl;

		    //set start values for theta
		    arma::vec& beta_lm = lm.Parameters();
		    cout << "The region MLEs are: " << endl;
		    std::cout <<"regionMLEs: "<<beta_lm(0)<<"\t"<<beta_lm(1)<<endl;

		    if(beta_lm(1)>0.01)
		    {
			//arma::vec temp = neighbFreq_lm - arma::trans( (arma::trans(beta_lm.subvec(1, beta_lm.n_elem - 1)) * fullExpMat_lm) + beta_lm(0));
			arma::vec temp = neighbFreq_lm - arma::trans((beta_lm(1) * fullExpMat_lm) + beta_lm(0));
			double sigmaSq = arma::dot(temp, temp) / (fullExpMat_size-2);//denominator: n-number of betas(parameters)

			//double var_err = lm.ComputeError(fullExpMat_lm, neighbFreq_lm);
			//cout<<"sigma squared: "<<sigmaSq<<endl;

			double sumX=0.0;
			double sumXsq=0.0;
			for(int i=0; i<fullExpMat_size; ++i)
			{
			    sumX+=fullExp[i];
			    sumXsq+=fullExp[i]*fullExp[i];
			}

			double beta_stderr = sqrt(fullExpMat_size*sigmaSq/(fullExpMat_size*sumXsq-sumX*sumX));
			cout<<"beta standard error: "<<beta_stderr<<endl;
			double fitPval=2 * gsl_cdf_gaussian_P(-fabs(beta_lm(1)/beta_stderr), 1.0);
			std::cout <<"regionMLEs: "<<beta_lm(0)<<"\t"<<beta_lm(1)<<"\t"<<fitPval<<"\t"<<colNum_mat<<endl;

			if (fitPval < PCUTOFF && colNum_mat > 1)
			{
			    testPosP<<domainIt<<"\t"<<regionIt<<"\t"<<testIt<<"\t"<<sitesMap[rowIt]<<"\t"<<sitesMap[colIt]<<"\t"<<beta_lm(1)<<"\t"<<fitPval<<"\n";

			    int distNum = colNum_mat;
			    vector<double> distVec_lasso(colNum_mat, 0.0);
			    vector< vector<double> > distMatrix_lasso;
			    for (int i=0; i<rowNum_mat; ++i)
			    {
				distMatrix_lasso.push_back(distVec_lasso);
			    }

			    int testIndex = 0;
			    int matColIndex = 0;
			    for (int neighbRowIt=(rowIt-NEIGHBDIS), squareRowIt = 0; squareRowIt < sideLen; ++neighbRowIt, ++squareRowIt)
			    {
				for (int neighbColIt=colIt-NEIGHBDIS, squareColIt = 0; squareColIt < sideLen; ++neighbColIt, ++squareColIt)
				{
				    int distMatrix_colIt=convtIndex(neighbRowIt, neighbColIt, siteNum);
				    //if (distMatrix_colIt==-1 || neighbColIt<=neighbRowIt || colIt<=neighbRowIt || neighbColIt<=rowIt || (neighbRowIt==rowIt && neighbColIt==colIt) )
				    if (distMatrix_colIt==-1 || neighbColIt<=neighbRowIt || colIt<=neighbRowIt || neighbColIt<=rowIt || freq_convert[distMatrix_colIt] <= TESTFREQTHRES)
				    //if (distMatrix_colIt==-1 || neighbColIt<=neighbRowIt || colIt<=neighbRowIt || neighbColIt<=rowIt || freq_convert[distMatrix_colIt] <= TESTFREQTHRES || abs(sitesMap[neighbRowIt]-sitesMap[neighbColIt])<20000)
				    {continue;}

				    if (neighbRowIt==rowIt && neighbColIt==colIt)
				    {
					testIndex=matColIndex;
					testPosP<<sitesMap[neighbRowIt]<<"\t"<<sitesMap[neighbColIt]<<"\t"<<testIndex<<"\n";
				    }
				    else
				    {
					testPosP<<sitesMap[neighbRowIt]<<"\t"<<sitesMap[neighbColIt]<<"\n";
				    }
				    ++matColIndex;
				    //distMatrix[make_pair(distMatrix_rowIt, distMatrix_colIt)]=exp(-MORANI*(abs(sitesMap[rowIt]-sitesMap[neighbRowIt])+abs(sitesMap[colIt]-sitesMap[neighbColIt])));
				}
			    }

			    regionIt++;
			    testIt++;

			    distMatrix_lasso.clear();
			    ++counter_fit;
			}
			else if(fitPval < PCUTOFF && colNum_mat==1)
			{
			    oneColTest<<sitesMap[rowIt]<<"\t"<<sitesMap[colIt]<<"\t"<<beta_lm(1)<<"\t"<<fitPval<<"\n";
			    //only 1 box has signal (freq>1), no need to use lasso to choose which one is real,or maybe could remove the ones due to the const terms, however glm should have taken care of that already
			}
		    }
                }
                else
                {
                    //#pragma omp critical
                    {++counter_skip;}
                    //#pragma omp end critical
                }
            }
        }
	oneColTest.close();
	testPosP.close();
	//outRegionData.close();
	//outDistMatrix.close();
    }
    //outDebug.close();
    
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
    return 1;
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


/*
nnlasso.normal<-function(x,y,lambda=NULL,intercept=TRUE,normalize=TRUE,tau=1,tol=1e-6,maxiter=1e5,nstep=100,min.lambda=1e-4,eps=1e-6,path=TRUE,SE=FALSE)
{
        np=dim(x)
        n=np[1]
        p=np[2]
        if (intercept)
        {
                meanx = colMeans(x)
                x = scale(x, meanx, FALSE)
                meany = sum(y)/n
                y = y - meany
        } else {
                meanx = rep(0, p)
                meany = 0
                }
        if (normalize)
        {
                normx = sqrt(colSums(x^2))
                x = scale(x, FALSE, normx)
        } else normx = rep(1, p)
        tx=t(x)
        xpy=tx%*%y
        max.lambda=max(abs(xpy))
        if (path)
        {
                stepsize=exp((log(min.lambda)-log(max.lambda))/nstep)
                lambdas=max.lambda*stepsize^((1:nstep)-1)
        } else {
                        nstep=2
                        lambdas=c(max.lambda,lambda)
                }
        coef=matrix(0,nstep,p)
        of.value=rep(0,nstep)
        coef[1,]=1e-2
        xbeta.old=x%*%coef[1,]
        if (n>p)
        {
                xpx=tx%*%x
                xpxbetaold=xpx%*%coef[1,]
        } else {
                        xpx=NULL
                        xpxbetaold=tx%*%xbeta.old
                }
        of.value[1]=-sum((y-xbeta.old)^2)/2-max.lambda*(tau*sum(coef[1,])+(1-tau)*sum(coef[1,]^2))
        lambda.iter=rep(0,nstep)
        g1=xpy
        for(iter in 2:nstep)
        {
                kkt=FALSE
                while(kkt==FALSE)
                {
                        res=nnlasso.normal.lambda(n,p,x,y,xpx,xpy,beta.old=coef[iter-1,],tau,lambda1=lambdas[iter],tol,maxiter,xbeta.old,eps,SE)
                        if (res$conv=="yes")
                        {
                                coef[iter,]=res$beta.new
                                xbeta=res$xbeta.new
                                if (n>p)
                                {
                                        xpxbetaold=xpx%*%res$beta.new
                                } else xpxbetaold=tx%*%xbeta
                                g1=xpy-xpxbetaold-lambdas[iter]*2*(1-tau)*coef[iter,]
                                indices=NULL
                                if (length(indices)==0)
                                {
                                        xbeta.old=res$xbeta.new
                                        lambda.iter[iter]=res$iter
                                        of.value[iter]=res$ofv.new
                                        kkt=TRUE
                                }
                        } else stop("The algorithm did not converge")
                }
        }
        coef=scale(coef,center=FALSE,scale=normx)
        if(SE)
        {
                vcov=res$vcov
                vcov=vcov/normx
                vcov=t(vcov)/normx
                se=sqrt(diag(vcov))
                if(intercept)
                {
                        se0=sqrt(t(meanx)%*%vcov%*%meanx)
                        se=c(se0,se)
                }
        } else se=NULL
        if (intercept) beta0=rep(meany,nstep)-coef%*%meanx else beta0=rep(0,nstep)
        L1norm=rowSums(abs(coef))
        norm.frac=L1norm/max(L1norm)
        obj=list(beta0=beta0,coef=coef,lambdas=lambdas,L1norm=L1norm,norm.frac=norm.frac,lambda.iter=lambda.iter,of.value=of.value,normx=normx,se=se)
        class(obj)='nnlasso'
        return(obj)
}

nnlasso.normal.lambda<-function(n,p,x,y,xpx,xpy,beta.old,tau,lambda1,tol,maxiter,xbeta.old,eps,SE=FALSE)
{
        epp=0.001
        if (n<=p) tx=t(x)
        ofv.old=-sum((y-xbeta.old)^2)/2-lambda1*(tau*sum(beta.old)+(1-tau)*sum(beta.old^2))
        for (iter in 1:maxiter)
        {
                if (n>p)
                {
                        xpxbetaold=xpx%*%beta.old
                } else xpxbetaold=tx%*%xbeta.old
                g1minus<-g1<-xpy-xpxbetaold
                g1minus[which(g1>0)]=0
                b=beta.old*(g1-lambda1*tau-2*lambda1*(1-tau)*beta.old)/(lambda1*tau+2*lambda1*(1-tau)*beta.old-g1minus+1e-16)
                beta.new=beta.old+b
                xb=x%*%b
                xbeta.new=xbeta.old+xb
                ofv.new=-sum((y-xbeta.new)^2)/2-lambda1*(tau*sum(beta.new)+(1-tau)*sum(beta.new^2))
                delta=1
                t1=epp*(sum(g1*b))
                while (ofv.new-delta*t1<ofv.old & delta>1e-5)
                {
                        delta=delta/2
                        beta.new=beta.old+delta*b
                        xbeta.new=xbeta.old+delta*xb
                        ofv.new=-sum((y-xbeta.new)^2)/2-lambda1*(tau*sum(beta.new)+(1-tau)*sum(beta.new^2))
                }
                if (ofv.new-delta*t1<ofv.old & delta<=1e-5)
                {
                        beta.new=beta.old
                        ofv.new=ofv.old
                        xbeta.new=xbeta.old
                        break
                }
                if(abs(ofv.old-ofv.new)<=tol) break
                beta.old=beta.new
                xbeta.old=xbeta.new
                ofv.old=ofv.new
        }
        if (iter<maxiter) conv="yes" else conv="no"
        if (SE==TRUE & conv=="yes")
        {
                grad=xpy-t(x)%*%xbeta.new-lambda1*tau-2*lambda1*(1-tau)*beta.new
                Finv<-F<-xpx+2*lambda1*(1-tau)
                index=which(beta.new<=eps & grad< -1e-2)
                if (length(index)>0)
                {
                        temp=F[-index,-index]
                        tempinv=solve(temp)
                        indexc=setdiff(1:p,index)
                        if (length(indexc)>0) Finv[indexc,indexc]=tempinv
                        Finv[index,index]=0
                } else Finv=solve(F)
                vc=Finv%*%F
                vcov=vc%*%t(Finv)
        } else vcov=NULL
        res=list(beta.new=beta.new,conv=conv,iter=iter,ofv.new=ofv.new,xbeta.new=xbeta.new,vcov=vcov)
        return(res)
}
*/

#endif
