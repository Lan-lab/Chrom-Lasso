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




int findIntraDomainInteraction(map<int, int>& empDisMap_fun, string jobChr_fun, vector< vector<int> >& domainSitesMap_fun, vector< vector< vector<int> > >& domainCSinterFreq_fun, map<int, int>& csInterChromTotalMap_fun, map<int, double>& disPolyCoef_local, map< int, map< pair<int, int>, double > >& pVal_fun, double corrThres_fun, double disIntcp_fun, double disCoeff_fun)
{
    string file_debug =  jobChr_fun + "_debug";
    ofstream outDebug(file_debug.c_str());

    int domainNum = domainSitesMap_fun.size();

    //#pragma omp parallel for
    for (int domainIt = 0; domainIt < domainNum; ++domainIt)
    {
        map< pair<int, int>, double>& domainPval = pVal_fun[domainIt];
        outDebug<<"domain "<<domainIt<<endl;
        vector<int>& sitesMap = domainSitesMap_fun[domainIt];
        vector< vector<int> >& csInterFreqMap = domainCSinterFreq_fun[domainIt];
        int siteNum = sitesMap.size();

        int size_convert = convtIndex(siteNum-2,siteNum-1,siteNum) + 1;//the convtIndex function's first two input is 0 based and the output is 0 based too
	vector<double> freq_convert;
        vector<double> offset;//distance effect
        vector<double> effVec;//efficiency of cutting site and mappability, measured by the total inter chromosome hybrid frags

	ofstream outDomainData("domainData");

        for (int rowIt=0; rowIt<siteNum; ++rowIt)
        {
            for (int colIt=rowIt+1; colIt<siteNum; ++colIt)
            {
                int index_convert=convtIndex(rowIt,colIt,siteNum);

		if(abs(sitesMap[rowIt]-sitesMap[colIt])>OPTDIST)//the index of freq_convert need to be changed because of this
		{
		    //log(y/ydist)=mu+betaeff*effvec
		    freq_convert.push_back(csInterFreqMap[rowIt][colIt]);
		    offset.push_back(calBgFrags(disPolyCoef_local, abs(sitesMap[rowIt]-sitesMap[colIt])));
		    //the digestion and ligation efficient of a cutting site should take log because, its log should proportional to log(freq)
		    effVec.push_back(log(sqrt(csInterChromTotalMap_fun[sitesMap[rowIt]])*sqrt(csInterChromTotalMap_fun[sitesMap[colIt]])+1));
		    outDomainData<<rowIt<<"\t"<<colIt<<"\t"<<csInterFreqMap[rowIt][colIt]<<"\t"<<csInterChromTotalMap_fun[sitesMap[rowIt]]<<"\t"<<csInterChromTotalMap_fun[sitesMap[colIt]]<<"\t"<<log(sqrt(csInterChromTotalMap_fun[sitesMap[rowIt]])*sqrt(csInterChromTotalMap_fun[sitesMap[colIt]])+1)<<"\t"<<sitesMap[rowIt]<<"\t"<<sitesMap[colIt]<<"\t"<<calBgFrags(disPolyCoef_local, abs(sitesMap[rowIt]-sitesMap[colIt])+OPTDIST)<<endl;
		}
            }
        }
	outDomainData.close();

	int size_convert_s = freq_convert.size();
	scythe::Matrix<double> freq_convert_s(size_convert_s, 1);
        scythe::Matrix<double> offset_s(size_convert_s,1);//distance effect
        scythe::Matrix<double> effVec_s(size_convert_s,2);//efficiency of cutting site and mappability, measured by the total inter chromosome hybrid frags

	ofstream outFreq("freq");
	ofstream outOffset("offset");
	ofstream outEffvec("effvec");
	for(int i=0;i<size_convert_s;++i)
	{
	    freq_convert_s[i,0]=freq_convert[i];
	    offset_s[i,0]=offset[i];
	    effVec_s[i,0]=1;
	    effVec_s[i,1]=effVec[i];
	    outFreq<<freq_convert_s[i,0]<<endl;
	    outOffset<<offset_s[i,0]<<endl;
	    outEffvec<<effVec_s[i,0]<<"\t"<<effVec_s[i,1]<<endl;
	}
	cout<<"construct matrix done."<<endl;
	outFreq.close();
	outOffset.close();
	outEffvec.close();

	//vector< vector<double> > covMatrix;
	//covMatrix.push_back(effVec);

	//IRLS irls("log-link");
	//bool quasi_lik = false;
	//irls.link->quasi = quasi_lik;
	//irls.load_data(freq_convert, covMatrix, offset);
	////outDebug<<"data loaded"<<endl;
	//irls.fit_model();
	////outDebug<<"model fitted"<<endl;
	//vector<double> coev = irls.get_coef();
	////outDebug<<"get coef done"<<endl;
	//vector<double> sev = irls.get_stderr();
	////outDebug<<"get stderr done"<<endl;

	//outDebug<<"Col\tEstimate\tStd.Error\tp-value"<<endl;
	//for(size_t i = 0; i < coev.size(); ++i)
	//{
	//    //printf("X%-9zu%12.9f%12.8f", i, coev[i], sev[i]);
	//    //X0 is intercept/mu
	//    outDebug<<"X"<<i<<"\t"<<coev[i]<<"\t"<<sev[i]<<"\t";
	//    if(! irls.link->quasi)
	//    {
	//	//printf("%15.6e\n", 2 * gsl_cdf_gaussian_P(-fabs(coev[i]/sev[i]), 1.0));
	//	outDebug<<2 * gsl_cdf_gaussian_P(-fabs(coev[i]/sev[i]), 1.0)<<endl;
	//    }
	//    else
	//    {
	//	//printf("%15.6e\n", 2 * gsl_cdf_tdist_P(-fabs(coev[i]/sev[i]), size_convert-irls.get_rank_X()));
	//	outDebug<<2 * gsl_cdf_tdist_P(-fabs(coev[i]/sev[i]), size_convert-irls.get_rank_X())<<endl;
	//    }
	//}
	//covMatrix.clear();
	//outDebug.close();

	PoissonModel poisson_model;
	poisson_model.y_ = freq_convert_s;
	poisson_model.X_ = effVec_s;
	poisson_model.offset_ = offset_s;
	scythe::mersenne myrng;
	cout<<"construt model done."<<endl;

	//set start values for theta
	scythe::Matrix<double> theta = invpd(crossprod(effVec_s+0.01)) * t(effVec_s+0.01) * scythe::log(freq_convert_s+1);
	//scythe::Matrix<double> theta(2,1);
	//theta(0,0)=1.0;
	//theta(1,0)=6.0;
	cout<<"theta initialization done, "<<theta<<endl;
	scythe::Matrix<double> beta_MLE = BFGS(poisson_model, theta, myrng, 100, 1e-5, true);
	cout<<"MLE done."<<endl;

	cout << "The MLEs are: " << endl;
	std::cout << t(beta_MLE) << "\n";


	double mu_eff = beta_MLE[0];
	//double alpha = coev[1];
	double beta_eff = beta_MLE[1];
	vector<double> constTerms(size_convert, 0.0);
        vector<double> effConsts(size_convert, 0.0);//efficiency constand mu+beta_eff*effVec[i]
	for (int i = 0; i < size_convert; ++i)
	{
	    //constTerms[i] = mu+dist_ij[i]*alpha+effVec[i]*beta;
	    effConsts[i] = mu_eff+effVec[i]*beta_eff;
	    constTerms[i] = exp(effConsts[i]+offset[i]);
	    //constTerms[i] = mu+effVec[i]*beta+offset[i];
	}

        int counter_fit=0;
        int counter_skip=0;

        //#pragma omp parallel for
        for (int rowIt=0; rowIt<siteNum; ++rowIt)
        {
            for (int colIt=rowIt+1; colIt<siteNum; ++colIt)
            {
                if(csInterFreqMap[rowIt][colIt]>TESTFREQTHRES && abs(sitesMap[rowIt]-sitesMap[colIt])>=20000)
                {
                    vector<double> neighbFreq;
                    vector<double> neighbConstTerms;
                    vector<double> fullExp;

                    int sideLen = 2*NEIGHBDIS + 1;
                    int colNum_mat = 0;
                    int rowNum_mat = 0;

		    ofstream outRegionData("regionData");

                    vector< pair<int, int> > oriCoor;
                    for (int neighbRowIt=(rowIt-NEIGHBDIS), squareRowIt = 0; squareRowIt < sideLen; ++neighbRowIt, ++squareRowIt)
                    {
                        for (int neighbColIt=colIt-NEIGHBDIS, squareColIt = 0; squareColIt < sideLen; ++neighbColIt, ++squareColIt)
                        {
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
			    outRegionData<<neighbRowIt<<"\t"<<neighbColIt<<"\t"<<distMatrix_colIt<<"\t"<<freq_convert[distMatrix_colIt]<<"\t"<<constTerms[distMatrix_colIt]<<"\t"<<mu_eff<<"\t"<<beta_eff<<"\t"<<effVec[distMatrix_colIt]<<"\t"<<sitesMap[rowIt]<<"\t"<<sitesMap[neighbRowIt]<<"\t"<<sitesMap[colIt]<<"\t"<<sitesMap[neighbColIt]<<"\t"<<disCoeff_fun<<"\t"<<disIntcp_fun<<"\t"<<calBgFrags(disPolyCoef_local, abs(sitesMap[rowIt]-sitesMap[neighbRowIt])+abs(sitesMap[colIt]-sitesMap[neighbColIt])+OPTDIST)<<endl;
                            //neighbConstTerms.push_back(constTerms[distMatrix_colIt]);
			    //
                            fullExp.push_back(exp(effConsts[distMatrix_colIt]+calBgFrags(disPolyCoef_local, abs(sitesMap[rowIt]-sitesMap[neighbRowIt])+abs(sitesMap[colIt]-sitesMap[neighbColIt])+OPTDIST)));
                            ++rowNum_mat;
                        }
                    }
		    outRegionData.close();
		    //exit(0);

                    //int neighbFreq_size = neighbFreq.size();
		    arma::vec neighbFreq_lm(neighbFreq);
 
		    int fullExpMat_size = fullExp.size();

		    arma::mat fullExpMat_lm(1,fullExpMat_size);

		    for(int i=0; i<fullExpMat_size; ++i)
		    {
			//fullExpMat_lm(i,0)=1;
			fullExpMat_lm(0,i)=fullExp[i];
		    }

		    mlpack::regression::LinearRegression lm(fullExpMat_lm, neighbFreq_lm,0.0,true);
		    //cout<<"construct region model done."<<endl;

		    //set start values for theta
		    arma::vec& beta_lm = lm.Parameters();
		    //cout << "The region MLEs are: " << endl;
		    //std::cout <<"regionMLEs: "<<beta_lm(0)<<"\t"<<beta_lm(1)<<endl;

		    if(beta_lm(1)>0.01){
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
		    //cout<<"beta standard error: "<<beta_stderr<<endl;
		    double fitPval=2 * gsl_cdf_gaussian_P(-fabs(beta_lm(1)/beta_stderr), 1.0);
		    std::cout <<"regionMLEs: "<<beta_lm(0)<<"\t"<<beta_lm(1)<<"\t"<<fitPval<<endl;

                    if (fitPval < PCUTOFF && colNum_mat > 1)
                    {
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
                                if (distMatrix_colIt==-1 || neighbColIt<=neighbRowIt || colIt<=neighbRowIt || neighbColIt<=rowIt || freq_convert[distMatrix_colIt] <= TESTFREQTHRES || abs(sitesMap[neighbRowIt]-sitesMap[neighbColIt])<20000)
                                {continue;}

                                int matRowIndex = 0;//inside loop is actually the row index for the final matrix
                                for (int neighbNeighbRowIt=(rowIt-NEIGHBDIS), squareNeighbRowIt = 0; squareNeighbRowIt < sideLen; ++neighbNeighbRowIt, ++squareNeighbRowIt)
                                {
                                    for (int neighbNeighbColIt=colIt-NEIGHBDIS, squareNeighbColIt = 0; squareNeighbColIt < sideLen; ++neighbNeighbColIt, ++squareNeighbColIt)
                                    {
                                        int distMatrix_rowIt=convtIndex(neighbNeighbRowIt, neighbNeighbColIt, siteNum);
                                        if (distMatrix_rowIt==-1 || neighbNeighbColIt<=neighbNeighbRowIt || colIt<=neighbNeighbRowIt || neighbNeighbColIt<=rowIt)
                                        {continue;}
                                        //distMatrix_lasso[matRowIndex][matColIndex]=exp(-MORANI*(abs(sitesMap[neighbNeighbRowIt]-sitesMap[neighbRowIt])+abs(sitesMap[neighbNeighbColIt]-sitesMap[neighbColIt])));
                                        distMatrix_lasso[matRowIndex][matColIndex]=exp(effConsts[distMatrix_rowIt]+calBgFrags(disPolyCoef_local, abs(sitesMap[neighbNeighbRowIt]-sitesMap[neighbRowIt])+abs(sitesMap[neighbNeighbColIt]-sitesMap[neighbColIt])+OPTDIST));
                                        ++matRowIndex;
                                    }
                                }
                                if (neighbRowIt==rowIt && neighbColIt==colIt)
                                {
                                    testIndex=matColIndex;
                                }
                                ++matColIndex;
                                //distMatrix[make_pair(distMatrix_rowIt, distMatrix_colIt)]=exp(-MORANI*(abs(sitesMap[rowIt]-sitesMap[neighbRowIt])+abs(sitesMap[colIt]-sitesMap[neighbColIt])));
                            }
                        }

                        int testFlag = 1;
                        bool flag_rmDep = false;

                        if (testFlag == 1)
                        {
                            //#pragma omp critical
                            {
                            //cout<<rowIt<<"\t"<<colIt<<"\ttest column number: "<<testIndex<<endl;

                            const int rowNum_r = distMatrix_lasso.size();
                            const int colNum_r = (distMatrix_lasso.begin())->size();
                            //cout<<rowNum_r<<"\t"<<colNum_r<<endl;
                //!!!!!!!!!!! if remove dependents is true, this need to change
		//
			    arma::mat fullExpMat_lars(colNum_r,rowNum_r);

			    ofstream outRegionData("distMatrix");
			    for(int i=0; i<colNum_r; ++i)
			    {
				for(int j=0; j<rowNum_r; ++j)
				{
				    fullExpMat_lars(i,j)=distMatrix_lasso[j][i];
				    outRegionData<<distMatrix_lasso[j][i]<<"\t";
				}
				outRegionData<<endl;
			    }
			    outRegionData.close();
			    exit(0);

			    distMatrix_lasso.clear();
			    //cout<<"construct lars exp matrix done."<<endl;

			    mlpack::regression::LARS lars_lm(true);
			    //cout<<"construct lars object done."<<endl;

			    arma::vec lars_betas;
			    lars_lm.Regress(fullExpMat_lars, neighbFreq_lm, lars_betas);
			    //cout<<"regression done."<<endl;
			    cout<<"lars betas"<<endl;
			    cout<<lars_betas<<endl;
			    cout<<"lars lambdas";
			    vector<double> lambdaP=lars_lm.LambdaPath();
			    for(int i=0;i<lambdaP.size();++i)
			    {
				cout<<"\t"<<lambdaP[i];
			    }
			    cout<<endl;

			    cout<<"lars beta path\n";
			    vector<arma::vec> betaP=lars_lm.BetaPath();
			    for(int i=0;i<betaP.size();++i)
			    {
				cout<<i<<"\n"<<betaP[i];
			    }
			    cout<<endl;

                            //#pragma omp critical
                            //domainPval[make_pair(sitesMap[rowIt], sitesMap[colIt])]=fitPval;
                            ++counter_fit;}//}
                            //#pragma omp end critical

                            //if (fitPval < PCUTOFF)
                            //{
                            //    
                            //}
                        }
                    }
                    else if(colNum_mat==1)
                    {
                        //only 1 box has signal (freq>1), no need to use lasso to choose which one is real,or maybe could remove the ones due to the const terms, however glm should have taken care of that already
                    }
		    }
                }
                else
                {
                    #pragma omp critical
                    {++counter_skip;}
                    //#pragma omp end critical
                }
            }
        }
        //freq_convert.clear();
        //distMatrix.clear();
        //dist_ij.clear();
        //offsets.clear();
        //effVec.clear();
        //outDebug<<"chr "<<jobChr_fun<<" domain "<<domainIt<<" siteNum "<<siteNum<<" convert size "<<size_convert<<"fit number "<<counter_fit<<" skip number "<<counter_skip<<endl;
    }
    outDebug.close();
    
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




#endif
