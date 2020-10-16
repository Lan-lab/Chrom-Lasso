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
int findIntraDomainInteraction(map<int, int>&, string, vector< vector<int> >&, vector< vector< vector<int> > >&, map<int, int>&, map< int, map< pair<int, int>, double > >&, double, double, double);
int rmDependentCol(vector< vector<double> >&, vector<double>&, vector< pair<int, int> >&, double, int&, vector< vector<double> >&);
int calStanScores(vector< vector<double> >&);
int calPearsonCorr(vector< vector<double> >&, vector< pair<int, int> >&, const double&, int, map<int, int>&);
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
    if (argc < 10)
    {
        usage();
        exit(1);
    }

    string empDisFile = argv[1];
    string ictmFile = argv[2];
    string dsmFile = argv[3];
    string difFile = argv[4];
    string chrNum = argv[5];
    double MORANI = atof(argv[6]);
    double corrThres = atof(argv[7]);
    double disIntcp = atof(argv[8]);
    double disCoeff = atof(argv[9]);


    map<int, int> empDis_map;
    map<int, int> ictm_map;
    vector< vector<int> > dsm_map;
    vector< vector< vector<int> > > dif_map;

    //reader_ictm(empDisFile, empDis_map);
    //cout<<"\nReading empDis file done."<<endl;
    reader_ictm(ictmFile, ictm_map);
    cout<<"\nReading ictm file done."<<endl;
    reader_dsm(dsmFile, dsm_map);
    cout<<"\nReading dsm file done."<<endl;
    reader_dif(difFile, dif_map);
    cout<<"\nReading dif file done."<<endl;

    map< int, map< pair<int, int>, double > > pVal_map;

    findIntraDomainInteraction(empDis_map, chrNum, dsm_map, dif_map, ictm_map, pVal_map, corrThres, disIntcp, disCoeff);

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




int findIntraDomainInteraction(map<int, int>& empDisMap_fun, string jobChr_fun, vector< vector<int> >& domainSitesMap_fun, vector< vector< vector<int> > >& domainCSinterFreq_fun, map<int, int>& csInterChromTotalMap_fun, map< int, map< pair<int, int>, double > >& pVal_fun, double corrThres_fun, double disIntcp_fun, double disCoeff_fun)
{
    string file_effBetas =  jobChr_fun + "_domainEffBetas";
    ofstream outEffBetas(file_effBetas.c_str());

    int domainNum = domainSitesMap_fun.size();

    //#pragma omp parallel for
    for (int domainIt = 0; domainIt < domainNum; ++domainIt)
    {
        map< pair<int, int>, double>& domainPval = pVal_fun[domainIt];
        outEffBetas<<"domain "<<domainIt;
        vector<int>& sitesMap = domainSitesMap_fun[domainIt];
        vector< vector<int> >& csInterFreqMap = domainCSinterFreq_fun[domainIt];
        int siteNum = sitesMap.size();

        int size_convert = convtIndex(siteNum-2,siteNum-1,siteNum) + 1;//the convtIndex function's first two input is 0 based and the output is 0 based too
	vector<double> freq_convert;
        vector<double> offset;//distance effect
        vector<double> effVec;//efficiency of cutting site and mappability, measured by the total inter chromosome hybrid frags

        for (int rowIt=0; rowIt<siteNum; ++rowIt)
        {
            for (int colIt=rowIt+1; colIt<siteNum; ++colIt)
            {
                int index_convert=convtIndex(rowIt,colIt,siteNum);

		//log(y/ydist)=mu+betaeff*effvec
                freq_convert.push_back(csInterFreqMap[rowIt][colIt]);
		offset.push_back(log(abs(sitesMap[rowIt]-sitesMap[colIt])+1)*disCoeff_fun+disIntcp_fun);
		//the digestion and ligation efficient of a cutting site should take log because, its log should proportional to log(freq)
                effVec.push_back(log(sqrt(csInterChromTotalMap_fun[sitesMap[rowIt]])*sqrt(csInterChromTotalMap_fun[sitesMap[colIt]])+1));
            }
        }
	cout<<"construt matrix done."<<endl;

	vector< vector<double> > covMatrix;
	covMatrix.push_back(effVec);

	IRLS irls("log-link");
	bool quasi_lik = false;
	irls.link->quasi = quasi_lik;
	irls.load_data(freq_convert, covMatrix, offset);
	//outEffBetas<<"data loaded"<<endl;
	irls.fit_model();
	//outEffBetas<<"model fitted"<<endl;
	vector<double> coev = irls.get_coef();
	//outEffBetas<<"get coef done"<<endl;
	vector<double> sev = irls.get_stderr();
	//outEffBetas<<"get stderr done"<<endl;

	//outEffBetas<<"Col\tEstimate\tStd.Error\tp-value"<<endl;
	for(size_t i = 0; i < coev.size(); ++i)
	{
	    //printf("X%-9zu%12.9f%12.8f", i, coev[i], sev[i]);
	    //X0 is intercept/mu
	    double p_eff=0.0;
	    if(! irls.link->quasi)
	    {
		//printf("%15.6e\n", 2 * gsl_cdf_gaussian_P(-fabs(coev[i]/sev[i]), 1.0));
		p_eff=2 * gsl_cdf_gaussian_P(-fabs(coev[i]/sev[i]), 1.0);
	    }
	    else
	    {
		//printf("%15.6e\n", 2 * gsl_cdf_tdist_P(-fabs(coev[i]/sev[i]), size_convert-irls.get_rank_X()));
		p_eff=2 * gsl_cdf_tdist_P(-fabs(coev[i]/sev[i]), size_convert-irls.get_rank_X());
	    }
	    outEffBetas<<"\t"<<coev[i]<<"\t"<<sev[i]<<"\t"<<p_eff;
	}
	outEffBetas<<endl;
	covMatrix.clear();
	//outEffBetas.close();
	continue;

	double mu = coev[0];
	//double alpha = coev[1];
	double beta = coev[1];
	vector<double> constTerms(size_convert, 0.0);
	for (int i = 0; i < size_convert; ++i)
	{
	    //constTerms[i] = mu+dist_ij[i]*alpha+effVec[i]*beta;
	    constTerms[i] = exp(mu+effVec[i]*beta+offset[i]);
	    //constTerms[i] = mu+effVec[i]*beta+offset[i];
	}

        int counter_fit=0;
        int counter_skip=0;

        //#pragma omp parallel for
        for (int rowIt=0; rowIt<siteNum; ++rowIt)
        {
            for (int colIt=rowIt+1; colIt<siteNum; ++colIt)
            {
                if(csInterFreqMap[rowIt][colIt]>TESTFREQTHRES)
                {
                    vector<double> neighbFreq;
                    vector<double> neighbConstTerms;
                    vector<double> distVec;

                    int sideLen = 2*NEIGHBDIS + 1;
                    int colNum_mat = 0;
                    int rowNum_mat = 0;

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
                            //neighbConstTerms.push_back(constTerms[distMatrix_colIt]);
                            distVec.push_back(log(abs(sitesMap[rowIt]-sitesMap[neighbRowIt])+abs(sitesMap[colIt]-sitesMap[neighbColIt])+1));
                            ++rowNum_mat;
                        }
                    }

		    IRLS irls("log-link");
		    bool quasi_lik = false;
		    irls.link->quasi = quasi_lik;
		    irls.load_data(freq_convert, covMatrix, offset);
		    //outEffBetas<<"data loaded"<<endl;
		    irls.fit_model();
		    //outEffBetas<<"model fitted"<<endl;
		    vector<double> coev = irls.get_coef();
		    //outEffBetas<<"get coef done"<<endl;
		    vector<double> sev = irls.get_stderr();
		    //outEffBetas<<"get stderr done"<<endl;

		    outEffBetas<<"Col\tEstimate\tStd.Error\tp-value"<<endl;
		    for(size_t i = 0; i < coev.size(); ++i)
		    {
			//printf("X%-9zu%12.9f%12.8f", i, coev[i], sev[i]);
			//X0 is intercept/mu
			outEffBetas<<"X"<<i<<"\t"<<coev[i]<<"\t"<<sev[i]<<"\t";
			if(! irls.link->quasi)
			{
			    //printf("%15.6e\n", 2 * gsl_cdf_gaussian_P(-fabs(coev[i]/sev[i]), 1.0));
			    outEffBetas<<2 * gsl_cdf_gaussian_P(-fabs(coev[i]/sev[i]), 1.0)<<endl;
			}
			else
			{
			    //printf("%15.6e\n", 2 * gsl_cdf_tdist_P(-fabs(coev[i]/sev[i]), size_convert-irls.get_rank_X()));
			    outEffBetas<<2 * gsl_cdf_tdist_P(-fabs(coev[i]/sev[i]), size_convert-irls.get_rank_X())<<endl;
			}
		    }
		    covMatrix.clear();
		    //outEffBetas.close();


                    //int neighbFreq_size = neighbFreq.size();
            arma::vec neighbFreq_lm(neighbFreq);
            //arma::vec neighbFreq_lm(neighbFreq_size);
            //for(int i=0; i<neighbFreq_size; ++i)
            //{
            //    neighbFreq_lm(i)=neighbFreq[i]-neighbConstTerms[i];
            //}

                    int distMat_size = distVec.size();
            arma::mat distMat_lm(distMat_size,2);
            for(int i=0; i<distMat_size; ++i)
            {
            distMat_lm(i,0)=1;
            distMat_lm(i,1)=distVec[i];
            }

            mlpack::regression::LinearRegression lm(distMat_lm, neighbFreq_lm);
            cout<<"construt region model done."<<endl;

            //set start values for theta
            arma::vec beta_lm = lm.Parameters();
            cout << "The region MLEs are: " << endl;
            std::cout <<beta_lm(0)<<"\t"<<beta_lm(1)<< "\n";

            //calculate std err using var(beta)=Var(error)*inv(t(X)X)
            //arma::mat<double> distMat_t = distMat_lm.t();
            //arma::mat<double> distMat_prod = distMat_t * distMat_lm;
            double var_err = lm.ComputeError(distMat_lm, neighbFreq_lm);
            arma::vec beta_stderr = sqrt(var_err * (distMat_lm.t()*distMat_lm).i());
            cout << "The stderr for betas are: " << endl;
            std::cout <<beta_stderr(0)<<"\t"<<beta_stderr(1)<< "\n";

            double fitPval=0.0;
            double beta_dis=0.0;
            exit(0);



                    if (fitPval < PCUTOFF && beta_dis < 0 && colNum_mat > 1)
                    {
                        int distNum = colNum_mat;
                        //colNum_mat++;//add const to the last column of the dist matrix
                        //colNum_mat;//add const to the last column of the dist matrix after rmIndependent cols is done
                        //vector<double> distVec_lasso(colNum_mat, 0.0);
                        vector<double> distVec_lasso((colNum_mat+1), 0.0);
                        vector< vector<double> > distMatrix_lasso;
                        for (int i=0; i<rowNum_mat; ++i)
                        {
                            distVec_lasso[colNum_mat] = neighbConstTerms[i];
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
                                if (distMatrix_colIt==-1 || neighbColIt<=neighbRowIt || colIt<=neighbRowIt || neighbColIt<=rowIt || freq_convert[distMatrix_colIt] < TESTFREQTHRES)
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
                                        distMatrix_lasso[matRowIndex][matColIndex]=log(abs(sitesMap[neighbNeighbRowIt]-sitesMap[neighbRowIt])+abs(sitesMap[neighbNeighbColIt]-sitesMap[neighbColIt])+1);
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

                        //vector< vector<double> > distMatrix_lasso_rmDep;

                        int testFlag = 1;
                        bool flag_rmDep = false;
                        //if (flag_rmDep)
                        //{
                        //    testFlag = rmDependentCol(distMatrix_lasso, neighbFreq, oriCoor, corrThres_fun, testIndex, distMatrix_lasso_rmDep);
                        //}
                        //else
                        //{
                        //    distMatrix_lasso_rmDep = distMatrix_lasso;
                        //}
            //distMatrix_lasso.clear();

                        //add const to the matrix
                        if (testFlag == 1)
                        {
                            //for (int i=0; i<rowNum_mat; ++i)
                            //{
                            //    distMatrix_lasso[i].push_back(neighbConstTerms[i]);
                            //}

/*Note Y need to substract offset before doing lasso fit, may use code from IRLS_glm*/

                            //#pragma omp critical
                            {
                            outEffBetas<<rowIt<<"\t"<<colIt<<"\ttest column number: "<<testIndex<<endl;

                            //if(rowIt==414 && colIt==415){
                            const int rowNum_r = distMatrix_lasso.size();
                            const int colNum_r = (distMatrix_lasso.begin())->size();
                //!!!!!!!!!!! if remove dependents is true, this need to change
                            //Rcpp::NumericVector freq_r(neighbFreq.begin(), neighbFreq.end());

                            outEffBetas<<rowNum_r<<"\t"<<colNum_r<<endl;

                            //Rcpp::NumericMatrix distMat_r(rowNum_r, colNum_r);
                            for (int i=0; i<rowNum_r; ++i)
                            {
                                for (int j=0; j<colNum_r; ++j)
                                {
                                    //double rndNum = preMaxRnd*rand()/(RAND_MAX + 1.0) + minRnd;
                                    //distMat_r(i,j)= floor((distMatrix_lasso[i][j]+rndNum) * PRECISION + 0.5) / PRECISION;
                                    //distMat_r(i,j)= floor((distMatrix_lasso[i][j]) * PRECISION + 0.5) / PRECISION;
                    outEffBetas<<i<<"\t"<<j<<endl;
                                    //distMat_r(i,j)= distMatrix_lasso[i][j];
                                    //distMat_r(i,j)= distMatrix_lasso[i][j] + rndNum;
                                }
                            }
                //distMatrix_lasso.clear();


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
        //outEffBetas<<"chr "<<jobChr_fun<<" domain "<<domainIt<<" siteNum "<<siteNum<<" convert size "<<size_convert<<"fit number "<<counter_fit<<" skip number "<<counter_skip<<endl;
    }
    outEffBetas.close();
    
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




#endif
