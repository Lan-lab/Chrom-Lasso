#ifndef HIC_MPLD_POSTDOC_H
#define HIC_MPLD_POSTDOC_H

#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
//#include <mpi.h>
#include <omp.h>
#include <vector>

#include <gsl/gsl_cdf.h>

#include <../IRLS_glm/IRLS.h>

#include <regionDetail2D_pairSite.hpp>
#include <bothEndsMappedFragInfo_withCuttingSite.hpp>

#include <RInside.h>

#define INTTAG 0
#define BOOLTAG 1
#define DIETAG 2

#define NEIGHBDIS 5
#define MORANI 0.001
#define PCUTOFF 0.004
#define PRECISION 100000

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


int postdoc(vector<double>, map<int, int>, const int, const string);
int receiveDomainMap(vector< pair<int, int> >&);
int write_domainFragInteractionMatrix(map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<fragInfo_bothEndsMapped_withCuttingSite> >&, map<int, int>&, vector< pair<int, int> >&, int, string&);
int getDomainCSinteractionMatrix(map< pair< pair<int, int>, pair<int, int> >, int >&, map< int, vector<int> >&, vector< pair<int, int> >&, int, vector< vector<int> >&, vector< vector< vector<int> > >&, map< int, map<int, int> >&);
inline int convtIndex(int, int, int);
int findIntraDomainInteraction(string, vector< vector<int> >&, vector< vector< vector<int> > >&, map<int, int>&, map< int, map< pair<int, int>, double > >&);
int outputDomainInfo(string, vector< vector<int> >&, vector< vector< vector<int> > >&, map<int, int>&, map< int, map< pair<int, int>, double > >&);
int findInteractingRegion(int, map< int, map< pair<int, int>, double> >&, vector< vector<int> >&, vector<regionDetail_pairSite_2D>&);
int searchNeighbSites(map< pair<int, int>, double >&, map<int, int>&, int, int, int, int, vector< pair< pair<int, int>, double > >&, int&, double&);
int writer_2Dregion(vector<regionDetail_pairSite_2D>&, const string&);
int glm_lasso_fit(vector<double>&, vector<double>&, vector< vector<double> >&);

#endif
