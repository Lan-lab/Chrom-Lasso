#ifndef HIC_MPLD_BOSS_H
#define HIC_MPLD_BOSS_H

//component number could be determined using peak scan + 3 (from 1 to hard threshold peak number + 3)
//set start locus at the prescaned peaks
//bayes' rule could be used to update the parameters, like in BALM (read papers about bayes' rule and see if it applicable)

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
//#include <mpi.h>
#include <omp.h>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <RInside.h>

//#include <boost/filesystem.hpp>

#include "peakDetail2D_pairFrag.hpp"
#include "bothEndsMappedFragInfo_withCuttingSite.hpp"

#define INTTAG 0
#define BOOLTAG 1
#define DIETAG 2

//using namespace std;
using std::basic_string;
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


//namespace bfs = boost::filesystem;
//
//using bfs::exists;
//using bfs::directory_iterator;
//using bfs::is_directory;
//using bfs::is_regular;
//using bfs::path;
//using bfs::system_complete;


bool comp_bothEndsMappedFrag_end1(const fragInfo_bothEndsMapped_withCuttingSite&, const fragInfo_bothEndsMapped_withCuttingSite&);
bool comp_bothEndsMappedFrag_end2(const fragInfo_bothEndsMapped_withCuttingSite&, const fragInfo_bothEndsMapped_withCuttingSite&);

int parser_enzymeCuttingMap(const string&, map< int, vector<int> >&);
void reader_interactingRegion(const string&, vector<fragInfo_bothEndsMapped_withCuttingSite>&);
//void findRealInteraction(vector<peakDetail_2D>&, vector<peakDetail_2D>&, vector<peakDetail_2D>&);
int parser_domainFile(const string&, map< int, vector< pair<int, int> > >&);
int extendCuttingSite(const map< int, vector<int> >&, map< int, map< bool, vector< pair<int, int> > > >&, const int);
pair <int, int> mapHybridFragToCuttingSite_end1(vector<fragInfo_bothEndsMapped_withCuttingSite>&, map< int, map< bool, vector< pair<int, int> > > >&, const int);
pair <int, int> mapHybridFragToCuttingSite_end2(vector<fragInfo_bothEndsMapped_withCuttingSite>&, map< int, map< bool, vector< pair<int, int> > > >&, const int);
int distributionDomainMap(string, vector<int>&);
int getCuttingSiteFrags(const map< int, vector<int> >&, const vector<fragInfo_bothEndsMapped_withCuttingSite>&, map< int, map< int, vector<fragInfo_bothEndsMapped_withCuttingSite> > >&, vector<int>&);
int genRanDis(const map< int, map< int, vector<fragInfo_bothEndsMapped_withCuttingSite> > >&, map< int, vector<int> >&, map<int, int>&, const int);
int freqToProb(const map<int, int>&, const int, map<int, double>&);
int pairFragsInteractionFrequency(const map< int, map< int, vector<fragInfo_bothEndsMapped_withCuttingSite> > >&, map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<fragInfo_bothEndsMapped_withCuttingSite> >&, map<int, int>&); //frag here is enzyme cutted frag, map<pair<chr1,pair<frag1Start, frag1End> >, pair<chr2,pair<frag2STart, frag2End> > >, pair<number of interaction, vector<interaction frag> > >
int transformSiteVecToFragMap(const map< int, vector<int> >&, map< int, map<int, int> >&, map<int, int>&);
int findInteractingSites(const map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<fragInfo_bothEndsMapped_withCuttingSite> >, map< int, map<int, int> >&, const double, vector<peakDetail_pairFrag_2D>&); //2Dpeakmap: Map< pair<pair<domain1 chr, domain2 chr>, vector(coordinates and score of each region in this peak)< pair< pair<domain1 position, domain2 position>,score> >, vector<peak property, p value, score etc> >
int searchNeighbourFrags(map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<fragInfo_bothEndsMapped_withCuttingSite> >&, map< int, map<int, int> >&, int, int, pair<int, int>, pair<int, int>, int, int, vector< pair< pair< pair<int, int>, pair<int, int> >, int > >&, double, int&, int&, vector<fragInfo_bothEndsMapped_withCuttingSite>&);
int getPeakSummit(vector<peakDetail_pairFrag_2D>&);
int getPvalCutoff(double*, int, double, double&);

int writer_hybridFrags(const vector<fragInfo_bothEndsMapped_withCuttingSite>&, const string&);
int writer_cuttingSiteFrags(const map< int, map< int, vector<fragInfo_bothEndsMapped_withCuttingSite> > >&, const string&);
int writer_disDistribution(const map<int, int>&, const string&);
int writer_disDistribution_prob(const map<int, double>&, const string&);
int write_fragInteractionMap(map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<fragInfo_bothEndsMapped_withCuttingSite> >&, const string&);
int writer_2Dpeak(vector<peakDetail_pairFrag_2D>&, const string&);
int write_fragInteractionMatrix(map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<fragInfo_bothEndsMapped_withCuttingSite> >&, map< int, map<int, int> >&, int, int, string&);

int boss(const string, const string, const int, const string, const int, map<int, int>, const double);


#endif




	//int interactionFrequencyDistribution(vector<fragInfo_bothEndsMapped_withCuttingSite>&, map< int, vector<int> >&, map< pair< pair<int, int>, pair<int, int> >, pair< int, vector< pair< pair<int, int>, pair<int, int> > > > >&);
	//int writer_interactionFreqMap(map< pair< pair<int, int>, pair<int, int> >, int >&, int, const string&);
	//map<double, double> calThres_2D(const map< pair< pair<int, int>, pair<int, int> >, int >&, const vector<double>, const int);
	//int findPeak_2D(map< pair< pair<int, int>, pair<int, int> >, int >, double, int, map< pair< pair<int, int>, vector< pair< pair<int, int>, int> > >, vector<int> >&);//2Dpeakmap: Map< pair<pair<domain1 chr, domain2 chr>, vector(coordinates and score of each region in this peak)< pair< pair<domain1 position, domain2 position>,score> >, vector<peak property, region number, total score etc> >
	////int findPeak_2D(map< pair< pair<int, int>, pair<int, int> >, pair< int, vector< pair< pair<int, int>, pair<int, int> > > > >, double, int, vector<peakDetail_2D>&);//2Dpeakmap: Map< pair<pair<domain1 chr, domain2 chr>, vector(coordinates and score of each region in this peak)< pair< pair<domain1 position, domain2 position>,score> >, vector<peak property, region number, total score etc> >
	//int searchNeighbour(map< pair< pair<int, int>, pair<int, int> >, pair< int, vector< pair< pair<int, int>, pair<int, int> > > > >&, int, int, int, int, int, int, vector< pair< pair<int, int>, int > >&, double, int&, int&, vector< pair< pair<int, int>, pair<int, int> > >&);
	//int searchNeighbour(map< pair< pair<int, int>, pair<int, int> >, int >&, int, int, int, int, int, int, vector< pair< pair<int, int>, int > >&, double, int&, int&);
	//int writer_2Dpeak(map< pair< pair<int, int>, vector< pair< pair<int, int>, int> > >, vector<int> >&, const string&, const int);
	//int findHighProbTwoEndsInteraction(map< pair< pair< pair<int, int>,pair<int, int> >, vector<int> >, int >&, map< pair< pair< pair<int, int>,pair<int, int> >, vector<int> >, int >&, map< pair< pair< pair<int, int>,pair<int, int> >, vector<int> >, int >&, int);
	//void writer_highProbTwoEndsMap_3C(const map< pair< pair< pair<int, int>,pair<int, int> >, vector<int> >, int >&, const string&);
	//void removeRepetitiveRegion(const map< pair< pair< pair<int, int>,pair<int, int> >, vector<int> >, int >&, const map< int, vector< pair<int, int> > >&, map< pair< pair< pair<int, int>,pair<int, int> >, vector<int> >, int >&);
	//void reader_repRegionFile(const string&, map< int, vector< pair<int, int> > >&);
	//int getReadsInfo(const string&, vector<readInfo_PE>&);
	//unsigned long long int loadFileToMemory(const char*, char**);
	//int reader_fastaFile(const string, vector<char>&);
	//void reader_oneEndMappedFrags(const string&, vector<fragInfo_oneEndMapped>&);
	//void reader_bothEndsMappedFrags(const string&, vector<fragInfo_bothEndsMapped_withCuttingSite>&);
	////int findBP(vector<peakDetail_2D>&, vector<fragInfo_oneEndMapped>&);
	////int mapReadToRegion(fragInfo_oneEndMapped&, peakDetail_2D&, ofstream&);


