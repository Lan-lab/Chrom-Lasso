#ifndef REGIONDETAILPAIRSITE2D_HPP
#define REGIONDETAILPAIRSITE2D_HPP

#include <string>
#include <vector>
#include <map>
#include <math.h>

//#include "bothEndsMappedFragInfo_withCuttingSite.hpp"

using std::pair;
using std::string;
using std::vector;
using std::map;
using std::make_pair;


//******concept of fragPair (a pair of enzyme digested frags) is different then frag (reads pair)********//
class regionDetail_pairSite_2D
{
public:
	int chrNo1;
	int chrNo2;
	int domain1;
	int domain2;
	vector< pair< pair<int, int>, double > > region_info;
	pair< pair<int, int>, double > summitPairSites;
	int sitePair_num;
	//int total_frag;
	//vector<fragInfo_bothEndsMapped_withCuttingSite> frag_info;
	pair< pair<int, int>, pair<int, int> > regionCoor;

	//int end1_regionStart;
	//int end1_regionEnd;

	//int end1_max;
	//int end1_min;
	//int end2_max;
	//int end2_min;

	double minPval;


	regionDetail_pairSite_2D(): chrNo1(0),chrNo2(0),domain1(0),domain2(0),sitePair_num(0),regionCoor(make_pair(make_pair(0,0), make_pair(0,0))),minPval(0.0){}

	int getRegionCoorSummit();

	//int getSumLogProb(const map<int, double>&, map<int, map< int, vector<fragInfo_bothEndsMapped_withCuttingSite> > >&, const double, const double, const int);
	//void findPeakRegion();
	//void getRegionSequence(map< int, vector<char> >&);
	//int getPeakRegion();
};


#endif
