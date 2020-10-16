#ifndef PEAKDETAILPAIRFRAG2D_HPP
#define PEAKDETAILPAIRFRAG2D_HPP

#include <string>
#include <vector>
#include <map>
#include <math.h>

#include "bothEndsMappedFragInfo_withCuttingSite.hpp"

using std::pair;
using std::string;
using std::vector;
using std::map;


//******concept of fragPair (a pair of enzyme digested frags) is different then frag (reads pair)********//
class peakDetail_pairFrag_2D
{
public:
	int chrNo_domain1;
	int chrNo_domain2;
	vector< pair< pair< pair<int, int>, pair<int, int> >, int > > region_info;
	int fragPair_num;
	int total_frag;
	vector<fragInfo_bothEndsMapped_withCuttingSite> frag_info;
	pair< pair<int, int>, pair<int, int> > peakRegion;

	bool end1_strand;
	bool end2_strand;
	
	int end1_regionStart;
	int end1_regionEnd;

	int end1_max;
	int end1_min;
	int end2_max;
	int end2_min;

	double sumLogProb;


	peakDetail_pairFrag_2D(): chrNo_domain1(0),chrNo_domain2(0),fragPair_num(0),total_frag(0),end1_max(0),end1_min(0),end2_max(0),end2_min(0),sumLogProb(0.0){}

	int getSumLogProb(const map<int, double>&, map<int, map< int, vector<fragInfo_bothEndsMapped_withCuttingSite> > >&, const double, const double, const int);
	//void findPeakRegion();
	//void getRegionSequence(map< int, vector<char> >&);
	//int getPeakRegion();
};


#endif
