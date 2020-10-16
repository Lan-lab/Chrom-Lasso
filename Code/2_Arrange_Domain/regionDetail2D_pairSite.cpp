#ifndef REGIONDETAILPAIRSITE2D_CPP
#define REGIONDETAILPAIRSITE2D_CPP


#include "regionDetail2D_pairSite.hpp"


//get the coor of the regions and find the summit of a region, if there are multiple paired enzyme cutted site, choose the one with lowest p value 
int regionDetail_pairSite_2D::getRegionCoorSummit()
{
	//int addPeakNum = 0;
	//sitePair_num = region_info.size();

	vector< pair< pair<int, int>, double > >::const_iterator pairSiteIt = region_info.begin(), pairSiteIt_end = region_info.end();
	regionCoor.first.second = pairSiteIt->first.first;
	regionCoor.first.first = pairSiteIt->first.first;
	regionCoor.second.second = pairSiteIt->first.second;
	regionCoor.second.first = pairSiteIt->first.second;
	if (pairSiteIt->second == minPval)
	{
		summitPairSites = *pairSiteIt;
	}
	++pairSiteIt;
	for ( ; pairSiteIt != pairSiteIt_end; ++pairSiteIt )
	{
		if (pairSiteIt->first.first > regionCoor.first.second)
		{
			regionCoor.first.second = pairSiteIt->first.first;
		}
		if (pairSiteIt->first.first < regionCoor.first.first)
		{
			regionCoor.first.first = pairSiteIt->first.first;
		}
		if (pairSiteIt->first.second > regionCoor.second.second)
		{
			regionCoor.second.second = pairSiteIt->first.second;
		}
		if (pairSiteIt->first.second < regionCoor.second.first)
		{
			regionCoor.second.first = pairSiteIt->first.second;
		}
		if (pairSiteIt->second == minPval)
		{
			summitPairSites = *pairSiteIt;
		}
	}

	return 0;
}


#endif
