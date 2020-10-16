#ifndef PEAKDETAILPAIRFRAG2D_CPP
#define PEAKDETAILPAIRFRAG2D_CPP


#include "peakDetail2D_pairFrag.hpp"


int peakDetail_pairFrag_2D::getSumLogProb(const map<int, double>& ranDisMap_local, map<int, map< int, vector<fragInfo_bothEndsMapped_withCuttingSite> > >& cuttingSiteFragMap_local, const double minProb_local, const double aveSiteFrag_local, const int binSize_local)
{
	vector<fragInfo_bothEndsMapped_withCuttingSite>::const_iterator fragIt = frag_info.begin(), fragIt_end = frag_info.end();
	for ( ; fragIt != fragIt_end; ++fragIt )
	{
		map< int, vector<fragInfo_bothEndsMapped_withCuttingSite> >::const_iterator end1Site = cuttingSiteFragMap_local[fragIt->end1_chr].find(fragIt->end1_cuttingSite);
		double correctionCoeff = (end1Site->second.size()) * 1.0 / aveSiteFrag_local;
		double minProb_corrected = minProb_local * correctionCoeff;
		double logMinProb = log(minProb_corrected);

		if (fragIt->end1_chr == fragIt->end2_chr)
		{
			int dis = abs(fragIt->end2_cuttingSite - fragIt->end1_cuttingSite);
			int binIndex = dis / binSize_local;

			map<int, double>::const_iterator disProbIt = ranDisMap_local.find(binIndex);
			if (disProbIt != ranDisMap_local.end())
			{
				sumLogProb += log((disProbIt->second)*correctionCoeff);
			}
			else
			{
				sumLogProb += logMinProb;
			}
		}
		else
		{
			sumLogProb += logMinProb;
		}
	}

	return 0;
}


#endif
