#ifndef HIC_MPLD_BOSS_CPP
#define HIC_MPLD_BOSS_CPP

//component number could be determined using peak scan + 3 (from 1 to hard threshold peak number + 3)
//set start locus at the prescaned peaks
//bayes' rule could be used to update the parameters, like in BALM (read papers about bayes' rule and see if it applicable)

#include "HiC_mixturePLD_boss.h"
#include "HiC_mixturePLD_postdoc.h"
#include "bothEndsMappedFragInfo_withCuttingSite.hpp"


bool comp_bothEndsMappedFrag_end1(const fragInfo_bothEndsMapped_withCuttingSite& lhs, const fragInfo_bothEndsMapped_withCuttingSite& rhs)
{
	if (lhs.end1_chr!=rhs.end1_chr)
	{
		return lhs.end1_chr<rhs.end1_chr;
	}
	else if (lhs.end1_pos!=rhs.end1_pos)
	{
		return lhs.end1_pos<rhs.end1_pos;
	}
	else if (lhs.end2_chr!=rhs.end2_chr)
	{
		return lhs.end2_chr<rhs.end2_chr;
	}
	else
	{
		return lhs.end2_pos<rhs.end2_pos;
	}
}




bool comp_bothEndsMappedFrag_end2(const fragInfo_bothEndsMapped_withCuttingSite& lhs, const fragInfo_bothEndsMapped_withCuttingSite& rhs)
{
	if (lhs.end2_chr!=rhs.end2_chr)
	{
		return lhs.end2_chr<rhs.end2_chr;
	}
	else if (lhs.end2_pos!=rhs.end2_pos)
	{
		return lhs.end2_pos<rhs.end2_pos;
	}
	else if (lhs.end1_chr!=rhs.end1_chr)
	{
		return lhs.end1_chr<rhs.end1_chr;
	}
	else
	{
		return lhs.end1_pos<rhs.end1_pos;
	}
}




int boss(const string file_cuttingSites_local, const string domainFile_local, const int cuttingSiteExtent_local, const string file_source_local, const int halfReadSize_local, map<int, int> chrLen, const double FDR)
{
	cout <<"Start processing....."<<endl;

	map< int, vector<int> > cuttingSitesMap; //map< chr, vector<pos> >
	int cuttingSiteTotal = parser_enzymeCuttingMap(file_cuttingSites_local, cuttingSitesMap);

	cout <<"Parse enzyme cutting map done."<<endl;

	map< int, map< bool, vector< pair<int, int> > > > cuttingSiteExtMap; //map< chr, map< strand, vector< pair<start, end> > > >
	extendCuttingSite(cuttingSitesMap, cuttingSiteExtMap, cuttingSiteExtent_local);

	cout <<"Extend cutting site region done."<<endl;

	map< int, vector< pair<int, int> > > domainMap;
	int domainNum = parser_domainFile(domainFile_local, domainMap);

	//vector<fragInfo_bothEndsMapped_withCuttingSite> interChrMappedVec;

	ifstream inputFile(file_source_local.c_str());
	if (!inputFile)
	{
		cout <<"\n"<< "Error opening " << file_source_local << "." << endl;
		exit(1);
	}

	fragInfo_bothEndsMapped_withCuttingSite prevFrag;

	map<int, int>::const_iterator chrIt = chrLen.begin(), chrIt_end = chrLen.end();
	int currChr=chrIt->first;
	int prevChr=0;
	//for ( ; chrIt!=chrIt_end; ++chrIt)
	//{
	char chrNo_end1[3];
	char chrNo_end2[3];
	char start_end1[15];
	char start_end2[15];
	char strand_end1 = '0';
	char strand_end2 = '0';
	char fragType_str[4];

	char* it_lineStr;
	char* it_token;

	char lineStr[512];
	bool endFileFlag=true; //false is end, true is not

	map< pair< pair<int, int>, pair<int, int> >, int > cuttingSiteInteractionMap_local;
	bool chrCompleteFlag=false;

	fragInfo_bothEndsMapped_withCuttingSite prevFragInfo;
	fragInfo_bothEndsMapped_withCuttingSite prevFragInfo_swap;

	while(endFileFlag)
	{
	    vector<fragInfo_bothEndsMapped_withCuttingSite> bothEndMappedVec;
	    if(prevFragInfo.end1_chr==currChr)
	    {
		bothEndMappedVec.push_back(prevFragInfo);
		bothEndMappedVec.push_back(prevFragInfo_swap);
		prevFragInfo.end1_chr=0;
	    }

	    for (int i=1;i<=10000000;++i)
	    {
		endFileFlag=inputFile.getline(lineStr,512);
		if(endFileFlag==false)
		{
		    break;
		}
		else
		{
		    if (lineStr != NULL)
		    {
			fragInfo_bothEndsMapped_withCuttingSite fragInfo;
			fragInfo_bothEndsMapped_withCuttingSite fragInfo_swap;
			it_lineStr = lineStr;

			it_token = chrNo_end1;
			while(*it_lineStr != '\t')
			{
				*it_token = *it_lineStr;
				++it_token;
				++it_lineStr;
			}
			*it_token = '\0';
			++it_lineStr;
			if (chrNo_end1[0] == 'X')
			{
				fragInfo.end1_chr = 23;
				fragInfo_swap.end2_chr = 23;
			}
			else if (chrNo_end1[0] == 'Y')
			{
				fragInfo.end1_chr = 24;
				fragInfo_swap.end2_chr = 24;
			}
			else if (chrNo_end1[0] == 'M')
			{
				continue;
			}
			else
			{
				fragInfo.end1_chr = atoi(chrNo_end1);
				fragInfo_swap.end2_chr = atoi(chrNo_end1);
			}

			it_token = start_end1;
			while(*it_lineStr != '\t')
			{
				*it_token = *it_lineStr;
				++it_token;
				++it_lineStr;
			}
			*it_token = '\0';
			++it_lineStr;
			fragInfo.end1_pos = atoi(start_end1);
			fragInfo_swap.end2_pos = atoi(start_end1);

			strand_end1 = *it_lineStr;
			if (strand_end1 == '0')
			{
				fragInfo.end1_strand = false;
				fragInfo_swap.end2_strand = false;
			}
			else
			{
				fragInfo.end1_strand = true;
				fragInfo_swap.end2_strand = true;
			}
			++it_lineStr;++it_lineStr;

			it_token = chrNo_end2;
			while(*it_lineStr != '\t')
			{
				*it_token = *it_lineStr;
				++it_token;
				++it_lineStr;
			}
			*it_token = '\0';
			++it_lineStr;
			if (chrNo_end2[0] == 'X')
			{
				fragInfo.end2_chr = 23;
				fragInfo_swap.end1_chr = 23;
			}
			else if (chrNo_end2[0] == 'Y')
			{
				fragInfo.end2_chr = 24;
				fragInfo_swap.end1_chr = 24;
			}
			else if (chrNo_end2[0] == 'M')
			{
				continue;
			}
			else
			{
				fragInfo.end2_chr = atoi(chrNo_end2);
				fragInfo_swap.end1_chr = atoi(chrNo_end2);
			}

			it_token = start_end2;
			while(*it_lineStr != '\t')
			{
				*it_token = *it_lineStr;
				++it_token;
				++it_lineStr;
			}
			*it_token = '\0';
			++it_lineStr;
			fragInfo.end2_pos = atoi(start_end2);
			fragInfo_swap.end1_pos = atoi(start_end2);

			strand_end2 = *it_lineStr;
			if (strand_end2 == '0')
			{
				fragInfo.end2_strand = false;
				fragInfo_swap.end1_strand = false;
			}
			else
			{
				fragInfo.end2_strand = true;
				fragInfo_swap.end1_strand = true;
			}
			++it_lineStr;++it_lineStr;

			it_token = fragType_str;
			while(*it_lineStr != '\t' && *it_lineStr != '\n' && *it_lineStr != '\0')
			{
				*it_token = *it_lineStr;
				++it_token;
				++it_lineStr;
			}
			*it_token = '\0';

			fragInfo.fragType = atoi(fragType_str);
			fragInfo_swap.fragType = atoi(fragType_str);

			if(fragInfo.end1_chr!=currChr)
			{
			    chrCompleteFlag=true;
			    currChr=fragInfo.end1_chr;
			    prevFragInfo=fragInfo;
			    prevFragInfo_swap=fragInfo_swap;
			    break;
			}
			else
			{
			    //if(fragInfo.end1_chr == fragInfo.end2_chr)
			    //{
				bothEndMappedVec.push_back(fragInfo);
				bothEndMappedVec.push_back(fragInfo_swap);
			    //}
			    //else
			    //{
			    //    interChrMappedVec.push_back(fragInfo);
			    //    interChrMappedVec.push_back(fragInfo_swap);
			    //}
			    prevChr=fragInfo.end1_chr;
			}
		    }
		}
	    }

	    cout<<"read both end mapped reads done."<<endl;
	    sort(bothEndMappedVec.begin(), bothEndMappedVec.end(), comp_bothEndsMappedFrag_end1);
	    cout<<"sort both end mapped reads done."<<endl;
	    //writer_hybridFrags(bothEndMappedVec, "interactingFrags");

	    vector<fragInfo_bothEndsMapped_withCuttingSite> interactingFrags;

	    vector<fragInfo_bothEndsMapped_withCuttingSite>::iterator bothEndMappedVecIt = bothEndMappedVec.begin(), bothEndMappedVecIt_end = bothEndMappedVec.end();
	    for ( ; bothEndMappedVecIt != bothEndMappedVecIt_end; ++bothEndMappedVecIt )
	    {
		    if ( ( (bothEndMappedVecIt->end1_chr != bothEndMappedVecIt->end2_chr) || (bothEndMappedVecIt->end1_pos - bothEndMappedVecIt->end2_pos > 1000) || (bothEndMappedVecIt->end2_pos - bothEndMappedVecIt->end1_pos > 1000) ) && (bothEndMappedVecIt->end1_chr != prevFrag.end1_chr || bothEndMappedVecIt->end1_pos != prevFrag.end1_pos || bothEndMappedVecIt->end2_chr != prevFrag.end2_chr || bothEndMappedVecIt->end2_pos != prevFrag.end2_pos) )
		    {
			    interactingFrags.push_back(*bothEndMappedVecIt);
		    }
		    prevFrag = *bothEndMappedVecIt;
	    }

	    bothEndMappedVec.clear();
	    cout<<"Remove duplicated tags done."<<endl;

	    //writer_hybridFrags(interactingFrags, "removeRedundant_done");

	    pair<int, int> end1MapInfo, end2MapInfo;
	    end1MapInfo = mapHybridFragToCuttingSite_end1(interactingFrags, cuttingSiteExtMap, halfReadSize_local);
	    cout<<"Map end 1 done."<<endl;
	    //writer_hybridFrags(interactingFrags, "end1_done");
	    sort(interactingFrags.begin(), interactingFrags.end(), comp_bothEndsMappedFrag_end2);
	    cout<<"Sort end 2 done."<<endl;
	    end2MapInfo = mapHybridFragToCuttingSite_end2(interactingFrags, cuttingSiteExtMap, halfReadSize_local);
	    cout<<"Map end 2 done."<<endl;

	    cout <<"End1:\t"<<end1MapInfo.first<<"\t"<<end1MapInfo.second<<endl;
	    cout <<"End2:\t"<<end2MapInfo.first<<"\t"<<end2MapInfo.second<<endl;

	    //1) filter out 0 or 1 end mapped to cutting site frags, and keep both end mapped. 2) link frag to the cutting site map 3) remove two ends mapped to the same cutting site or adjacent cutting site if the two ends following the self loop characteristic (end mapped to the precedent cutting site is -, and the other end is +)
	    //!!!!IMPORTANT selfloop can also cross multiple enzyme digested fragments.
	    //NOTE new program will take the Ren's group HiC summary file, so frag need to be duplicated since one frag only have one record in the input
	    //map< int, map< int, vector<fragInfo_bothEndsMapped_withCuttingSite> > > cuttingSiteFragsMap; //map< chr, map< pos, vector<frags linked to the cutting site> > >, because one frag is already duplicated in the fragMap (each end assign end1 and end2 respectively), in this step each frag in the fragMap will only assigned to one cutting site.
	    //
	    //interchromosome interactions do not need to be saved here (cuttingSiteInteractionMap_local)
	    //only need to count the number of inter chr interactions for each cutting sites.
	    vector<fragInfo_bothEndsMapped_withCuttingSite>::const_iterator fragIt=interactingFrags.begin(), fragIt_end=interactingFrags.end();
	    for ( ; fragIt!=fragIt_end; ++fragIt)
	    {
		cuttingSiteInteractionMap_local[make_pair(make_pair(fragIt->end1_chr, fragIt->end1_cuttingSite), make_pair(fragIt->end2_chr, fragIt->end2_cuttingSite))]++;
	    }
	    interactingFrags.clear();

	    if(chrCompleteFlag==true || endFileFlag==false)
	    {
		vector< vector<int> > domainSitesMap_local;
		vector< vector < vector<int> > > domainCSinterFreq_local;
		map<int, map<int, int> > csInterChromTotalMap_local;
		if(domainMap.find(prevChr)!=domainMap.end())
		{
		    getDomainCSinteractionMatrix(cuttingSiteInteractionMap_local, cuttingSitesMap, domainMap[prevChr], prevChr, domainSitesMap_local, domainCSinterFreq_local, csInterChromTotalMap_local);


		    char prevChrName[6];
		    if (prevChr<23)
		    {
			    sprintf(prevChrName, "chr%d", prevChr);
		    }
		    else if (prevChr==23)
		    {
			    sprintf(prevChrName, "chrX");
		    }
		    else if (prevChr==24)
		    {
			    sprintf(prevChrName, "chrY");
		    }
		    else if (prevChr==25)
		    {
			    sprintf(prevChrName, "inter");
		    }
		    else
		    {
			    cout<<"Warning! Chromosome not found in Human or Mouse."<<endl;
		    }

		    map< int, map< pair<int, int>, double > > pVal_local;
		    outputDomainInfo(prevChrName, domainSitesMap_local, domainCSinterFreq_local, csInterChromTotalMap_local[prevChr], pVal_local);
		}
		else
		{
		    cerr<<"Warning! Chromosome "<<prevChr<<" not found!"<<endl;
		}

		chrCompleteFlag=false;
		cuttingSiteInteractionMap_local.clear();
	    }
	}

	//double aveSiteFrag = totalHybridFrag * 1.0 / cuttingSiteTotal;
	//cout<<"total hybrid frags "<<totalHybridFrag<<" average site frag "<<aveSiteFrag<<"\nconstruct cutting site frag map done."<<endl;
	//
	return 0;
}




int extendCuttingSite(const map< int, vector<int> >& cuttingSiteMap_local, map< int, map< bool, vector< pair<int, int> > > >& cuttingSiteExtMap_local, const int cuttingSiteExtent_local) //leave the read size alone, it will be taken cared of during mapping all the reads to the cutting site, in which reads should be shitfed half the read size
{
	map< int, vector<int> >::const_iterator chrIt = cuttingSiteMap_local.begin(), chrIt_end = cuttingSiteMap_local.end();
	for ( ; chrIt != chrIt_end; ++chrIt )
	{
		vector<int>::const_iterator siteIt = chrIt->second.begin(), siteIt_end = chrIt->second.end();

		//if true strand, site should be second, if false then should be first
		//first site
		if (*siteIt <= cuttingSiteExtent_local)
		{
			cuttingSiteExtMap_local[chrIt->first][true].push_back(make_pair(1, *siteIt));
		}
		else
		{
			cuttingSiteExtMap_local[chrIt->first][true].push_back(make_pair((*siteIt)-cuttingSiteExtent_local, *siteIt));
		}
		int preSite = *siteIt;
		++siteIt;

		for ( ; siteIt != siteIt_end; ++siteIt )
		{
			if (*siteIt-preSite > cuttingSiteExtent_local)
			{
				cuttingSiteExtMap_local[chrIt->first][false].push_back(make_pair(preSite, preSite+cuttingSiteExtent_local));
				cuttingSiteExtMap_local[chrIt->first][true].push_back(make_pair(*siteIt-cuttingSiteExtent_local, *siteIt));
			}
			else
			{
				cuttingSiteExtMap_local[chrIt->first][false].push_back(make_pair(preSite, *siteIt-1));
				cuttingSiteExtMap_local[chrIt->first][true].push_back(make_pair(preSite+1, *siteIt));
			}

			preSite = *siteIt;
		}

		//last site
		--siteIt;
		cuttingSiteExtMap_local[chrIt->first][false].push_back(make_pair(*siteIt, *siteIt+cuttingSiteExtent_local)); //it seems that it would not be a problem for the last region to exceed the chromosomal len
	}

	return 0;
}




//map hybrid reads to cutting sites
pair <int, int> mapHybridFragToCuttingSite_end1(vector<fragInfo_bothEndsMapped_withCuttingSite>& hybridFragVec_local, map< int, map< bool, vector< pair<int, int> > > >& cuttingSiteExtMap_local, const int halfReadSize_local)
{
	int regionStart = 0;
	int regionEnd = 0;
	int totalFragsInRegions = 0;
	int totalRegionNum = 0;
	double totalAveDensity = 0.0;

	map< int, map< bool, vector< pair<int, int> > > >::iterator chrIt_regionMap = cuttingSiteExtMap_local.begin(), chrIt_regionMap_end = cuttingSiteExtMap_local.end();
	vector<fragInfo_bothEndsMapped_withCuttingSite>::iterator fragIt_begin = hybridFragVec_local.begin(), fragIt = hybridFragVec_local.begin(), fragIt_chrStart = hybridFragVec_local.begin(), fragIt_end = hybridFragVec_local.end();

	//true strand first, then false strand
	for ( ; chrIt_regionMap != chrIt_regionMap_end; ++chrIt_regionMap )
	{
		//cout <<"\nentered"<<i<<endl;
		while (fragIt->end1_chr < chrIt_regionMap->first && fragIt != fragIt_end)
		{
			++fragIt;
		}

		if (fragIt == fragIt_end)
		{
			break;
		}
		else if (fragIt->end1_chr > chrIt_regionMap->first)
		{
			continue;
		}
		else if (fragIt->end1_chr == chrIt_regionMap->first)
		{
			fragIt_chrStart = fragIt;
			//while (fragIt != fragIt_begin && fragIt->end1_chr == chrIt_regionMap->first)
			//{
			//	--fragIt;
			//}
			//if (fragIt != fragIt_end && fragIt->end1_chr < chrIt_regionMap->first)
			//{
			//	++fragIt;
			//}

			vector< pair<int, int> >::const_iterator regionIt = chrIt_regionMap->second[true].begin(), regionIt_end = chrIt_regionMap->second[true].end();
			while (regionIt != regionIt_end && fragIt->end1_chr == chrIt_regionMap->first)
			{
				//cout <<"entered while"<<endl;
				regionStart = regionIt->first;
				regionEnd = regionIt->second;

				while ( fragIt != fragIt_end && fragIt->end1_chr == chrIt_regionMap->first && (fragIt->end1_pos + halfReadSize_local) < regionStart )
				{
					++fragIt;
					//cout <<"ET1W1\t";
				}

				bool moveBackFlag = false;
				if ( fragIt != fragIt_end && fragIt->end1_chr == chrIt_regionMap->first && fragIt != fragIt_chrStart && (fragIt->end1_pos + halfReadSize_local) >= regionStart ) //--fragIt is to counter the effect of the following ++fragIt, when it would not go into the following while loop
				{
					--fragIt;
					moveBackFlag = true;
				}

				while ( fragIt != fragIt_end && fragIt->end1_chr == chrIt_regionMap->first && fragIt != fragIt_chrStart && (fragIt->end1_pos + halfReadSize_local) >= regionStart )
				{
					--fragIt;
				}
				//cout <<"OT1W2\t";

				if ( fragIt != fragIt_end && fragIt->end1_chr == chrIt_regionMap->first && moveBackFlag )
				{
					++fragIt;
				}
				else if ( fragIt == fragIt_end || fragIt->end1_chr != chrIt_regionMap->first )
				{
					break;
				}


				//cout <<"OT1W2\t";
				int fragsCounter = 0;
				//int startingFrag = 0;

				for ( ; fragIt != fragIt_end && fragIt->end1_chr == chrIt_regionMap->first && (fragIt->end1_pos + halfReadSize_local) <= regionEnd ; ++fragIt )
				//for (startingRead_forward = *fragIt_forward; *fragIt_forward < fragRangeEnd_forward; ++fragIt_forward)
				{
					//cout <<"ET1F\t";
					if (fragIt->end1_strand == true) //this made the totalRegionNum meaningless. Even the for loop is entered, there might be no frags in this region since the frag might only on the other strand
					{
						fragIt->end1_cuttingSite = regionEnd; //true strand = regionEnd and false strand = regionStart
						fragsCounter += 1;
					}
				}

				++totalRegionNum;
				totalFragsInRegions += fragsCounter;

				++regionIt;
			}
		}
	}



	//false strand
	chrIt_regionMap = cuttingSiteExtMap_local.begin(), chrIt_regionMap_end = cuttingSiteExtMap_local.end();
	fragIt_begin = hybridFragVec_local.begin(), fragIt = hybridFragVec_local.begin(), fragIt_chrStart = hybridFragVec_local.begin(), fragIt_end = hybridFragVec_local.end();

	for ( ; chrIt_regionMap != chrIt_regionMap_end; ++chrIt_regionMap )
	{
		//cout <<"\nentered"<<i<<endl;
		while (fragIt->end1_chr < chrIt_regionMap->first && fragIt != fragIt_end)
		{
			++fragIt;
		}

		if (fragIt == fragIt_end)
		{
			break;
		}
		else if (fragIt->end1_chr > chrIt_regionMap->first)
		{
			continue;
		}
		else if (fragIt->end1_chr == chrIt_regionMap->first)
		{
			fragIt_chrStart = fragIt;
			//while (fragIt != fragIt_begin && fragIt->end1_chr == chrIt_regionMap->first)
			//{
			//	--fragIt;
			//}
			//if (fragIt != fragIt_end && fragIt->end1_chr < chrIt_regionMap->first)
			//{
			//	++fragIt;
			//}

			vector< pair<int, int> >::const_iterator regionIt = chrIt_regionMap->second[false].begin(), regionIt_end = chrIt_regionMap->second[false].end();
			while (regionIt != regionIt_end && fragIt->end1_chr == chrIt_regionMap->first)
			{
				//cout <<"entered while"<<endl;
				regionStart = regionIt->first;
				regionEnd = regionIt->second;

				while ( fragIt != fragIt_end && fragIt->end1_chr == chrIt_regionMap->first && (fragIt->end1_pos + halfReadSize_local) < regionStart ) //--fragIt is to counter the effect of the following ++fragIt, when it would not go into the following while loop
				{
					++fragIt;
					//cout <<"ET1W1\t";
				}

				bool moveBackFlag = false;
				if ( fragIt != fragIt_end && fragIt->end1_chr == chrIt_regionMap->first && fragIt != fragIt_chrStart && (fragIt->end1_pos + halfReadSize_local) >= regionStart ) //--fragIt is to counter the effect of the following ++fragIt, when it would not go into the following while loop
				{
					--fragIt;
					moveBackFlag = true;
				}

				while ( fragIt != fragIt_end && fragIt->end1_chr == chrIt_regionMap->first && fragIt != fragIt_chrStart && (fragIt->end1_pos + halfReadSize_local) >= regionStart )
				{
					--fragIt;
				}
				//cout <<"OT1W2\t";

				if ( fragIt != fragIt_end && fragIt->end1_chr == chrIt_regionMap->first && moveBackFlag )
				{
					++fragIt;
				}
				else if ( fragIt == fragIt_end || fragIt->end1_chr != chrIt_regionMap->first )
				{
					break;
				}


				//cout <<"OT1W2\t";
				int fragsCounter = 0;
				//int startingFrag = 0;

				for ( ; fragIt != fragIt_end && fragIt->end1_chr == chrIt_regionMap->first && (fragIt->end1_pos + halfReadSize_local) <= regionEnd ; ++fragIt )
				//for (startingRead_forward = *fragIt_forward; *fragIt_forward < fragRangeEnd_forward; ++fragIt_forward)
				{
					//cout <<"ET1F\t";
					if (fragIt->end1_strand == false) //this made the totalRegionNum meaningless. Even the for loop is entered, there might be no frags in this region since the frag might only on the other strand
					{
						fragIt->end1_cuttingSite = regionStart; //true strand = regionEnd and false strand = regionStart
						fragsCounter += 1;
					}
				}

				++totalRegionNum;
				totalFragsInRegions += fragsCounter;

				++regionIt;
			}
		}
	}

	return make_pair(totalRegionNum, totalFragsInRegions);
	//return totalFragsInPeaks;
}




//map hybrid reads to cutting sites
pair <int, int> mapHybridFragToCuttingSite_end2(vector<fragInfo_bothEndsMapped_withCuttingSite>& hybridFragVec_local, map< int, map< bool, vector< pair<int, int> > > >& cuttingSiteExtMap_local, const int halfReadSize_local)
{
	int regionStart = 0;
	int regionEnd = 0;
	int totalFragsInRegions = 0;
	int totalRegionNum = 0;
	double totalAveDensity = 0.0;

	map< int, map< bool, vector< pair<int, int> > > >::iterator chrIt_regionMap = cuttingSiteExtMap_local.begin(), chrIt_regionMap_end = cuttingSiteExtMap_local.end();
	vector<fragInfo_bothEndsMapped_withCuttingSite>::iterator fragIt_begin = hybridFragVec_local.begin(), fragIt = hybridFragVec_local.begin(), fragIt_chrStart = hybridFragVec_local.begin(), fragIt_end = hybridFragVec_local.end();

	//true strand first, then false strand
	for ( ; chrIt_regionMap != chrIt_regionMap_end; ++chrIt_regionMap )
	{
		//cout <<"\nentered"<<i<<endl;
		while (fragIt->end2_chr < chrIt_regionMap->first && fragIt != fragIt_end)
		{
			++fragIt;
		}

		if (fragIt == fragIt_end)
		{
			break;
		}
		else if (fragIt->end2_chr > chrIt_regionMap->first)
		{
			continue;
		}
		else if (fragIt->end2_chr == chrIt_regionMap->first)
		{
			fragIt_chrStart = fragIt;
			//while (fragIt != fragIt_begin && fragIt->end2_chr == chrIt_regionMap->first)
			//{
			//	--fragIt;
			//}
			//if (fragIt != fragIt_end && fragIt->end2_chr < chrIt_regionMap->first)
			//{
			//	++fragIt;
			//}

			vector< pair<int, int> >::const_iterator regionIt = chrIt_regionMap->second[true].begin(), regionIt_end = chrIt_regionMap->second[true].end();
			while (regionIt != regionIt_end && fragIt->end2_chr == chrIt_regionMap->first)
			{
				//cout <<"entered while"<<endl;
				regionStart = regionIt->first;
				regionEnd = regionIt->second;

				while ( fragIt != fragIt_end && fragIt->end2_chr == chrIt_regionMap->first && (fragIt->end2_pos + halfReadSize_local) < regionStart ) //--fragIt is to counter the effect of the following ++fragIt, when it would not go into the following while loop
				{
					++fragIt;
					//cout <<"ET1W1\t";
				}

				bool moveBackFlag = false;
				if ( fragIt != fragIt_end && fragIt->end2_chr == chrIt_regionMap->first && fragIt != fragIt_chrStart && (fragIt->end2_pos + halfReadSize_local) >= regionStart )
				{
					--fragIt;
					moveBackFlag = true;
				}

				while ( fragIt != fragIt_end && fragIt->end2_chr == chrIt_regionMap->first && fragIt != fragIt_chrStart && (fragIt->end2_pos + halfReadSize_local) >= regionStart )
				{
					--fragIt;
				}
				//cout <<"OT1W2\t";

				if ( fragIt != fragIt_end && fragIt->end2_chr == chrIt_regionMap->first && moveBackFlag )
				{
					++fragIt;
				}
				else if ( fragIt == fragIt_end || fragIt->end2_chr != chrIt_regionMap->first )
				{
					break;
				}


				//cout <<"OT1W2\t";
				int fragsCounter = 0;
				//int startingFrag = 0;

				for ( ; fragIt != fragIt_end && fragIt->end2_chr == chrIt_regionMap->first && (fragIt->end2_pos + halfReadSize_local) <= regionEnd ; ++fragIt )
				//for (startingRead_forward = *fragIt_forward; *fragIt_forward < fragRangeEnd_forward; ++fragIt_forward)
				{
					//cout <<"ET1F\t";
					if (fragIt->end2_strand == true) //this made the totalRegionNum meaningless. Even the for loop is entered, there might be no frags in this region since the frag might only on the other strand
					{
						fragIt->end2_cuttingSite = regionEnd; //true strand = regionEnd and false strand = regionStart
						fragsCounter += 1;
					}
				}

				++totalRegionNum;
				totalFragsInRegions += fragsCounter;

				++regionIt;
			}
		}
	}



	//false strand
	chrIt_regionMap = cuttingSiteExtMap_local.begin(), chrIt_regionMap_end = cuttingSiteExtMap_local.end();
	fragIt_begin = hybridFragVec_local.begin(), fragIt = hybridFragVec_local.begin(), fragIt_chrStart = hybridFragVec_local.begin(), fragIt_end = hybridFragVec_local.end();

	for ( ; chrIt_regionMap != chrIt_regionMap_end; ++chrIt_regionMap )
	{
		//cout <<"\nentered"<<i<<endl;
		while (fragIt->end2_chr < chrIt_regionMap->first && fragIt != fragIt_end)
		{
			++fragIt;
		}

		if (fragIt == fragIt_end)
		{
			break;
		}
		else if (fragIt->end2_chr > chrIt_regionMap->first)
		{
			continue;
		}
		else if (fragIt->end2_chr == chrIt_regionMap->first)
		{
			fragIt_chrStart = fragIt;
			//while (fragIt != fragIt_begin && fragIt->end2_chr == chrIt_regionMap->first)
			//{
			//	--fragIt;
			//}
			//if (fragIt != fragIt_end && fragIt->end2_chr < chrIt_regionMap->first)
			//{
			//	++fragIt;
			//}

			vector< pair<int, int> >::const_iterator regionIt = chrIt_regionMap->second[false].begin(), regionIt_end = chrIt_regionMap->second[false].end();
			while (regionIt != regionIt_end && fragIt->end2_chr == chrIt_regionMap->first)
			{
				//cout <<"entered while"<<endl;
				regionStart = regionIt->first;
				regionEnd = regionIt->second;

				while ( fragIt != fragIt_end && fragIt->end2_chr == chrIt_regionMap->first && (fragIt->end2_pos + halfReadSize_local) < regionStart ) //--fragIt is to counter the effect of the following ++fragIt, when it would not go into the following while loop
				{
					++fragIt;
					//cout <<"ET1W1\t";
				}

				bool moveBackFlag = false;
				if ( fragIt != fragIt_end && fragIt->end2_chr == chrIt_regionMap->first && fragIt != fragIt_chrStart && (fragIt->end2_pos + halfReadSize_local) >= regionStart )
				{
					--fragIt;
					moveBackFlag = true;
				}

				while ( fragIt != fragIt_end && fragIt->end2_chr == chrIt_regionMap->first && fragIt != fragIt_chrStart && (fragIt->end2_pos + halfReadSize_local) >= regionStart )
				{
					--fragIt;
				}
				//cout <<"OT1W2\t";

				if ( fragIt != fragIt_end && fragIt->end2_chr == chrIt_regionMap->first && moveBackFlag )
				{
					++fragIt;
				}
				else if ( fragIt == fragIt_end || fragIt->end2_chr != chrIt_regionMap->first )
				{
					break;
				}

				//cout <<"OT1W2\t";
				int fragsCounter = 0;
				//int startingFrag = 0;

				for ( ; fragIt != fragIt_end && fragIt->end2_chr == chrIt_regionMap->first && (fragIt->end2_pos + halfReadSize_local) <= regionEnd ; ++fragIt )
				//for (startingRead_forward = *fragIt_forward; *fragIt_forward < fragRangeEnd_forward; ++fragIt_forward)
				{
					//cout <<"ET1F\t";
					if (fragIt->end2_strand == false) //this made the totalRegionNum meaningless. Even the for loop is entered, there might be no frags in this region since the frag might only on the other strand
					{
						fragIt->end2_cuttingSite = regionStart; //true strand = regionEnd and false strand = regionStart
						fragsCounter += 1;
					}
				}

				++totalRegionNum;
				totalFragsInRegions += fragsCounter;

				++regionIt;
			}
		}
	}

	return make_pair(totalRegionNum, totalFragsInRegions);
	//return totalFragsInPeaks;
}



//1) filter out 0 or 1 end mapped to cutting site frags, and keep both end mapped. 2) link frag to the cutting site map 3) remove two ends mapped to the same cutting site or adjacent cutting site if the two ends following the self loop characteristic (end mapped to the precedent cutting site is -, and the other end is +) 4), frag needed to be duplicated since only 1 record of each frag is in the input file (wait till next step to reduce communication between nodes)
int getCuttingSiteFrags(const map< int, vector<int> >& cuttingSitesMap_local, const vector<fragInfo_bothEndsMapped_withCuttingSite>& interactingFrags_local, map< int, map< int, vector<fragInfo_bothEndsMapped_withCuttingSite> > >& cuttingSiteFragsMap_local, vector<int>& mastersIds_local)
{
    return 0;
}


//construct distribution of random interaction genomic distance
int genRanDis(const map< int, map< int, vector<fragInfo_bothEndsMapped_withCuttingSite> > >& cuttingsitefragsmap_local, map< int, vector<int> >& cuttingsitesmap_local, map<int, int>& randismap_local, int binsize_local)
{
	return 0;
}




int freqToProb(const map<int, int>& ranDisMap_local, const int intraChrFrag_local, map<int, double>& ranDisMap_prob_local)
{
	map<int, int>::const_iterator disIt = ranDisMap_local.begin(), disIt_end = ranDisMap_local.end();
	for ( ; disIt != disIt_end; ++disIt )
	{
		ranDisMap_prob_local[disIt->first] = disIt->second * 1.0 / intraChrFrag_local;
	}

	return 0;
}




//find interaction frequency between each pair of ENZYME DIGESTED FRAGMENTS, then search for freq that above the threshold (including the neighbouring frag)
int pairFragsInteractionFrequency(const map< int, map< int, vector<fragInfo_bothEndsMapped_withCuttingSite> > >& cuttingSiteFragsMap_local, map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<fragInfo_bothEndsMapped_withCuttingSite> >& freqMaptoWrite, map<int, int>& chrLen_local) //map<pair<chr1,pair<frag1Start, frag1End> >, pair<chr2,pair<frag2STart, frag2End> > >, pair<number of interaction, vector<interaction frag> > >
{
	freqMaptoWrite.clear();
	map< int, map< int, vector<fragInfo_bothEndsMapped_withCuttingSite> > >::const_iterator chrIt = cuttingSiteFragsMap_local.begin(), chrIt_end = cuttingSiteFragsMap_local.end();
	for ( ; chrIt != chrIt_end; ++chrIt )
	{
		map< int, vector<fragInfo_bothEndsMapped_withCuttingSite> >::const_iterator siteIt_begin = chrIt->second.begin(), siteIt_frag1End = chrIt->second.begin(), siteIt = chrIt->second.begin(), siteIt_end = chrIt->second.end();
		for ( ; siteIt != siteIt_end; ++siteIt )
		{
			vector<fragInfo_bothEndsMapped_withCuttingSite>::const_iterator fragIt = siteIt->second.begin(), fragIt_end = siteIt->second.end();
			for ( ; fragIt != fragIt_end; ++fragIt )
			{
				int frag1_start = 0; //a pair of interacting frag, frag1 is end1, frag2 is end2
				int frag1_end = 0;
				siteIt_frag1End = siteIt;
				if (fragIt->end1_strand)
				{
					frag1_end = siteIt->first;
					if (siteIt != siteIt_begin)
					{
						siteIt_frag1End--;
						frag1_start = siteIt_frag1End->first;
					}
					else
					{
						frag1_start = 1;
					}
				}
				else
				{
					frag1_start = siteIt->first;

					siteIt_frag1End++;
					if (siteIt_frag1End != siteIt_end)
					{
						frag1_end = siteIt_frag1End->first;
					}
					else
					{
						frag1_end = chrLen_local[fragIt->end1_chr];
					}
				}

				map< int, map< int, vector<fragInfo_bothEndsMapped_withCuttingSite> > >::const_iterator chrIt_end2 = cuttingSiteFragsMap_local.find(fragIt->end2_chr);
				if (chrIt_end2 != chrIt_end)
				{
					map< int, vector<fragInfo_bothEndsMapped_withCuttingSite> >::const_iterator siteIt_end2_begin = chrIt_end2->second.begin(), siteIt_frag2End = chrIt_end2->second.begin(), siteIt_end2 = chrIt_end2->second.find(fragIt->end2_cuttingSite), siteIt_end2_end = chrIt_end2->second.end();
					int frag2_start = 0;
					int frag2_end = 0;
					siteIt_frag2End = siteIt_end2;
					if (fragIt->end2_strand)
					{
						frag2_end = siteIt_end2->first;
						if (siteIt_end2 != siteIt_end2_begin)
						{
							siteIt_frag2End--;
							frag2_start = siteIt_frag2End->first;
						}
						else
						{
							frag2_start = 1;
						}
					}
					else
					{
						frag2_start = siteIt_end2->first;

						siteIt_frag2End++;
						if (siteIt_frag2End != siteIt_end2_end)
						{
							frag2_end = siteIt_frag2End->first;
						}
						else
						{
							frag2_end = chrLen_local[fragIt->end2_chr];
						}
					}
					freqMaptoWrite[make_pair(make_pair(chrIt->first, make_pair(frag1_start, frag1_end)), make_pair(chrIt_end2->first, make_pair(frag2_start, frag2_end)))].push_back(*fragIt);
				}
			}
		}
	}

	return 0;
}




int transformSiteVecToFragMap(const map< int, vector<int> >& cuttingSitesMap_local, map< int, map<int, int> >& digestedFragMap_local, map<int, int>& chrLen_local)
{
	map< int, vector<int> >::const_iterator chrIt = cuttingSitesMap_local.begin(), chrIt_end = cuttingSitesMap_local.end();
	for ( ; chrIt != chrIt_end; ++chrIt )
	{
		int prevSite = 1;
		vector<int>::const_iterator siteIt = chrIt->second.begin(), siteIt_end = chrIt->second.end();
		for ( ; siteIt != siteIt_end; ++siteIt )
		{
			digestedFragMap_local[chrIt->first][prevSite] = *siteIt;
			prevSite = *siteIt;
		}
		digestedFragMap_local[chrIt->first][prevSite] = chrLen_local[chrIt->first];
	}

	return 0;
}




//int getInteractingSites(map< int, map< int, vector<fragInfo_bothEndsMapped_withCuttingSite> > >& cuttingSiteFragMap_local, 
int findInteractingSites(map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<fragInfo_bothEndsMapped_withCuttingSite> > fragInteractionMap_local, map< int, map<int, int> >& digestedFragMap_local, const double thres, vector<peakDetail_pairFrag_2D>& peakMaptoWrite_2D)//2Dpeakmap: Map< pair<pair<domain1 chr, domain2 chr>, vector(coordinates and score of each region in this peak)< pair< pair<domain1 position, domain2 position>,score> >, vector<peak property, p value, score etc> >
{
	//cout<<"freq map duplicated"<<endl;
	peakMaptoWrite_2D.clear();

	int countPeaknum = 0;
	const double doubleThres = 2 * thres;

	map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<fragInfo_bothEndsMapped_withCuttingSite> >::iterator fragPairIt = fragInteractionMap_local.begin(), fragPairIt_end = fragInteractionMap_local.end();

	//recurrsively search for neighbour region that has a higher score than the threshold
	while(fragPairIt != fragPairIt_end)
	{
		//bool flag_peak = (regionIt->second >= doubleThres);//first scored bin

		int chr_domain1 = fragPairIt->first.first.first;
		int chr_domain2 = fragPairIt->first.second.first;
		pair<int, int> position_frag1 = fragPairIt->first.first.second;
		pair<int, int> position_frag2 = fragPairIt->first.second.second;
		//int oldPosition_domain1 = regionIt->first.first.second;
		//int oldPosition_domain2 = regionIt->first.second.second;
		int fragPairScore = fragPairIt->second.size();
		int totalScore = fragPairIt->second.size();
		//double old_totalScore = regionIt->second;
		//double aveScore = regionIt->second;
		//double old_aveScore = regionIt->second;

		int fragPairNum = 1;//my $binNo = 0;

		vector< pair< pair< pair<int, int>, pair<int, int> >, int > > peakRegionsCoordinates;
		peakDetail_pairFrag_2D peakInfo;
		vector<fragInfo_bothEndsMapped_withCuttingSite> peakFragInfo;
		peakFragInfo = fragPairIt->second;

		peakRegionsCoordinates.push_back( make_pair(make_pair(position_frag1,position_frag2), fragPairScore) );
		fragInteractionMap_local.erase(fragPairIt);

		if (fragPairScore >= thres)
		{
			searchNeighbourFrags(fragInteractionMap_local, digestedFragMap_local, chr_domain1, chr_domain2, position_frag1, position_frag2, 0, 0, peakRegionsCoordinates, thres, fragPairNum, totalScore, peakFragInfo);
		}

		if (totalScore >= doubleThres)//if regionIt itself or the sum of all scores in the peak is larger than double threshold -> peak
		{
			peakInfo.chrNo_domain1 = chr_domain1;
			peakInfo.chrNo_domain2 = chr_domain2;
			peakInfo.fragPair_num = fragPairNum;
			peakInfo.total_frag = totalScore;
			peakInfo.frag_info = peakFragInfo;
			peakInfo.region_info = peakRegionsCoordinates;
			//peakInfo.getPeakSummit();

			peakMaptoWrite_2D.push_back(peakInfo);

			++countPeaknum;
		}

		fragPairIt = fragInteractionMap_local.begin();
	}

	return countPeaknum;
}




int searchNeighbourFrags(map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<fragInfo_bothEndsMapped_withCuttingSite> >& freqMaptoSearch_local, map< int, map<int, int> >& digestedFragMap_search, int chr_domain1_local, int chr_domain2_local, pair<int, int> position_domain1_local, pair<int, int> position_domain2_local, int parentSetoff_x, int parentSetoff_y, vector< pair< pair< pair<int, int>, pair<int, int> >, int > >& peakRegionsCoordinates_local, double thres_local, int& pairFragNum_local, int& totalScore_local, vector<fragInfo_bothEndsMapped_withCuttingSite>& peakFragInfo_local)
{
	//cout <<chr_domain1_local<<"\t"<<chr_domain2_local<<"\t"<<position_domain1_local.first<<"\t"<<position_domain2_local.first<<"\t"<<parentSetoff_x<<"\t"<<parentSetoff_y<<"\t"<<thres_local<<"\t"<<pairFragNum_local<<"\t"<<totalScore_local<<endl;
	for (int x=-1; x<=1; ++x)
	{
		//cout<<"enter x."<<endl;
		map<int, int>::const_iterator targetFrag_x = digestedFragMap_search[chr_domain1_local].find(position_domain1_local.first);

		if ( x==-1 && targetFrag_x != digestedFragMap_search[chr_domain1_local].begin() && targetFrag_x != digestedFragMap_search[chr_domain1_local].end() ) //search for the previous one
		{
			//cout<<x<<endl;
			targetFrag_x--;
		}
		else if ( x==1 && targetFrag_x != digestedFragMap_search[chr_domain1_local].end() )
		{
			//cout<<x<<endl;
			targetFrag_x++;
			if (targetFrag_x == digestedFragMap_search[chr_domain1_local].end())
			{
				//cout<<"if continue"<<endl;
				continue;
			}
		}
		else if ( x==0 && targetFrag_x != digestedFragMap_search[chr_domain1_local].end() )
		{
			//cout<<x<<endl;
		}
		else
		{
			//cout<<"else continue"<<endl;
			continue;
		}

		for (int y=-1; y<=1; ++y)
		{
			//cout<<"enter y."<<endl;
			map<int, int>::const_iterator targetFrag_y = digestedFragMap_search[chr_domain2_local].find(position_domain2_local.first);

			if ( y==-1 && targetFrag_y != digestedFragMap_search[chr_domain2_local].begin() && targetFrag_y != digestedFragMap_search[chr_domain2_local].end() ) //search for the previous one
			{
				targetFrag_y--;
			}
			else if ( y==1 && targetFrag_y != digestedFragMap_search[chr_domain2_local].end() )
			{
				targetFrag_y++;
				if (targetFrag_y == digestedFragMap_search[chr_domain2_local].end())
				{
					continue;
				}
			}
			else if ( y==0 && targetFrag_y != digestedFragMap_search[chr_domain2_local].end() ) {}
			else
			{
				continue;
			}

			if ( !( (x==0 && y==0) || ( x==(0-parentSetoff_x) && y==(0-parentSetoff_y) ) ) )
			{
				//cout<<"enter xy."<<endl;
				map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<fragInfo_bothEndsMapped_withCuttingSite> >::iterator neighbourPairFrag = freqMaptoSearch_local.find( make_pair( make_pair(chr_domain1_local, *targetFrag_x), make_pair(chr_domain2_local, *targetFrag_y) ) );

				if (neighbourPairFrag != freqMaptoSearch_local.end())
				{
					//cout<<"enter recurrsive."<<endl;
					int neighbourScore = neighbourPairFrag->second.size();
					pair<int, int> neighbourPosition_domain1 = neighbourPairFrag->first.first.second, neighbourPosition_domain2 = neighbourPairFrag->first.second.second;
					vector<fragInfo_bothEndsMapped_withCuttingSite> neighbourPairFragInfo = neighbourPairFrag->second;
					freqMaptoSearch_local.erase(neighbourPairFrag);//delete the element before pass down to recurrsive checking, if not, for example, in 3*3 array, [0,0][0,1][1,1],has value > thres, [0,0]->[0,1]->[1,1]->[0,0]->.....

					if (neighbourScore >= thres_local)
					{
						totalScore_local += neighbourScore;
						++pairFragNum_local;
						peakRegionsCoordinates_local.push_back( make_pair( make_pair( neighbourPosition_domain1, neighbourPosition_domain2 ), neighbourScore ) );
						vector<fragInfo_bothEndsMapped_withCuttingSite>::iterator fragIt = neighbourPairFragInfo.begin(), fragIt_end = neighbourPairFragInfo.end();
						for ( ; fragIt != fragIt_end; ++fragIt)
						{
							peakFragInfo_local.push_back(*fragIt);
						}
						searchNeighbourFrags(freqMaptoSearch_local, digestedFragMap_search, chr_domain1_local, chr_domain2_local, neighbourPosition_domain1, neighbourPosition_domain2, x, y, peakRegionsCoordinates_local, thres_local, pairFragNum_local, totalScore_local, peakFragInfo_local);
					}
				}
			}
		}
	}
	return 0;
}




//find the summit of a peak, if there are multiple paired enzyme cutted frags, choose the one with highest score, if there multiple with highest score then: connect them if continuous, make two peak if seperated
int getPeakSummit(vector<peakDetail_pairFrag_2D>& peakMap_2D_local)
{
	int addPeakNum = 0;

	vector<peakDetail_pairFrag_2D> tempPeakVec;
	vector<peakDetail_pairFrag_2D>::iterator peakIt = peakMap_2D_local.begin(), peakIt_end = peakMap_2D_local.end();
	for ( ; peakIt != peakIt_end; ++peakIt )
	{
		vector< pair< pair< pair<int, int>, pair<int, int> >, int > >::const_iterator regionIt = peakIt->region_info.begin(), regionIt_end = peakIt->region_info.end();
		int maxScore =0;
		for ( ; regionIt != regionIt_end; ++regionIt )
		{
			if (regionIt->second > maxScore)
			{
				maxScore = regionIt->second;
			}
		}

		vector< pair< pair< pair<int, int>, pair<int, int> >, int > > maxScorePairFrags;
		regionIt = peakIt->region_info.begin();
		for ( ; regionIt != regionIt_end; ++regionIt )
		{
			if ( regionIt->second == maxScore )
			{
				maxScorePairFrags.push_back(*regionIt);
			}
		}

		regionIt = maxScorePairFrags.begin();
		regionIt_end = maxScorePairFrags.end();

		peakIt->peakRegion.first.first = regionIt->first.first.first;
		peakIt->peakRegion.first.second = regionIt->first.first.second;
		peakIt->peakRegion.second.first = regionIt->first.second.first;
		peakIt->peakRegion.second.second = regionIt->first.second.second;

		++regionIt;
		for ( ; regionIt != regionIt_end; ++regionIt )
		{
			if ( peakIt->peakRegion.first.first > regionIt->first.first.second || peakIt->peakRegion.first.second < regionIt->first.first.first || peakIt->peakRegion.second.first > regionIt->first.second.second || peakIt->peakRegion.second.second < regionIt->first.second.first )
			{
				if (tempPeakVec.empty())
				{
					peakDetail_pairFrag_2D tempPeakInfo;
					tempPeakInfo.chrNo_domain1 = peakIt->chrNo_domain1;
					tempPeakInfo.chrNo_domain2 = peakIt->chrNo_domain2;
					tempPeakInfo.fragPair_num = 1;
					tempPeakInfo.total_frag = regionIt->second;
					tempPeakInfo.frag_info = peakIt->frag_info; //added peak has no frag information
					tempPeakInfo.region_info.push_back(*regionIt);

					tempPeakInfo.peakRegion.first.first = regionIt->first.first.first;
					tempPeakInfo.peakRegion.first.second = regionIt->first.first.second;
					tempPeakInfo.peakRegion.second.first = regionIt->first.second.first;
					tempPeakInfo.peakRegion.second.second = regionIt->first.second.second;

					tempPeakInfo.sumLogProb = peakIt->sumLogProb;
					//peakInfo.getPeakSummit();

					tempPeakVec.push_back(tempPeakInfo);

					++addPeakNum;
				}
				else
				{
					vector<peakDetail_pairFrag_2D>::reverse_iterator lastAddedPeak = tempPeakVec.rbegin();
					if ( lastAddedPeak->peakRegion.first.first > regionIt->first.first.second || lastAddedPeak->peakRegion.first.second < regionIt->first.first.first || lastAddedPeak->peakRegion.second.first > regionIt->first.second.second || lastAddedPeak->peakRegion.second.second < regionIt->first.second.first )
					{
						peakDetail_pairFrag_2D tempPeakInfo;
						tempPeakInfo.chrNo_domain1 = peakIt->chrNo_domain1;
						tempPeakInfo.chrNo_domain2 = peakIt->chrNo_domain2;
						tempPeakInfo.fragPair_num = 1;
						tempPeakInfo.total_frag = regionIt->second;
						tempPeakInfo.frag_info = peakIt->frag_info; //added peak has no frag information
						tempPeakInfo.region_info.push_back(*regionIt);

						tempPeakInfo.peakRegion.first.first = regionIt->first.first.first;
						tempPeakInfo.peakRegion.first.second = regionIt->first.first.second;
						tempPeakInfo.peakRegion.second.first = regionIt->first.second.first;
						tempPeakInfo.peakRegion.second.second = regionIt->first.second.second;

						tempPeakInfo.sumLogProb = peakIt->sumLogProb;
						//peakInfo.getPeakSummit();

						tempPeakVec.push_back(tempPeakInfo);

						++addPeakNum;
					}
					else
					{
						if (regionIt->first.first.first < lastAddedPeak->peakRegion.first.first)
						{
							lastAddedPeak->peakRegion.first.first = regionIt->first.first.first;
						}
						if (regionIt->first.first.second > lastAddedPeak->peakRegion.first.second)
						{
							lastAddedPeak->peakRegion.first.second = regionIt->first.first.second;
						}
						if (regionIt->first.second.first < lastAddedPeak->peakRegion.second.first)
						{
							lastAddedPeak->peakRegion.second.first = regionIt->first.second.first;
						}
						if (regionIt->first.second.second > lastAddedPeak->peakRegion.second.second)
						{
							lastAddedPeak->peakRegion.second.second = regionIt->first.second.second;
						}
					}

				}

			}
			else
			{
				if (regionIt->first.first.first < peakIt->peakRegion.first.first)
				{
					peakIt->peakRegion.first.first = regionIt->first.first.first;
				}
				if (regionIt->first.first.second > peakIt->peakRegion.first.second)
				{
					peakIt->peakRegion.first.second = regionIt->first.first.second;
				}
				if (regionIt->first.second.first < peakIt->peakRegion.second.first)
				{
					peakIt->peakRegion.second.first = regionIt->first.second.first;
				}
				if (regionIt->first.second.second > peakIt->peakRegion.second.second)
				{
					peakIt->peakRegion.second.second = regionIt->first.second.second;
				}
			}
		}
	}
	peakMap_2D_local.insert(peakMap_2D_local.end(), tempPeakVec.begin(), tempPeakVec.end());

	return addPeakNum;
}




void reader_interactingRegion(const string& fileToread, vector<fragInfo_bothEndsMapped_withCuttingSite>& interactingRegionMap_local)
{
	ifstream inputFile(fileToread.c_str());
	if (!inputFile)
	{
		cout <<"\n"<< "Error opening " << fileToread << "." << endl;
		exit(1);
	}

	char chrNo_end1[3];
	char chrNo_end2[3];
	char start_end1[15];
	char start_end2[15];
	char strand_end1 = '0';
	char strand_end2 = '0';
	char fragType_str[4];

	char* it_lineStr;
	char* it_token;

	char lineStr[512];
	while(inputFile.getline(lineStr,512))
	{
		if (lineStr != NULL)
		{
			fragInfo_bothEndsMapped_withCuttingSite fragInfo;
			it_lineStr = lineStr;

			it_token = chrNo_end1;
			while(*it_lineStr != '\t')
			{
				*it_token = *it_lineStr;
				++it_token;
				++it_lineStr;
			}
			*it_token = '\0';
			++it_lineStr;
			if (chrNo_end1[0] == 'X')
			{
				fragInfo.end1_chr = 23;
			}
			else if (chrNo_end1[0] == 'Y')
			{
				fragInfo.end1_chr = 24;
			}
			else if (chrNo_end1[0] == 'M')
			{
				continue;
			}
			else
			{
				fragInfo.end1_chr = atoi(chrNo_end1);
			}

			it_token = start_end1;
			while(*it_lineStr != '\t')
			{
				*it_token = *it_lineStr;
				++it_token;
				++it_lineStr;
			}
			*it_token = '\0';
			++it_lineStr;
			fragInfo.end1_pos = atoi(start_end1);

			strand_end1 = *it_lineStr;
			if (strand_end1 == '0')
			{
				fragInfo.end1_strand = false;
			}
			else
			{
				fragInfo.end1_strand = true;
			}
			++it_lineStr;++it_lineStr;

			it_token = chrNo_end2;
			while(*it_lineStr != '\t')
			{
				*it_token = *it_lineStr;
				++it_token;
				++it_lineStr;
			}
			*it_token = '\0';
			++it_lineStr;
			if (chrNo_end2[0] == 'X')
			{
				fragInfo.end2_chr = 23;
			}
			else if (chrNo_end2[0] == 'Y')
			{
				fragInfo.end2_chr = 24;
			}
			else if (chrNo_end2[0] == 'M')
			{
				continue;
			}
			else
			{
				fragInfo.end2_chr = atoi(chrNo_end2);
			}

			it_token = start_end2;
			while(*it_lineStr != '\t')
			{
				*it_token = *it_lineStr;
				++it_token;
				++it_lineStr;
			}
			*it_token = '\0';
			++it_lineStr;
			fragInfo.end2_pos = atoi(start_end2);

			strand_end2 = *it_lineStr;
			if (strand_end2 == '0')
			{
				fragInfo.end2_strand = false;
			}
			else
			{
				fragInfo.end2_strand = true;
			}
			++it_lineStr;++it_lineStr;

			it_token = fragType_str;
			while(*it_lineStr != '\t' && *it_lineStr != '\n' && *it_lineStr != '\0')
			{
				*it_token = *it_lineStr;
				++it_token;
				++it_lineStr;
			}
			*it_token = '\0';

			fragInfo.fragType = atoi(fragType_str);

			interactingRegionMap_local.push_back(fragInfo);
		}
	}
}




int parser_domainFile(const string& fileToread, map< int, vector< pair<int, int> > >& peakMaptoWrite)
{
	ifstream inputFile(fileToread.c_str());
	if (!inputFile)
	{
		//char a;
		//cin >>a;
		cout <<"\n"<< "Error opening " << fileToread << "." << endl;
		exit(1);
	}

	char lineStr[20480];
	char chrNo_str[3] = "";
	int chrNo_int = 0;
	char start_str[15];
	int start_int = 0;
	char end_str[15];
	int end_int = 0;
	int col_chrNo = 0;
	int col_start = 1;
	int col_end = 2;

	int column = 0;

	char* it_lineStr;
	char* it_token;

	int peakNum = 0;

	//copy data to a map, which chrNo is the key, the value is a vector of ints(read's starting point)
	//ofstream testBadFile("test");
	while(inputFile.getline(lineStr,20480))
	{
		it_lineStr = lineStr;

		column = 0;

		if (lineStr != NULL)
		{
			//testBadFile<<lineStr<<endl;

			//pass head white space
			while(*it_lineStr == ' '|| *it_lineStr == '\t')
			{
				++it_lineStr;
			}

			while(*it_lineStr != '\0')
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
					if (column == col_chrNo)
					{
						it_token = chrNo_str;
						++it_lineStr;++it_lineStr;++it_lineStr;
						while(*it_lineStr != ' ' && *it_lineStr != '\t' && *it_lineStr != '\0')
						{
							*it_token = *it_lineStr;
							++it_token;
							++it_lineStr;
							//++chrNo_len;
						}
						*it_token = '\0';
					}
					else if (column == col_start)
					{
						it_token = start_str;
						while(*it_lineStr != ' ' && *it_lineStr != '\t' && *it_lineStr != '\0' && *it_lineStr != '\n')
						{
							*it_token = *it_lineStr;
							++it_token;
							++it_lineStr;
						}
						*it_token = '\0';
					}
					else if (column == col_end)
					{
						it_token = end_str;
						while(*it_lineStr != ' ' && *it_lineStr != '\t' && *it_lineStr != '\0' && *it_lineStr != '\n')
						{
							*it_token = *it_lineStr;
							++it_token;
							++it_lineStr;
						}
						*it_token = '\0';
					}
					else
					{
						while(*it_lineStr != ' ' && *it_lineStr != '\t' && *it_lineStr != '\0' && *it_lineStr != '\n')
						{
							++it_lineStr;
						}
					}
				}
			}

			if (chrNo_str[0]!='M')
			{
				if (chrNo_str[0] == 'X')
				{
					chrNo_int = 23;
				}
				else if (chrNo_str[0] == 'Y')
				{
					chrNo_int = 24;
				}
				else
				{
					chrNo_int = atoi(chrNo_str);
				}

				start_int = atoi(start_str);
				end_int = atoi(end_str);

				peakMaptoWrite[chrNo_int].push_back(make_pair(start_int, end_int));

				++peakNum;
			}
		}
	}

	inputFile.close();

	return peakNum;

}




//read enzyme cutting site map file, like HindIII_hg18.bed
int parser_enzymeCuttingMap(const string& inputMapFile, map< int, vector<int> >& cuttingMaptoWrite)
{
	//cout <<"parser_eland started!"<<endl;
	cuttingMaptoWrite.clear();
	int cuttingSiteLen = 0;
	int halfCuttingSiteLen = 0;
	int sitesNum = 0;

	ifstream inputFile(inputMapFile.c_str());
	if (!inputFile)
	{
		cout <<"\n"<< "Error opening " <<inputMapFile<< "." <<endl;
		exit(1);
	}

	char lineStr[512];
	const int col_chrNo = 0;
	const int col_start = 1;
	const int col_end = 2;
	char chrNo[6];
	int chrNo_int = 0;
	char start_str[15];
	int start_int = 0;
	char end_str[15];
	int end_int = 0;
	const int max_column = 1;
	int column = 0;

	bool format_cuttingSiteLen = 0;
	bool format_start = 0;
	bool format_end = 0;
	bool format_chrNo = 0;

	int errorCode = 0;

	//istringstream buf_lineStr;
	//string dump = "";
	//string dump1 = "";
	//string dump2 = "";
	//string dump3 = "";
	//string chrNo = "";

	char* it_lineStr;
	char* it_token;

	//bool random_flag = false;


	//determine if the file format is appropriate and also determine the read length
	inputFile.getline(lineStr, 512);
	if (lineStr != NULL)
	{
		//cout<<lineStr<<endl;
		it_lineStr = lineStr;
		//cout<<*it_lineStr<<endl;
		//pass head white space
		while(*it_lineStr == ' '|| *it_lineStr == '\t')
		{
			++it_lineStr;
			//cout<<"1"<<endl;
		}

		while(*it_lineStr != '\0' && column <= 2)//need to test end position which is at column 2, as well
		{
			if (*it_lineStr == ' ' || *it_lineStr == '\t')
			{
				//remove white space
				while(*it_lineStr == ' ' || *it_lineStr == '\t')
				{
					++it_lineStr;
					//cout<<"2"<<endl;
				}

				++column;
				//cout<<"3"<<endl;
			}
			else
			{
				//cout<<"4"<<endl;
				switch (column)
				{
					case col_chrNo:
						if (*it_lineStr != 'c')
						{
							format_chrNo = 1;
							errorCode += 100;
							cout <<"Format error detected in first line at field chromosome ID: "<<*it_lineStr<<endl;
							break;
						}
						++it_lineStr;
						if (*it_lineStr != 'h')
						{
							format_chrNo = 1;
							errorCode += 1000;
							cout <<"Format error detected in first line at field chromosome ID: "<<*it_lineStr<<endl;
							break;
						}
						++it_lineStr;
						if (*it_lineStr != 'r')
						{
							format_chrNo = 1;
							errorCode += 10000;
							cout <<"Format error detected in first line at field chromosome ID: "<<*it_lineStr<<endl;
							break;
						}
						while(*it_lineStr != ' ' && *it_lineStr != '\t' && *it_lineStr != '\n')
						{
							++it_lineStr;
						}
						break;

					case col_start:
						it_token = start_str;//read in start_str
						while(*it_lineStr != ' ' && *it_lineStr != '\t' && *it_lineStr != '\n')
						{
							if (*it_lineStr < '0' || *it_lineStr > '9')
							{
								format_start = 1;
								errorCode += 100000;
								cout <<"Format error detected in first line at field start position: "<<*it_lineStr<<endl;
								break;
							}
							else
							{
								*it_token = *it_lineStr;
								++it_token;
							}
							++it_lineStr;
						}
						*it_token = '\0';
						start_int = atoi(start_str);
						break;

					case col_end:
						it_token = end_str;
						while(*it_lineStr != ' ' && *it_lineStr != '\t' && *it_lineStr != '\n')
						{
							if (*it_lineStr < '0' || *it_lineStr > '9')
							{
								format_end = 1;
								errorCode += 100000;
								cout <<"Format error detected in first line at field start position: "<<*it_lineStr<<endl;
								break;
							}
							else
							{
								*it_token = *it_lineStr;
								++it_token;
							}
							++it_lineStr;
						}
						*it_token = '\0';
						end_int = atoi(end_str);
						break;

					default:
						//cout<<"default"<<endl;
						while(*it_lineStr != ' ' && *it_lineStr != '\t' && *it_lineStr != '\n')
						{
							++it_lineStr;
						}
						break;
				}

			}
		}

		//cout<<"chrNoSp_len: "<<chrNoSp_len<<"\treadSize: "<<readSize<<endl;

		cuttingSiteLen = end_int - start_int;
		//cout<<"end_int\t"<<end_int<<endl;
		//cout<<"start_int\t"<<start_int<<endl;
		//cout<<"cuttingSiteLen\t"<<cuttingSiteLen<<endl;
		halfCuttingSiteLen = cuttingSiteLen/2;
		if (cuttingSiteLen > 20 || format_start || format_end || format_chrNo)
		{
			cout <<"\n\n\t"<<inputMapFile<<": File format not recognized.\n\n\
	Acceptable aligned file format: TGAGTGAGGTGTGGGCTCCACACCC 12500 1 chr19:4124128 F TGAGTGAGGTGTGGGCTCCACACCC 9359\n\
	Acceptable extended eland file format: >HWI-EAS435:4:1:2:1706#0/1      GGATGGAGTGCAGTGCTGCAATCATGGTTCACTGAA    0:1:86  chr12.fa:121208577R15G20\n\
	Acceptable sorted file format: Seq_4383668     ChIPSeq	1       +       559765  559799\n\
	error code: "<<errorCode<<endl;
			exit(1);
		}
	}
	//cout << readSize << endl;
	
	inputFile.clear();
	inputFile.seekg(0);

	int readsNum = 0;

	//copy data to a map, which chrNo is the key, the value is a vector of ints(read's starting point)
	while(inputFile.getline(lineStr,512))
	{
		it_lineStr = lineStr;
		//cout<<*it_lineStr<<endl;
		column = 0;

		if (lineStr != NULL && *it_lineStr != '#')
		{
			++readsNum;

			//pass head white space
			while(*it_lineStr == ' '|| *it_lineStr == '\t')
			{
				++it_lineStr;
				//cout<<"1"<<endl;
			}

			while(*it_lineStr != '\0' && column <= max_column)
			{
				if (*it_lineStr == ' ' || *it_lineStr == '\t')
				{
					//remove white space
					while(*it_lineStr == ' ' || *it_lineStr == '\t')
					{
						++it_lineStr;
						//cout<<"2"<<endl;
					}

					++column;
					//cout<<"3"<<endl;
				}
				else
				{
					//cout<<"4"<<endl;
					switch (column)
					{
						case col_chrNo:
							//cout<<"6"<<endl;
							it_token = chrNo;//read in chrNo
							++it_lineStr;++it_lineStr;++it_lineStr;//skip "chr"
							while(*it_lineStr != ' ' && *it_lineStr != '\t' && *it_lineStr != '\n')
							{
								//cout<<"before"<<endl;
								*it_token = *it_lineStr;
								++it_token;
								++it_lineStr;
							}
							*it_token = '\0';
							break;

						case col_start:
							it_token = start_str;//read in start_str
							while(*it_lineStr != ' ' && *it_lineStr != '\t' && *it_lineStr != '\n')
							{
								//cout<<"after"<<endl;
								*it_token = *it_lineStr;
								++it_token;
								++it_lineStr;
							}
							*it_token = '\0';
							break;

						case col_end:
							it_token = end_str;//read in start_str
							while(*it_lineStr != ' ' && *it_lineStr != '\t' && *it_lineStr != '\n')
							{
								//cout<<"after"<<endl;
								*it_token = *it_lineStr;
								++it_token;
								++it_lineStr;
							}
							*it_token = '\0';
							break;

						default:
							//cout<<"default"<<endl;
							while(*it_lineStr != ' ' && *it_lineStr != '\t' && *it_lineStr != '\n')
							{
								++it_lineStr;
							}
							break;
					}
				}

			}

			if (chrNo[0] != 'M')
			{
				if (chrNo[0] == 'X')
				{
					chrNo_int = 23;
				}
				else if (chrNo[0] == 'Y')
				{
					chrNo_int = 24;
				}
				else
				{
					chrNo_int = atoi(chrNo);
				}


				start_int = atoi(start_str);
				start_int += halfCuttingSiteLen;

				cuttingMaptoWrite[chrNo_int].push_back(start_int);

				//MPI_Barrier(MPI_COMM_WORLD);
				//MPI_Bcast(&start_int, 1, MPI_INT, mastersIds_local[0], MPI_COMM_WORLD);
				//MPI_Barrier(MPI_COMM_WORLD);

				++sitesNum;
			}
		}
	}
	//cout << "parser_eland done!"<<endl;
	inputFile.close();

	return sitesNum;
}




int writer_interactionFreqMap(map< pair< pair<int, int>, pair<int, int> >, int >& freqMaptoRead, int binSize_local, const string& freqFiletoWrite)
{
	ofstream freqFile(freqFiletoWrite.c_str());

	map< pair< pair<int, int>, pair<int, int> >, int >::const_iterator binIt = freqMaptoRead.begin(), binIt_end = freqMaptoRead.end();

	for ( ; binIt != binIt_end; ++binIt )
	{
		freqFile <<binIt->first.first.first<<"\t"<<binIt->first.first.second*binSize_local+1<<"\t"<<binIt->first.second.first<<"\t"<<binIt->first.second.second*binSize_local+1<<"\t"<<binIt->second<<endl;
	}

	freqFile.close();

	return 0;
}




int writer_2Dpeak(map< pair< pair<int, int>, vector< pair< pair<int, int>, int> > >, vector<int> >& peakMaptoRead_2D, const string& peakFiletoWrite_2D, const int binSize_local)
{
	ofstream output2DpeakFile(peakFiletoWrite_2D.c_str());
	map< pair< pair<int, int>, vector< pair< pair<int, int>, int> > >, vector<int> >::const_iterator peakIt = peakMaptoRead_2D.begin(), peakIt_end = peakMaptoRead_2D.end();
	for ( ; peakIt != peakIt_end; ++peakIt )
	{
		output2DpeakFile <<peakIt->first.first.first<<"\t"<<peakIt->first.first.second<<"\t"<<peakIt->second[0]<<"\t"<<peakIt->second[1];
		vector< pair< pair<int, int>, int > >::const_iterator regionIt = peakIt->first.second.begin(), regionIt_end = peakIt->first.second.end();
		for ( ; regionIt != regionIt_end; ++regionIt )
		{
			output2DpeakFile <<"\t"<<(regionIt->first.first)*binSize_local+1<<"\t"<<(regionIt->first.second)*binSize_local+1<<"\t"<<regionIt->second;
		}
		output2DpeakFile<<endl;
	}
	output2DpeakFile.close();

	return 0;
}




int writer_2Dpeak(vector<peakDetail_pairFrag_2D>& peakMaptoRead_2D, const string& peakFiletoWrite_2D)
{
	ofstream output2DpeakFile(peakFiletoWrite_2D.c_str());
	vector<peakDetail_pairFrag_2D>::const_iterator peakIt = peakMaptoRead_2D.begin(), peakIt_end = peakMaptoRead_2D.end();
	for ( ; peakIt != peakIt_end; ++peakIt )
	{
		output2DpeakFile <<">\t"<<peakIt->chrNo_domain1<<"\t"<<peakIt->chrNo_domain2<<"\t"<<peakIt->fragPair_num<<"\t"<<peakIt->total_frag<<"\t"<<(peakIt->peakRegion.first.first)<<"\t"<<(peakIt->peakRegion.first.second)<<"\t"<<peakIt->peakRegion.second.first<<"\t"<<peakIt->peakRegion.second.second<<"\t"<<peakIt->sumLogProb;
		vector< pair< pair< pair<int, int>, pair<int, int> >, int > >::const_iterator regionIt = peakIt->region_info.begin(), regionIt_end = peakIt->region_info.end();
		for ( ; regionIt != regionIt_end; ++regionIt )
		{
			output2DpeakFile <<"\t"<<(regionIt->first.first.first)<<"\t"<<(regionIt->first.first.second)<<"\t"<<(regionIt->first.second.first)<<"\t"<<(regionIt->first.second.second)<<"\t"<<regionIt->second;
		}
		output2DpeakFile<<endl;

		vector<fragInfo_bothEndsMapped_withCuttingSite>::const_iterator fragIt = peakIt->frag_info.begin(), fragIt_end = peakIt->frag_info.end();
		for ( ; fragIt != fragIt_end; ++fragIt)
		{
			output2DpeakFile <<fragIt->end1_chr<<"\t"<<fragIt->end1_pos<<"\t"<<fragIt->end1_strand<<"\t"<<fragIt->end1_cuttingSite<<"\t"<<fragIt->end2_chr<<"\t"<<fragIt->end2_pos<<"\t"<<fragIt->end2_strand<<"\t"<<fragIt->end2_cuttingSite<<"\t"<<fragIt->fragType<<endl;
		}
	}
	output2DpeakFile.close();

	return 0;
}




int writer_hybridFrags(const vector<fragInfo_bothEndsMapped_withCuttingSite>& fragMaptoRead, const string& fragFiletoWrite)
{
	ofstream outputFragFile(fragFiletoWrite.c_str());
	vector<fragInfo_bothEndsMapped_withCuttingSite>::const_iterator fragIt = fragMaptoRead.begin(), fragIt_end = fragMaptoRead.end();
	for ( ; fragIt != fragIt_end; ++fragIt )
	{
		outputFragFile <<fragIt->end1_chr<<"\t"<<fragIt->end1_pos<<"\t"<<fragIt->end1_strand<<"\t"<<fragIt->end1_cuttingSite<<"\t"<<fragIt->end2_chr<<"\t"<<fragIt->end2_pos<<"\t"<<fragIt->end2_strand<<"\t"<<fragIt->end2_cuttingSite<<"\t"<<fragIt->fragType<<endl;
	}
	outputFragFile.close();

	return 0;
}




int writer_cuttingSiteFrags(const map< int, map< int, vector<fragInfo_bothEndsMapped_withCuttingSite> > >& cuttingSiteFragsMap_local, const string& siteFragFiletoWrite)
{
	ofstream outputFragFile(siteFragFiletoWrite.c_str());
	map< int, map< int, vector<fragInfo_bothEndsMapped_withCuttingSite> > >::const_iterator chrIt = cuttingSiteFragsMap_local.begin(), chrIt_end = cuttingSiteFragsMap_local.end();
	for ( ; chrIt != chrIt_end; ++chrIt )
	{
		map< int, vector<fragInfo_bothEndsMapped_withCuttingSite> >::const_iterator siteIt = chrIt->second.begin(), siteIt_end = chrIt->second.end();
		for ( ; siteIt != siteIt_end; ++siteIt )
		{
			outputFragFile <<">\t"<<chrIt->first<<"\t"<<siteIt->first<<endl;
			vector<fragInfo_bothEndsMapped_withCuttingSite>::const_iterator fragIt = siteIt->second.begin(), fragIt_end = siteIt->second.end();
			for ( ; fragIt != fragIt_end; ++fragIt )
			{
				outputFragFile <<fragIt->end1_chr<<"\t"<<fragIt->end1_pos<<"\t"<<fragIt->end1_strand<<"\t"<<fragIt->end1_cuttingSite<<"\t"<<fragIt->end2_chr<<"\t"<<fragIt->end2_pos<<"\t"<<fragIt->end2_strand<<"\t"<<fragIt->end2_cuttingSite<<"\t"<<fragIt->fragType<<endl;
			}
		}
	}
	outputFragFile.close();

	return 0;
}




int write_fragInteractionMap(map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<fragInfo_bothEndsMapped_withCuttingSite> >& fragInteractionMap_local, const string& fragInteractionFiletoWrite)
{
	ofstream outputFragInteractionFile(fragInteractionFiletoWrite.c_str());
	map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<fragInfo_bothEndsMapped_withCuttingSite> >::const_iterator fragPairIt = fragInteractionMap_local.begin(), fragPairIt_end = fragInteractionMap_local.end();
	for ( ; fragPairIt != fragPairIt_end; ++fragPairIt )
	{
		outputFragInteractionFile <<">\t"<<fragPairIt->first.first.first<<"\t"<<fragPairIt->first.first.second.first<<"\t"<<fragPairIt->first.first.second.second<<"\t"<<fragPairIt->first.second.first<<"\t"<<fragPairIt->first.second.second.first<<"\t"<<fragPairIt->first.second.second.second<<endl;
		vector<fragInfo_bothEndsMapped_withCuttingSite>::const_iterator fragIt = fragPairIt->second.begin(), fragIt_end = fragPairIt->second.end();
		for ( ; fragIt != fragIt_end; ++fragIt )
		{
			outputFragInteractionFile <<fragIt->end1_chr<<"\t"<<fragIt->end1_pos<<"\t"<<fragIt->end1_strand<<"\t"<<fragIt->end1_cuttingSite<<"\t"<<fragIt->end2_chr<<"\t"<<fragIt->end2_pos<<"\t"<<fragIt->end2_strand<<"\t"<<fragIt->end2_cuttingSite<<"\t"<<fragIt->fragType<<endl;
		}
	}
	outputFragInteractionFile.close();

	return 0;
}




int writer_disDistribution(const map<int, int>& distributionMap_local, const string& distributionFiletoWrite)
{
	ofstream outputFile(distributionFiletoWrite.c_str());

	map<int, int>::const_iterator binIt = distributionMap_local.begin(), binIt_end = distributionMap_local.end();
	for ( ; binIt != binIt_end; ++binIt )
	{
		outputFile<<binIt->first<<"\t"<<binIt->second<<endl;
	}

	outputFile.close();

	return 0;
}




int writer_disDistribution_prob(const map<int, double>& distributionMap_local, const string& distributionFiletoWrite)
{
	ofstream outputFile(distributionFiletoWrite.c_str());

	map<int, double>::const_iterator binIt = distributionMap_local.begin(), binIt_end = distributionMap_local.end();
	for ( ; binIt != binIt_end; ++binIt )
	{
		outputFile<<binIt->first<<"\t"<<binIt->second<<endl;
	}

	outputFile.close();

	return 0;
}




int write_fragInteractionMatrix(map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<fragInfo_bothEndsMapped_withCuttingSite> >& fragInteractionFreqMap_fun, map< int, map<int, int> >& cuttingSiteFragsMap_fun, int start_fun, int end_fun, string& file_fragContactMap_fun)
{
	ofstream outputFile(file_fragContactMap_fun.c_str());

	map< int, map<int, int> >::const_iterator chrIt = cuttingSiteFragsMap_fun.begin(), chrIt_end = cuttingSiteFragsMap_fun.end();
	for ( ; chrIt != chrIt_end; ++chrIt )
	{
		if ( chrIt->first == 11 )
		{
		map<int, int>::const_iterator fragIt_1 = chrIt->second.begin(), fragIt_1_end = chrIt->second.end();
		for ( ; fragIt_1 != fragIt_1_end; ++fragIt_1 )
		{
			if (fragIt_1->first >= start_fun && fragIt_1->second < end_fun)
			{
				pair< int, pair<int, int> > frag1 = make_pair(chrIt->first, make_pair(fragIt_1->first, fragIt_1->second));
				map<int, int>::const_iterator fragIt_2 = chrIt->second.begin(), fragIt_2_end = chrIt->second.end();
				for ( ; fragIt_2 != fragIt_2_end; ++fragIt_2 )
				{
					if (fragIt_2->first >= start_fun && fragIt_2->second < end_fun)
					{
						pair< int, pair<int, int> > frag2 = make_pair(chrIt->first, make_pair(fragIt_2->first, fragIt_2->second));
						map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<fragInfo_bothEndsMapped_withCuttingSite> >::const_iterator interactionIt = fragInteractionFreqMap_fun.find(make_pair(frag1, frag2));
						if (interactionIt != fragInteractionFreqMap_fun.end())
						{
							outputFile<<interactionIt->second.size()<<"\t";
						}
						else
						{
							outputFile<<"0\t";
						}
					}
				}
				outputFile<<endl;
			}
		}
		}
	}

	outputFile.close();

	return 0;
}




//int distributionDomainMap(string domainFile_fun, vector<int>& mastersIds_local)
//{
//	map< int, map< pair<int, int>, int > > domainMap;
//	int domainNum_local = parser_domainFile(domainFile_fun, domainMap);
//
//	for ( map< int, map< pair<int, int>, int > >::const_iterator chrIt = domainMap.begin(), chrIt_end = domainMap.end(); chrIt != chrIt_end; ++chrIt )
//	{
//		for ( map< pair<int, int>, int >::const_iterator domainIt = chrIt->second.begin(), domainIt_end = chrIt->second.end(); domainIt != domainIt_end; ++domainIt )
//		{
//			int domainInfo[2];
//			domainInfo[0] = domainIt->first.first;
//			domainInfo[1] = domainIt->first.second;
//			MPI_Send(domainInfo, 2, MPI_INT, mastersIds_local[chrIt->first], INTTAG, MPI_COMM_WORLD);
//		}
//	}
//	vector<int>::iterator nodeIt = mastersIds_local.begin(), nodeIt_end = mastersIds_local.end();
//	++nodeIt;
//	for ( ; nodeIt != nodeIt_end; ++nodeIt )
//	{
//		MPI_Send(0, 0, MPI_INT, *nodeIt, DIETAG, MPI_COMM_WORLD);
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	return 0;
//}



int getPvalCutoff(double* pvalVec_local, int totalSize_local, double FDR_local, double& pvalCutoff_local)
{
	try
	{
		int contructR_argc = 1;
		char* constructR_argv[1];
		RInside R(contructR_argc, constructR_argv);          // create an embedded R instance, without arguments product compilation error "parseEvalQ is of non-class type"

		string rcommand = "suppressMessages(library(qvalue))";
		R.parseEvalQ(rcommand);              // load library, no return value

		Rcpp::NumericVector pval_r(pvalVec_local, &pvalVec_local[totalSize_local]);
		R["pval_R"] = pval_r;
		cout<<"pvalues are converted for R."<<endl;
		rcommand = "suppressMessages(qObj <- qvalue(as.matrix(pval_R)))";
		R.parseEvalQ(rcommand);
		rcommand = "suppressMessages(pdf(\"test.pdf\"))";
		R.parseEvalQ(rcommand);
		rcommand = "suppressMessages(qplot(qObj))";
		R.parseEvalQ(rcommand);
		rcommand = "suppressMessages(dev.off())";
		R.parseEvalQ(rcommand);
		rcommand = "suppressMessages(qwrite(qObj, filename=\"test.txt\"))";
		R.parseEvalQ(rcommand);

		ostringstream FDRbuffer;
		FDRbuffer<<FDR_local;
		rcommand = "suppressMessages(max(qObj$pvalues[qObj$qvalues <= " + FDRbuffer.str() + "]))";
		FDRbuffer.clear();
		pvalCutoff_local = Rcpp::as< double >(R.parseEval(rcommand));

	} catch(std::exception& ex) {
		std::cerr << "RInside exception caught: " << ex.what() << std::endl;
	} catch(...) {
		std::cerr << "RInside unknown exception caught" << std::endl;
	}

	return 0;
}



#endif

