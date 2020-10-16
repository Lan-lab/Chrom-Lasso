#ifndef HIC_MPLD_POSTDOC_CPP
#define HIC_MPLD_POSTDOC_CPP

//component number could be determined using peak scan + 3 (from 1 to hard threshold peak number + 3)
//set start locus at the prescaned peaks
//bayes' rule could be used to update the parameters, like in BALM (read papers about bayes' rule and see if it applicable)

#include "HiC_mixturePLD_postdoc.h"
#include "HiC_mixturePLD_boss.h"


bool sort2Dregions(const regionDetail_pairSite_2D& lhs, const regionDetail_pairSite_2D& rhs)
{
	if (lhs.chrNo1!=rhs.chrNo1)
	{
		return lhs.chrNo1<rhs.chrNo1;
	}
	else if (lhs.domain1!=rhs.domain1)
	{
		return lhs.domain1<rhs.domain1;
	}
	else if (lhs.summitPairSites.first!=rhs.summitPairSites.first)
	{
		return lhs.summitPairSites.first<rhs.summitPairSites.first;
	}
	else if (lhs.chrNo2!=rhs.chrNo2)
	{
		return lhs.chrNo2<rhs.chrNo2;
	}
	else if (lhs.domain2!=rhs.domain2)
	{
		return lhs.domain2<rhs.domain2;
	}
	else
	{
		return lhs.summitPairSites.second<rhs.summitPairSites.second;
	}
}




int postdoc(vector<fragInfo_bothEndsMapped_withCuttingSite> interactingFrags, map<int, int> chrlenMap_local, const int modelBinSize_local, const string file_source_local)
{
	return 0;
}



int receiveDomainMap(vector< pair<int, int> >& domainMap_fun)
{
	return 0;
}





int write_domainFragInteractionMatrix(map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<fragInfo_bothEndsMapped_withCuttingSite> >& fragInteractionFreqMap_fun, map<int, int>& digestedFragsMap_fun, vector< pair<int, int> >& domainMap_fun, int jobChr_fun, string& file_domainFragContactMap_fun)
{
	ofstream outputFile(file_domainFragContactMap_fun.c_str());

	int counterIntra = 0;
	int counterDomain = 1;
	outputFile<<">domain "<<counterDomain<<endl;
	int domainFragNum = 0;
	int intraDomainCounter = 0;

	map<int, int>::const_iterator fragIt = digestedFragsMap_fun.begin(), fragIt_end = digestedFragsMap_fun.end();
	vector< pair<int, int> >::iterator domainIt = domainMap_fun.begin(), domainIt_end = domainMap_fun.end();

	map<int, int> domainInteractionFreqMatrix;
	while (domainIt != domainIt_end && fragIt != fragIt_end)
	{
		//cout<<"11"<<endl;
		int domainStart = domainIt->first;
		int domainEnd = domainIt->second;
		int fragStart = fragIt->first;
		int fragEnd = fragIt->second;

		if (fragStart > domainEnd)
		{
			outputFile<<">domain "<<counterDomain<<" start "<<domainStart<<" end "<<domainEnd<<" fragNum "<<domainFragNum<<" totalFreq "<<intraDomainCounter<<endl;
			counterDomain++;
			domainFragNum = 0;
			counterIntra += intraDomainCounter;
			intraDomainCounter = 0;
			//cout<<"10"<<endl;
			++domainIt;
			if (domainIt != domainIt_end)
			{
				outputFile<<">domain "<<counterDomain<<endl;
			}
			continue;
		}
		//check distance
		else if ( (fragStart >= domainStart && fragEnd <= domainEnd) )
		{
			//cout<<"5"<<endl;
			map<int, int>::const_iterator domainStartFrag = fragIt;
			while (fragIt != fragIt_end && fragIt->second <= domainEnd)
			{
				++fragIt;
				++domainFragNum;
			}
			for (map<int, int>::const_iterator domainFragIt1 = domainStartFrag; domainFragIt1 != fragIt; ++domainFragIt1)
			{
				pair< int, pair<int, int> > frag1 = make_pair(jobChr_fun, make_pair(domainFragIt1->first, domainFragIt1->second));
				for (map<int, int>::const_iterator domainFragIt2 = domainStartFrag; domainFragIt2 != fragIt; ++domainFragIt2)
				{
					pair< int, pair<int, int> > frag2 = make_pair(jobChr_fun, make_pair(domainFragIt2->first, domainFragIt2->second));
					map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<fragInfo_bothEndsMapped_withCuttingSite> >::const_iterator interactionIt = fragInteractionFreqMap_fun.find(make_pair(frag1, frag2));
					if (interactionIt != fragInteractionFreqMap_fun.end())
					{
						int freq = interactionIt->second.size();
						intraDomainCounter += freq;
						outputFile<<freq<<"\t";
					}
					else
					{
						//outputFile<<fragIt->first<<"\t"<<fragIt->second<<"\t";
						outputFile<<"0\t";
					}
				}
				outputFile<<endl;
			}
		}
		if (fragIt != fragIt_end)
		{
			++fragIt;
		}
		else
		{
			outputFile<<">domain "<<counterDomain<<" start "<<domainStart<<" end "<<domainEnd<<" fragNum "<<domainFragNum<<" totalFreq "<<intraDomainCounter<<endl;
			counterIntra += intraDomainCounter;
		}
	}

	cout<<"Total intra domain interaction: "<<counterIntra<<endl;

	outputFile.close();

	return counterIntra;
}



//one chromosome at a time
int getDomainCSinteractionMatrix(map< pair< pair<int, int>, pair<int, int> >, int >& cuttingSiteInteractionMap_fun, map< int, vector<int> >& cuttingSiteMap_fun, vector< pair<int, int> >& domainMap_fun, int jobChr_fun, vector< vector<int> >& domainSitesMap_fun, vector< vector< vector<int> > >& domainCSinterFreq_fun, map< int, map<int, int> >& CStotalInterMap_fun)
{
	//ofstream outputFile(file_domainCScontactMap_fun.c_str());

	int counterIntra = 0;
	int counterDomain = 1;
	int domainSiteNum = 0;
	int intraDomainCounter = 0;

	for ( map< int, vector<int> >::const_iterator chrIt = cuttingSiteMap_fun.begin(), chrIt_end = cuttingSiteMap_fun.end(); chrIt != chrIt_end; ++chrIt )
	{
		vector<int>::const_iterator siteIt = chrIt->second.begin(), siteIt_end = chrIt->second.end();
		for ( ; siteIt != siteIt_end; ++siteIt )
		{
			CStotalInterMap_fun[chrIt->first][*siteIt] = 0;
		}
	}

	//outputFile<<">domain "<<counterDomain<<endl;

	vector<int>::const_iterator siteIt = cuttingSiteMap_fun[jobChr_fun].begin(), siteIt_end = cuttingSiteMap_fun[jobChr_fun].end();
	vector< pair<int, int> >::iterator domainIt = domainMap_fun.begin(), domainIt_end = domainMap_fun.end();

	while (domainIt != domainIt_end && siteIt != siteIt_end)
	{
		//cout<<"11"<<endl;
		int domainStart = domainIt->first;
		int domainEnd = domainIt->second;

		if (*siteIt > domainEnd)
		{
			//outputFile<<">domain "<<counterDomain<<" start "<<domainStart<<" end "<<domainEnd<<" fragNum "<<domainSiteNum<<" totalFreq "<<intraDomainCounter<<endl;
			counterDomain++;
			domainSiteNum = 0;
			counterIntra += intraDomainCounter;
			intraDomainCounter = 0;
			//cout<<"10"<<endl;
			++domainIt;
			//if (domainIt != domainIt_end)
			//{
			//	outputFile<<">domain "<<counterDomain<<endl;
			//}
			continue;
		}
		//check distance
		else if (*siteIt < domainStart)
		{
			while (siteIt!=siteIt_end && *siteIt<domainStart)
			{
				++siteIt;
			}
		}
		else if ( (*siteIt >= domainStart && *siteIt <= domainEnd) )
		{
			//cout<<"5"<<endl;
			vector<int>::const_iterator domainStartSite = siteIt;
			while (siteIt != siteIt_end && *siteIt <= domainEnd)
			{
				++siteIt;
				++domainSiteNum;
			}
			vector<int> tempDomainSites;
			vector< vector<int> > tempDomainInterFreq;
			for (vector<int>::const_iterator domainSiteIt1 = domainStartSite; domainSiteIt1 != siteIt; ++domainSiteIt1)
			{
				tempDomainSites.push_back(*domainSiteIt1);
				vector<int> tempSiteInterFreq;
				pair<int, int> site1 = make_pair(jobChr_fun, *domainSiteIt1);
				for (vector<int>::const_iterator domainSiteIt2 = domainStartSite; domainSiteIt2 != siteIt; ++domainSiteIt2)
				{
					pair<int, int> site2 = make_pair(jobChr_fun, *domainSiteIt2);
					map< pair< pair<int, int>, pair<int, int> >, int >::const_iterator interactionIt = cuttingSiteInteractionMap_fun.find(make_pair(site1, site2));
					if (interactionIt != cuttingSiteInteractionMap_fun.end())
					{
						int freq = interactionIt->second;
						intraDomainCounter += freq;
						tempSiteInterFreq.push_back(freq);
						//outputFile<<freq<<"\t";
					}
					else
					{
						//outputFile<<siteIt->first<<"\t"<<siteIt->second<<"\t";
						tempSiteInterFreq.push_back(0);
						//outputFile<<"0\t";
					}
				}
				tempDomainInterFreq.push_back(tempSiteInterFreq);
				//outputFile<<endl;
			}
			domainSitesMap_fun.push_back(tempDomainSites);
			domainCSinterFreq_fun.push_back(tempDomainInterFreq);
		}
		if (siteIt == siteIt_end)
		{
			//outputFile<<">domain "<<counterDomain<<" start "<<domainStart<<" end "<<domainEnd<<" siteNum "<<domainSiteNum<<" totalFreq "<<intraDomainCounter<<endl;
			counterIntra += intraDomainCounter;
		}
	}

	map<int, int> CStotalIntraMap;
	siteIt = cuttingSiteMap_fun[jobChr_fun].begin();
	for ( ; siteIt != siteIt_end; ++siteIt )
	{
		CStotalIntraMap[*siteIt] = 0;
	}
	map< pair< pair<int, int>, pair<int, int> >, int >::const_iterator interactionIt = cuttingSiteInteractionMap_fun.begin(), interactionIt_end = cuttingSiteInteractionMap_fun.end();
	for ( ; interactionIt != interactionIt_end; ++interactionIt )
	{
		//since one pairend reads only have one record in the inputfile, both end should be checked when summing the total interaction frequency for the cutting sites.
		CStotalIntraMap[interactionIt->first.first.second] += interactionIt->second;
		CStotalIntraMap[interactionIt->first.second.second] += interactionIt->second;
	}
	//string file_CStotalIntra = file_domainCScontactMap_fun + "total";
	//ofstream outputCS(file_CStotalIntra.c_str());
	//for ( map<int, int>::const_iterator csIt = CStotalIntraMap.begin(), csIt_end = CStotalIntraMap.end(); csIt != csIt_end; ++csIt )
	//{
	//	outputCS<<jobChr_fun<<"\t"<<csIt->first<<"\t"<<csIt->second<<endl;
	//}
	//outputCS.close();

	map< pair< pair<int, int>, pair<int, int> >, int >::const_iterator sitePairIt = cuttingSiteInteractionMap_fun.begin(), sitePairIt_end = cuttingSiteInteractionMap_fun.end();
	for ( ; sitePairIt != sitePairIt_end; ++sitePairIt )
	{
		CStotalInterMap_fun[sitePairIt->first.first.first][sitePairIt->first.first.second] += sitePairIt->second;
		CStotalInterMap_fun[sitePairIt->first.second.first][sitePairIt->first.second.second] += sitePairIt->second;
	}
	//for (map< int, map<int, int> >::const_iterator chrIt = CStotalInterMap.begin(), chrIt_end = CStotalInterMap.end(); chrIt != chrIt_end; ++chrIt)
	//{
	//	for (map<int, int>::const_iterator siteIt = chrIt->second.begin(), siteIt_end = chrIt->second.end(); siteIt != siteIt_end; ++siteIt)
	//	{
	//		outputFile<<chrIt->first<<"\t"<<siteIt->first<<"\t"<<siteIt->second<<endl;
	//	}
	//}
	cout<<"Total intra domain interaction: "<<counterIntra<<endl;

	//outputFile.close();

	return counterIntra;
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



int findIntraDomainInteraction(string jobChr_fun, vector< vector<int> >& domainSitesMap_fun, vector< vector< vector<int> > >& domainCSinterFreq_fun, map<int, int>& csInterChromTotalMap_fun, map< int, map< pair<int, int>, double > >& pVal_fun)
{
	string file_debug =  jobChr_fun + "_debug";
	ofstream outDebug(file_debug.c_str());

	srand((unsigned)time(0));
	double preMaxRnd = 0.9/PRECISION;
	double minRnd = 0.1/PRECISION;
	//for (map<int, int>::const_iterator tempIt = csInterChromTotalMap_fun.begin(), tempIt_end = csInterChromTotalMap_fun.end(); tempIt != tempIt_end; ++tempIt)
	//{
	//	cout<<jobChr_fun<<"\t"<<tempIt->first<<"\t"<<tempIt->second<<endl;
	//}
	int contructR_argc = 1;
	char* constructR_argv[1];
	RInside R(contructR_argc, constructR_argv);// create an embedded R instance, without arguments product compilation error "parseEvalQ is of non-class type"
	//cout<<"R instance created"<<endl;
	string rcommand = "suppressMessages(library(glmnet))";
	R.parseEvalQ(rcommand);
	rcommand = "suppressMessages(library(lars))";
	R.parseEvalQ(rcommand);
	rcommand = "suppressMessages(library(glmpath))";
	R.parseEvalQ(rcommand);
	rcommand = "suppressMessages(source(\"~/program/HiC_mixPLD/covTest/R/funcs.R\"))";
	R.parseEvalQ(rcommand);
	//cout<<"Library loaded"<<endl;

	int domainNum = domainSitesMap_fun.size();

	//#pragma omp parallel for
	for (int domainIt = 0; domainIt < domainNum; ++domainIt)
	{
		map< pair<int, int>, double>& domainPval = pVal_fun[domainIt];
		outDebug<<"domain "<<domainIt<<endl;
		vector<int>& sitesMap = domainSitesMap_fun[domainIt];
		vector< vector<int> >& csInterFreqMap = domainCSinterFreq_fun[domainIt];
		int siteNum = sitesMap.size();
		//outDebug<<"domain "<<domainIt<<" siteNum "<<siteNum<<endl;

		int size_convert = convtIndex(siteNum-2,siteNum-1,siteNum) + 1;//the convtIndex function's first two input is 0 based and the output is 0 based too
		vector<double> freq_convert(size_convert, 0.0);
		//map< pair<int, int>, double > distMatrix;
		vector<double> dist_ij(size_convert, 0.0);
		//vector<double> tempDist(size_convert, 0);
		//tempDist[0]=1;
		//for ( int i = 0; i < size_convert; ++i )
		//{
		//	distMatrix.push_back(tempDist);
		//}
		vector<double> effVec(size_convert, 0.0);//efficiency of cutting site and mappability, measured by the total inter chromosome hybrid frags

		for (int rowIt=0; rowIt<siteNum; ++rowIt)
		{
			for (int colIt=rowIt+1; colIt<siteNum; ++colIt)
			{
				//if (colIt<=rowIt) {continue;}
				int index_convert=convtIndex(rowIt,colIt,siteNum);
				//outDebug<<"distMatrix_rowIt "<<distMatrix_rowIt<<endl;
				dist_ij[index_convert]=log(abs(sitesMap[rowIt]-sitesMap[colIt]));

				freq_convert[index_convert]=csInterFreqMap[rowIt][colIt];
				effVec[index_convert]=sqrt(csInterChromTotalMap_fun[sitesMap[rowIt]])*sqrt(csInterChromTotalMap_fun[sitesMap[colIt]]);
				//outDebug<<sitesMap[rowIt]<<"\t"<<sitesMap[colIt]<<"\t"<<effVec[index_convert]<<endl;
				//distMatrix[make_pair(distMatrix_rowIt,1)]=log(dist);
				//y2[n]=A2[i,j]
				//d[n]=abs(i-j)

				//calculate the Weight matrix
			}
		}
		//distMatrix: col1:1, col2:log(dist), col3....distMatrix
		vector<double> offsets(size_convert,0.0);
		//for (int i = 0; i < size_convert; ++i)
		//{
		//	distMatrix[make_pair(i, i)]=1;
		//	offsets.push_back(0.0);
		//}

		vector< vector<double> > covMatrix;
		vector<double> tempIntercept(size_convert, 1.0); 
		covMatrix.push_back(tempIntercept);
		covMatrix.push_back(dist_ij);
		covMatrix.push_back(effVec);

		IRLS irls("log-link");
		bool quasi_lik = false;
		irls.link->quasi = quasi_lik;
		//outDebug<<"freq_convert "<<freq_convert.size()<<" matrix raw "<<covMatrix.size()<<" matrix col "<<covMatrix[0].size()<<endl;
		irls.load_data(freq_convert, covMatrix, offsets);
		//outDebug<<"data loaded"<<endl;
		irls.fit_model();
		//outDebug<<"model fitted"<<endl;
		vector<double> coev = irls.get_coef();
		//outDebug<<"get coef done"<<endl;
		vector<double> sev = irls.get_stderr();
		//outDebug<<"get stderr done"<<endl;

		// print the results
		//outDebug<<"domain: "<<domainIt+1<<endl;
		//printf("dispersion=%.4f\n", irls.get_dispersion());
		//outDebug<<"dispersion "<<irls.get_dispersion()<<endl;
		//printf("%10s%12s%12s%15s\n", "", "Estimate", "Std.Error", "p-value");
		outDebug<<"Col\tEstimate\tStd.Error\tp-value"<<endl;
		for(size_t i = 0; i < coev.size(); ++i)
		{
			//printf("X%-9zu%12.9f%12.8f", i, coev[i], sev[i]);
			outDebug<<"X"<<i<<"\t"<<coev[i]<<"\t"<<sev[i]<<"\t";
			if(! irls.link->quasi)
			{
				//printf("%15.6e\n", 2 * gsl_cdf_gaussian_P(-fabs(coev[i]/sev[i]), 1.0));
				outDebug<<2 * gsl_cdf_gaussian_P(-fabs(coev[i]/sev[i]), 1.0)<<endl;
			}
			else
			{
				//printf("%15.6e\n", 2 * gsl_cdf_tdist_P(-fabs(coev[i]/sev[i]), size_convert-irls.get_rank_X()));
				outDebug<<2 * gsl_cdf_tdist_P(-fabs(coev[i]/sev[i]), size_convert-irls.get_rank_X())<<endl;
			}
		}
		covMatrix.clear();
		//outDebug.close();

		double mu = coev[0];
		double alpha = coev[1];
		double beta = coev[2];
		vector<double> constTerms(size_convert, 0.0);
		for (int i = 0; i < size_convert; ++i)
		{
			constTerms[i] = mu+dist_ij[i]*alpha+effVec[i]*beta;
		}

		int counter_fit=0;
		int counter_skip=0;

		//#pragma omp parallel for
		for (int rowIt=0; rowIt<siteNum; ++rowIt)
		{
			for (int colIt=rowIt+1; colIt<siteNum; ++colIt)
			{
				//if (colIt<=rowIt) {continue;}
				//int distMatrix_rowIt=convtIndex(rowIt,colIt,siteNum);
				if(csInterFreqMap[rowIt][colIt]>1)
				{
					vector<double> neighbFreq;
					vector<double> neighbConstTerms;
					vector<double> distVec;
					int sideLen = 2*NEIGHBDIS + 1;
					int colNum_mat = 0;
					int rowNum_mat = 0;

					for (int neighbRowIt=(rowIt-NEIGHBDIS), squareRowIt = 0; squareRowIt < sideLen; ++neighbRowIt, ++squareRowIt)
					{
						for (int neighbColIt=colIt-NEIGHBDIS, squareColIt = 0; squareColIt < sideLen; ++neighbColIt, ++squareColIt)
						{
							int distMatrix_colIt=convtIndex(neighbRowIt, neighbColIt, siteNum);
							//if (distMatrix_colIt==-1 || neighbColIt<=neighbRowIt || colIt<=neighbRowIt || neighbColIt<=rowIt || (neighbRowIt==rowIt && neighbColIt==colIt) )
							if (distMatrix_colIt==-1 || neighbColIt<=neighbRowIt || colIt<=neighbRowIt || neighbColIt<=rowIt)
							{continue;}
							neighbFreq.push_back(freq_convert[distMatrix_colIt]);
							neighbConstTerms.push_back(constTerms[distMatrix_colIt]);
							distVec.push_back(exp(-MORANI*(abs(sitesMap[rowIt]-sitesMap[neighbRowIt])+abs(sitesMap[colIt]-sitesMap[neighbColIt]))));
							++rowNum_mat;
							if (csInterFreqMap[neighbRowIt][neighbColIt] != 0)
							{
								++colNum_mat;
							}
						}
					}
					vector< vector<double> > distMatrix;
					distMatrix.push_back(distVec);

					IRLS irls("log-link");
					bool quasi_lik = false;
					irls.link->quasi = quasi_lik;
					irls.load_data(neighbFreq, distMatrix, neighbConstTerms);
					irls.fit_model();
					vector<double> coev_fitFreq = irls.get_coef();
					vector<double> sev_fitFreq = irls.get_stderr();

					double fitPval = 0.0;
					if(! irls.link->quasi)
					{
						fitPval = 2 * gsl_cdf_gaussian_P(-fabs(coev_fitFreq[0]/sev_fitFreq[0]), 1.0);
					}
					else
					{
						fitPval = 2 * gsl_cdf_tdist_P(-fabs(coev_fitFreq[0]/sev_fitFreq[0]), size_convert-irls.get_rank_X());
					}

					if (fitPval < PCUTOFF && colNum_mat > 1)
					{
						int distNum = colNum_mat;
						colNum_mat++;//add const to the last column of the dist matrix
						vector<double> distVec_lasso(colNum_mat, 0.0);
						vector< vector<double> > distMatrix_lasso;
						for (int i=0; i<rowNum_mat; ++i)
						{
							distVec_lasso[distNum] = neighbConstTerms[i];
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
								if (distMatrix_colIt==-1 || neighbColIt<=neighbRowIt || colIt<=neighbRowIt || neighbColIt<=rowIt || csInterFreqMap[neighbRowIt][neighbColIt] == 0)
								{continue;}

								int matRowIndex = 0;//inside loop is actually the row index for the final matrix
								for (int neighbNeighbRowIt=(rowIt-NEIGHBDIS), squareNeighbRowIt = 0; squareNeighbRowIt < sideLen; ++neighbNeighbRowIt, ++squareNeighbRowIt)
								{
									for (int neighbNeighbColIt=colIt-NEIGHBDIS, squareNeighbColIt = 0; squareNeighbColIt < sideLen; ++neighbNeighbColIt, ++squareNeighbColIt)
									{
										int distMatrix_rowIt=convtIndex(neighbNeighbRowIt, neighbNeighbColIt, siteNum);
										if (distMatrix_rowIt==-1 || neighbNeighbColIt<=neighbNeighbRowIt || colIt<=neighbNeighbRowIt || neighbNeighbColIt<=rowIt)
										{continue;}
										distMatrix_lasso[matRowIndex][matColIndex]=exp(-MORANI*(abs(sitesMap[neighbNeighbRowIt]-sitesMap[neighbRowIt])+abs(sitesMap[neighbNeighbColIt]-sitesMap[neighbColIt])));
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

/*Note Y need to substract offset before doing lasso fit, may use code from IRLS_glm*/

					//outDebug<<"neighbFreq"<<endl;
					//for (vector<double>::const_iterator vecIt = neighbFreq.begin(), vecIt_end = neighbFreq.end(); vecIt != vecIt_end; ++vecIt)
					//{
					//	outDebug<<*vecIt<<endl;
					//}
					//outDebug<<"neighbConstTerms"<<endl;
					//for (vector<double>::const_iterator vecIt = neighbConstTerms.begin(), vecIt_end = neighbConstTerms.end(); vecIt != vecIt_end; ++vecIt)
					//{
					//	outDebug<<*vecIt<<endl;
					//}
					//outDebug<<"distMatrix"<<endl;
					//for (vector< vector<double> >::const_iterator matIt = distMatrix_lasso.begin(), matIt_end = distMatrix_lasso.end(); matIt != matIt_end; ++matIt)
					//{
					//	vector<double>::const_iterator vecIt = (*matIt).begin(), vecIt_end = (*matIt).end();
					//	outDebug<<*vecIt;
					//	++vecIt;
					//	for ( ; vecIt != vecIt_end; ++vecIt)
					//	{
					//		outDebug<<"\t"<<*vecIt;
					//	}
					//	outDebug<<endl;
					//}
						#pragma omp critical
						{
						outDebug<<rowIt<<"\t"<<colIt<<endl;
						Rcpp::NumericVector freq_r(neighbFreq.begin(), neighbFreq.end());
						Rcpp::NumericMatrix distMat_r(rowNum_mat, colNum_mat);
						for (int i=0; i<rowNum_mat; ++i)
						{
							for (int j=0; j<colNum_mat; ++j)
							{
								double rndNum = preMaxRnd*rand()/(RAND_MAX + 1.0) + minRnd;
								//distMat_r(i,j)= floor((distMatrix_lasso[i][j]+rndNum) * PRECISION + 0.5) / PRECISION;
								//distMat_r(i,j)= distMatrix_lasso[i][j];
								distMat_r(i,j)= distMatrix_lasso[i][j] + rndNum;
							}
						}
						ostringstream rowNum_str;
						rowNum_str<<rowNum_mat;

						R["freq_R"] = freq_r;
						rcommand = "write.table(freq_R, \"" + jobChr_fun + "_freq\", col.names = FALSE, row.names = FALSE, sep = \"\\t\", quote = FALSE)";
						R.parseEvalQ(rcommand);
						R["distMat_R"] = distMat_r;
						rcommand = "write.table(distMat_R, \"" + jobChr_fun + "_distMat\", col.names = FALSE, row.names = FALSE, sep = \"\\t\", quote = FALSE)";
						R.parseEvalQ(rcommand);
						//cout<<"matrix and vector are converted for R."<<endl;

						//rcommand = "round(distMat_R, 8)";
						rcommand = "suppressMessages(lars.glm(distMat_R, freq_R, family=\"poisson\", status=NULL, standardize=TRUE))";
						Rcpp::List larsGlm_fitObj = R.parseEval(rcommand);
						outDebug<<"lars.glm done."<<endl;
						R["larsGlm_fitObj_R"] = larsGlm_fitObj;
						rcommand = "suppressMessages(covTest(larsGlm_fitObj_R, distMat_R, freq_R, sigma.est=\"full\",status=NULL, maxp=" + rowNum_str.str() + "))";
						Rcpp::List lasso_p = R.parseEval(rcommand);
						rowNum_str.clear();
						//show(glmnetFit);

						//cout<<"glmnet fit done."<<endl;

						//glm_lasso_fit(neighbFreq, neighbConstTerms, distMatrix);

						//IRLS irls("log-link");
						//bool quasi_lik = false;
						//irls.link->quasi = quasi_lik;
						//irls.load_data(neighbFreq, distMatrix, neighbConstTerms);
						//irls.fit_model();
						//vector<double> coev_fitFreq = irls.get_coef();
						//vector<double> sev_fitFreq = irls.get_stderr();

						//double fitPval = 0.0;
						//if(! irls.link->quasi)
						//{
						//	fitPval = 2 * gsl_cdf_gaussian_P(-fabs(coev_fitFreq[0]/sev_fitFreq[0]), 1.0);
						//}
						//else
						//{
						//	fitPval = 2 * gsl_cdf_tdist_P(-fabs(coev_fitFreq[0]/sev_fitFreq[0]), size_convert-irls.get_rank_X());
						//}
						//distMatrix.clear();
						//neighbFreq.clear();
						//neighbConstTerms.clear();
						//distVec.clear();

						//#pragma omp critical
						//domainPval[make_pair(sitesMap[rowIt], sitesMap[colIt])]=fitPval;
						++counter_fit;}
						//#pragma omp end critical

						//if (fitPval < PCUTOFF)
						//{
						//	
						//}
					}
					else if(colNum_mat==1)
					{
						//only 1 box has signal (freq>1), no need to use lasso to choose which one is real
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
		freq_convert.clear();
		//distMatrix.clear();
		dist_ij.clear();
		offsets.clear();
		effVec.clear();
		//outDebug<<"chr "<<jobChr_fun<<" domain "<<domainIt<<" siteNum "<<siteNum<<" convert size "<<size_convert<<"fit number "<<counter_fit<<" skip number "<<counter_skip<<endl;
	}
	outDebug.close();
	
	return 0;
}




int outputDomainInfo(string jobChr_fun, vector< vector<int> >& domainSitesMap_fun, vector< vector< vector<int> > >& domainCSinterFreq_fun, map<int, int>& csInterChromTotalMap_fun, map< int, map< pair<int, int>, double > >& pVal_fun)
{
	string file_dsm =  jobChr_fun + "_domainSitesMap";
	string file_dif =  jobChr_fun + "_domainCSinterFreq";
	string file_ictm =  jobChr_fun + "_csInterChromTotalMap";
	string file_p =  jobChr_fun + "_pVal";
	ofstream outDsm(file_dsm.c_str());
	ofstream outDif(file_dif.c_str());
	ofstream outIctm(file_ictm.c_str());
	//ofstream outP(file_p.c_str());

	for (vector< vector<int> >::const_iterator domnIt = domainSitesMap_fun.begin(), domnIt_end = domainSitesMap_fun.end(); domnIt != domnIt_end; ++domnIt)
	{
		vector<int>::const_iterator siteIt = domnIt->begin(), siteIt_end = domnIt->end();
		outDsm<<*siteIt;
		++siteIt;
		for ( ; siteIt != siteIt_end; ++siteIt)
		{
			outDsm<<"\t"<<*siteIt;
		}
		outDsm<<endl;
	}
	outDsm.close();

	for (vector< vector< vector<int> > >::const_iterator domnIt = domainCSinterFreq_fun.begin(), domnIt_end = domainCSinterFreq_fun.end(); domnIt != domnIt_end; ++domnIt)
	{
		outDif<<"domain "<<domnIt-domainCSinterFreq_fun.begin()<<endl;
		for (vector< vector<int> >::const_iterator rowIt = domnIt->begin(), rowIt_end = domnIt->end(); rowIt != rowIt_end; ++rowIt)
		{
			vector<int>::const_iterator colIt = rowIt->begin(), colIt_end = rowIt->end();
			outDif<<*colIt;
			++colIt;
			for ( ; colIt != colIt_end; ++colIt)
			{
				outDif<<"\t"<<*colIt;
			}
			outDif<<endl;
		}
	}
	outDif.close();

	for (map<int, int>::const_iterator csIt = csInterChromTotalMap_fun.begin(), csIt_end = csInterChromTotalMap_fun.end(); csIt != csIt_end; ++csIt)
	{
		outIctm<<csIt->first<<"\t"<<csIt->second<<endl;
	}
	outIctm.close();

	//for (map< int, map< pair<int, int>, double > >::const_iterator domnIt = pVal_fun.begin(), domnIt_end = pVal_fun.end(); domnIt != domnIt_end; ++domnIt)
	//{
	//	for (map< pair<int, int>, double >::const_iterator csIt = domnIt->second.begin(), csIt_end = domnIt->second.end(); csIt != csIt_end; ++csIt)
	//	{
	//		outP<<domnIt->first<<"\t"<<csIt->first.first<<"\t"<<csIt->first.second<<"\t"<<csIt->second<<endl;
	//	}
	//}
	//outP.close();

	return 0;
}




int findInteractingRegion(int jobChr_fun, map< int, map< pair<int, int>, double> >& sigInteractionMap_fun, vector< vector<int> >& domainSitesMap_fun, vector<regionDetail_pairSite_2D>& regionMaptoWrite_2D)
{
	//cout<<"freq map duplicated"<<endl;
	regionMaptoWrite_2D.clear();

	int domainNum = domainSitesMap_fun.size();

	int countRegionnum = 0;

	#pragma omp parallel for
	for (int domnIt = 0; domnIt < domainNum; ++domnIt)
	{
		map< pair<int, int>, double > domainSigInteraction = sigInteractionMap_fun[domnIt];
		map<int, int> domainCS;
		for (vector<int>::const_iterator siteIt = domainSitesMap_fun[domnIt].begin(), siteIt_end = domainSitesMap_fun[domnIt].end(); siteIt != siteIt_end; ++siteIt)
		{
			domainCS[*siteIt] = 0;
		}

		map< pair<int, int>, double >::iterator csPairIt = domainSigInteraction.begin(), csPairIt_end = domainSigInteraction.end();
		//recurrsively search for neighbour region that has a higher score than the threshold
		while(csPairIt != csPairIt_end)
		{
			int pos_cs1 = csPairIt->first.first;
			int pos_cs2 = csPairIt->first.second;

			//int csPval = csPairIt->second;
			double minPval = csPairIt->second;

			int csPairNum = 1;//my $binNo = 0;

			vector< pair< pair<int, int>, double > > regionsCoordinates;
			regionDetail_pairSite_2D regionInfo;
			//vector<csInfo_bothEndsMapped_withCuttingSite> regionFragInfo;
			//regionFragInfo = csPairIt->second;

			regionsCoordinates.push_back(*csPairIt);
			domainSigInteraction.erase(csPairIt);

			//if (csPairScore >= thres)
			//{
				searchNeighbSites(domainSigInteraction, domainCS, pos_cs1, pos_cs2, 0, 0, regionsCoordinates, csPairNum, minPval);
			//}

			//if (totalScore >= doubleThres)//if regionIt itself or the sum of all scores in the region is larger than double threshold -> region
			//{
				regionInfo.chrNo1 = jobChr_fun;
				regionInfo.chrNo2 = jobChr_fun;
				regionInfo.domain1 = domnIt;
				regionInfo.domain2 = domnIt;
				regionInfo.sitePair_num = csPairNum;
				regionInfo.minPval = minPval;
				//regionInfo.cs_info = regionFragInfo;
				regionInfo.region_info = regionsCoordinates;
				//regionInfo.getPeakSummit();

				#pragma omp critical
				{
				regionMaptoWrite_2D.push_back(regionInfo);
				++countRegionnum;
				}
				//#pragma omp end critical
			//}

			csPairIt = domainSigInteraction.begin();
		}
	}

	return countRegionnum;
}




int searchNeighbSites(map< pair<int, int>, double >& csInteraction_search, map<int, int>& csMap_search, int pos_domain1_local, int pos_domain2_local, int parentSetoff_x, int parentSetoff_y, vector< pair< pair<int, int>, double > >& regionsCoordinates_local, int& pairSiteNum_local, double& minPval_local)
{
	//cout <<chr_domain1_local<<"\t"<<chr_domain2_local<<"\t"<<pos_domain1_local.first<<"\t"<<pos_domain2_local.first<<"\t"<<parentSetoff_x<<"\t"<<parentSetoff_y<<"\t"<<thres_local<<"\t"<<pairSiteNum_local<<"\t"<<totalScore_local<<endl;
	for (int x=-1; x<=1; ++x)
	{
		//cout<<"enter x."<<endl;
		map<int, int>::const_iterator targetFrag_x = csMap_search.find(pos_domain1_local);

		if ( x==-1 && targetFrag_x != csMap_search.begin() && targetFrag_x != csMap_search.end() ) //search for the previous one
		{
			//cout<<x<<endl;
			targetFrag_x--;
		}
		else if ( x==1 && targetFrag_x != csMap_search.end() )
		{
			//cout<<x<<endl;
			targetFrag_x++;
			if (targetFrag_x == csMap_search.end())
			{
				//cout<<"if continue"<<endl;
				continue;
			}
		}
		else if ( x==0 && targetFrag_x != csMap_search.end() )
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
			map<int, int>::const_iterator targetFrag_y = csMap_search.find(pos_domain2_local);

			if ( y==-1 && targetFrag_y != csMap_search.begin() && targetFrag_y != csMap_search.end() ) //search for the previous one
			{
				targetFrag_y--;
			}
			else if ( y==1 && targetFrag_y != csMap_search.end() )
			{
				targetFrag_y++;
				if (targetFrag_y == csMap_search.end())
				{
					continue;
				}
			}
			else if ( y==0 && targetFrag_y != csMap_search.end() ) {}
			else
			{
				continue;
			}

			if ( !( (x==0 && y==0) || ( x==(0-parentSetoff_x) && y==(0-parentSetoff_y) ) ) )
			{
				//cout<<"enter xy."<<endl;
				map< pair<int, int>, double >::iterator neighbPairSites = csInteraction_search.find(make_pair(targetFrag_x->first, targetFrag_y->first));

				if (neighbPairSites != csInteraction_search.end())
				{
					//cout<<"enter recurrsive."<<endl;
					double neighbPval = neighbPairSites->second;
					int neighbPos_domain1 = neighbPairSites->first.first, neighbPos_domain2 = neighbPairSites->first.second;
					//vector<csInfo_bothEndsMapped_withCuttingSite> neighbPairSitesInfo = neighbPairSites->second;

					if (neighbPval < minPval_local)
					{
						minPval_local = neighbPval;
					}
					++pairSiteNum_local;
					regionsCoordinates_local.push_back(*neighbPairSites);
					//vector<csInfo_bothEndsMapped_withCuttingSite>::iterator csIt = neighbPairSitesInfo.begin(), csIt_end = neighbPairSitesInfo.end();
					//for ( ; csIt != csIt_end; ++csIt)
					//{
					//	regionFragInfo_local.push_back(*csIt);
					//}

					csInteraction_search.erase(neighbPairSites);//delete the element before pass down to recurrsive checking, if not, for example, in 3*3 array, [0,0][0,1][1,1],has value > thres, [0,0]->[0,1]->[1,1]->[0,0]->.....
					searchNeighbSites(csInteraction_search, csMap_search, neighbPos_domain1, neighbPos_domain2, x, y, regionsCoordinates_local, pairSiteNum_local, minPval_local);
					//}
				}
			}
		}
	}
	return 0;
}



int writer_2Dregion(vector<regionDetail_pairSite_2D>& regionMaptoRead_2D, const string& regionFiletoWrite_2D)
{
	ofstream output2DregionFile(regionFiletoWrite_2D.c_str());
	vector<regionDetail_pairSite_2D>::const_iterator regionIt = regionMaptoRead_2D.begin(), regionIt_end = regionMaptoRead_2D.end();
	for ( ; regionIt != regionIt_end; ++regionIt )
	{
		output2DregionFile <<regionIt->chrNo1<<"\t"<<regionIt->chrNo2<<"\t"<<regionIt->domain1<<"\t"<<regionIt->domain2<<"\t"<<regionIt->sitePair_num<<"\t"<<(regionIt->regionCoor.first.first)<<"\t"<<(regionIt->regionCoor.first.second)<<"\t"<<regionIt->regionCoor.second.first<<"\t"<<regionIt->regionCoor.second.second<<"\t"<<regionIt->summitPairSites.first.first<<"\t"<<regionIt->summitPairSites.first.second<<"\t"<<regionIt->minPval;
		vector< pair< pair<int, int>, double > >::const_iterator siteIt = regionIt->region_info.begin(), siteIt_end = regionIt->region_info.end();
		for ( ; siteIt != siteIt_end; ++siteIt )
		{
			output2DregionFile <<"\t"<<(siteIt->first.first)<<"\t"<<(siteIt->first.second)<<"\t"<<siteIt->second;
		}
		output2DregionFile<<endl;

		//vector<fragInfo_bothEndsMapped_withCuttingSite>::const_iterator fragIt = regionIt->frag_info.begin(), fragIt_end = regionIt->frag_info.end();
		//for ( ; fragIt != fragIt_end; ++fragIt)
		//{
		//	output2DregionFile <<fragIt->end1_chr<<"\t"<<fragIt->end1_pos<<"\t"<<fragIt->end1_strand<<"\t"<<fragIt->end1_cuttingSite<<"\t"<<fragIt->end2_chr<<"\t"<<fragIt->end2_pos<<"\t"<<fragIt->end2_strand<<"\t"<<fragIt->end2_cuttingSite<<"\t"<<fragIt->fragType<<endl;
		//}
	}
	output2DregionFile.close();

	return 0;
}



int glm_lasso_fit(vector<double>& freq_fun, vector<double>& consts_fun, vector< vector<double> >& distMat_fun)
{
	try
	{
		//cout<<"entered try"<<endl;
		int rowNum = distMat_fun.size();
		int colNum = (distMat_fun.begin())->size();
		int contructR_argc = 1;
		char* constructR_argv[1];
		RInside R(contructR_argc, constructR_argv);// create an embedded R instance, without arguments product compilation error "parseEvalQ is of non-class type"
		//cout<<"R instance created"<<endl;
		//string rcommand = "suppressMessages(library(Matrix))";
		//R.parseEvalQ(rcommand);// load library, no return value
		//rcommand = "suppressMessages(library(lattice))";
		//R.parseEvalQ(rcommand);// load library, no return value
		string rcommand = "suppressMessages(library(glmnet))";
		R.parseEvalQ(rcommand);// load library, no return value
		//cout<<"Library loaded"<<endl;

		Rcpp::NumericVector freq_r(freq_fun.begin(), freq_fun.end());
		R["freq_R"] = freq_r;
		Rcpp::NumericMatrix distMat_r(rowNum, colNum);
		for (int i=0; i<rowNum; ++i)
		{
			for (int j=0; j<colNum; ++j)
			{
				distMat_r(i,j)=distMat_fun[i][j];
			}
		}
		R["distMat_R"] = distMat_r;

		cout<<"matrix and vector are converted for R."<<endl;

		//rcommand = "glmnet.fit=glmnet(distMat_R, freq_R, family=\"poisson\")";
		rcommand = "glmnet(distMat_R, freq_R, family=\"poisson\")";
		Rcpp::List glmnetFit = R.parseEval(rcommand);
		//show(glmnetFit);

		cout<<"glmnet fit done."<<endl;

		//ostringstream FDRbuffer;
		//FDRbuffer<<FDR_local;
		//rcommand = "suppressMessages(max(qObj$pvalues[qObj$qvalues <= " + FDRbuffer.str() + "]))";
		//FDRbuffer.clear();
		//pvalCutoff_local = Rcpp::as< double >(R.parseEval(rcommand));

	} catch(std::exception& ex) {
		std::cerr << "RInside exception caught: " << ex.what() << std::endl;
	} catch(...) {
		std::cerr << "RInside unknown exception caught" << std::endl;
	}
	cout << "All done, past catch()"<<endl;

	return 0;
}

#endif
