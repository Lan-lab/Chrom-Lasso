#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <omp.h>
//#include <mpi.h>

//#include <boost/filesystem.hpp>

#include "HiC_mixturePLD_boss.h"
#include "HiC_mixturePLD_postdoc.h"

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


//namespace bfs = boost::filesystem;
//
//using bfs::exists;
//using bfs::directory_iterator;
//using bfs::is_directory;
//using bfs::is_regular;
//using bfs::path;
//using bfs::system_complete;


//argument should contain nodes number and how many cores per node
//chromosome X and Y, are represented as 23, 24 in original script, however in Ren's data is X Y
//MPI_Send(void* buf, int count, MPI_Datatype datatype,int dest, int tag, MPI_Comm comm)
//MPI_Recv(void* buf, int count, MPI_Datatype datatype,int source, int tag, MPI_Comm comm,MPI_Status *status)

//get detailed interacting region info from hybrid fragments, example, HiC_PLD -w 76(read size) 500(dis to cutting site threshold) 100000(model bin size) 26(number of nodes), 8(cores per nodes)

//string getPrefix(const string&);
void usage();

inline void usage()
{
	cout <<"Usage: HiCanalysis [OPTIONs] [PARAMETERS] INPUT ELAND OR BED FILE"<<endl;
}


int main(int argc, char* argv[])
{
	//parse arguments
	bool flag_diffSexChr = false;
	string dir_genomeFasta = "";
	string file_cuttingSites = "";
	string file_domain = "";
	string outputDir = "";

	double fdr = 0.01;
	map<int, int> chrlenMap;
	vector<int> intArgVec;
	vector<double> dblArgVec;
	vector<string> fileNames_forCount;

	if (argc < 2)
	{
		usage();
		exit(1);
	}

	int optionArg = 0;//number of option arguments (arguments that not sourcefile name)
	optionArg = argc - 2;

	int intArg = 0;
	double dblArg = 0.0;

	//options -w(bin size)[20,1000], -p(threshold percentage)[0.5,1.0), -s(input file is sorted), -f(calculate FDR);
	int argIt = 1;
	while ( argIt <= optionArg )
	{
		if (argv[argIt][0] == '-')
		{
			switch (argv[argIt][1])
			{
				case 'g':
					if ( strlen(argv[argIt]) != 2)
					{
						cout <<"Error: -g option can not be combined with other options."<<endl;
						exit(1);
					}

					++argIt;
					if ( argv[argIt][0] == '-' || argIt > optionArg )
					{
						cout <<"Specify the directory contains genome sequence information following -g option."<<endl;
						exit(1);
					}
					else
					{
						dir_genomeFasta = argv[argIt];
						++argIt;
						if ( argIt <= optionArg && argv[argIt][0] != '-' )
						{
							cout <<"Only one genome could be choosen for one data."<<endl;
							exit(1);
						}
					}
					break;

				case 'w':
					if ( strlen(argv[argIt]) != 2)
					{
						cout <<"Error: -w option can not be combined with other options."<<endl;
						exit(1);
					}

					++argIt;
					if ( argv[argIt][0] == '-' || argIt > optionArg )
					{
						cout <<"Error: no bin size specified following -w option."<<endl;
						exit(1);
					}
					while ( argv[argIt][0] != '-' && argIt <= optionArg )
					{
						intArgVec.push_back(atoi(argv[argIt]));
						++argIt;
					}
					break;

				case 'f'://HiCanalysis -f 0.01 -w 50..... FDR
					if ( strlen(argv[argIt]) != 2)
					{
						cout <<"Error: -f option can not be combined with other options."<<endl;
						exit(1);
					}

					++argIt;
					if ( argv[argIt][0] == '-' || argIt > optionArg )
					{
						cout <<"Error: no FDR specified following -f option."<<endl;
						exit(1);
					}
					else
					{
						fdr = atof(argv[argIt]);
						++argIt;
						if ( argIt <= optionArg && argv[argIt][0] != '-' )
						{
							cout <<"Only one fdr could be choosen for the process."<<endl;
							exit(1);
						}
					}
					break;

				case 'T'://HiCanalysis -T 4.5 6.4 -w 50.....
					if ( strlen(argv[argIt]) != 2)
					{
						cout <<"Error: -T option can not be combined with other options."<<endl;
						exit(1);
					}

					++argIt;
					if ( argv[argIt][0] == '-' || argIt > optionArg )
					{
						cout <<"Error: no threshold specified following -T option."<<endl;
						exit(1);
					}
					while ( argv[argIt][0] != '-' && argIt <= optionArg )
					{
						dblArgVec.push_back(atof(argv[argIt]));
						++argIt;
					}
					break;

				case 'q'://use -a option at the beginning to specify the raw data file names that need to be counted
					if ( strlen(argv[argIt]) != 2)
					{
						cout <<"Error: -a option can not be combined with other options."<<endl;
						exit(1);
					}

					++argIt;
					if ( argv[argIt][0] == '-' || argIt > optionArg )
					{
						cout <<"Error: no file specified following -a option."<<endl;
						exit(1);
					}
					while ( argv[argIt][0] != '-' && argIt <= optionArg )
					{
						fileNames_forCount.push_back(argv[argIt]);
						++argIt;
					}
					break;

				case 'd':
					if ( strlen(argv[argIt]) != 2)
					{
						cout <<"Error: -d option can not be combined with other options."<<endl;
						exit(1);
					}

					++argIt;
					if ( argv[argIt][0] == '-' || argIt > optionArg )
					{
						cout <<"Specify the name of the domain file following -d option."<<endl;
						exit(1);
					}
					else
					{
						file_domain = argv[argIt];
						++argIt;
						if ( argIt <= optionArg && argv[argIt][0] != '-' )
						{
							cout <<"Only one domain file could be used for one data."<<endl;
							exit(1);
						}
					}
					break;


				case 'c':
					if ( strlen(argv[argIt]) != 2)
					{
						cout <<"Error: -c option can not be combined with other options."<<endl;
						exit(1);
					}

					++argIt;
					if ( argv[argIt][0] == '-' || argIt > optionArg )
					{
						cout <<"Specify the name of the enzyme cutting site map file following -c option."<<endl;
						exit(1);
					}
					else
					{
						file_cuttingSites = argv[argIt];
						++argIt;
						if ( argIt <= optionArg && argv[argIt][0] != '-' )
						{
							cout <<"Only one enzyme cutting site map file could be used for one data."<<endl;
							exit(1);
						}
					}
					break;

				case 'h':
					if ( strlen(argv[argIt]) != 2)
					{
						cout <<"Error: -h option can not be combined with other options."<<endl;
						exit(1);
					}

					++argIt;
					if ( argv[argIt][0] == '-' || argIt > optionArg )
					{
						cout <<"Specify the output directory following -h option, default output directory is the current working directory"<<endl;
						exit(1);
					}
					else
					{
						outputDir = argv[argIt];
						++argIt;
						if ( argIt <= optionArg && argv[argIt][0] != '-' )
						{
							cout <<"Only one output directory could be used."<<endl;
							exit(1);
						}
					}
					break;

				default:
					for (int inargIt = 1; inargIt != strlen(argv[argIt]); ++inargIt)
					{
						switch (argv[argIt][inargIt])
						{
							case 'y':
								flag_diffSexChr = true;
								break;
							default:
								cout <<"Error: Undefined option: -"<<argv[argIt][inargIt]<<endl;
								usage();
								exit(1);
						}
					}
					++argIt;
			}
		}
		else
		{
			cout <<"Error: Undefined argument: "<<argv[argIt]<<endl;
			usage();
			exit(1);
		}
	}

	string file_source = argv[argc-1];

	if (dir_genomeFasta.size() == 0)
	{
		dir_genomeFasta = "hg18";
		//cout <<"genomeSequencedSize: "<<genomeSequencedSize<<endl;
	}

	//for easy converting chrX Y to number, all X will be 23 and all Y will be 24, unless there are more than 22 pairs of autochromosome
	if (dir_genomeFasta == "hg18")
	{
		chrlenMap[1] = 247249719;
		chrlenMap[2] = 242951149;
		chrlenMap[3] = 199501827;
		chrlenMap[4] = 191273063;
		chrlenMap[5] = 180857866;
		chrlenMap[6] = 170899992;
		chrlenMap[7] = 158821424;
		chrlenMap[8] = 146274826;
		chrlenMap[9] = 140273252;
		chrlenMap[10] = 135374737;
		chrlenMap[11] = 134452384;
		chrlenMap[12] = 132349534;
		chrlenMap[13] = 114142980;
		chrlenMap[14] = 106368585;
		chrlenMap[15] = 100338915;
		chrlenMap[16] = 88827254;
		chrlenMap[17] = 78774742;
		chrlenMap[18] = 76117153;
		chrlenMap[19] = 63811651;
		chrlenMap[20] = 62435964;
		chrlenMap[21] = 46944323;
		chrlenMap[22] = 49691432;
		chrlenMap[23] = 154913754;
		chrlenMap[24] = 57772954;
		chrlenMap[25] = 16571;
	}
	else if (dir_genomeFasta == "hg19")
	{
		chrlenMap[1] = 249250621;
		chrlenMap[2] = 243199373;
		chrlenMap[3] = 198022430;
		chrlenMap[4] = 191154276;
		chrlenMap[5] = 180915260;
		chrlenMap[6] = 171115067;
		chrlenMap[7] = 159138663;
		chrlenMap[8] = 146364022;
		chrlenMap[9] = 141213431;
		chrlenMap[10] = 135534747;
		chrlenMap[11] = 135006516;
		chrlenMap[12] = 133851895;
		chrlenMap[13] = 115169878;
		chrlenMap[14] = 107349540;
		chrlenMap[15] = 102531392;
		chrlenMap[16] = 90354753;
		chrlenMap[17] = 81195210;
		chrlenMap[18] = 78077248;
		chrlenMap[19] = 59128983;
		chrlenMap[20] = 63025520;
		chrlenMap[21] = 48129895;
		chrlenMap[22] = 51304566;
		chrlenMap[23] = 155270560;
		chrlenMap[24] = 59373566;
	}
	else if (dir_genomeFasta == "mm8")
	{
		chrlenMap[1] =  197069962;
		chrlenMap[2] =  181976762;
		chrlenMap[3] =  159872112;
		chrlenMap[4] =  155029701;
		chrlenMap[5] =  152003063;
		chrlenMap[6] =  149525685;
		chrlenMap[7] =  145134094;
		chrlenMap[8] =  132085098;
		chrlenMap[9] =  124000669;
		chrlenMap[10] =  129959148;
		chrlenMap[11] =  121798632;
		chrlenMap[12] =  120463159;
		chrlenMap[13] =  120614378;
		chrlenMap[14] =  123978870;
		chrlenMap[15] =  103492577;
		chrlenMap[16] =  98252459;
		chrlenMap[17] =  95177420;
		chrlenMap[18] =  90736837;
		chrlenMap[19] =  61321190;
		chrlenMap[23] =  165556469;
		chrlenMap[24] =  16029404;
	}
	else if (dir_genomeFasta == "mm9")
	{
		chrlenMap[1] = 197195432;
		chrlenMap[2] = 181748087;
		chrlenMap[3] = 159599783;
		chrlenMap[4] = 155630120;
		chrlenMap[5] = 152537259;
		chrlenMap[6] = 149517037;
		chrlenMap[7] = 152524553;
		chrlenMap[8] = 131738871;
		chrlenMap[9] = 124076172;
		chrlenMap[10] = 129993255;
		chrlenMap[11] = 121843856;
		chrlenMap[12] = 121257530;
		chrlenMap[13] = 120284312;
		chrlenMap[14] = 125194864;
		chrlenMap[15] = 103494974;
		chrlenMap[16] = 98319150;
		chrlenMap[17] = 95272651;
		chrlenMap[18] = 90772031;
		chrlenMap[19] = 61342430;
		chrlenMap[23] = 166650296;
		chrlenMap[24] = 15902555;
	}
	else if (dir_genomeFasta == "mm10")
	{
		chrlenMap[1] = 195471971;
		chrlenMap[2] = 182113224;
		chrlenMap[3] = 160039680;
		chrlenMap[4] = 156508116;
		chrlenMap[5] = 151834684;
		chrlenMap[6] = 149736546;
		chrlenMap[7] = 145441459;
		chrlenMap[8] = 129401213;
		chrlenMap[9] = 124595110;
		chrlenMap[10] = 130694993;
		chrlenMap[11] = 122082543;
		chrlenMap[12] = 120129022;
		chrlenMap[13] = 120421639;
		chrlenMap[14] = 124902244;
		chrlenMap[15] = 104043685;
		chrlenMap[16] = 98207768;
		chrlenMap[17] = 94987271;
		chrlenMap[18] = 90702639;
		chrlenMap[19] = 61431566;
		chrlenMap[23] = 171031299;
		chrlenMap[24] = 91744699;
	}

	//start process
	int readSize = intArgVec[0];
	int halfReadSize = readSize/2;
	int cuttingSiteExtent = intArgVec[1];//1000?
	int modelBinSize = intArgVec[2];
	int nodesNum = intArgVec[3];
	int procsPerNode = intArgVec[4];
	int bossId = 0;

	boss(file_cuttingSites, file_domain, cuttingSiteExtent, file_source, halfReadSize, chrlenMap, fdr);

	return 0;
}


