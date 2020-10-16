#ifndef BOTHENDSMAPPEDFRAGINFO
#define BOTHENDSMAPPEDFRAGINFO

class fragInfo_bothEndsMapped_withCuttingSite
{
public:
	int end1_chr;
	int end2_chr;
	int end1_pos;
	int end1_cuttingSite;
	int end2_pos;
	int end2_cuttingSite;
	bool end1_strand;
	bool end2_strand;

	int fragType;

	fragInfo_bothEndsMapped_withCuttingSite():end1_chr(0),end1_pos(0),end1_strand(false),end1_cuttingSite(0),end2_chr(0),end2_pos(0),end2_strand(false),end2_cuttingSite(0),fragType(0){}
};

#endif
