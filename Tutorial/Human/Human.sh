#####Prepare Input File#####
Chrom-Lasso needs 2 Input files prepared according to the Hi-C experimental design,
One is cutting site file, which presents the cutting sites of restriction enzyme used
in the Hi-C experiment on genome, and the other one is domain file, which describes the 
compartmentalization of TAD on genome. 
For the cutting site file, we directly take advantage of Oligomatch, an open-source tool to find cutting sites in all genome for a specific restriction endonuclease,
to genereta this file(the details are in /Prepare_Input_File/Cutting_Site_File/Prepare_CuttingSite_File_By_Oligomatch), 
and for domain file, we upload domain files for mouse mm10 genome, and for human hg19 genome,
and researchers can transform the genome version via LiftOver.
More details can be found in folder "Prepare_Input_File"
################################################################################################
#####Preprocess#####
The preprocessing of Hi-C data mainly contains two parts, mappping and filtering reads.
Here we directly use the preprocess protocol of JUICER to handle the raw sequencing file of Hi-C library.
In the preprocess part, the raw "fastq" file will be transformed into "sortChr" file.

###1. Mapping and sorting.
refSeq=/"The path to your reference genome for mapping"/
name1=/"The path to your paired-end fastq file 1"/
name2=/"The path to your paired-end fastq file 2"/
bwa mem -t 16 $refSeq $name1 > Human_1.sam
bwa mem -t 16 $refSeq $name2 > Human_2.sam 
sort -k1,1f Human_1.sam > Human_1_sort.sam
sort -k1,1f Human_2.sam > Human_2_sort.sam
awk 'BEGIN{OFS="\t"}NF >= 11{$1 = $1"/1";print}' Human_1_sort.sam > Human_1_sort1.sam
awk 'BEGIN{OFS="\t"}NF >= 11{$1 = $1"/2";print}' Human_2_sort.sam > Human_2_sort1.sam
sort -k1,1f -m  Human_1_sort1.sam Human_2_sort1.sam > Human.sam

###2. Filtering reads.
awk -v "fname1"=Human_norm.txt -v "fname2"=Human_abnorm.sam -v "fname3"=Human_unmapped.sam -f /Code/1_Preprocess/2_Filter_Reads/chimeric_blacklist.awk  Human.sam
awk '{printf("%s %s %s %d %s %s %s %d", $1, $2, $3, 0, $4, $5, $6, 1); for (i=7; i<=NF; i++) {printf(" %s",$i);}printf("\n");}' Human_norm.txt > Human.frag.txt
sort -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n Human.frag.txt > Human.sort.txt
awk -v "name"=Human  -f /Code/1_Preprocess/2_Filter_Reads/dups.awk  Human.sort.txt
awk '{if($1==0){strand1=1;}else if($1==16) {strand1=0;} if($5==0){strand2=1;}else if($5==16) {strand2=0;} print $2"\t"$3"\t"strand1"\t"$6"\t"$7"\t"strand2"\t0";}' Human_merged_nodups.txt > Human.formatted
for chr in {1..22} X;do awk '{if($1=="'$chr'") print $0;}' Human.formatted;done > Human.sortChr
################################################################################################
#####Identify Chromatin Interactions by Chrom-Lasso#####
After the preprocessing of Hi-C data, it comes to the interaction calling by Chrom-Lasso, 
In this file we illustrate the analysis of mouse Hi-C data step by step using Chrom-Lasso,
starting from the "sortChr" file derived from preprocessing step. 

###1. Rearrange reads according to cutting site file and domain file.###
Parameter: 
-g: set the genome version
-d: path to domain file
-c: path to cutting site file

Order:
/Code/2_Arrange_Domain/HiC_mixturePLD_singleThread -g hg19 -w 51 500 10000 26 1 -d /Prepare_Input_File/Domain_File/total.domain.hg19  -c /Prepare_Input_File/Cutting_Site_File/MboI.hg19.bed Human.sortChr > "Human.sortChr_summary"

Results:
For each chromosome, it generates 3 files, "chr_csInterChromTotalMap", "chr_domainCSinterFreq", "chr_domainSitesMap".

###2. Model the background distribution via polynomial regression for Y(ligation frequency) and X(genomic distance).###
for chr in chr{1..22} chrX;do awk 'function abs(x){return ((x < 0.0) ? -x : x)}BEGIN{i=0;}{if($2!=0){map[i]=$1;++i;}}END{for(j=0;j<i;++j){for(k=j+1;k<i;++k){dis=abs(map[j]-map[k]);if(dis<100000){dist[int(dis/100)]++;}else{break;}}}for(i=0;i<1000;++i) print i"\t"dist[i];}' "$chr"_csInterChromTotalMap;done | awk '{map[$1]+=$2;}END{for(i=0;i<1000;++i) print i"\t"map[i];}' > csDistanceDistr
cat Human.sortChr | awk 'function abs(x){return ((x < 0.0) ? -x : x)}{dis=abs($2-$5);if($1==$4 && dis<10000000){map[int(dis/1000)]++;}}END{for(i=0;i<10000;++i) print i"\t"map[i];}' > Human.empericalDist
cat Human.sortChr | awk 'function abs(x){return ((x < 0.0) ? -x : x)}{dis=abs($2-$5);if($1==$4 && dis<1000000){map[int(dis/100)]++;}}END{for(i=0;i<10000;++i) print i"\t"map[i];}' > Human.empericalDist.bin100				
cat Human.sortChr | awk 'function abs(x){return ((x < 0.0) ? -x : x)}BEGIN{nonloop=0;}{if($3!=$6 && (($2<$5 && $3==1) || ($2>$5 && $3==0))){nonloop++;}else{dis=abs($2-$5);if($1==$4 && dis<1000000){map[int(dis/100)]++;}}}END{for(i=0;i<10000;++i) {print i"\t"map[i];} print i"\t"nonloop}' > Human.bin100.withNonloop
Rscript /Code/3_Model_Distribution/empericalDist.r

Parameter:
You shoule revise the input files for variable "s" and "x", and also the name of output file ".polyCoef" in the R script.

Results:
"Mouse.polyCoef"

###3. Identify possible intradomain interactions.###
Order:
for chr in {1..22} X;
do 
   mkdir chr$chr
   cd chr$chr
   /Code/4_Find_IntraDomain_Interaction/findIntraDomainInteraction  /Output_path/chr$chr"_csInterChromTotalMap" /Output_path/chr$chr"_domainSitesMap" /Output_path/chr$chr"_domainCSinterFreq"  chr$chr"_out"  /Output_path/Human.polyCoef \-13.95 1.854 > chr$chr"_out_summary"
   cd ../
done;

Results:
For each chromosome, it generates a independent folder containing many "regionData" and "distMatrix" files, a "chr_out_debug" file and a "chr_out_summary" file.

###4. Lasso regression to select true interaction center.###
Order:
for chr in {1..22} X;
do 
   cd chr$chr
   domainNum=`tail -n1 chr"$chr"_out_debug | awk '{print $2}'`
   for (( c=0; c<=$domainNum; c++ ));
   do
   Rscript /Code/5_Lasso_Determine_Center/nnlasso.HiC.r regionData_$c distMatrix_$c > nnlassoOut_$c
   done;
cd ../
done;

Results:
In each chromosome folder, it generates many "nnlassoOut" files.

###5. Test reads distribution of possible interactions.###
Order:
for chr in {1..22} X;
do 
   cd chr$chr
   /Code/6_Test_Distribution/outputTestPosPvalue  /Output_path/chr$chr"_csInterChromTotalMap" /Output_path/chr$chr"_domainSitesMap" /Output_path/chr$chr"_domainCSinterFreq" chr$chr"_posP" /Output_path/Human.polyCoef \-13.95 1.854
   cd ../
done;

Results:
For each chromosome, it generates many "oneCol" and "testPosP" files.

###6. Collect all interactions.###
for chr in {1..22} X
do
    cd chr$chr
    domainNum=`tail -n1 chr"$chr"_out_debug | awk '{print $2}'`
    for (( c=0; c<=$domainNum; c++ ))
    do
        testNum1=`cat nnlassoOut_"$c" | tail -n2 | head -1 | awk '{print $2}'`;testNum2=`cat regionData_"$c" | tail -n2 | head -1 | awk '{print $2}'`;testNum3=`wc -l regionData_"$c" | awk '{print $1;}'`
        if [ "$testNum1" = "$testNum2" ] || [ "$testNum3" == 0 ]
        then
            echo $chr.$c.match
        else
            echo $chr.$c.notMatch
        fi
    done
    cd ../
done > domainNnlassoResults
for chr in {1..22} X
do
    echo $chr
    cd chr$chr
    domainNum=`tail -n1 chr"$chr"_out_debug | awk '{print $2}'`
    for (( c=0; c<=$domainNum; c++ ))
    do
        testNum1=`cat nnlassoOut_"$c" | tail -n2 | head -1 | awk '{print $2}'`;testNum2=`cat testPosP_"$c" | awk '{if(NF==7) {m=$3;}}END{print m;}'`;testNum3=`wc -l regionData_"$c" | awk '{print $1;}'`
        if [ "$testNum1" = "$testNum2" ] || [ "$testNum3" == 0 ]
        then
            echo $chr.$c.match
        else
            echo $chr.$c.notMatch
        fi
    done
    cd ../
done > domainTestPos
for chr in {1..22} X
do
    cd chr$chr
    domainNum=`tail -n1 chr"$chr"_out_debug | awk '{print $2}'`
    for (( c=0; c<=$domainNum; c++ ))
    do
	cat testPosP_"$c" | awk '{if(NF==7) {printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.10f\t%.10f\n",$1,$2,$3,$4,$5,$6,$7,$6,$7);}else{print $0;}}' > temp1
	cat nnlassoOut_"$c" | awk '{if($1!~/^\[/) {printf("%.10f",$1);for(i=2;i<=NF;i++){printf("\t%.10f",$i);}printf("\n");}else{print $0;}}' > temp2
	awk 'FNR==NR{if(NF==9){i=$2;j=1;beta[i]=$6;pval[i]=$7}else{if(NF==3){map[i,j]=$0;}j++;} next;}{if($1~/^\[/){i=$2;}else{if(i in pval){for(j=1;j<=NF;j++) {if(((i,j) in map) && $j!="Error!"){print map[i,j]"\t"beta[i]"\t"pval[i]"\t"$j;}}}}}' temp1  temp2 | sort -k1,1n -k2,2n -u | gzip > sigTest_"$c".all.gz
    done
    cd ../
done
for chr in {1..22} X
do
    echo $chr
    cd chr$chr
    domainNum=`tail -n1 chr"$chr"_out_debug | awk '{print $2}'`
    for (( c=0; c<=$domainNum; c++ ))
    do
	testNum1=`cat nnlassoOut_"$c" | tail -n2 | head -1 | awk '{print $2+1}'`;testNum2=`zcat sigTest_"$c".all.gz | wc -l | awk '{print $1;}'`;testNum3=`wc -l regionData_"$c" | awk '{print $1;}'`
	if [ "$testNum1" = "$testNum2" ] || [ "$testNum3" == 0 ]
	then
	    echo $chr.$c.match
	else
	    echo $chr.$c.notMatch
	fi
    done
    cd ../
done > domainSigTest.all
for chr in {1..22} X;do zcat chr$chr/sigTest_*.all.gz | awk '{print "'$chr'\t"$0;}';cat chr$chr/oneCol_* | awk '{print "'$chr'\t"$1"\t"$2"\t0\t"$3"\t"$4"\t"$3;}'; done | sort -k1,1 -k2,2n -k3,3n > Human_all_interactions

###7. Model random significance level.###
Order:
for chr in {1..22} X;
do 
    cd chr$chr
   /Code/7_Present_Significance/randomSamplingForFDR /Output_path/chr$chr"_csInterChromTotalMap" /Output_path/chr$chr"_domainSitesMap" /Output_path/chr$chr"_domainCSinterFreq" chr$chr /Output_path/Human.polyCoef \-13.95 1.854
   cd ../
done;

Results:
For each chromosome, it generates a "randomSamples_chr" file.

###8. Present significance level.###
Order:
cat chr*/randomSamples* > randomSamples.combined
Rscript /Code/7_Present_Significance/fdrFromRandomSamples.r

Results:
It generates a "randomSamples.combined.fdr" file.

Get BH_0.05 for FDR cutoff:
awk '{if($5<0.01 && $5>0.0099) print $0;}' randomSamples.combined.fdr | sort -k5,5n | tail -n1 | awk '{print $4}' > BH_0.01 (R:p.adjust function selects method"fdr")
awk '{if($6<0.01 && $6>0.0099) print $0;}' randomSamples.combined.fdr | sort -k6,6n | tail -n1 | awk '{print $4}' > BY_0.01 (R:p.adjust function selects method"BY")
#####BH_0.01=0.00000587

Get interactions below the FDR cutoff:
awk '{if($6<=0.00000587) print $0;}' Human_all_interactions > Human_fdr0.01
################################################################################################













 






 
 