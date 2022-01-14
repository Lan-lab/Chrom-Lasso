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