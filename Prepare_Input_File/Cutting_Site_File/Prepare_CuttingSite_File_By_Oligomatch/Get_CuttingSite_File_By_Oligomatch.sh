oligoMatch /Prepare_CuttingSite_File/MboI.fa  /path_to_reference_genome/Mouse_mm10.genome.fa MboI.mm10.bed
sort -k1,1 MboI.mm10.bed > MboI.mm10.sort.bed
awk 'BEGIN{OFS="\t"}{$1="chr"$1;print}' MboI.mm10.sort.bed > MboI.mm10.final.bed