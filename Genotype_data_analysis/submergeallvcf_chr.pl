open(OUTPUT,"+>submergeallvcf_chr.sh");
for($i=1;$i<=22;$i++){
print OUTPUT "qsub -l nodes=1:ppn=4,vmem=50gb -m ae -o mergeallvcf_chr$i.o -e mergeallvcf_chr$i.e -N mergeallvcf_chr$i mergeallvcf_chr$i.sh\n";
}