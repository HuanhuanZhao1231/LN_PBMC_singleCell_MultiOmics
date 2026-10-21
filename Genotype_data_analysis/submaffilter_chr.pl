open(OUTPUT,"+>submaffilter_chr.sh");
for($i=1;$i<=22;$i++){
print OUTPUT "qsub -l nodes=1:ppn=4,vmem=50gb -m ae -o maffilter_chr$i.o -e maffilter_chr$i.e -N maffilter_chr$i maffilter_chr$i.sh\n";
}