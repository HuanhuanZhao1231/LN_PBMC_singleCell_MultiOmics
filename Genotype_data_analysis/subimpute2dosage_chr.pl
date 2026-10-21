open(OUTPUT,"+>subimpute2dosage_chr.sh");
for($i=1;$i<=22;$i++){
print OUTPUT "qsub -l nodes=1:ppn=4,vmem=50gb -m ae -o impute2dosage_chr$i.o -e impute2dosage_chr$i.e -N impute2dosage_chr$i impute2dosage_chr$i.sh\n";
}