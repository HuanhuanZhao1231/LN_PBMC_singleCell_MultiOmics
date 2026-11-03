for($i=1;$i<=100;$i++){
open(OUTPUT,"+>fastqtl_chunk$i.sh");
print OUTPUT "cd ~/eQTL/mydata/eQTL/99sample/fastQTL\n~/software/fastqtl-6p/bin/fastQTL.static --vcf ~/eQTL/mydata/eQTL/99sample/AllDosage.vcf.gz --bed phenotype.bed.gz --permute 10000 --out permutations.chunk$i.txt.gz --window 1e6 --cov peer10_nokgp.txt.gz --chunk $i 100\n";
};