#!/usr/bin/perl
open(OUTPUT,"+>catVCF.sh");
print OUTPUT "cd ~/eQTL/mydata/eQTL/99sample\ncat snptitle.txt ";
for($i=1;$i<=22;$i++){
print OUTPUT "chr$i.Dosage.vcf ";
}
print OUTPUT "> AllDosage.vcf";

