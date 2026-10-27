#!/usr/bin/perl
open(INPUT,"TWB2_BBJ_eGFR_hg19_METAL_FUMA_noNA");
open(OUTPUT,"+>inputforswissGWAS.txt");
print OUTPUT "CHR\tPOS\tREF\tALT\tMARKER_ID\tPVALUE\tN\n";
while(<INPUT>){
$row++;
if($row>=2){
@a=split;
#@b=split/\:/,$a[0];
print OUTPUT "$a[8]\t$a[9]\t";
$ref=uc($a[1]);
$alt=uc($a[2]);

print OUTPUT "$ref\t$alt\t$a[8]:$a[9]_$ref/$alt\t$a[5]\t244952\n";
}
}
