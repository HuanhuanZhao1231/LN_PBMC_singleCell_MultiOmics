#!/usr/bin/perl
open(INPUT1,"/public/database/share/11_LN_Chen/DataShare/ToXin2/eqtl_zhh/AllDosage.vcf");
open(OUTPUT,"+>eqtlSNPsID.txt");
while(<INPUT1>){
$row++;
if($row>=3){
@a=split;
#@ran=split/r/,$a[0];

$maker="$a[0]:$a[1]_$a[3]/$a[4]";
print OUTPUT "$a[2]\t$maker\n";
}
}
