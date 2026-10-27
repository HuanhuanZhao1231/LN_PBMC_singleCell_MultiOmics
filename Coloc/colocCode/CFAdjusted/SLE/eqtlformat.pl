#!/usr/bin/perl
open(INPUT,"eqtlCare.txt");
open(INPUT2,"/public/database/share/11_LN_Chen/DataShare/ToXin2/eqtl_zhh/matrixeqtl/snpsloc.txt");
open(INPUT3,"../MAF.txt");
open(OUTPUT,"+>eqtlFormated.txt");
while(<INPUT3>){
$row2++;
if($row2>=1){
@c=split;
$maf{$c[0]}=$c[1];
}
}
while(<INPUT2>){
@b=split;
$row++;
if($row>=2){
$pos{$b[0]}=$b[2];
$chr{$b[0]}=$b[1];
}
}
print OUTPUT "SNP\tCHR\tPOS\tPVAL\tBETA\tSTD\tMARKLOC\tMAF\tN\teGene\n";
while(<INPUT>){
@a=split;
#beta/tstat=std
$std=$a[6]/$a[3];
print OUTPUT "$a[1]\t$chr{$a[1]}\t$pos{$a[1]}\t$a[4]\t$a[6]\t$std\t$chr{$a[1]}:$pos{$a[1]}\t$maf{$a[1]}\t303\t$a[2]\n";
}
