#!/usr/bin/perl
open(INPUT1,"SNPsinGWAS.txt");
open(INPUT2,"../eqtlSNPsID.txt");
open(INPUT3,"GCST90011866_Asian.txt");
open(OUTPUT,"+>GWAS.txt");
while(<INPUT2>){
@c=split;
$row++;
if($row>=1){
@s1=split/\_/,$c[1];
$pos{$c[0]}="chr$s1[0]";
$raw{$pos{$c[0]}}=$c[1];
}
}
while(<INPUT1>){
@a=split;
@b=split/\,/,$a[0];
$ge=@b;
for($i=0;$i<$ge;$i++){
$snps{$pos{$b[$i]}}=1;#记录下所有需要的snp
}
}
##没有maf###暂时不显示呢##eGFR的allele1好像实际上是ref,2是alt#
print OUTPUT "SNP\tCHR\tPOS\tPVAL\tBETA\tSTD\tMARKLOC\tN\tGWASID\tmqtlID\n";
while(<INPUT3>){
@d=split;
$hang++;
if($hang>=2){
$mark="chr$d[1]:$d[2]";
$gwasid=$d[1].":".$d[2]."_".uc($d[3])."/".uc($d[4]);
print OUTPUT "$d[0]\tchr$d[1]\t$d[2]\t$d[8]\t$d[6]\t$d[7]\t$mark\t$d[10]\t";
print OUTPUT "$gwasid\t$raw{$mark}\n";
}
}