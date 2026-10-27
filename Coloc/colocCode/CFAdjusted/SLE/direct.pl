#!/usr/bin/perl
open(INPUT1,"GWAS.txt");
open(OUTPUT,"+>GWAS2.txt");
while(<INPUT1>){
$row++;
@a=split;
if($row>=2){

if($a[9] eq $a[10]){print OUTPUT "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\n";}
else{
$new=-$a[4];
print OUTPUT "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$new\t$a[5]\t$a[6]\t$a[7]\t$a[8]\n";}
}else{
print OUTPUT "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\n";
}
}
