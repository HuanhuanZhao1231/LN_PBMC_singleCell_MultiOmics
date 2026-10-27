#!/usr/bin/perl
open(INPUT1,"myresults.clump");
open(OUTPUT,"+>gwas100kb.bed");
while(<INPUT1>){
$row++;
if($row>=2){
@a=split/\t/,$_;
$chr="chr".$a[0];
$loc=$a[1];
$start=$loc-100000-1;
$end=$loc+100000;
print OUTPUT "$chr\t$start\t$end\t$a[4]\n";
}
}
