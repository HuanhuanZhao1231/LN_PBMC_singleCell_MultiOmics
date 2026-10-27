#!/usr/bin/perl
open(INPUT2,"eQTLCarePairs.txt");
open(INPUT,"/public/home/shengxin/coloc/Allcispart1_peer10_nokgp.txt");
open(OUTPUT,"+>eqtlCare.txt");
while(<INPUT2>){
@b=split;
$mark="$b[1],$b[0]";
$ji{$mark}=1;
}
while(<INPUT>){
@c=split;
$biao="$c[1],$c[2]";
if(exists $ji{$biao}){
print OUTPUT "$_";
}


}
