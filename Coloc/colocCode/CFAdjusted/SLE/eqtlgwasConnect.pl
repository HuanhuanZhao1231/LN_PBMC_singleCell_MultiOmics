#!/usr/bin/perl
open(INPUT1,"SNP.in.gwasCare.bed");#所有落在最终选定的gwas区域的snps
open(INPUT2,"egenes.in.gwaseGFR.bed");
open(INPUT3,"gwasOverCKDGen.bed");#final version of used gwas
open(OUTPUT,"+>eQTLCarePairs.txt");
open(OUTPUT2,"+>eGenesOfEGWAS.txt");
open(OUTPUT3,"+>SNPsinGWAS.txt");
while(<INPUT1>){#无重复的记录每个gwas区域所拥有的snp
@a=split;
$mark="$a[7],$a[3]";
if(exists $snpIngwas{$a[7]}){
if(exists $ji{$mark}){$skip=1;}else{
$snpIngwas{$a[7]}="$snpIngwas{$a[7]},$a[3]";}
}else{
$snpIngwas{$a[7]}="$a[3]";
$ji{$mark}=1;
}
}
while(<INPUT2>){#无重复的记录每个gwas区域的egene
@b=split;
$biao="$b[7],$b[3]";
if(exists $mqtlIngwas{$b[7]}){
if(exists $you{$biao}){$skip=1;}else{
$mqtlIngwas{$b[7]}="$mqtlIngwas{$b[7]},$b[3]";
}}else{
$mqtlIngwas{$b[7]}="$b[3]";
$you{$biao}=1;
}
}
while(<INPUT3>){#通过gwas区域将每个配对结合起来并输出
@c=split;
@s=split/\,/,$snpIngwas{$c[3]};
$ge=@s;
@q=split/\,/,$mqtlIngwas{$c[3]};
$l=@q;
for($j=0;$j<$l;$j++){
#print OUTPUT "$q[$j]\t";
for($i=0;$i<$ge;$i++){
print OUTPUT "$q[$j]\t$s[$i]\t$c[0]\t$c[1]\t$c[2]\t$c[3]\n";

}
}
print OUTPUT2 "$mqtlIngwas{$c[3]}\t$c[3]\n";
print OUTPUT3 "$snpIngwas{$c[3]}\t$c[3]\n";
}
