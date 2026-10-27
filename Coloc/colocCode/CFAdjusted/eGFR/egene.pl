#!/usr/bin/perl
open(INPUT1,"eGenesOfEGWAS.txt");
open(INPUT3,"gwasOverCKDGen.bed");
open(OUTPUT,"+>eqtl.in.finalgwas.txt");
while(<INPUT1>){
@a=split;
$egene{$a[1]}=$a[0];
}
while(<INPUT3>){
@c=split;
@g=split/\,/,$egene{$c[3]};
$eall=@g;
for($i=0;$i<$eall;$i++){
print OUTPUT "$g[$i]\t$c[3]\n";

}
}
