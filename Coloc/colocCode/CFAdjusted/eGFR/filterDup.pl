#!/usr/bin/perl
open(INPUT,"eqtlFormated.txt");
open(INPUT2,"dup.txt");
open(OUTPUT,"+>eqtlFormated2.txt");
while(<INPUT2>){
@a=split;
$ji{$a[1]}=1;
}
while(<INPUT>){
@b=split;
$row++;
if($row==1){print OUTPUT "$_";}
elsif(exists $ji{$b[0]}){
$skip=1;
}
else{print OUTPUT "$_";}
}