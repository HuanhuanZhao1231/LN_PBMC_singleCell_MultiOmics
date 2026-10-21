#!/usr/bin/perl
open(INPUT,"GCST90011866_Asian.txt");
open(OUTPUT,"+>inputforswissGWAS.txt");
print OUTPUT "CHR\tPOS\tREF\tALT\tMARKER_ID\tPVALUE\tN\n";
while(<INPUT>){
$row++;
if($row>=2){
@a=split;
#@b=split/\:/,$a[0];
print OUTPUT "$a[1]\t$a[2]\t";
$ref=uc($a[3]);
$alt=uc($a[4]);
print OUTPUT "$ref\t$alt\t$a[1]:$a[2]_$ref/$alt\t$a[8]\t$a[10]\n";
}
}
