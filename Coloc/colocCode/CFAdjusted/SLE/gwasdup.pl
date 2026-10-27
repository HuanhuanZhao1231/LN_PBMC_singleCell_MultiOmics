#!/usr/bin/perl
open(INPUT,"GWAS2.txt");
open(OUTPUT,"+>GWAS3.txt");
while(<INPUT>){
$row++;@a=split;
if($row==1){print OUTPUT "$_";}
else{
if(exists $ji{$a[6]}){$skip=1;}
else{
print OUTPUT "$_";
$ji{$a[6]}=1;
}
}
}
