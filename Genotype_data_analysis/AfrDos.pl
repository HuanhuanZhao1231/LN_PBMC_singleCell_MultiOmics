#!/usr/bin/perl
open(INPUT,"AllDosage.vcf");
open(OUTPUT,"+>Dosage2.txt");
open(OUTPUT2,"+>snpsloc.txt");
print OUTPUT2 "SNP\tChr\tPos\n";
while(<INPUT>){
$row++;
@a=split;


if($row==2){

print OUTPUT "$a[2]\t";
$geshu=@a;
$last=$geshu-1;
for($i=9;$i<$last;$i++){
print OUTPUT "$a[$i]\t";
}print OUTPUT "$a[$last]\n";
}
elsif($row>2){
print OUTPUT2 "$a[2]\t$a[0]\t$a[1]\n";
print OUTPUT "$a[2]\t";
$geshu=@a;
$last=$geshu-1;
for($i=9;$i<$last;$i++){
@c=split/\:/,$a[$i];

if($c[0] eq './.'){$dos=$c[1];}
elsif($c[0] eq '1/1'){
$pos=$c[1]/2;
if($pos<0.9){$dos=$c[1];}
else{$dos="2.00";}
}
elsif($c[0] eq '0/0'){
$dos=0;
}
elsif($c[0] eq '0/1'){
if($c[1]<0.9){
$dos=$c[1];}
else{$dos='1.00';}
}
else{$dos=$c[1];}

print OUTPUT "$dos\t";
}
@c=split/\:/,$a[$last];
if($c[0] eq './.'){$dos=$c[1];}
elsif($c[0] eq '1/1'){
$pos=$c[1]/2;
if($pos<0.9){$dos=$c[1];}
else{$dos="2.00";}
}
elsif($c[0] eq '0/0'){
$dos=0;
}
elsif($c[0] eq '0/1'){
if($c[1]<0.9){
$dos=$c[1];}
else{$dos="1.00";}
}
else{$dos=$c[1];}
print OUTPUT "$dos\n";
}
}
