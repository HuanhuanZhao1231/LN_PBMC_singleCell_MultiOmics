#!/usr/bin/perl

use Getopt::Long;
use POSIX;
my $version=1.00;

my %opts;
GetOptions(\%opts,"i=s");           #i:NUM
if (!(defined $opts{i})) { #necessary arguments
&usage;
}
my $input=$opts{'i'};

sub usage{
print <<"USAGE";
Version $version
Usage:
$0 -i 
options:
-i chromosome 
USAGE
exit(1);
}
open(INPUT,"meth_".$input.".gen");
open(INPUT2,"meth.sample");
while(<INPUT2>){
$row8++;
}
open(OUTPUT,"+>".$input.".dosage.vcf");
#open(OUTPUT2,"+>".$input.".Onlydosage.vcf");
print OUTPUT "##fileformat=VCFv4.1\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype calls\">\n##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage\">\n";
print OUTPUT "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t\n";
open(INPUT2,"meth.sample");
while(<INPUT2>){
$row5++;
@ks=split;
if($row5>=2){
if($row5!=$row8){
#print OUTPUT "$ks[0]\t";
}
else{
#print OUTPUT "$ks[0]\n";
}
}
}
@ran=split/r/,$input;

while(<INPUT>){
my $j;
my $pop;
@a=split;
$geshu=@a;
#$copy{$a[1]}=0;

for($i=5;$i<$geshu;$i++){
$j++;
if($j%3==1){$pop++;$mark="$pop,$a[1]";$copy{$mark}=0*$a[$i];
if($a[$i]>=0.9){$genotype{$mark}="0/0";}
}
if($j%3==2){
#if($a[$i]>0)
{$copy{$mark}=$copy{$mark}+1*$a[$i];
}
if($a[$i]>=0.9){$genotype{$mark}="0/1";}
}
if($j%3==0){
#if($a[$i]>0)
{$copy{$mark}=$copy{$mark}+2*$a[$i];}
if($a[$i]>=0.9){$genotype{$mark}="1/1";}
}
#print "$mark\t$copy{$mark},$a[$i],$j,$i\n";

}
$renshu=($geshu-5)/3;

#@ff=split/\:/,$a[1];
print OUTPUT "$ran[1]\t$a[2]\t$a[1]\t$a[3]\t$a[4]\t.\t.\t.\tGT:DS\t";
#print OUTPUT2 "$ran[1]\t$a[1]\t";
$person=$renshu-1;

for($s=1;$s<=$renshu;$s++){
	$jilu="$s,$a[1]";
if(exists $genotype{$jilu}){$skip=1;}
else{$genotype{$jilu}="./.";}
#print"$person,$renshu;,$s\n";
if($s==$renshu){

print OUTPUT "$genotype{$jilu}:$copy{$jilu}\n";
#print OUTPUT2 "$copy{$jilu}\n";

}else{
	
	print OUTPUT "$genotype{$jilu}:$copy{$jilu}\t";
	#print OUTPUT2 "$copy{$jilu}\t";
	}


}}


