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
open(OUTPUT,"+>$input.Dosage.vcf");
open(INPUT2,"$input.SNPs.needtobeFiltered.lowMAF.txt");
while(<INPUT2>){
@b=split;
$ji{$b[0]}=1;
}
open(INPUT,"$input.dosage.vcf");
my $row;
while(<INPUT>){
$row++;
if(($row<=4)&&($i==1)){print OUTPUT "$_";}
else{
if($row>=5){
@c=split;
if(exists $ji{$c[2]}){$skip=1;}
else{
print OUTPUT "$_";}}
}

}

