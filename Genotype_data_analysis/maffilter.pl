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
open(OUTPUT,"+>$input.SNPs.dosage.txt");
open(OUTPUT2,"+>$input.SNPs.needtobeFiltered.lowMAF.txt");
print OUTPUT "SNP\tMAF\tAF_AlT\n";
open(INPUT,"$input.dosage.vcf");
my $row;

while(<INPUT>){
$row++;
my $pop;
if($row>=5){
@a=split;
$geshu=@a;
#print $geshu;
for($s=9;$s<$geshu;$s++){
@b=split/\:/,$a[$s];
#print "ff$b[0]\n";
$pop++;
$snps{$a[2]}=1;
$marker="$a[2],$pop";
if($b[0]=~/\.\/\./){$dosage{$marker}=0;}
elsif($b[0]=~/0\/0/){$dosage{$marker}=0;}
elsif($b[0]=~/0\/1/){$dosage{$marker}=1;}
elsif($b[0]=~/1\/1/){$dosage{$marker}=2;}
#print "$b[0]\t$marker\t$dosage{$marker}\n";
}
$num=$pop;#sample total counts
#到这里，对每行的snp每个人都dosage已经得到了。
#print "$marker\t$dosage{$marker}\n";
for($f=1;$f<=$num;$f++){
$mark2="$a[2],$f";
#下面这步得到所有人中的剂量
$add{$a[2]}=$add{$a[2]}+$dosage{$mark2};
#print "$f\t$num\t$a[2]\t$add{$a[2]}\t$dosage{$mark2}\t$mark2\n";
}

$af{$a[2]}=$add{$a[2]}/($num*2);#这步得到alt ref的等位频率
if($af{$a[2]}<0.5){
$maf{$a[2]}=$af{$a[2]};}
else{$maf{$a[2]}=1-$af{$a[2]};}

print OUTPUT "$a[2]\t$maf{$a[2]}\t$af{$a[2]}\n";
if($maf{$a[2]}<0.05){print OUTPUT2 "$a[2]\t$maf{$a[2]}\t$af{$a[2]}\n";}
}
}



