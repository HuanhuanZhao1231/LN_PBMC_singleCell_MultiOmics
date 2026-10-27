#!/usr/bin/perl
open(INPUT,"SNPsinGWAS.txt");
open(OUTPUT,"+>gwasSNPs.txt");
while(<INPUT>){
@a=split;
@b=split/\,/,$a[0];
$ge=@b;
for($i=0;$i<$ge;$i++){
print OUTPUT "$b[$i]\t$a[1]\n";
}

}