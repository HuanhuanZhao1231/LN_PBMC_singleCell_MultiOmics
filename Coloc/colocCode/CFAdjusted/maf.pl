#!/usr/bin/perl
open(INPUT,"/public/database/share/11_LN_Chen/DataShare/ToXin2/eqtl_zhh/matrixeqtl/Dosage2.txt");
open(OUTPUT,"+>MAF.txt");
while(<INPUT>){
$row++;
if($row>1){
@a=split;
$ge=@a;
$ren=($ge-1)*2;
print "$ren\n";
for($i=1;$i<$ge;$i++){
$count{$a[0]}=$count{$a[0]}+$a[$i];

}
$af{$a[0]}=$count{$a[0]}/$ren;
if($af{$a[0]}>0.5){
$maf{$a[0]}=1-$af{$a[0]};}
else{$maf{$a[0]}=$af{$a[0]};}
print OUTPUT "$a[0]\t$maf{$a[0]}\n";
}
}
