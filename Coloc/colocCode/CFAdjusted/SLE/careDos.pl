open(INPUT,"/public/database/share/11_LN_Chen/DataShare/ToXin2/eqtl_zhh/matrixeqtl/Dosage2.txt");
open(INPUT2,"SNP.in.gwasCare.bed");
open(OUTPUT,"+>CareDosage.txt");
while(<INPUT2>){
@a=split;
#$mark="$a[0]\t$a[1]\t$a[2]\t$a[3]\n";
$ji{$a[3]}=1;
}
while(<INPUT>){
$row++;
if($row==1){print OUTPUT "$_";}
else{
@b=split;
if(exists $ji{$b[0]}){print OUTPUT "$_";}
}
}