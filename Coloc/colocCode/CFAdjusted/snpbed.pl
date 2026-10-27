open(INPUT,"/public/database/share/11_LN_Chen/DataShare/ToXin2/eqtl_zhh/matrixeqtl/snpsloc.txt");
open(OUTPUT,"+>allSNPs.bed");
while(<INPUT>){
$row++;
if($row>=2){
@a=split;
$start=$a[2]-1;
print OUTPUT "$a[1]\t$start\t$a[2]\t$a[0]\n";
}
}

