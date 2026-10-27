open(INPUT1,"GenCode.v19.PCGgeneSymbolList.bed");
open(INPUT2,"/public/database/share/11_LN_Chen/DataShare/ToXin2/forColoc/final_eqtl_10peer_nokgp_formated.txt");
open(OUTPUT,"+>eGenePCG.bed");
while(<INPUT2>){
@a=split;
$ji{$a[1]}=1;
}
while(<INPUT1>){
@b=split;
if(exists $ji{$b[4]}){print OUTPUT "$b[0]\t$b[1]\t$b[2]\t$b[4]\n";}
}