open(INPUT1,"egenes.in.gwaseGFR.bed");
open(OUTPUT,"+>gwasOverCKDGen.bed");
while(<INPUT1>){
@a=split;
$mark="$a[4]\t$a[5]\t$a[6]\t$a[7]\n";
$ji{$mark}=1;
}
foreach $r (sort keys % ji){
print OUTPUT "$r";

};
