#!/usr/bin/perl
open(INPUT,"eqtlFormated.txt");
#open(OUTPUT1,"+>eqtlFormated2.txt");
open(OUTPUT2,"+>multiSNPs.txt");
while(<INPUT>){
@a=split;
$row++;
if($row>=2){
$mark="$a[6],$a[0]";
if(exists $you{$mark}){
$skip=1;
}else{
if(exists $ji{$a[6]}){$ji{$a[6]}=$ji{$a[6]}.",".$a[0];
$you{$mark}=1;
}else{
$ji{$a[6]}=$a[0];
$you{$mark}=1;
}
}
}
}
foreach $r (keys % ji){
@m=split/\,/,$ji{$r};
$ge=@m;
if($ge>1){
print OUTPUT2 "$r\t$ji{$r}\t$ge\n";}

}