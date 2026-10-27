#!/usr/bin/perl
open(INPUT1,"multiSNPs.txt");
open(INPUT2,"inputforswissGWAS.txt");
open(OUTPUT,"+>dup.txt");
while(<INPUT1>){
@a=split;
$ji{$a[0]}=$a[1];
}
while(<INPUT2>){
@b=split;
$row++;
if($row>=2){
$mark="chr$b[0]:$b[1]";
if(exists $ji{$mark}){
@s=split/\,/,$ji{$mark};
$ge=@s;
for($i=0;$i<$ge;$i++){
@m=split/\:/,$s[$i];
$cob1="$b[0]:$b[1]_$m[2]/$m[3]";
$cob2="$b[0]:$b[1]_$m[3]/$m[2]";
if(($b[4] eq $cob1)||($b[4] eq $cob2)){
$skip=1;
}else{
print OUTPUT "$mark\t$s[$i]\n";
}

}
#$count{$mark}++;
#$record{$mark}=$_;
}
}
}
#foreach $r (sort keys% record) {print OUTPUT "$count{$r}\t$record{$r}";}
