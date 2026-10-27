#!/usr/bin/perl
`mkdir split`;
use POSIX;
open(INPUT1,"eqtl.in.finalgwas.txt");
while(<INPUT1>){
$total++;
}
close INPUT1;
#$total=8199;
$per=POSIX::floor($total/40);
for($i=1;$i<=40;$i++){
$start=($i-1)*$per+1;
$end=($i)*$per;
if($i==40){$end=$total;
}
print "$start\t$end\n";

open(OUTPUT,"+>split/eqtl.in.finalgwasPart$i.txt");
my $row;
open(INPUT1,"eqtl.in.finalgwas.txt");
while(<INPUT1>){
$row++;
if($row>=$start&&$row<=$end){
print OUTPUT "$_";
}

}close INPUT1;
}
