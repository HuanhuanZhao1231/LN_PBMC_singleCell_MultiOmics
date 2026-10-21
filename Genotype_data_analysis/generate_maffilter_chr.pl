for($i=1;$i<=22;$i++){
open(OUTPUT,"+>maffilter_chr$i.sh");
print OUTPUT "cd /public/home/zhaohuanhuan/eQTL/mydata/eQTL/99sample\nperl maffilter.pl -i chr$i\n";
};