for($i=1;$i<=22;$i++){
open(OUTPUT,"+>mergeallvcf_chr$i.sh");
print OUTPUT "cd ~/eQTL/mydata/eQTL/99sample\nperl mergeallvcf.pl -i chr$i\n";
};