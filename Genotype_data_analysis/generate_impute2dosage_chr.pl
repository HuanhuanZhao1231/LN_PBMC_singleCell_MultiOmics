for($i=1;$i<=22;$i++){
open(OUTPUT,"+>impute2dosage_chr$i.sh");
print OUTPUT "cd ~/eQTL/mydata/eQTL/99sample\nperl impute2dosage.pl -i chr$i\n";
};