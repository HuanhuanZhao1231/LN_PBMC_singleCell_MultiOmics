open(OUTPUT1,"+>getcode.sh");
open(OUTPUT2,"+>split/subcode.sh");
for($i=1;$i<=40;$i++){
open(OUTPUT3,"+>split/runcoloc$i.sh");
print OUTPUT1 "cat coloc.R > split/coloc$i.R\nsed -i s#Part1#Part$i#g split/coloc$i.R\n";
print OUTPUT3 "Rscript coloc$i.R\n";
print OUTPUT2 "nohup sh runcoloc$i.sh &\n";
}
