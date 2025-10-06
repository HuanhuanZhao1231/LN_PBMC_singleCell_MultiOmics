#!/usr/bin/perl

# Define an array of cell types
my @cell_types = ("BIN","BMem","ABC","Plasma","CD4NC","CD4ET","CD8NC","CD8ET","NK","NKR","MonoC","MonoNC-I","MonoNC","DC","Neu");
#"BIN",这个已经运行过
# Loop through different cell types
foreach my $cell_type (@cell_types) {
    # Loop through different -t values
    for my $t (0..5) {
        my $pbs_script = "./${cell_type}/${cell_type}_t_${t}.pbs";
        
        open my $fh, '>', $pbs_script or die "Cannot open $pbs_script: $!";
        
        print $fh <<EOF;
# set the queue
#PBS -q batch

# set the path of wrong file
#PBS -e ${cell_type}${t}.e

# set the path of right file
#PBS -o ${cell_type}${t}.o

# nodes: number of nodes requested by job
# ppn: the number of processors per node requested by job
##############################################################
cd /public/home/zhaohuanhuan/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/gkmSVM/${cell_type}
~/software/lsgkm-0.1.1/src/gkmtrain -t $t -x 10 -m 16384 -T 4 ./posSet_${cell_type}.fa ./negSet_${cell_type}.fa ./kernal_${cell_type}_t_$t
EOF
        
        close $fh;
    }
}
