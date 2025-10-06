# the name of job
#PBS -N H2

# set the queue
#PBS -q batch

#set the path of wrong file
#PBS -e H2.e

#set the path of  right file 
#PBS -o H2.o

# nodes: number of nodes requested by job
# ppn: the number of processors per node requested by job
#PBS -l nodes=2:ppn=24
#PBS -l walltime=10000:00:00
############################################
wd="~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/results/GWAS/ldsc"
sumstat_dir="~/GWAS/LDSC/sumstats"
baseline_dir=${wd}/1000G_Phase3_baselineLD_v2.2
# Weights:
weights_dir=${wd}/1000G_Phase3_weights_hm3_no_MHC
#frqfile=1000G_Phase3_frq.tgz
frq_dir=${wd}/1000G_Phase3_frq
# Compute LD scores with annot files
############################################

# Activate conda environment
source activate ldsc
cat ~/snATAC/B/ArchR/snATAC_ArchR_basedHairCode/scripts/LDSC/cellList.txt | while read c
do
cts_result_dir="${wd}/cts_results/${c}_cts_results"
ldscore_dir="${wd}/ldscore/${c}_ldscores"
# Directory for storing partitioned heritability results
h2_result_dir="${wd}/h2_results/${c}_h2_results"
for sumstat_f in ${sumstat_dir}/*.sumstats.gz
do
    # Get the prefix for output from sumstats file
    trait=$(basename ${sumstat_f} | sed 's/.sumstats.gz//')
    ~/software/ldsc/ldsc.py \
    --h2 ${sumstat_f} \
    --ref-ld-chr ${ldscore_dir}/${c}_,${baseline_dir}/baselineLD. \
    --w-ld-chr ${weights_dir}/weights.hm3_noMHC. \
    --frqfile-chr ${frq_dir}/1000G.EUR.QC. \
    --out ${h2_result_dir}/${c}.${trait} \
    --overlap-annot
done
done