######0.QC#######
trim_galore --fastqc --retain_unpaired --paired /data/center/11_LN_Chen/5_rnaseq/01.RawData/LR1_1.fq.gz /data/center/11_LN_Chen/5_rnaseq/01.RawData/LR1_2.fq.gz -o /public/home/zhaohuanhuan/RNAseq/test/chuli
###1.Index######
cd ~/lunix_lesson/rnaseq/raw_data
mkdir -p /n/scratch2/username/chr1_hg38_index
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir /public/home/zhaohuanhuan/RNAseq/reference/index/STAR/hg19/ --genomeFastaFiles /public/home/zhaohuanhuan/RNAseq/reference/index/STAR/GRCh37.p13.genome.fa --sjdbGTFfile /public/home/zhaohuanhuan/RNAseq/reference/index/STAR/gencode.v19.annotation.gtf --sjdbOverhang 149 ####sjdbOverhang测序长度最大值减1
###2.Align####
mkdir ../results/STAR #
STAR --runMode alignReads --runThreadN 16 \  
--genomeDir /public/home/zhaohuanhuan/RNAseq/reference/index/STAR/hg19 \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outFilterScoreMinOverLread 0.33 \
--outFilterMatchNminOverLread 0.33 \
--readFilesIn /public/home/zhaohuanhuan/RNAseq/test/chuli/LR1_1_val_1.fq.gz /public/home/zhaohuanhuan/RNAseq/test/chuli/LR1_2_val_2.fq.gz \
--readFilesCommand zcat \
--outFileNamePrefix LN1 \
--alignSoftClipAtReferenceEnds Yes \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outSAMunmapped Within KeepPairs \
--chimSegmentMin 15 \
--chimJunctionOverhangMin 15 \
--chimOutType WithinBAM SoftClip \
--outSAMattributes NH HI AS nM NM MD jM jI XS \
--outSAMattrRGline ID:rg1 SM:sm1
##3.MarkDuplicates##
cat /public/home/zhaohuanhuan/RNAseq/alignOut/RSEM/samplelist.txt | while read i 
do
python3 -u /public/home/zhaohuanhuan/eQTL/gtex-pipeline-master/rnaseq/src/run_MarkDuplicates.py /public/home/zhaohuanhuan/RNAseq/alignOut/STAR/${i}Aligned.sortedByCoord.out.bam ${i}
done
#####4.RSEM###########
mkdir rsem
cd rsem
#####4.1RSEM-Index###
rsem-prepare-reference --gtf /public/home/zhaohuanhuan/RNAseq/reference/index/STAR/gencode.v19.annotation.gtf \  ####gtf文件
--star \  ###STAR
-p 16 \
/public/home/zhaohuanhuan/RNAseq/reference/index/STAR/GRCh37.p13.genome.fa \  
/public/home/zhaohuanhuan/RNAseq/reference/RSEM/hg19 
####4.2transcript quantification###
cd /public/home/zhaohuanhuan/RNAseq/test/RSEM
rsem-calculate-expression --num-threads 16 \
--fragment-length-max 1000 \
--no-bam-output \
--paired-end \
--estimate-rspd \
--forward-prob 0 \
--bam /public/home/zhaohuanhuan/RNAseq/test/chuli/LN1Aligned.toTranscriptome.out.bam \
/public/home/zhaohuanhuan/RNAseq/reference/RSEM/hg19 \ 
LN1  ####
##4.3###conduct counts-matrix using RSEM results
rsem-generate-data-matrix *.genes.results > /public/home/zhaohuanhuan/RNAseq/alignOut/deseq_out/output.matrix
###
a=`ls *.genes.results | tr "\n" " "`
paste $a > /public/home/zhaohuanhuan/RNAseq/alignOut/deseq_out/alloutput.txt