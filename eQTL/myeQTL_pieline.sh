######1.RNAseq数据提供---expression文件#####
expressionINT <- read.table("~/RNAseq/alignOut/99sample/expression/expressionINT.txt",header = TRUE,row.names =1,sep = "\t")
mapping <- read.table("~/eQTL/mydata/eQTL/99sample/LR_LN.txt",header = TRUE,sep = "\t")
###将expressionINT转置成行名为样本
t_expressionINT <- t(expressionINT)
###将rownames添加为列
row_names <- rownames(t_expressionINT)
t_expressionINT <- cbind(ID = row_names, t_expressionINT)
t_expressionINT[,"ID"]
rownames(mapping) <- mapping$LN
merged_matrix <- merge(t_expressionINT, mapping, by.x = "ID", by.y = "LR", all = FALSE)
rownames(merged_matrix) <- merged_matrix$LN
merged_matrix <- merged_matrix [, !(colnames(merged_matrix) %in% c("LN", "ID"))]
#####按covirates和snp文件的样本顺序排序
# 提取列名中的数字部分
ln_numbers <- as.numeric(sub("LN(\\d+)", "\\1", rownames(merged_matrix)))
# 根据ln_numbers对merged_matrix进行排序
merged_matrix <- merged_matrix [order(ln_numbers),]
#merged_matrix <- merged_matrix [, !(colnames(merged_matrix) %in% c("LN", "ID"))]
write.table(merged_matrix,file="~/eQTL/mydata/eQTL/99sample/expression/expressionINT.txt" , sep ="\t", row.names =TRUE,col.names =TRUE, quote =FALSE)
#expression <- read.table("~/eQTL/mydata/eQTL/99sample/expression/expressionINT.txt",header=TRUE)
merged_matrix <- t(merged_matrix)
write.table(merged_matrix,file="~/eQTL/mydata/eQTL/99sample/expression/expressionINTforeQTL.txt",sep="\t",quote=F)
##########################22222222222222genelocation看gene.pl和pheoyidy.pl################
Rscript PEER10_nokgp.R
Rscript prepare_peer.R
Rscript plotPeer.R
Rscript plotPeer.R

##############2.Genotype信息###
####99 samples#####
for ((i=1;i<=22;i++))
do
gtool -S \
--g ~/eQTL/mydata/eQTL/meth_chr${i}.gen \
--s ~/eQTL/mydata/eQTL/meth.sample \
--og ~/eQTL/mydata/eQTL/99sample/meth_chr${i}.gen \
--sample_id ~/eQTL/mydata/eQTL/filter.sample.id.txt
done

perl generate_impute2dosage_chr.pl
perl subimpute2dosage_chr.pl
sh subimpute2dosage_chr.sh
perl generate_maffilter_chr.pl
perl submaffilter_chr.pl
sh submaffilter_chr.sh
perl generate_mergeallvcf_chr.pl
perl submergeallvcf_chr.pl
sh submergeallvcf_chr.sh
perl snptitle.pl 
perl catvcf.pl
qsub -l nodes=1:ppn=4,vmem=50gb -o catVCF.o -e catVCF.e -N catVCF catVCF.sh
qsub -l nodes=1:ppn=4,vmem=50gb -o vcf2dos.o -e vcf2dos.e -N vcf2dos vcf2dos.sh
#将vcf文件压缩并索引
#Linux
condaac bioinfo
vcf-sort AllDosage.vcf | bgzip -c > AllDosage.vcf.gz
tabix -p vcf AllDosage.vcf.gz
###############3.确认不同文件样本顺序一致##########
Rscript check_sample_order.R
Rscript cpg_locationfiles_prepare.R 
#################4.运行MatrixeQTL#########
qsub -l nodes=1:ppn=4,vmem=64gb -o MatrixEQTL_30peer.o -e MatrixEQTL_30peer.e -N MatrixEQTL_30peer MatrixEQTL_30peer.sh
qsub -l nodes=1:ppn=4,vmem=64gb -o MatrixEQTL_25peer.o -e MatrixEQTL_25peer.e -N MatrixEQTL_25peer MatrixEQTL_25peer.sh
qsub -l nodes=1:ppn=4,vmem=64gb -o MatrixEQTL_20peer.o -e MatrixEQTL_20peer.e -N MatrixEQTL_20peer MatrixEQTL_20peer.sh
qsub -l nodes=1:ppn=4,vmem=64gb -o MatrixEQTL_15peer.o -e MatrixEQTL_15peer.e -N MatrixEQTL_15peer MatrixEQTL_15peer.sh
qsub -l nodes=1:ppn=4,vmem=64gb -o MatrixEQTL_10peer.o -e MatrixEQTL_10peer.e -N MatrixEQTL_10peer MatrixEQTL_10peer.sh
qsub -l nodes=1:ppn=4,vmem=64gb -o MatrixEQTL_5peer.o -e MatrixEQTL_5peer.e -N MatrixEQTL_5peer MatrixEQTL_5peer.sh
qsub -l nodes=1:ppn=4,vmem=64gb -o MatrixEQTL_no_peer.o -e MatrixEQTL_no_peer.e -N MatrixEQTL_no_peer MatrixEQTL_no_peer.sh
qsub -l nodes=1:ppn=4,vmem=64gb -o wc.o -e wc.e -N wc wc.sh
qsub -l nodes=1:ppn=4,vmem=256gb -o choose_peer_no.o -e choose_peer_no.e -N choose_peer_no choose_peer_no.sh
#5.运行fastqtl
qsub -l nodes=1:ppn=4,vmem=50gb -o phenotype_file_prepare.o -e phenotype_file_prepare.e -N phenotype_file_prepare phenotype_file_prepare.sh
head -1 title.txt > phenotype_title.txt
sort -k1,1 -k2,2n phenotype.unsorted.bed > phenotype.sorted.bed
cat phenotype_title.txt phenotype.sorted.bed > phenotype.bed
qsub -l nodes=1:ppn=4,vmem=50gb -o gzipbed.o -e gzipbed.e -N gzipbed gzipbed.sh
qsub -l nodes=1:ppn=4,vmem=50gb -o gzipvcf.o -e gzipvcf.e -N gzipvcf gzipvcf.sh
Rscript check_rownames_for_fastqtl.R
bgzip -k covariate_for_MatrixEQTL_10peer.txt
perl generate_fastqtl.pl
perl subfastqtl.pl
###########
sh subfastqtl.sh
################6.筛选significant eGene, eSNPs
qsub -l nodes=1:ppn=4,vmem=50gb -o zcat.o -e zcat.e -N zcat zcat.sh
qsub -l nodes=1:ppn=4,vmem=50gb -o find_cutoff.o -e find_cutoff.e -N find_cutoff find_cutoff.sh
qsub -l nodes=1:ppn=4,vmem=50gb -o find_final_eqtl.o -e find_final_eqtl.e -N find_final_eqtl find_final_eqtl.sh

#7.计算pi值
qsub -l nodes=1:ppn=4,vmem=50gb -o format_snpfile.o -e format_snpfile.e -N format_snpfile format_snpfile.sh
qsub -l nodes=1:ppn=4,vmem=50gb -o getpi.o -e getpi.e -N getpi getpi.sh











