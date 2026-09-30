library(MatrixEQTL)
library(data.table)
useModel = modelLINEAR;
SNP_file_name = paste("Dosage2.txt", sep="");
snps_location_file_name = paste("snpsloc.txt", sep="");
expression_file_name = paste("exp.txt", sep="");
gene_location_file_name = paste("gene.Pos.txt", sep="");
covariates_file_name = "peer10_nokgp_noCF.txt"
output_file_name_cis = tempfile("cisTransQQ","tmp");
output_file_name_tra = tempfile("transTransQQ","tmp");
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 0;
errorCovariance = numeric();
cisDist = 1e6;
snps = SlicedData$new();
snps$fileDelimiter = "\t";
snps$fileOmitCharacters="NA";
snps$fileSkipRows=1;
snps$fileSkipColumns=1;
snps$fileSliceSize = 2000;
snps$LoadFile(SNP_file_name);
gene = SlicedData$new();
gene$fileDelimiter = "\t";      
gene$fileOmitCharacters = "NA"; 
gene$fileSkipRows = 1;          
gene$fileSkipColumns = 1;       
gene$fileSliceSize = 2000;      
gene$LoadFile(expression_file_name);
#save(gene,file="MvalueTrans.RData")
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$LoadFile(covariates_file_name)
cvrt
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE,sep="\t");
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE,sep="\t");
me = Matrix_eQTL_main(
snps = snps, 
gene = gene, 
cvrt = cvrt,
output_file_name     = output_file_name_tra,
pvOutputThreshold     = pvOutputThreshold_tra,
useModel = useModel, 
errorCovariance = errorCovariance, 
verbose = TRUE, 
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos, 
genepos = genepos,
cisDist = cisDist,
pvalue.hist = "qqplot",
min.pv.by.genesnp = TRUE,
noFDRsaveMemory = FALSE);
me$cis$eqtl->cisres
write.table(cisres,file="Allcispart1_peer10_nokgp_noCF.txt",quote=F,sep='\t')