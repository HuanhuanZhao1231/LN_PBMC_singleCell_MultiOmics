#zcat permutations.chunk*.txt.gz | gzip -c > permutations.all.chunks.txt.gz
opt_fdr<-0.05
opt_n<-99

library(data.table)
D = fread("permutations.all.chunks.txt.gz")
D <- data.frame(D)
head(D)
#D = read.table("perm.txt", hea=F, stringsAsFactors=F)
#rite.table(head(D),"~/test1.txt",sep="\t",quote=F,row.names=F,col.names = T)
colnames(D) = c("pid", "nvar", "shape1", "shape2", "dummy","sid", "dist","X2" ,"X3","X4","X5","npval", "slope","slope_se","ppval", "bpval")
###去掉低表达量的基因####
MASK=!is.na(D[,16])
D[is.na(D$bpval),]
Dnas=D[!MASK,]
D = D[MASK,]
library(qvalue)
D$st = qvalue(D$bpval)$qvalues
head(D,10)
dim(D[D$st <= 0.05, ])
#write.table(head(D,10), "~/test2.txt", sep="\t",quote=F, row.names=F, col.names=T)
write.table(D[D$st <= 0.05, ], "permutations.all.chunks.storeyFiltered.txt", quote=F, row.names=F, col.names=T)
Q = qvalue(D[,16]);
D$qval = NA;
head(D)
D$qval= Q$qvalue;
head(D)
set0 = D[which(D$qval <= opt_fdr),]
set1 = D[which(D$qval > opt_fdr),]
sort(set0$bpval)
pthreshold = (sort(set1$bpval)[1] - sort(-1.0 * set0$bpval)[1]) / 2
cat("  * Corrected p-value threshold = ", pthreshold, "\n")
pval0 = qbeta(pthreshold, D[,3], D[,4], ncp = 0, lower.tail = TRUE, log.p = FALSE)
test0 = qf(pval0, 1, D[,5], ncp = 0, lower.tail = FALSE, log.p = FALSE)
corr0 = sqrt(test0 / (D[,5] + test0))
test1 = opt_n * corr0 * corr0 / (1 - corr0 * corr0)
pval1 = pf(test1, 1, opt_n, ncp = 0, lower.tail = FALSE, log.p = FALSE)
cat("  * pval0 = ", mean(pval0), " +/- ", sd(pval0), "\n")
cat("  * test0 = ", mean(test0), " +/- ", sd(test0), "\n")
cat("  * corr0 = ", mean(corr0), " +/- ", sd(corr0), "\n")
cat("  * test1 = ", mean(test1), " +/- ", sd(test1), "\n")
cat("  * pval1 = ", mean(pval1), " +/- ", sd(pval1), "\n")
head(D)
D$nthresholds = pval1
D1=D[, c(1, 19)]
head(Dnas)
D2=Dnas[, c(1,16)]
names(D2)=names(D1)
D3=rbind(D1, D2)  
write.table(D3, "Cutoff_eVariantFiltered.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)