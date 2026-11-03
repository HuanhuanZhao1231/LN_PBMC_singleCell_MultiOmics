library(data.table)
cutoff_wo_dummy<-fread("Cutoff_eVariantFiltered.txt")
data_wo_dummy<-fread("~/eQTL/mydata/eQTL/99sample/matrixeqtl/Allcispart1_peer10_nokgp.txt")
head(cutoff_wo_dummy)
names(cutoff_wo_dummy)<-c("gene","cutoff")
sum(duplicated(cutoff_wo_dummy$gene))
cutoff_wo_dummy1 <- cutoff_wo_dummy[!duplicated(cutoff_wo_dummy$gene), ]####去除重复
head(data_wo_dummy)
eqtl<-merge(data_wo_dummy,cutoff_wo_dummy1,by="gene",sort=F)
head(eqtl)
eqtl<-eqtl[eqtl$`pvalue`<eqtl$cutoff,]
write.table(eqtl,"final_eqtl.txt",sep="\t",quote=F,row.names=F,col.names = T)
eqtl<-fread("final_eqtl.txt")
head(eqtl)
eqtl<-eqtl[!duplicated(eqtl$gene),]
nrow(eqtl)
write.table(eqtl,"final_egene.txt",sep="\t",quote=F,row.names=F,col.names = T)
