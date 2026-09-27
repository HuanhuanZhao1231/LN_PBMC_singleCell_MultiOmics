#读入gwas位点和snp位点的关联文件##SLE-GWAS,eGFR_GWAS
library(coloc)
read.table("../gwasSNPs.txt",header=F)->gwassnps
#读入所有eqtl
read.table("../eqtlFormated2.txt",header=T)->eqtl
#读入所有gwas
read.table("../GWAS3.txt",header=T)->gwas
gwas$N<-as.numeric(as.character(gwas$N))
#读入egene和gwas的连接文件
read.table("eqtl.in.finalgwasPart40.txt",header=F)->allconnect
#数据整理
colnames(allconnect)<-c("egene","gwas")
allconnect$egene<-as.character(allconnect$egene)
allconnect$gwas<-as.character(allconnect$gwas)
colnames(gwassnps)<-c("snp","gwas")
gwassnps$snp<-as.character(gwassnps$snp)
gwassnps$gwas<-as.character(gwassnps$gwas)
eqtl$SNP<-as.character(eqtl$SNP)
eqtl$MARKLOC<-as.character(eqtl$MARKLOC)
eqtl$eGene<-as.character(eqtl$eGene)
gwas$MARKLOC<-as.character(gwas$MARKLOC)
#先提取allconnect的第一个gwas区域
pp4<-NULL
num<-NULL
snpPmax<-NULL
bestSNP<-NULL
for(j in 1:dim(allconnect)[1]){
tpmgwas<-allconnect[j,]$gwas
tpmgene<-allconnect[j,]$egene
#取得第一个区域中需要研究的snps (gwas,eqtl,mqtl)
gwassnps[gwassnps$gwas %in% tpmgwas,]$snp->gwasNeed
#提取该区域中第一个关注的egene
if(length(gwasNeed)>0){
eqtl[eqtl$eGene %in% tpmgene,]->qtlInfo
gwasmark<-qtlInfo[qtlInfo$SNP %in% gwasNeed,]$MARKLOC
gwas[gwas$MARKLOC %in% gwasmark,]->gwasInfo
#求交集
intersect(gwasInfo$MARKLOC,qtlInfo$MARKLOC)->over2
#loc2是最终的gwas,mqtl,eqtl都取了交集的位点
qtlInfo[qtlInfo$MARKLOC %in% over2,]->eqtlF
gwas[gwas$MARKLOC %in% over2,]->gwas1
rownames(gwas1[grep("\\.",gwas1$SNP),])->del
gwas1[!(rownames(gwas1) %in% del),]->gwasF
rownames(gwasF)<-gwasF$MARKLOC
rownames(eqtlF)<-eqtlF$MARKLOC
rownames(eqtlF)->xu
gwasF[xu,]->gwasF
eqtlF[xu,]->eqtlF
all(rownames(eqtlF)==rownames(gwasF))
gwasF$SNP<-gwasF$MARKLOC
eqtlF$SNP<-eqtlF$MARKLOC
gwasF$Var<-gwasF$STD^2
eqtlF$Var<-eqtlF$STD^2
my.res<-coloc.abf(dataset1=list(snp=gwasF$MARKLOC,beta=gwasF$BETA,varbeta=gwasF$Var,N=gwasF$N,type="cc"),dataset2=list(snp=eqtlF$MARKLOC, beta=eqtlF$BETA,varbeta=eqtlF$Var,N=99,type="quant"),MAF=eqtlF$MAF)#coloc 
my.res$summary[6]->colocValue
SNPNum<-dim(eqtlF)[1]
t(as.matrix(my.res[[1]]))->abf
my.res[[2]]->snpsPP
max(snpsPP$SNP.PP.H4)->snpPPmax
snpsPP[snpsPP$SNP.PP.H4==snpPPmax,]$snp->snpid
strsplit(snpid,".",fixed=TRUE)[[1]][2]->ci
ci<-as.numeric(as.character(ci))
eqtlF$SNP[ci]->topSNP
pp4<-rbind(pp4,colocValue)
num<-rbind(num,SNPNum)
snpPmax<-rbind(snpPmax,snpPPmax)
bestSNP<-rbind(bestSNP,topSNP)
}
else{
pp4<-rbind(pp4,"NA")
num<-rbind(num,"NA")
snpPmax<-rbind(snpPmax,"NA")
bestSNP<-rbind(bestSNP,"NA")
}
}
allconnect$mark<-paste0(allconnect$gwas,",",allconnect$egene)
rownames(pp4)<-allconnect$mark
rownames(num)<-allconnect$mark
cbind(num,pp4,snpPmax,bestSNP)->res
colnames(res)<-c("NumSNPs","PP4","bestcolocSNP_PP4","bestcolocSNP")
write.table(res,file="AllcolocResGWASNCPart40.txt",quote=F,sep='\t')
