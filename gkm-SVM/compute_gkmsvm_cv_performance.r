library(data.table)
library(ROCR)
cell_type <- rep(c("BIN","BMem","ABC","Plasma","CD4NC","CD4ET","CD8NC","CD8ET","NK","NKR","MonoC","MonoNC-I","MonoNC","DC","Neu"),each=6)
t <- rep(0:5,15)
kernal_function <- rep(c("gapped-kmer",
"estimated 1-mer with full filter",
"estimated 1-mer with truncated filter",
"gkm + RBF", 
"gkm + center weighted",
"gkm + center weighted + RBF"),15)
a <- NULL
for(i in c("BIN","BMem","ABC","Plasma","CD4NC","CD4ET","CD8NC","CD8ET","NK","NKR","MonoC","MonoNC-I","MonoNC","DC","Neu")){
    for(j in 0:5){
        data <- fread(paste0("kernel_",i,"_t_",j,".cvpred.txt"))
        pred <- prediction(data$V2,data$V3)
        auc <- performance(pred,"auc")
        print(paste(i,"t",j,"auc=",round(auc@y.values[[1]],3),sep=""))
        a <- c(a,round(auc@y.values[[1]],3))
    }
}
df <- data.frame(cell_type=cell_type,t=t,kernel_function=kernel_function,auc=a)
head(df)
library(openxlsx)
write.xlsx(df,"auc.xlsx")