library(data.table)
library(openxlsx)
load("~/eQTL/mydata/eQTL/99sample/biaoxing33.RData")
pc<-fread("~/eQTL/mydata/eQTL/99sample/tenPCs.txt")
clinical <- my.biao[,c("SampleNo","Age","Gender")]
colnames(clinical) <- c("ME_id","Age","Sex")
index<-read.xlsx("~/eQTL/mydata/eQTL/99sample/LN_ME_index.xlsx")
names(index) <- c("ME_id","LN_id")
clinical <- merge(clinical, index, by = "ME_id", sort = FALSE)
clinical$Sex <- ifelse(clinical$Sex == "F", 0,
                       ifelse(clinical$Sex == "M", 1, NA))
CF <- read.table("~/scRNA/B/forCibersort/CellFraction_fromCibersort.txt")
index2 <- read.table("~/eQTL/mydata/eQTL/99sample/LR_LN.txt",header=T)
setnames(CF, 1, "LR")
CF <- merge(CF, index2[, .(LR, LN)],by = "LR", all = FALSE)
names(CF)[1] <- "LR"
CF <- merge(CF, index2, by = "LR", sort = FALSE)
names(CF)[names(CF) == "LN"] <- "LN_id"
CF$LR <- NULL
pc <- pc[, c("LN_id", "PC1", "PC2", "PC3")]
df <- merge(CF, pc, by = "LN_id", sort = FALSE)
df <- merge(df, clinical[, c("LN_id", "Age", "Sex")],
            by = "LN_id", sort = FALSE)
# 7. 保存文件
write.table(df, "~/eQTL/mydata/eQTL/99sample/covariates/covariates_nokgp_nopeer.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
####后续与peer整合生成完整covarite文件
