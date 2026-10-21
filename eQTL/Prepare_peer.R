setwd("/public/home/zhaohuanhuan/eQTL/mydata/eQTL/99sample/peer")
load("/public/home/zhaohuanhuan/eQTL/mydata/eQTL/99sample/peer/PEER50_kgp.RData")
factors1 <- t(factors)
factors2 <- factors1[-12,]
peer1_nokgp <- factors2[1:12,]
peer2_nokgp <- factors2[1:13,]
peer3_nokgp <- factors2[1:14,]
peer4_nokgp <- factors2[1:15,]
peer5_nokgp <- factors2[1:16,]
peer10_nokgp <- factors2[1:21,]
peer15_nokgp <- factors2[1:26,]
peer20_nokgp <- factors2[1:31,]
peer25_nokgp <- factors2[1:36,]
peer30_nokgp <- factors2[1:41,]
peer35_nokgp <- factors2[1:46,]
peer40_nokgp <- factors2[1:51,]
peer45_nokgp <- factors2[1:56,]
peer50_nokgp <- factors2[1:61,]
cov <- read.table("~/eQTL/mydata/eQTL/99sample/covariates/covariates_nokgp_nopeer.txt",row.name = 1, header=TRUE)
row_1peer <- c(rownames(cov),"X1")
row_2peer <- c(row_1peer,"X2")
row_3peer <- c(row_2peer,"X3")
row_4peer <- c(row_3peer,"X4")
row_5peer <- c(row_4peer,"X5")
row_10peer <- c(row_5peer, "X6","X7","X8","X9","X10")
row_15peer <- c(row_10peer, "X11","X12","X13","X14","X15")
row_20peer <- c(row_15peer, "X16","X17","X18","X19","X20")
row_25peer <- c(row_20peer, "X21","X22","X23","X24","X25")
row_30peer <- c(row_25peer, "X26","X27","X28","X29","X30")
row_35peer <- c(row_30peer, "X31","X32","X33","X34","X35")
row_40peer <- c(row_35peer, "X36","X37","X38","X39","X40")
row_45peer <- c(row_40peer, "X41","X42","X43","X44","X45")
row_50peer <- c(row_45peer, "X46","X47","X48","X49","X50")
rownames(peer1_nokgp) <- row_1peer
colnames(peer1_nokgp) <- colnames(cov)
write.table(peer1_nokgp,file="peer1_nokgp.txt",sep="\t",quote=F)
rownames(peer2_nokgp) <- row_2peer
colnames(peer2_nokgp) <- colnames(cov)
write.table(peer2_nokgp,file="peer2_nokgp.txt",sep="\t",quote=F)
rownames(peer3_nokgp) <- row_3peer
colnames(peer3_nokgp) <- colnames(cov)
write.table(peer3_nokgp,file="peer3_nokgp.txt",sep="\t",quote=F)
rownames(peer4_nokgp) <- row_4peer
colnames(peer4_nokgp) <- colnames(cov)
write.table(peer4_nokgp,file="peer4_nokgp.txt",sep="\t",quote=F)
rownames(peer5_nokgp) <- row_5peer
colnames(peer5_nokgp) <- colnames(cov)
write.table(peer5_nokgp,file="peer5_nokgp.txt",sep="\t",quote=F)
rownames(peer10_nokgp) <- row_10peer
colnames(peer10_nokgp) <- colnames(cov)
write.table(peer10_nokgp,file="peer10_nokgp.txt",sep="\t",quote=F)
rownames(peer15_nokgp) <- row_15peer
colnames(peer15_nokgp) <- colnames(cov)
write.table(peer15_nokgp,file="peer15_nokgp.txt",sep="\t",quote=F)
rownames(peer20_nokgp) <- row_20peer
colnames(peer20_nokgp) <- colnames(cov)
write.table(peer20_nokgp,file="peer20_nokgp.txt",sep="\t",quote=F)
rownames(peer25_nokgp) <- row_25peer
colnames(peer25_nokgp) <- colnames(cov)
write.table(peer25_nokgp,file="peer25_nokgp.txt",sep="\t",quote=F)
rownames(peer30_nokgp) <- row_30peer
colnames(peer30_nokgp) <- colnames(cov)
write.table(peer30_nokgp,file="peer30_nokgp.txt",sep="\t",quote=F)
rownames(peer35_nokgp) <- row_35peer
colnames(peer35_nokgp) <- colnames(cov)
write.table(peer35_nokgp,file="peer35_nokgp.txt",sep="\t",quote=F)
rownames(peer40_nokgp) <- row_40peer
colnames(peer40_nokgp) <- colnames(cov)
write.table(peer40_nokgp,file="peer40_nokgp.txt",sep="\t",quote=F)
rownames(peer45_nokgp) <- row_45peer
colnames(peer45_nokgp) <- colnames(cov)
write.table(peer45_nokgp,file="peer45_nokgp.txt",sep="\t",quote=F)
rownames(peer50_nokgp) <- row_50peer
colnames(peer50_nokgp) <- colnames(cov)
write.table(peer50_nokgp,file="peer50_nokgp.txt",sep="\t",quote=F)