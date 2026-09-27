setwd("~/eQTL/mydata/eQTL/99sample/matrixeqtl")
library(data.table)
data0<-fread("Allcispart1_nopeer_nokgp.txt")
data1 <- fread("Allcispart1_peer1_nokgp.txt")
data2 <- fread("Allcispart1_peer2_nokgp.txt")
data3 <- fread("Allcispart1_peer3_nokgp.txt")
data4 <- fread("Allcispart1_peer4_nokgp.txt")
data5 <- fread("Allcispart1_peer5_nokgp.txt")
data10 <- fread("Allcispart1_peer10_nokgp.txt")
data15 <- fread("Allcispart1_peer15_nokgp.txt")
data20 <- fread("Allcispart1_peer20_nokgp.txt")
data25 <- fread("Allcispart1_peer25_nokgp.txt")
data30 <- fread("Allcispart1_peer30_nokgp.txt")



data0 <- data0[data0$FDR<0.05,]
data0 <- data0[!duplicated(data0$gene),]
data1 <- data1[data1$FDR<0.05,]
data1 <- data1[!duplicated(data1$gene),]
data2 <- data2[data2$FDR<0.05,]
data2 <- data2[!duplicated(data2$gene),]
data3 <- data3[data3$FDR<0.05,]
data3 <- data3[!duplicated(data3$gene),]
data4 <- data4[data4$FDR<0.05,]
data4 <- data4[!duplicated(data4$gene),]
data5 <- data5[data5$FDR<0.05,]
data5 <- data5[!duplicated(data5$gene),]
data10 <- data10[data10$FDR<0.05,]
data10 <- data10[!duplicated(data10$gene),]
data15 <- data15[data15$FDR<0.05,]
data15 <- data15[!duplicated(data15$gene),]
data20 <- data20[data20$FDR<0.05,]
data20 <- data20[!duplicated(data20$gene),]
datda25 <- data25[data25$FDR<0.05,]
data25 <- data25[!duplicated(data25$gene),]
data30 <- data30[data30$FDR<0.05,]
data30 <- data30[!duplicated(data30$gene),]

print (paste("0 peer",nrow(data0)))
print (paste("1 peer",nrow(data1)))
print (paste("2 peer",nrow(data2)))
print (paste("3 peer",nrow(data3)))
print (paste("4 peer",nrow(data4)))
print (paste("5 peer",nrow(data5)))
print (paste("10 peer",nrow(data10)))
print (paste("15 peer",nrow(data15)))
print (paste("20 peer",nrow(data20)))
print (paste("25 peer",nrow(data25)))
print (paste("30 peer",nrow(data30)))
