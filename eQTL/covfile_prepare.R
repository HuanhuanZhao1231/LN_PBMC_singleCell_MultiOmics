library(data.table)
library(openxlsx)
load("/public/database/share/11_LN_Chen/DataShare/ToHH/biaoxing33.RData")
pc<-fread("/public/database/share/11_LN_Chen/DataShare/ToHH/tenPCs.txt")
data<-read.csv("epic.csv")
head(data)
data$ID1<-paste(data$Sentrix_ID,data$Sentrix_Position,sep="_")
data$ID1 <- paste0("X",data$ID1)
data$person <- sub("-1","",data$person)
index<-read.xlsx("LN-ME编号对应.xlsx")
names(index) <- c("person","ID2")
data<-merge(data,index,by="person",sort=F)
data <- data[,7:8]
rownames(my.biao)<-paste0("X",rownames(my.biao))
my.biao<-my.biao[,c(2,3,8,11:17)]
head(my.biao)
names(my.biao)
my.biao$ID1 <- row.names(my.biao)
df <- merge(data,my.biao,by="ID1",sort=F)
head(df)
nrow(df)
dim(df)
load("NoDummyPEER_MvaluePC3New.RData")
factors[1:7,]
factors[,1]==my.biao$Age
head(my.biao)
df0 <- data.frame(factors[,11:40])
row.names(df0) <- row.names(my.biao)
df0$ID1 <- row.names(df0)
df <- merge(df0,df,by="ID1",sort=F)
names(pc)[1]<-"ID2"
df1 <- merge(df,pc,by="ID2")
dim(df1)
names(df1)
df1 <- df1[,c(1:2,33:53,3:32)]
df1 <- df1[,c(-2,-13)]
names(df1)[1] <-"id" 
row.names(df1) <- df1$id
df1$Slide <- as.character(df1$Slide)
table(df1$Slide)
head(df1)
table(df1$Slide)
test <- model.matrix(~.,data=df1[,-1])[,-1]
colnames(test)
df1 <- test[,c(-16,-21)]
df2<-t(df1)
covariate <- data.frame(df2)
covariate$id <- row.names(covariate)
colnames(covariate)
dim(covariate)
covariate <- covariate[,c(100,1:99)]
rownames(covariate)
dim(covariate)
# save(covariate,file="covariate.RData")
write.table(covariate,"covariate_for_MatrixEQTL_30peer.txt",sep="\t",quote=F,col.names=T,row.names = F)
print("30 done!")
covariate <- covariate[1:54,]
write.table(covariate,"covariate_for_MatrixEQTL_25peer.txt",sep="\t",quote=F,col.names=T,row.names = F)
print("25 done!")

covariate <- covariate[1:49,]
write.table(covariate,"covariate_for_MatrixEQTL_20peer.txt",sep="\t",quote=F,col.names=T,row.names = F)
print("20 done!")

covariate <- covariate[1:44,]
write.table(covariate,"covariate_for_MatrixEQTL_15peer.txt",sep="\t",quote=F,col.names=T,row.names = F)
print("15 done!")

covariate <- covariate[1:39,]
write.table(covariate,"covariate_for_MatrixEQTL_10peer.txt",sep="\t",quote=F,col.names=T,row.names = F)
print("10 done!")

covariate <- covariate[1:34,]
write.table(covariate,"covariate_for_MatrixEQTL_5peer.txt",sep="\t",quote=F,col.names=T,row.names = F)
print("5 done!")

covariate <- covariate[1:29,]
write.table(covariate,"covariate_for_MatrixEQTL_no_peer.txt",sep="\t",quote=F,col.names=T,row.names = F)
print("0 done!")


