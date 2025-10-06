library(peer)
#load phenotype data and TPM data matrix (INT_Transformed)
load("$QTLDATAFILENAME.RData")
$INT_TRANSFORMED_TPM_DAT_AMATRIX->dat.int
exp <- dat.int
covs<-$PHENOTYPE_DATA_MATRIX
expr <- t(exp)
model = PEER()
PEER_setPhenoMean(model, as.matrix(expr))
PEER_setNk(model,10)
as.matrix(covs)->finalcov
PEER_setCovariates(model, finalcov)
PEER_setAdd_mean(model, TRUE)
PEER_update(model)
factors = PEER_getX(model)
residuals = PEER_getResiduals(model)
dat <- t(residuals)
dat.int <- matrix(,nrow(dat),ncol(dat))
for (i in 1:nrow(dat)){
  dat.int[i,] <- qqnorm(dat[i,],plot.it=F)$x
}
dat.int<-data.frame(dat.int)
rownames(dat.int)<-rownames(exp)
colnames(dat.int)<-colnames(exp)
dat.int->$TPMEXP
