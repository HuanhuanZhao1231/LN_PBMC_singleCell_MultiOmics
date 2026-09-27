
#library(patchwork)
library(Seurat)
library(harmony)
library(dplyr)
setwd("/public/home/zhaohuanhuan/scRNA/B/10HC11LN")
#hc2.data<-Read10X(data.dir="/public/home/zhaohuanhuan/scRNA/pbmcControl/HC2/outs/filtered_feature_bc_matrix")
#hc2<-CreateSeuratObject(counts=hc2.data, project="hc2", min.cells=5, min.features=400)
#hc2[["percent.mt"]] <-PercentageFeatureSet(hc2, pattern="^MT-")
#hc2<-subset(hc2, subset=nFeature_RNA>400&nFeature_RNA<2500&percent.mt<20)
#hc3.data<-Read10X(data.dir="/public/home/zhaohuanhuan/scRNA/pbmcControl/HC3/outs/filtered_feature_bc_matrix")
#hc3<-CreateSeuratObject(counts=hc3.data, project="hc3", min.cells=5, min.features=400)
#hc3[["percent.mt"]] <-PercentageFeatureSet(hc3, pattern="^MT-")
#hc3<-subset(hc3, subset=nFeature_RNA>400&nFeature_RNA<2500&percent.mt<20)
#hc4.data<-Read10X(data.dir="/public/home/zhaohuanhuan/scRNA/pbmcControl/HC4/outs/filtered_feature_bc_matrix")
#hc4<-CreateSeuratObject(counts=hc4.data, project="hc4", min.cells=5, min.features=400)
#hc4[["percent.mt"]] <-PercentageFeatureSet(hc4, pattern="^MT-")
#hc4<-subset(hc4, subset=nFeature_RNA>400&nFeature_RNA<2500&percent.mt<20)
#hc5.data<-Read10X(data.dir="/public/home/zhaohuanhuan/scRNA/pbmcControl/HC5/outs/filtered_feature_bc_matrix")
#hc5<-CreateSeuratObject(counts=hc5.data, project="hc5", min.cells=5, min.features=400)
#hc5[["percent.mt"]] <-PercentageFeatureSet(hc5, pattern="^MT-")
#hc5<-subset(hc5, subset=nFeature_RNA>400&nFeature_RNA<2500&percent.mt<20)
#hc6.data<-Read10X(data.dir="/public/home/zhaohuanhuan/scRNA/pbmcControl/HC6/outs/filtered_feature_bc_matrix")
#hc6<-CreateSeuratObject(counts=hc6.data, project="hc6", min.cells=5, min.features=400)
#hc6[["percent.mt"]] <-PercentageFeatureSet(hc6, pattern="^MT-")
#hc6<-subset(hc6, subset=nFeature_RNA>400&nFeature_RNA<2500&percent.mt<20)
#hc7.data<-Read10X(data.dir="/public/home/zhaohuanhuan/scRNA/pbmcControl/HC7/outs/filtered_feature_bc_matrix")
#hc7<-CreateSeuratObject(counts=hc7.data, project="hc7", min.cells=5, min.features=400)
#hc7[["percent.mt"]] <-PercentageFeatureSet(hc7, pattern="^MT-")
#hc7<-subset(hc7, subset=nFeature_RNA>400&nFeature_RNA<2500&percent.mt<20)
#hc8.data<-Read10X(data.dir="/public/home/zhaohuanhuan/scRNA/pbmcControl/HC8/outs/filtered_feature_bc_matrix")
#hc8<-CreateSeuratObject(counts=hc8.data, project="hc8", min.cells=5, min.features=400)
#hc8[["percent.mt"]] <-PercentageFeatureSet(hc8, pattern="^MT-")
#hc8<-subset(hc8, subset=nFeature_RNA>400&nFeature_RNA<2500&percent.mt<20)
#hc9.data<-Read10X(data.dir="/public/home/zhaohuanhuan/scRNA/pbmcControl/HC9/outs/filtered_feature_bc_matrix")
#hc9<-CreateSeuratObject(counts=hc9.data, project="hc9", min.cells=5, min.features=400)
#hc9[["percent.mt"]] <-PercentageFeatureSet(hc9, pattern="^MT-")
#hc9<-subset(hc9, subset=nFeature_RNA>400&nFeature_RNA<2500&percent.mt<20)
#hc10.data<-Read10X(data.dir="/public/home/zhaohuanhuan/scRNA/pbmcControl/HC10/outs/filtered_feature_bc_matrix")
#hc10<-CreateSeuratObject(counts=hc10.data, project="hc10", min.cells=5, min.features=400)
#hc10[["percent.mt"]] <-PercentageFeatureSet(hc10, pattern="^MT-")
#hc10<-subset(hc10, subset=nFeature_RNA>400&nFeature_RNA<2500&percent.mt<20)
#hc11.data<-Read10X(data.dir="/public/home/zhaohuanhuan/scRNA/pbmcControl/HC11/outs/filtered_feature_bc_matrix")
#hc11<-CreateSeuratObject(counts=hc11.data, project="hc11", min.cells=5, min.features=400)
#hc11[["percent.mt"]] <-PercentageFeatureSet(hc11, pattern="^MT-")
#hc11<-subset(hc11, subset=nFeature_RNA>400&nFeature_RNA<2500&percent.mt<20)
#pbmc1.data<-Read10X(data.dir="/public/home/zhaohuanhuan/scRNA/B/cellranger/code/LN1_B/outs/filtered_feature_bc_matrix")
#pbmc1<-CreateSeuratObject(counts=pbmc1.data, project="pbmc1k", min.cells=3, min.features=20)
#pbmc1[["percent.mt"]] <-PercentageFeatureSet(pbmc1, pattern="^MT-")
#pbmc1<-subset(pbmc1, subset=nFeature_RNA>200&nFeature_RNA<4000&percent.mt<10)
#pbmc2.data<-Read10X(data.dir="/public/home/zhaohuanhuan/scRNA/B/cellranger/code/LN2_B/outs/filtered_feature_bc_matrix")
#pbmc2<-CreateSeuratObject(counts=pbmc2.data, project="pbmc2k", min.cells=3, min.features=200)
#pbmc2[["percent.mt"]] <-PercentageFeatureSet(pbmc2, pattern="^MT-")
#pbmc2<-subset(pbmc2, subset=nFeature_RNA>200&nFeature_RNA<4500&percent.mt<10)
#pbmc3.data<-Read10X(data.dir="/public/home/zhaohuanhuan/scRNA/B/cellranger/code/LN3_B/outs/filtered_feature_bc_matrix")
#pbmc3<-CreateSeuratObject(counts=pbmc3.data, project="pbmc3k", min.cells=3, min.features=200)
#pbmc3[["percent.mt"]] <-PercentageFeatureSet(pbmc3, pattern="^MT-")
#pbmc3<-subset(pbmc3, subset=nFeature_RNA>200&nFeature_RNA<5000&percent.mt<10)
#pbmc4.data<-Read10X(data.dir="/public/home/zhaohuanhuan/scRNA/B/cellranger/code/LN4_B/outs/filtered_feature_bc_matrix")
#pbmc4<-CreateSeuratObject(counts=pbmc4.data, project="pbmc4k", min.cells=3, min.features=200)
#pbmc4[["percent.mt"]] <-PercentageFeatureSet(pbmc4, pattern="^MT-")
#pbmc4<-subset(pbmc4, subset=nFeature_RNA>200&nFeature_RNA<4000&percent.mt<10)
#pbmc5.data<-Read10X(data.dir="/public/home/zhaohuanhuan/scRNA/B/cellranger/code/LN5_B/outs/filtered_feature_bc_matrix")
#pbmc5<-CreateSeuratObject(counts=pbmc5.data, project="pbmc5k", min.cells=3, min.features=20)
#pbmc5[["percent.mt"]] <-PercentageFeatureSet(pbmc5, pattern="^MT-")
#pbmc5<-subset(pbmc5, subset=nFeature_RNA>200&nFeature_RNA<4500&percent.mt<10)
#pbmc6.data<-Read10X(data.dir="/public/home/zhaohuanhuan/scRNA/B/cellranger/code/LN6_B/outs/filtered_feature_bc_matrix")
#pbmc6<-CreateSeuratObject(counts=pbmc6.data, project="pbmc6k", min.cells=3, min.features=200)
#pbmc6[["percent.mt"]] <-PercentageFeatureSet(pbmc6, pattern="^MT-")
#pbmc6<-subset(pbmc6, subset=nFeature_RNA>200&nFeature_RNA<4500&percent.mt<10)
#pbmc7.data<-Read10X(data.dir="/public/home/zhaohuanhuan/scRNA/B/cellranger/code/LN7_B/outs/filtered_feature_bc_matrix")
#pbmc7<-CreateSeuratObject(counts=pbmc7.data, project="pbmc7k", min.cells=3, min.features=200)
#pbmc7[["percent.mt"]] <-PercentageFeatureSet(pbmc7, pattern="^MT-")
#pbmc7<-subset(pbmc7, subset=nFeature_RNA>200&nFeature_RNA<5000&percent.mt<10)
#pbmc8.data<-Read10X(data.dir="/public/home/zhaohuanhuan/scRNA/B/cellranger/code/LN8_B/outs/filtered_feature_bc_matrix")
#pbmc8<-CreateSeuratObject(counts=pbmc8.data, project="pbmc8k", min.cells=3, min.features=20)
#pbmc8[["percent.mt"]] <-PercentageFeatureSet(pbmc8, pattern="^MT-")
#pbmc8<-subset(pbmc8, subset=nFeature_RNA>200&nFeature_RNA<4500&percent.mt<10)
#pbmc9.data<-Read10X(data.dir="/public/home/zhaohuanhuan/scRNA/B/cellranger/code/LN9_B/outs/filtered_feature_bc_matrix")
#pbmc9<-CreateSeuratObject(counts=pbmc9.data, project="pbmc9k", min.cells=3, min.features=200)
#pbmc9[["percent.mt"]] <-PercentageFeatureSet(pbmc9, pattern="^MT-")
#pbmc9<-subset(pbmc9, subset=nFeature_RNA>200&nFeature_RNA<4500&percent.mt<10)
#pbmc10.data<-Read10X(data.dir="/public/home/zhaohuanhuan/scRNA/B/cellranger/code/LN10_B/outs/filtered_feature_bc_matrix")
#pbmc10<-CreateSeuratObject(counts=pbmc10.data, project="pbmc10k", min.cells=3, min.features=200)
#pbmc10[["percent.mt"]] <-PercentageFeatureSet(pbmc10, pattern="^MT-")
#pbmc10<-subset(pbmc10, subset=nFeature_RNA>200&nFeature_RNA<5000&percent.mt<15)
#pbmc11.data<-Read10X(data.dir="/public/home/zhaohuanhuan/scRNA/B/cellranger/code/LN11_B/outs/filtered_feature_bc_matrix")
#pbmc11<-CreateSeuratObject(counts=pbmc11.data, project="pbmc11k", min.cells=3, min.features=200)
#pbmc11[["percent.mt"]] <-PercentageFeatureSet(pbmc11, pattern="^MT-")
#pbmc11<-subset(pbmc11, subset=nFeature_RNA>200&nFeature_RNA<3500&percent.mt<20)
#pbmc.list<-list(hc2,hc3,hc4,hc5,hc6,hc7,hc8,hc9,hc10,hc11,pbmc1,pbmc2,pbmc3,pbmc4,pbmc5,pbmc6,pbmc7,pbmc8,pbmc9,pbmc10,pbmc11)
# normalize and identify variable features for each dataset independently
#pbmc.list <- lapply(X = pbmc.list, FUN = function(x) {
#    x <- NormalizeData(x)
#    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
#})
# select features that are repeatedly variable across datasets for integration
#features <- SelectIntegrationFeatures(object.list = pbmc.list)
###1快速整合用rpca
#pbmc.list<-lapply(X =pbmc.list, FUN =function(x){
#x<-ScaleData(x, features =features, verbose =FALSE)
#x<-RunPCA(x, features =features, verbose =FALSE)
#})
#immune.anchors <- FindIntegrationAnchors(object.list = pbmc.list, anchor.features = features, reduction = "rpca")
# this command creates an 'integrated' data assay
#immune.combined <- IntegrateData(anchorset = immune.anchors)
#save(immune.combined,file="integrated_immune.combined.Rdata")
#immune.combined <- ScaleData(immune.combined, verbose = FALSE)###必须先缩放再PCA
#immune.combined <- RunPCA(immune.combined, features = VariableFeatures(immune.combined), npcs = 30)
#run harmony
#immune.combined <- RunHarmony(immune.combined, group.by.vars = 'orig.ident',project.dim = F)######不加project.dim = F 会报错Error in data.use %*% cell.embeddings : non-conformable arguments
#save(immune.combined,file="afterHarmony_immune.combined.Rdata")
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
#immune.combined <- JackStraw(immune.combined, num.replicate = 100)
#immune.combined <- ScoreJackStraw(immune.combined, dims = 1:20)
#pdf("elbowplot.pdf")
#ElbowPlot(immune.combined)
#dev.off()
load("afterHarmony_immune.combined.Rdata")
DefaultAssay(immune.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
#for (i in c(6,7,8,9,10,11,12,13,14)) {
#    for (k in c(20,40,60)) {
#immune.combined <-immune.combined %>% 
#RunUMAP(reduction = "harmony", n.neighbors=22,min.dist = 0.2, dims = 1:i) %>% 
#FindNeighbors(reduction = "harmony", dims = 1:i,k.param =k) ####默认k.param =20，根据需要可改为其他
#immune.combined <- FindClusters(immune.combined,resolution =0.8,graph.name = "RNA_snn")###不加graph.name = "RNA_snn"报错 Provided graph.name not present in Seurat object
# Visualization
#pdfname <- paste0("integrationmin0.2","dim",i,"k",k,".pdf")
#p1 <- DimPlot(immune.combined, reduction = "umap",group.by="orig.ident",raster=FALSE)
#p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE,raster=FALSE)
#p <- p1 + p2
#ggsave(pdfname,plot=p,width=13,height=5)
#}
#    }
immune.combined <-immune.combined %>% 
RunUMAP(reduction = "harmony", n.neighbors=22,min.dist = 0.2, dims = 1:15) %>% 
FindNeighbors(reduction = "harmony", dims = 1:15,k.param =20) ####默认k.param =20，根据需要可改为其他
immune.combined <- FindClusters(immune.combined,resolution = 0.8,graph.name = "RNA_snn")###不加graph.name = "RNA_snn"报错 Provided graph.name not present in Seurat object
# Visualization
pdf("integration_15.8K20umap.pdf",height=5,width=13)
p1 <- DimPlot(immune.combined, reduction = "umap",group.by="orig.ident",raster=FALSE)
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE,raster=FALSE)
p <- p1+p2
p
dev.off()
DefaultAssay(immune.combined) <- "RNA"
#VlnPlot(immune.combined, features = c("FCGR3B","CSF3R","CMTM2","NCF1","SRGN","FPR1","TREM1"),pt.size=0,raster=FALSE)#LDG
#VlnPlot(immune.combined, features = c("C1QA","C1QB","C1QC","CD68"),pt.size=0,raster=FALSE)##C1Q+monocyte
pdf(file="umap15.8_minD0.2_lymph_Vlnplot.pdf",height=30,width=30)
VlnPlot(immune.combined, features = c("CD34","SOX4",#Progenitor
"STMN1","TOP2A",#Prolif
"CD3E",
"CD4","CD40LG","IL7R","TNFRSF4","RTKN2","FOXP3",
"CD8A",#CD8B没有CD8A好用
"CCR7",
"PRF1","GZMH","GZMK",
"NCAM1","GZMA","GZMB","KLRB1","GNLY","NKG7","FCGR3A",
"TCL1A","BANK1","MZB1",
"FCRL5",
"CD19",
"MS4A1",
"TCL1A","FCER2","IL4R",#B-immature,
"CD27","TNFRSF13B",#B-memory
"TNFRSF17","IGJ",##plasma
"FUT4","FCGR3B",#LDG
"ITGAM","RETN","FCGR2A","CSF3R","CMTM2","NCF1","SRGN","FPR1","TREM1"
),pt.size=0,raster=FALSE)
pdf(file="umap15.8_minD0.2_Monosubset_Vlnplot2.pdf",height=15,width=30)
VlnPlot(immune.combined, features = c("CD14","LYZ","C1QA","C1QB","CD68","ISG15","MX1","MX2","SOD2","CYBA","VIM","HLA-DRA","HLA-DQA1","HLA-DPB1","FCGR3A","MS4A7"
),pt.size=0,raster=FALSE)##monocyte亚型
pdf(file="umap15.8_minD0.2_DC_NeuVlnplot2.pdf",height=10,width=30)
VlnPlot(immune.combined, features = c("FCER1A","CST3","CD68",
"LILRA4","IL3RA",#pDC
"CD1C","ITGAX",#cDC
"CLEC10A","CD74",
"HBB","PPBP",
"SSC","ITGAM","RETN","FCGR2A","S100A8","IL1B","CSF3R","CMTM2","NCF1","SRGN","FPR1","TREM1",#Neu
"FCGR3B","FUT4"#LDG
),pt.size=0,raster=FALSE)#DC
save(immune.combined,file="unlabledImmune.combined15.8minD0.2.RData")
immune.combined<-RenameIdents(immune.combined, `0` ="CD8 ET", `1` ="Mono C", `2` ="CD4 NC",
    `3` ="Mono NC-I", `4` ="CD8 NC", `5` ="CD4 ET", `6` ="NK", `7` ="CD8 ET", `8` ="Mono C", `9` ="B IN",
    `10` ="CD4 NC", `11` ="B Mem",`12` ="Neu", `13` ="CD8 ET", `14` ="Mono C", `15` ="CD8 NC", `16` ="CD8 NC", `17` ="CD8 ET",
    `18` ="NKR", `19` ="Mega",`20` ="Neu", `21` ="Mono NC-I",`22` ="LDG", `23` ="B IN", `24` ="Prolif", `25` ="Mono NC", `26` ="CD8 ET",`27` ="Plasma",
    `28` ="CD8 ET", `29` ="Mono C",`30` ="cDC", `31` ="CD4 NC",`32` ="Mono C", `33` ="Mono C", `34` ="B IN", `35` ="B IN", `36` ="pDC",`37` ="CD8 NC"
    )###Mono3种

pdf("labled_integration_15.8umapmin0.2.pdf",height=5,width=13)
p1 <- DimPlot(immune.combined, reduction = "umap",group.by="orig.ident",raster=FALSE)
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE,raster=FALSE)
p1+p2
dev.off()
###美化Umap图
library(ggplot2)
library(dplyr)
library(scales)
library(ggrepel)
umap = immune.combined@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = Idents(immune.combined)) # 注释后的label信息 
head(umap)
#cellorder <- c("CD4T","Plasma","CD8 NC","CD8 ET","Mono C", "Mono NC","Mono T","Mono B","NKT","NK","B IN","B Mem","Prolif", "Progen","Neu","LDG","cDC","Mega")
#cellorder = c("CD8 ET","CD8 NC","CD4T","B IN", "B Mem","Neu","Mono T","Mono C","NK","Mono Int","NKT","Mono T","Prolif", "LDG","Mega","Mono NC","cDC","Plasma","Mono B","Progen")#目前的默认顺序
#umap$cell_type <- factor(umap$cell_type,levels = cellorder)
#umap <- umap[order(umap$cell_type), ]

pdf(file="MonoInt_umap_15.8_2.pdf",,width=5.5,height=5.5)
p <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  geom_point(size = 0.03 , alpha =1 )  +  scale_color_manual(values = allcolour)
#######调整umap图 - theme，去掉网格线，坐标轴和背景色即可
p2 <- p  +
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))
####3.2 调整umap图 - legend，legeng部分去掉legend.title后，调整标签大小，标签点的大小以及 标签之间的距离
p3 <- p2 +         
        theme(
          legend.title = element_blank(), #去掉legend.title 
          legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=20), #设置legend标签的大小
        legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5))) #设置legend中 点的大小 
####3.3 调整umap图 - annotation，坐标轴放到左下角可以通过ggplot2添加箭头和文本实现。
p4 <- p3 + 
  geom_segment(aes(x = min(umap$UMAP_1) , y = min(umap$UMAP_2) ,
                   xend = min(umap$UMAP_1) +3, yend = min(umap$UMAP_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(umap$UMAP_1)  , y = min(umap$UMAP_2)  ,
                   xend = min(umap$UMAP_1) , yend = min(umap$UMAP_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap$UMAP_1) +1.5, y = min(umap$UMAP_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(umap$UMAP_1) -1, y = min(umap$UMAP_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90) 
#3.4 调整umap图 - repel - labels
#1）计算每个cluster的median 坐标位置        
cell_type_med <- umap %>%
  group_by(cell_type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
##2）geom_label_repel 添加注释
##使用ggrepel包的repel函数可以使注释的标签不重叠。
###3）去掉legend
p4 +
geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med,
                   point.padding=unit(0.5, "lines")) +
  theme(legend.position = "none")
dev.off()
####添加分组信息
#
group <- immune.combined@meta.data$orig.ident
for(i_group in c("hc2",'hc3',"hc4","hc5","hc6",'hc7',"hc8","hc9","hc10","hc11") ){
    group[which(group==i_group)] <- 'Healthy'
}
for(i_group in c("pbmc4k",'pbmc10k',"pbmc11k","pbmc1k","pbmc8k",'pbmc6k','pbmc7k','pbmc5k','pbmc3k','pbmc9k','pbmc2k')){
    group[which(group==i_group)] <- 'LN'}
unique(group)##可以看group的值
immune.combined@meta.data[['Group']] <- group
#AIgroup
group <- immune.combined@meta.data$orig.ident
for(i_group in c("hc2",'hc3',"hc4","hc5","hc6",'hc7',"hc8","hc9","hc10","hc11") ){
    group[which(group==i_group)] <- 'Healthy'
}
for(i_group in c("pbmc4k",'pbmc10k',"pbmc11k","pbmc1k","pbmc8k",'pbmc6k','pbmc7k')){
    group[which(group==i_group)] <- 'Moderate'
}
for(i_group in c('pbmc5k','pbmc3k','pbmc9k','pbmc2k')){
    group[which(group==i_group)] <- 'Severe'
}####按顺序生成一个group列表
unique(group)##可以看group的值
immune.combined@meta.data[['AI2Group']] <- group
#----------------------DAI分组------------------------
group <- immune.combined@meta.data$orig.ident
for(i_group in c("hc2",'hc3',"hc4","hc5","hc6",'hc7',"hc8","hc9","hc10","hc11") ){
    group[which(group==i_group)] <- 'Healthy'
}
for(i_group in c('pbmc10k',"pbmc11k","pbmc5k","pbmc2k")){
    group[which(group==i_group)] <- 'Mild'
}
for(i_group in c('pbmc4k',"pbmc1k","pbmc8k")){
    group[which(group==i_group)] <- 'Moderate'
}
for(i_group in c('pbmc3k','pbmc6k','pbmc9k','pbmc7k')){
    group[which(group==i_group)] <- 'Severe'
}####按顺序生成一个group列表
unique(group)##可以看group的值
immune.combined@meta.data[['DAIGroup']] <- group
save(immune.combined,file="labled15.8immune.combinedCelltypeGroup.RData")
###分群比例
table(Idents(immune.combined),immune.combined$orig.ident)
Cellratio <- prop.table(table(Idents(immune.combined),immune.combined$orig.ident), margin = 2)
Cellratio
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("cluster","sample","proportion")
sample_order <- c("hc2",'hc3',"hc4","hc5","hc6",'hc7',"hc8","hc9","hc10","hc11","pbmc4k",'pbmc10k',"pbmc11k","pbmc1k","pbmc8k",'pbmc6k','pbmc7k','pbmc5k','pbmc3k','pbmc9k','pbmc2k')##AI排序
#sample_order2 <- c("hc1k",'hc2k',"hc3k","hc4k",'pbmc10k',"pbmc11k","pbmc5k","pbmc2k",'pbmc4k',"pbmc1k","pbmc8k",'pbmc3k','pbmc6k','pbmc9k','pbmc7k')##DAI排序
custom_order <- c("CD8 NC","CD8 ET","CD4 NC","CD4 ET","NK","NKR","Prolif","B IN","B Mem","Plasma","Mono C","Mono NC-I","Mono NC","Neu","LDG","cDC","pDC","Mega")

Cellratio$cluster <- factor(Cellratio$cluster, levels = custom_order, ordered = TRUE)
Cellratio <- Cellratio[order(Cellratio$cluster),]
Cellratio$sample <- factor(Cellratio$sample, levels = sample_order, ordered = TRUE)
Cellratio <- Cellratio[order(Cellratio$sample),]
library(ggplot2)
pdf(file="15.8cell_ratio_AI_Mono3.pdf",width=10,height=7)
ggplot(Cellratio) + 
  geom_bar(aes(x =sample, y= proportion, fill = cluster),stat = "identity",position="fill",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = cellratioColour)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
#-------------------------------------------------------------------------
table(immune.combined$orig.ident)#查看各组细胞数
prop.table(table(Idents(immune.combined)))
table(Idents(immune.combined), immune.combined$orig.ident)#各组不同细胞群细胞数
Cellratio <- prop.table((table(Idents(immune.combined), immune.combined$orig.ident)), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- data.frame(Cellratio)
library(reshape2)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]
###添加分组信息
group <- as.character(rownames(cellper))
for(i_group in c("hc2",'hc3',"hc4","hc5","hc6",'hc7',"hc8","hc9","hc10","hc11") ){
    group[which(group==i_group)] <- 'Healthy'
}
for(i_group in c("pbmc4k",'pbmc10k',"pbmc11k","pbmc1k","pbmc8k",'pbmc6k','pbmc7k','pbmc5k','pbmc3k','pbmc9k','pbmc2k')){
    group[which(group==i_group)] <- 'LN'
}
####按顺序生成一个group列表
unique(group)
cellper$group <- group
cellper$sample <- rownames(cellper)
###作图展示
pplist = list()
sce_groups = c("CD8 ET","CD8 NC","CD4 NC","CD4 ET","NK","NKR","B IN", "B Mem","Plasma","Neu","LDG","Mono C","Mono NC","Mono NC-I","Prolif", "cDC","pDC","Mega")
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
for(group_ in sce_groups){
  cellper_  = cellper %>% select(one_of(c('sample','group',group_)))
  colnames(cellper_) = c('sample','group','proportion')
  cellper_$proportion = as.numeric(cellper_$proportion)
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(proportion, 0.75), 
                                                      lower = quantile(proportion, 0.25),
                                                      mean = mean(proportion),
                                                      median = median(proportion))
  print(group_)
  print(cellper_$median)
  pp1 = ggplot(cellper_,aes(x=group,y=proportion)) + 
    geom_jitter(shape = 21,aes(fill=group),width = 0.25) + 
    stat_summary(fun=mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Percentage') +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
  
  ###组间t检验分析
  labely = max(cellper_$proportion)
  compare_means(proportion ~ group,  data = cellper_)
  my_comparisons <- list(c("Healthy","LN"))
  #my_comparisons <- list( c("Healthy", "Moderate"), c("Healthy", "Severe"), c("Moderate", "Severe") )#AIgroup
#my_comparisons <- list( c("Healthy", "Mild"), c("Healthy", "Moderate"), c("Healthy", "Severe"), c("Moderate", "Severe"), c("Mild","Moderate"),c("Mild","Severe") )#DAIgroup
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}

library(cowplot)
sce_groups = c("CD8 ET","CD8 NC","CD4 NC","CD4 ET","NK","NKR","Prolif","B IN", "B Mem","Plasma","Neu","LDG","Mono C","Mono NC","Mono NC-I", "cDC","pDC","Mega")
pdf("cellrationMonoC_NC_HCvsLN.pdf",width=12,height=15)
plot_grid(pplist[['CD8 ET']],
          pplist[['CD8 NC']],
          pplist[['CD4 NC']],
          pplist[['CD4 ET']],
          pplist[['B IN']],
          pplist[['B Mem']],
          pplist[['Plasma']],
          pplist[['Neu']],
          pplist[['LDG']],
          pplist[['Mono C']],
          pplist[['Mono NC']],
          pplist[['Mono NC-I']],
          pplist[['NK']],
          pplist[['NKR']],
          pplist[['Prolif']],
          pplist[['cDC']],
          pplist[['pDC']],
          pplist[['Mega']]
          )
##-------------------AUCell-------------------------
####准备矩阵###
library(GSEABase)
library(AUCell)
library(DelayedArray)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)
exprMatrix <- immune.combined@assays$RNA@data
###########准备GeneSets##########
gene_list <- read.table("~/scRNA/B/integration/5sample/data/ISG_Nature_immunity_heatmap_100gene.txt",header=FALSE)
isg_genes <- gene_list$V1
gene_list2 <- read.table("~/scRNA/B/integration/5sample/data/type1IFN.txt",header=FALSE)
ifn1_genes <- gene_list2$V1
gene_list3 <- read.table("~/scRNA/B/integration/5sample/data/Inflammatory.txt",header=FALSE)
infla_genes <- gene_list3$V1
gene_list4 <- read.table("~/scRNA/B/integration/5sample/data/cytokines.txt",header=FALSE)
cyto_genes <- gene_list4$V1
geneSets <- GeneSet(unique(cyto_genes), setName="geneSet4")
geneSets
## setName: geneSet1 
## geneIds: gene1, gene2, gene3 (total: 3)
## geneIdType: Null
## collectionType: Null 
## details: use 'details(object)'
cells_rankings <- AUCell_buildRankings(exprMatrix)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, 
                            aucMaxRank=nrow(cells_rankings)*0.05)###0.05试一下或者0.1
AUCell_auc <- as.numeric(getAUC(cells_AUC))
immune.combined@meta.data[['cytoscore']] <- AUCell_auc
#############################分样本画#######################
pdf("ISGscore_AUCell_DAIgroup_p.pdf",width=4,height=5)
#custom_order <- c("LN10K", "LN11K", "LN5K","LN2K","LN4K","LN1K","LN8K","LN6K","LN3K","LN9K","LN7K")
data <- FetchData(immune.combined,vars = c("ISGscore","DAIGroup"))
#data$orig.ident <- factor(data$orig.ident, levels = custom_order, ordered = TRUE)
#data <- data[order(data$orig.ident), ]
#my_comparisons <- list( c("Healthy", "Moderate"), c("Healthy", "Severe"), c("Moderate", "Severe"))
my_comparisons <- list( c("Healthy", "Mild"), c("Healthy", "Moderate"), c("Healthy", "Severe"), c("Mild", "Severe"), c("Mild", "Moderate"), c("Moderate", "Severe"))
#可以查看计算的p值stat.test <- data %>% t_test(ISGscore ~ DAIGroup) %>% adjust_pvalue(method = "bonferroni") %>% add_significance("p.adj")
p <- ggplot(data = data,aes(DAIGroup,ISGscore,fill = DAIGroup))
p + geom_boxplot(outlier.size=0.2) + theme_bw() + RotatedAxis() + labs(title = "ISGscore",y = "Score") + 
theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size = 10,face = "bold"),axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12)) +
coord_cartesian(ylim = c(0, 0.4)) +
stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test",label = "p.signif",label.y = c(0.34,0.36,0.38))
###p.signif是星，p.format是数值，geom_boxplot(outlier.size=0.2)设置点的大小
###
pdf("Inflascore_cell_DAIGroup.pdf",width=20,height=7)
data <- FetchData(immune.combined,vars = c("celltype","DAIGroup","Inflascore"))
#custom_order <- c("Endo", "PT", "PT2","LOH","DCT","PC","IC","MC","Tcell","Bcell","Macro","Neutro")
#data$celltype <- factor(data$celltype, levels = custom_order, ordered = TRUE)  # 将celltype转换为有序因子，并指定顺序  # 将celltype变量转换为因子
#data <- data[order(data$celltype), ]
group_order <- c("Healthy","Mild","Moderate","Severe")
#group_order <- c("Healthy","LN")
data$Group <- factor(data$DAIGroup, levels = group_order, ordered = TRUE) #将celltype转换为有序因子，并指定顺序  # 将celltype变量转换为因子
data <- data[order(data$DAIGroup), ]
p <- ggplot(data = data,aes(celltype,Inflascore,fill = DAIGroup))
p + geom_boxplot() + 
theme_bw() + RotatedAxis() + 
labs(title = "Inflascore",y = "Score") + 
theme(plot.title = element_text(hjust = 0.5),
axis.text = element_text(size = 10,face = "bold"),
axis.title.x = element_text(size = 12),
axis.title.y = element_text(size = 12))
dev.off()
#-------------------------------AUCell-FeaturePlot---------------------------
library(ggraph)
pdf("Cytoscore_Featureplot.pdf",width=10,height=10)
ggplot(data.frame(immune.combined@meta.data, immune.combined@reductions$umap@cell.embeddings), aes(UMAP_1, UMAP_2, color=cytoscore)
) + geom_point( size=0.1
) + scale_color_viridis(option="A")  + theme_light(base_size = 15)+labs(title = "cytoscore")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

