
setwd("D:/yzw")


library(Seurat)
library(tidyverse)
library(harmony)
library(ggdist)
library(RColorBrewer)
library(scales)
library(patchwork)



library(ggplot2)
sample_info <- read.csv("./group.csv")
sce_list <- list()
setwd("D:/yzw/新建文件夹")
list <- list.files()
list
for (i in list) {
  # file = "GF1"
  path <- i
  Data10x <- read.csv(i)
  row.names(Data10x) <- Data10x$X
  Data10x <- Data10x[,-1]
  sce <- CreateSeuratObject(Data10x,min.cells = 3,min.features = 200,project = i)
  sce$sample <- i
  sce$sample <- substr(sce$sample,1,10)
  sce$group <- sample_info$group[which(sample_info$sample == i)]
  sce$background <- sample_info$background[which(sample_info$sample == i)]
  
  sce$age <- sample_info$age[which(sample_info$sample == i)] 
  sce$group1 <- sample_info$group1[which(sample_info$sample == i)]
  sce$Day <- sample_info$Day[which(sample_info$sample == i)]
  sce_list <- c(sce_list,sce)
}
sce_all <- merge(sce_list[[1]],sce_list[2:length(sce_list)],add.cell.ids = sample_info$sample)
table(sce_all$orig.ident)
View(sce_all@meta.data)

sce_all$percent.mt <- PercentageFeatureSet(sce_all, pattern = "^mt-")
sce_all[["ribo.mt"]] <- PercentageFeatureSet(sce_all, pattern = "^Rp")

col <- CellChat::scPalette(length(unique(sce_all$group)))
VlnPlot(sce_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","ribo.mt"),pt.size = 0,group.by = "group",ncol = 2)&scale_fill_manual(values = c("#E41A1C","#377EB8","#4DAF4A"))&theme(axis.title.x.bottom = element_blank())     

sce_all <- subset(sce_all,nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 20 & ribo.mt >1)
VlnPlot(sce_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","ribo.mt"),pt.size = 0,group.by = "group",ncol = 2)&scale_fill_manual(values = c("#E41A1C","#377EB8","#4DAF4A"))&theme(axis.title.x.bottom = element_blank())     
sce_split <- SplitObject(sce_all,split.by = "sample")
class(sce_split)

sce_split  <- lapply(X = sce_split, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = sce_split)
immune.anchors <- FindIntegrationAnchors(object.list = sce_split, anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- RunTSNE(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters

saveRDS(immune.combined,file = 'immune.combined.RDS')

DimPlot(immune.combined, reduction = "umap")


a=c("Tmem119","P2ry12","Cx3cr1",#Microglial 
    "Mrc1","F13a1","Cd68",#Macrophage
    "Ly6c2","Ccr2",#Monocyte
    "Snap25",#Neuron
    "Rgs5","Acta2",#VSMC
    "Cd19","Cd79a",#B
    "Cd44",#T
    "Flt1",#Endothelial
    "Mbp",#Oligodendrocyte
    "Gfap","Aqp4",#Astrocyte
    "Ttr",# schwann cell
    "Dcn",#Fibroblast
    "Abcc9",
    "Pdgfrb",
    "Camp","Cd69",
    "Cd52","Nkg7","Trbc2",
    "Cbr2","Cd163",
    "Lyz2","S100a8",
    "Kit","Tpsb2"
)
DotPlot(immune.combined, features = a, assay = "RNA") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),axis.title.x = element_blank()) +
  scale_color_gradientn(colours = c('#acd3e0','white','#e38f7a')) +
  labs(y = 'Subclusters')

#Microglial  0,2,12
#Monocyte  1
#Macrophage 3 11
#VSMC  23
#B  7
#T  4,8
#Lymphocyte 5

###Idents函数
Idents(immune.combined)##查看当前应用的metadata列
#Idents(immune.combined) <- "orig.ident "
Idents(immune.combined) <- "seurat_clusters" 

DefaultAssay(immune.combined)="integrated"
DefaultAssay(immune.combined) = "RNA"
#marks = FindMarkers(immune.combined,ident.1="2",logfc.threshold = 0.5,only.pos = T)
#ident.1##实验组
#ident.2##对照组  如不设置，则跟其余所有进行对比

#marks$a=marks$pct.1-marks$pct.2
#immune.combined
#DefaultAssay(immune.combined)

#marks = FindMarkers(immune.combined,ident.1="4",logfc.threshold = 0.5,only.pos = T)
#marks$a=marks$pct.1-marks$pct.2

#marks = FindMarkers(immune.combined,ident.1="5",logfc.threshold = 0.5,only.pos = T)
#marks$a=marks$pct.1-marks$pct.2

#marks = FindMarkers(immune.combined,ident.1="6",logfc.threshold = 0.5,only.pos = T)
#marks$a=marks$pct.1-marks$pct.2

#marks = FindMarkers(immune.combined,ident.1="4",logfc.threshold = 0.5,only.pos = T)
#marks$a=marks$pct.1-marks$pct.2

#marks = FindMarkers(immune.combined,ident.1="4",logfc.threshold = 0.5,only.pos = T)
#marks$a=marks$pct.1-marks$pct.2

immune.combined.markers <- FindAllMarkers(immune.combined, only.pos = TRUE)
immune.combined.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)


immune.combined.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n=3) %>%
  ungroup() -> top3
DoHeatmap(immune.combined, features = top3$gene) + NoLegend()

Idents(immune.combined)
a=subset(immune.combined,downsample=100)
table(a$seurat_clusters)


immune.combined.markers <- FindAllMarkers(a, only.pos = TRUE)
immune.combined.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

immune.combined.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n=3) %>%
  ungroup() -> top3
DoHeatmap(a, features = top3$gene) + NoLegend()














FeaturePlot(immune.combined, features = c("P2ry12",#Microglial 
                                          "Flt1","Cldn5","Esam","Ly6c1",#Endothelial
                                          "Mbp","Fa2h","Olig1",#Oligodendrocyte
                                          "Lyz2","Mrc1","Cd38","H2-Aa",#BAM
                                          "Gfap","Aqp4","Mfge8","Aldh1l1",#Astrocyte
))


celltype = data.frame(ClusterID = 0:25,celltype = 'NA')
#celltype[celltype$ClusterID %in% c( 0,5 ),2] = "Monocyte" 
#celltype[celltype$ClusterID %in% c( 1,9 ),2] = 'Macrophage'
#celltype[celltype$ClusterID %in% c( 2,3,7 ),2] = 'VSMC'
celltype[celltype$ClusterID %in% c( 0,2,11,12,15 ),2] = 'Microglia'
#celltype[celltype$ClusterID %in% c( 8,14 ),2] = 'Pericyte' 
#celltype[celltype$ClusterID %in% c( 10 ),2] = "Astrocyte" 
#celltype[celltype$ClusterID %in% c( 13 ),2] = "T" 
#celltype[celltype$ClusterID %in% c( 15 ),2] = 'BMEC'
#celltype[celltype$ClusterID %in% c( 16 ),2] = 'Oligo'
#celltype[celltype$ClusterID %in% c( 17 ),2] = 'DC' 
#celltype[celltype$ClusterID %in% c( 18,21 ),2] = 'OPC' 
#celltype[celltype$ClusterID %in% c( 19 ),2] = 'B' 
#celltype[celltype$ClusterID %in% c( 20 ),2] = 'Neuron' 
#celltype[celltype$ClusterID %in% c( 22 ),2] = 'Ependymal'  


head(celltype)
sce=immune.combined
sce$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$RNA_snn_res.0.4 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce@meta.data$celltype) 

saveRDS(immune.combined,file = "../immune.combine.2.rds")

immune.combined$group2=paste(immune.combined$group,immune.combined$group1,sep="_")
