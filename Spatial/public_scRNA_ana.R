library(dplyr)
library(Seurat)
library(ggplot2)
library(reticulate)
library(reshape2)
library(cowplot)
library(patchwork)

set.seed(1234)

WHO_GEX.list <- list()
dir_list <- list.files("./", pattern="TD")

#for (prj in dir_list){
for ( i in c(1:length(dir_list))) {
	prj <- dir_list[i]
	SeuratOBJ.data <- Read10X(data.dir = prj)
	SeuratOBJ <- CreateSeuratObject(counts = SeuratOBJ.data, project = prj)
	SeuratOBJ[["percent.mt"]] <- PercentageFeatureSet(SeuratOBJ, pattern = "^MT-")
	SeuratOBJ[["hemoglobin"]] <- PercentageFeatureSet(SeuratOBJ, pattern = "^HB[AB]")
	subset_list <- SeuratOBJ$nFeature_RNA > 200 &
					SeuratOBJ$nFeature_RNA < 10000 &
					SeuratOBJ$nCount_RNA > 500 &
					SeuratOBJ$percent.mt < 10 &
					SeuratOBJ$hemoglobin < 10 

	SeuratOBJ <- SeuratOBJ[, subset_list]
	SeuratOBJ <- NormalizeData(SeuratOBJ)
	SeuratOBJ <- FindVariableFeatures(SeuratOBJ, selection.method = "vst", nfeatures =2000)
	SeuratOBJ <- RenameCells(SeuratOBJ, prj)
	WHO_GEX.list[i] <- SeuratOBJ
}


names(WHO_GEX.list) <- dir_list

features <- SelectIntegrationFeatures(object.list = WHO_GEX.list, nfeatures =2000)
Merge.anchors <- FindIntegrationAnchors(object.list = WHO_GEX.list ,  anchor.features = features)
filterSeuratOBJ <- IntegrateData(anchorset = Merge.anchors)
saveRDS(filterSeuratOBJ,"filterSeuratOBJ.rds" )


rm(WHO_GEX.list)
gc();gc()

DefaultAssay(filterSeuratOBJ) <- "integrated"
filterSeuratOBJ@meta.data %>% group_by(orig.ident) %>% summarise(cell=n(), medianUMI = median(nCount_RNA), medianGene = median(nFeature_RNA),)

##=----------------------------------
#  orig.ident  cell medianUMI medianGene
#  <chr>      <int>     <dbl>      <dbl>
#1 TD1        13741      2272        955
#2 TD2        17193      2455        940
#3 TD3        11565      3082       1377
#4 TD4        11249      3469       1328
#5 TD5        18593      2386        997
#6 TD6         8095      3149       1194
#7 TD7         3370      2366        982
#8 TD8        15181      2873       1179
#9 TD9        16260      2575       1053
##=----------------------------------


filterSeuratOBJ <- ScaleData(filterSeuratOBJ)
filterSeuratOBJ <- RunPCA(filterSeuratOBJ, features = VariableFeatures(object = filterSeuratOBJ))

pcadim=c(1:30)
resol = 0.5

filterSeuratOBJ <- FindNeighbors(filterSeuratOBJ, reduction = "pca", dims = pcadim)
filterSeuratOBJ <- FindClusters(filterSeuratOBJ, resolution = resol)
filterSeuratOBJ <- RunUMAP(filterSeuratOBJ, reduction = "pca", dims = 1:50)

saveRDS(filterSeuratOBJ,"filterSeuratOBJ_clustered.rds")

freq_table <- prop.table(x = table(Idents(filterSeuratOBJ), filterSeuratOBJ@meta.data[, "orig.ident"]), margin = 2)
freq_mmod<-melt(t(freq_table),id=c("sample","cluster","val"))
colnames(freq_mmod) <- c("sample","cluster","val")
freq_mmod$cluster <- as.factor( freq_mmod$cluster )

sink("frek_table_RNA.txt")
table(Idents(filterSeuratOBJ), filterSeuratOBJ@meta.data[, "orig.ident"])
sink()

sink("freq_tableprop_RNA.txt")
freq_table
sink()

filePath="./umap_RNA_seuraCluster.pdf"
p1 <- DimPlot(object = filterSeuratOBJ, group.by = "seurat_clusters", reduction = "umap", pt.size = 0.2,
        label = T) + ggtitle(label = "UMAP Seurat cluster")
p2 <- DimPlot(object = filterSeuratOBJ , group.by = "orig.ident", reduction = "umap", pt.size = 0.2,
        label = F) + ggtitle(label = "UMAP Sample")
p <- p1 + p2
ggsave(file = filePath, plot = p, dpi = 100, width = 12, height = 5)

filePath = "./umap_ClusterOnly.pdf"
p <- DimPlot(object = filterSeuratOBJ , group.by="seurat_clusters", reduction = "umap", pt.size = 0.2,
        label = T) + ggtitle(label = "UMAP Seurat cluster")
ggsave(file = filePath, plot = p, dpi = 100, width = 6, height = 5)

filePath <- paste0("./ClusterQC_Vinplot_.pdf")
p <- VlnPlot(filterSeuratOBJ, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = "seurat_clusters", ncol = 3, pt.size = 0) + NoLegend()
ggsave(file = filePath, plot = p, dpi = 100, width = 12, height = 5)


filePath = "./umap_ClusterOnly_Split_sample.pdf"
p <- DimPlot(object = filterSeuratOBJ, split.by = "orig.ident", reduction = "umap", pt.size = 0.2,
        label = T) + ggtitle(label = "UMAP Seurat cluster")
ggsave(file = filePath, plot = p, dpi = 100, width = 20, height = 5)


filterSeuratOBJ@meta.data %>% group_by(seurat_clusters) %>% summarise(cell=n(), medianUMI = median(nCount_RNA), medianGene = median(nFeature_RNA),)

##=------------------------------------------
#   seurat_clusters  cell medianUMI medianGene
#1                0 19625    2373.0      917.0
#2                1 15222    2041.0      925.0
#3                2 10028    2377.0      931.0
#4                3  7695    3174.0     1224.0
#5                4  6554    4479.5     1520.5
#6                5  6022    1779.0      877.0
#7                6  5732    2282.0     1088.0
#8                7  4999   14203.0     3300.0
#9                8  4760    2419.5     1063.0
#10               9  4565    1357.0      676.0
#11              10  4414    6040.5     1927.0
#12              11  4007    2725.0     1423.0
#13              12  3751    5550.0     1916.0
#14              13  3237    5083.0     2016.0
#15              14  3151    6036.0     1917.0
#16              15  2734    3651.5      847.5
#17              16  2516    9461.5     2663.0
#18              17  2179     899.0      456.0
#19              18  1396    6146.5     2122.5
#20              19   929    2606.0     1292.0
#21              20   667    6704.0     2898.0
#22              21   344    6946.5     1900.0
#23              22   241    5187.0     1929.0
#24              23   190    1814.0      805.0
#25              24   146    6585.5     2040.0
#26              25   143   22950.0     5285.0
##=------------------------------------------


DefaultAssay(filterSeuratOBJ) <- "RNA"

MerkerGene <- list( 
	"Epithelial"  = c("CAPS","SCGB1A1","WFDC2","KRT8","KRT19","SCGB3A1","SCGB3A2","KRT18","EPCAM","MUC1","SFTPC","SFTPA1","SFTPA2","SFTPB","PGC","AGER","CAV1","CYP4B1") ,
	"T_NK"        = c("GZMK","GZMB","NKG7","GNLY","KLRD1","KLRC1","TIGIT","LTB","TNFRSF4","TNFRSF18","IL2RA","GZMH","CD3D","FOXP3") ,
	"Myeloid"     = c("LYZ","CTSD","CD14","HLA-DRA","HLA-DRB1","TNF","FCGR3A","CD63","CD163","FCGR2A","TPSB2","TPSAB1","CPA3","HPGDS","MS4A2","AIF1","CD68","CCL3","FABP4","SPP1","CXCL9") ,
	"B"           = c("JCHAIN","IGHM","IGHG1","IGHA1","IGHG4","IGHA2","IGHG3","IGHG2","MZB1","IGHD","CD79A","MS4A1"),
	"Endothelial" = c("CLDN5","RAMP2","CAV1","VWF","PECAM1","CDH5","EMCN","CD34","CD31","THBD"),
	"Fibroblasts" = c("DCN","LUM","COL1A2","COL3A1","COL1A1","FN1","COL6A2","COL6A3 ACTA2","COL6A1","COL5A2","COL4A2")
)

MerkerGene2 <- list( 
	"Epithelial"  = c("SCGB1A1""SCGB3A1","EPCAM","MUC1","SFTPC","SFTPB") ,
	"T_NK"        = c("GZMK","GZMB","NKG7","GNLY","CD3D") ,
	"Myeloid"     = c("LYZ",) ,
	"Myeloid"     = c("LYZ","CTSD","CCL18","CD68","CD14","HLA-DRA","HLA-DRB1","TNF","FCGR3A","CD63","CD163","FCGR2A","TPSB2","TPSAB1","CPA3","HPGDS","MS4A2","AIF1") ,
	"B"           = c("CD79A","MS4A1"),
	"Immuno"      = c("JCHAIN","IGHM","IGHG1","IGHA1","IGHG4","IGHA2","IGHG3","IGHG2","MZB1","IGHD"),
	"Endothelial" = c("CLDN5","PECAM1"),
	"Fibroblasts" = c("DCN","LUM","COL1A2","COL3A1","COL1A1","FN1","COL6A2","COL6A3 ACTA2","COL6A1","COL5A2","COL4A2")
)


for (i in c(1:length(MerkerGene))){
	p <- DotPlot(filterSeuratOBJ, features = MerkerGene[i], cols = c("blue", "red"), dot.scale = 8, group.by="seurat_clusters") 
	filePath = paste0("./DotPlot_EXP_", names(MerkerGene)[i], ".pdf")
	ggsave(file = filePath, plot = p, dpi = 100, width = 15, height = 8)
}

CellType1 <- c("T","T","B","MAST","Epithelial","other","NK","Myeloid","T","Myeloid","Myeloid","Endothelial","Myeloid","Fibroblast","Myeloid","Immnuno_cell","Epithelial","other","other",
	"Epithelial","Myeloid","other","B","T","other","other")
filterSeuratOBJ[["CellType1"]] <- "NA"
for (i in c(1:length(CellType1))){
	filterSeuratOBJ@meta.data$CellType1[filterSeuratOBJ@meta.data$seurat_clusters == (i-1)] <- CellType1[i]
}

filePath = "./umap_ClusterOnly_annotation_tmp.pdf"
p<- DimPlot(object = filterSeuratOBJ , group.by = "CellType1", reduction = "umap", pt.size = 0.2,
        label = T) + ggtitle(label = "UMAP Seurat cluster")
ggsave(file = filePath , plot = p, dpi = 100, width = 6, height = 5)


## Subcluter
Ep_obj <- filterSeuratOBJ[, rownames(filterSeuratOBJ@meta.data)[filterSeuratOBJ@meta.data$CellType1 == "Epithelial"]]

DefaultAssay(Ep_obj) <- "integrated"
Ep_obj <- ScaleData(Ep_obj)
Ep_obj <- RunPCA(Ep_obj, features = VariableFeatures(object = Ep_obj))

pcadim = c(1:30)
resol = 0.5

Ep_obj <- FindNeighbors(Ep_obj, reduction = "pca", dims = pcadim)
Ep_obj <- FindClusters(Ep_obj, resolution = resol)
Ep_obj <- RunUMAP(Ep_obj, reduction = "pca", dims = 1:50)

filePath = "./Ep_obj_umap_RNA_seuraCluster.pdf"
p1 <- DimPlot(object = Ep_obj , group.by = "seurat_clusters", reduction = "umap", pt.size = 0.2,
        label = T) + ggtitle(label = "UMAP Seurat cluster")
p2 <- DimPlot(object = Ep_obj, group.by = "orig.ident", reduction = "umap", pt.size = 0.2,
        label = F) + ggtitle(label = "UMAP Sample")
p<- p1 + p2
ggsave(file = filePath, plot = p, dpi = 100, width = 12, height = 5)

table(Ep_obj@,etadata)

EpMarkers <- FindAllMarkers(Ep_obj, assay="RNA", group.by ="seurat_clusters", min.pct = 0.3, logfc.threshold = 0.5, only.pos =T)
fname <- paste0("Type_Marker_Clster", "_EP", ".txt")
write.table(EpMarkers, file=fname, quote=F, col.names=NA, sep="\t")

DefaultAssay(Ep_obj) <- "RNA"
Gene_list <- c("TM4SF1","CRABP2","UBE2C","PTPRC","CD8A","SFTPC","SFTPB","SFTPA1","SFTPA2","SCGB3A1","LAMP1","NAPSA","AGER","CAV1","BPIFB1","BPIFA1","SAA1")
p <- DotPlot( Ep_obj , features = Gene_list, cols = c("blue", "red"), dot.scale = 8, group.by = "seurat_clusters")
filePath = paste0("./Ep_DotPlot_EXP.pdf")
ggsave(file = filePath, plot = p, dpi=100, width=15 , height=8)

Ep_CellType1 <- c("Epithelial-1","Epithelial-2","Epithelial-2","Epithelial-1","Epithelial-1","AT2","AT1","Epithelial-2","Epithelial-1","Other","Club","Other","Other","Other","Epithelial-2")
for ( i in c(1:length(Myo_CellType1))){
    Ep_obj@meta.data$MyCellType[ Ep_obj@meta.data$seurat_clusters == (i-1)] <- Ep_CellType1[i]
}

filePath = "./Ep_obj_annotatedUmap.pdf"
p1 <- DimPlot(object = Ep_obj, group.by="seurat_clusters", reduction = "umap", pt.size = 0.2,
        label = T) + ggtitle(label = "UMAP Seurat cluster")
p2 <- DimPlot(object = Ep_obj, group.by = "MyCellType", reduction = "umap", pt.size = 0.2,
        label = F) + ggtitle(label = "UMAP Sample")
p <- p1 + p2
ggsave(file = filePath, plot = p, dpi = 100, width = 12, height = 5)

DefaultAssay(Ep_obj) <- "RNA"
p<- DotPlot(Ep_obj, features = Gene_list, cols = c("blue", "red"), dot.scale = 8, group.by = "MyCellType")
filePath = paste0("./Ep_DotPlot_EXPanno.pdf")
ggsave(file = filePath, plot = p, dpi = 100, width = 15, height = 8)

saveRDS(Ep_obj, "sub_Epithelial.obj")


# Subcluter
Myo_obj <- filterSeuratOBJ[, rownames(filterSeuratOBJ@meta.data)[ filterSeuratOBJ@meta.data$CellType1 %in% c("MAST","Myeloid")]]
DefaultAssay(Myo_obj) <- "integrated"

Myo_obj <- ScaleData(Myo_obj)
Myo_obj <- RunPCA(Myo_obj, features = VariableFeatures(object = Myo_obj))

pcadim = c(1:30)
resol = 0.5

Myo_obj <- FindNeighbors(Myo_obj, reduction = "pca", dims = pcadim)
Myo_obj <- FindClusters(Myo_obj, resolution = resol)
Myo_obj <- RunUMAP(Myo_obj, reduction = "pca", dims = 1:30)

freq_table <- prop.table(x = table(Idents(Myo_obj), Myo_obj@meta.data[, "orig.ident"]), margin = 2)
freq_mmod <- melt(t(freq_table), id = c("sample","cluster","val"))
colnames(freq_mmod) <- c("sample","cluster","val")
freq_mmod$cluster <- as.factor(freq_mmod$cluster)

filePath = "./Myo_obj_umap_RNA_seuraCluster.pdf"
p1 <- DimPlot(object = Myo_obj, group.by = "seurat_clusters", reduction = "umap", pt.size = 0.2,
        label = T) + ggtitle(label = "UMAP Seurat cluster")
p2 <- DimPlot(object = Myo_obj, group.by = "orig.ident", reduction = "umap", pt.size = 0.2,
        label = F) + ggtitle(label = "UMAP Sample")
p <- p1 + p2
ggsave(file = filePath, plot = p, dpi = 100, width = 12, height = 5)

DefaultAssay(Myo_obj) <- "RNA"
Gene_list <- c("TPSB2","CPA3","CD68","CCL17","FABP4","SPP1","S100B","CLEC9A","S100A9","FCN1","LYZ")
p <- DotPlot( Myo_obj, features = Gene_list, cols = c("blue", "red"), dot.scale = 8, group.by = "seurat_clusters")
filePath = paste0("./Myo_DotPlot_EXP.pdf")
ggsave(file = filePath, plot = p, dpi = 100, width = 15, height = 8)

MyoMarkers <- FindAllMarkers(Myo_obj, assay = "RNA", group.by = "seurat_clusters", min.pct = 0.1, logfc.threshold = 0.25)
fname <- paste0("Type_Marker_Clster","_Myo",".txt")
write.table(MyoMarkers, file=fname, quote=F, col.names=NA, sep="\t")


Myo_CellType1 <- c("Mastocyte","FABP4+ macrophage","CCL17+ macrophage","Monocyte","Mastocyte","Macrophage","Monocyte","SPP1+ Macrophage","Monocyte","DC","Other","Other","Other_myeloid","Mastocyte","DC","Mastocyte","Mastocyte")

for (i in c(1:length(Myo_CellType1))){
    Myo_obj@meta.data$MyCellType[ Myo_obj@meta.data$seurat_clusters == (i-1)] <- Myo_CellType1[i]
}

filePath = "./Myo_obj_annotatedUmap.pdf"
p1 <- DimPlot(object = Myo_obj, group.by = "seurat_clusters", reduction = "umap", pt.size = 0.2,
        label = T) + ggtitle(label = "UMAP Seurat cluster")
p2 <- DimPlot(object = Myo_obj, group.by = "MyCellType", reduction = "umap", pt.size = 0.2,
        label = F) + ggtitle(label = "UMAP Sample")
p <- p1 + p2
ggsave(file = filePath, plot = p, dpi = 100, width = 12, height = 5)

DefaultAssay(Myo_obj) <- "RNA"
Gene_list <- c("TPSB2","CPA3","CD68","CCL17","FABP4","SPP1","S100B","CLEC9A","S100A9","FCN1","LYZ")
p<- DotPlot(Myo_obj, features = Gene_list , cols = c("blue", "red"), dot.scale = 8, group.by = "MyCellType")
filePath = paste0("./Myo_DotPlot_EXPanno.pdf")
ggsave(file = filePath, plot = p, dpi = 100, width = 15, height = 8)

saveRDS(Myo_obj, "sub_myeloid.obj")


# Subcluter
Lym_obj <- filterSeuratOBJ[, rownames(filterSeuratOBJ@meta.data)[ filterSeuratOBJ@meta.data$CellType1 %in% c("T","NK","B")]]
DefaultAssay(Lym_obj) <- "integrated"
Lym_obj <- ScaleData(Lym_obj )
Lym_obj <- RunPCA(Lym_obj, features = VariableFeatures(object = Lym_obj))

pcadim = c(1:30)
resol = 0.5

Lym_obj <- FindNeighbors(Lym_obj, reduction = "pca", dims = pcadim)
Lym_obj <- FindClusters(Lym_obj, resolution = resol)
Lym_obj <- RunUMAP(Lym_obj, reduction = "pca", dims = 1:30)

filePath = "./Lym_obj_umap_RNA_seuraCluster.pdf"
p1 <- DimPlot(object = Lym_obj , group.by = "seurat_clusters", reduction = "umap", pt.size = 0.2,
        label = T) + ggtitle(label = "UMAP Seurat cluster")
p2 <- DimPlot(object = Lym_obj, group.by = "orig.ident", reduction = "umap", pt.size = 0.2,
        label = F)+ ggtitle(label = "UMAP Sample")
p <- p1 + p2
ggsave(file = filePath, plot = p, dpi = 100, width = 12, height = 5)


DefaultAssay(Lym_obj) <- "RNA"
Gene_list <- c("CD4","CD8A","FOXP3","IL7R","CCR7","LAG3","NKG7","NCAM1","MS4A1")
p <- DotPlot( Lym_obj, features = Gene_list , cols = c("blue", "red"), dot.scale = 8, group.by = "seurat_clusters")
filePath = paste0("./Lym_DotPlot_EXP.pdf")
ggsave(file = filePath, plot = p, dpi = 100, width = 15, height = 8)

LymMarkers <- FindAllMarkers(Myo_obj, assay="RNA", group.by ="seurat_clusters", min.pct = 0.1, logfc.threshold = 0.25)
fname <- paste0("Type_Marker_Clster","_Lym",".txt")
write.table(LymMarkers, file = fname, quote = F, col.names = NA, sep = "\t")

Lym_CellType1 <-c("CD4+ T","B","NK","CD8+ T","CD8+ T","CD4+ Treg","other T-1","CD8+ T","NKT","other T-2","other T-3","B")
for (i in c(1:length(Lym_CellType1))){
    Lym_obj@meta.data$MyCellType[Lym_obj@meta.data$seurat_clusters == (i-1)] <- Lym_CellType1[i]
}

filePath = "./Lym_obj_annotatedUmap.pdf"
p1 <- DimPlot(object = Lym_obj, group.by = "seurat_clusters", reduction = "umap", pt.size = 0.2,
        label = T) + ggtitle(label = "UMAP Seurat cluster")
p2 <- DimPlot(object = Lym_obj, group.by = "MyCellType", reduction = "umap", pt.size = 0.2,
        label = F) + ggtitle(label = "UMAP Sample")
p <- p1 + p2
ggsave(file = filePath, plot = p, dpi = 100, width = 12, height = 5)

DefaultAssay(Lym_obj) <- "RNA"
Gene_list <- c("CD4","CD8A","FOXP3","IL7R","CCR7","LAG3","NKG7","NCAM1","MS4A1")
p <- DotPlot(Lym_obj , features = Gene_list, cols = c("blue", "red"), dot.scale = 8, group.by = "MyCellType")
filePath = paste0("./Lym_DotPlot_EXPanno.pdf")
ggsave(file = filePath, plot = p, dpi = 100, width = 15, height = 8)

saveRDS(Lym_obj, "sub_lymphoid.obj")



filterSeuratOBJ@meta.data$CellType1[filterSeuratOBJ@meta.data$seurat_clusters == (i-1)] <- CellType1[i]
filterSeuratOBJ[["MyAnn"]] <- filterSeuratOBJ@meta.data$CellType1
filterSeuratOBJ@meta.data[rownames(Ep_obj@meta.data), "MyAnn"] <- Ep_obj@meta.data$MyCellType
filterSeuratOBJ@meta.data[rownames(Myo_obj@meta.data), "MyAnn"] <- Myo_obj@meta.data$MyCellType
filterSeuratOBJ@meta.data[rownames(Lym_obj@meta.data), "MyAnn"] <- Lym_obj@meta.data$MyCellType


filePath="./umap_RNA_seuraCluster_Myannotaton_v2.pdf"
p1 <- DimPlot(object = filterSeuratOBJ, group.by = "MyAnn",reduction = "umap", pt.size = 0.2,
        label = F) + ggtitle(label = "UMAP Seurat cluster")
p2 <- DimPlot(object = filterSeuratOBJ, group.by = "orig.ident", reduction = "umap", pt.size = 0.2,
        label = F) + ggtitle(label = "UMAP Sample")
p <- p1 + p2
ggsave(file = filePath, plot = p, dpi = 100, width = 12, height = 5)

saveRDS(filterSeuratOBJ, "filterSeuratOBJ_clustered_MyAnno.rds")
