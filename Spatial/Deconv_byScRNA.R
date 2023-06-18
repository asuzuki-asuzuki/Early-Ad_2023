library(dplyr)
library(tidyverse)
library(pheatmap)
library(Seurat)
library(ggplot2)
library(patchwork)
library(future)

plan("multisession", workers = 10)
options('future.globals.maxSize' = 500000*1024^2)

library(spacexr)
set.seed(1234)

Fun_lotate <- function(loc, do){
    c = do *pi/180
    x = loc[1]
    y = loc[2]
    x2 = x*cos(c) - y * sin(c)
    y2 = x*sin(c) + y * cos(c)
    return(c(x2, y2))
}


# Annotated scRNA-seq reference dataset
scRef.obj <- readRDS("filterSeuratOBJ_clustered_MyAnno.rds")

Idents(scRef.obj) <- "MyAnn"
unique(Idents(scRef.obj))[-grep("other", unique(Idents(scRef.obj)))]
scRef.obj <- subset(scRef.obj, ident= unique(Idents(scRef.obj))[-grep("other", unique(Idents(scRef.obj)))])

refCount <- as.matrix(GetAssayData(object = scRef.obj,assay="RNA", slot = "counts"))

cell_types <- scRef.obj@meta.data$MyAnn
names(cell_types) <- colnames(scRef.obj)
cell_types <- as.factor(cell_types)
nUMI <- scRef.obj@meta.data$nCount_RNA
names(nUMI) <- colnames(scRef.obj)
reference <- Reference(refCount, cell_types, nUMI, require_int = F, n_max_cells = 5000, min_UMI = 100)


# Visium dataset
Visium.obj <- readRDS("lung_visium.rds")

Visium_count <- as.matrix(GetAssayData(object = Visium.obj, assay = "Spatial", slot = "counts"))
nUMI_Visium  <- Visium.obj@meta.data$nCount_Spatial
names(nUMI_Visium) <- colnames(Visium.obj)
Visium_data  <- SpatialRNA(counts = Visium_count, nUMI = nUMI_Visium, use_fake_coords=T)

myRCTD <- create.RCTD(Visium_data, reference, max_cores = 10, CELL_MIN_INSTANCE = 20)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

saveRDS(myRCTD, 'myRCTD_Visium_full.rds')
#
myRCTD <- readRDS("myRCTD_Visium_full.rds")

weights <- myRCTD@results$weights
norm_weights <- normalize_weights(weights)
Visium.obj <- AddMetaData(Visium.obj, metadata = norm_weights)
SpatialFeaturePlot(Visium.obj, features = c("RCTDpredict_c9", alpha = c(0.1, 1)))


library(pheatmap)
visium_cluster = Visium.obj@meta.data[rownames(norm_weights), "seurat_clusters"]
cluste_inf <- as.data.frame(cbind(Visium_cluster = paste0(rep("C", n = length(visium_cluster)), visium_cluster)))
rownames(cluste_inf) <- rownames(norm_weights)
cluster_order <- order((cluste_inf$Visium_cluster))

p <- pheatmap(norm_weights[cluster_order,], show_rownames = F, show_colnames = T, annotation_row = cluste_inf, cluster_row = T)
filePath = "./Deconv_heatmap.pdf"
ggsave(file = filePath, plot = p, dpi = 100, width = 10, height = 8)
write.table(norm_weights, file = "RCTD_semibulk.txt", quote = F, col.names = NA, sep = "\t")


library(Polychrome)
source("./VIS_topic.R")
theta = norm_weights
Pos = Visium.obj@images$slice1@coordinates[,c(4,5)]
new_coord <-  t(apply(Pos, 1, Fun_lotate, -90))

colnames(new_coord) = c("x","y")
topicCol <- DiscretePalette( ncol(theta), palette = "polychrome", shuffle = FALSE)
require(scales)
topicCol <- hue_pal()(length(colnames(norm_weights) ))


p <- vizAllTopics(theta, new_coord, topicOrder = seq(ncol(theta)), topicCols = topicCol , showLegend = TRUE, r=60, lwd = 0)
filePath = "./Decon_VisTopic.pdf"
ggsave(file = filePath, plot = p, dpi = 100, width = 12, height = 8)

p1_lists <- list()

for (clus in colnames(norm_weights)) {
	p1 <- SpatialFeaturePlot(Visium.obj, features = c(clus, alpha = c(0.1, 1)))	
	p1_lists <- c(p1_lists, list(p1))
}

filePath <- paste0("./Decov_pct.pdf")
ggsave(file = filePath, plot = wrap_plots(p1_lists), dpi = 100, width = 30, height = 30)

filePath <- paste0("./Decov_pct.svg")
ggsave(file = filePath, plot = wrap_plots(p1_lists), dpi = 100, width = 30, height = 30)
