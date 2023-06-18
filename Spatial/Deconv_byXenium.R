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
    return(c(x2,y2))
}

#Load Xenium and Visium dataset
xenium.obj <- readRDS("lung_xenium.rds")
Visium.obj <- readRDS("lung_visium.rds")


# Make Xenium Reference data
refCount <- as.matrix(GetAssayData(object = xenium.obj, assay="Xenium", slot = "counts"))
DimPlot(xenium.obj)
cell_types <- xenium.obj@meta.data$seurat_clusters
names(cell_types) <- colnames(xenium.obj)
cell_types <- as.factor(cell_types)
nUMI <- xenium.obj@meta.data$nCount_Xenium
names(nUMI) <- colnames(xenium.obj)
reference <- Reference(refCount, cell_types, nUMI, require_int = F, n_max_cells = 100000, min_UMI = 10)

# Visium
Visium_count <- as.matrix(GetAssayData(object = Visium.obj, assay = "Spatial", slot = "counts"))
nUMI_Visium <- Visium.obj@meta.data$nCount_Spatial
names(nUMI_Visium) <- colnames(Visium.obj)
Visium_data  <- SpatialRNA(counts = Visium_count, nUMI = nUMI_Visium, use_fake_coords = T)

# Run deconvolution
myRCTD <- create.RCTD(Visium_data, reference, max_cores = 10, CELL_MIN_INSTANCE = 20)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

saveRDS(myRCTD, 'myRCTD_Visium_full.rds')

weights <- myRCTD@results$weights
norm_weights <- normalize_weights(weights)

colnames(norm_weights) <- paste0(rep("RCTDpredict_c", ncol(norm_weights)), colnames(norm_weights))
Visium.obj <- AddMetaData(Visium.obj, metadata = norm_weights)
SpatialFeaturePlot(Visium.obj, features = c("RCTDpredict_c9", alpha = c(0.1, 1)))


library(pheatmap)
visium_cluster = Visium.obj@meta.data[rownames(norm_weights), "seurat_clusters"]
cluste_inf <- as.data.frame(cbind(Visium_cluster = paste0(rep("C", n = length(visium_cluster)), visium_cluster)))
rownames(cluste_inf) <- rownames(norm_weights)
cluster_order <- order((cluste_inf$Visium_cluster))

p < -pheatmap(norm_weights[cluster_order,], show_rownames = F, show_colnames = T, annotation_row = cluste_inf, cluster_row = T)
filePath = "./Deconv_heatmap.pdf"
ggsave(file = filePath, plot = p, dpi = 100, width = 10, height = 8)
write.table(norm_weights, file="RCTD_semibulk.txt", quote = F, col.names = NA, sep = "\t")


library(Polychrome)
source("./VIS_topic.R")
theta = norm_weights
Pos = Visium.obj@images$slice1@coordinates[, c(4,5)]
new_coord <- t(apply(Pos, 1, Fun_lotate, -90))

colnames(new_coord) = c("x", "y")
topicCol <- DiscretePalette(ncol(theta), palette = "polychrome", shuffle = FALSE)

p <- vizAllTopics(theta, new_coord, topicOrder = seq(ncol(theta)), topicCols = topicCol, showLegend = TRUE, r = 60)
filePath = "./Decon_VisTopic.pdf"
ggsave(file = filePath, plot = p, dpi = 100, width = 12, height = 8)

p1_lists <- list()

for (clus in colnames(norm_weights)) {
	p1 < -SpatialFeaturePlot(Visium.obj, features = c(clus, alpha = c(0.1, 1)))	
	p1_lists <- c(p1_lists, list(p1))
}

filePath <- paste0("./Decov_pct.pdf")
ggsave(file = filePath, plot = wrap_plots(p1_lists), dpi = 100, width = 20, height = 16)
