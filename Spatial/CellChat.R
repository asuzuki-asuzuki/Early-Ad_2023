library(dplyr)
library(tidyverse)
library(pheatmap)
library(CellChat)
library(Seurat)
library(ggplot2)
library(patchwork)

library(future)
plan("multisession", workers = 10)
options('future.globals.maxSize' = 800000*1024^2)

set.seed(1234)

Fun_lotate <- function( loc , do){
    c = do *pi/180
    x = loc[1]
    y = loc[2]
    x2 = x*cos(c) - y * sin(c)
    y2 = x*sin(c) + y * cos(c)
    return(c(x2,y2))
}


# Read Visium Seurat obj
Visium.obj <- readRDS("lung_visium.rds")

# Modify cluster names
levels_data <- paste0(rep("Cluster", n = length(unique(Visium.obj@meta.data$ClusterID))), c(0:(length(unique(Visium.obj@meta.data$ClusterID))-1)))

Visium.obj@meta.data$ClusterID <- paste0(rep("Cluster", n = ncol(Visium.obj)), Visium.obj@meta.data$seurat_clusters) %>% factor(levels = levels_data)

Idents(Visium.obj)<- "ClusterID"
data.input = GetAssayData(Visium.obj, slot = "data", assay = "SCT") # normalized data matrix
meta = data.frame(labels = Idents(Visium.obj), row.names = names(Idents(Visium.obj)))


spatial.locs = GetTissueCoordinates(Visium.obj, scale = NULL, cols = c("imagerow", "imagecol")) 
scale.factors = jsonlite::fromJSON(txt = 'scalefactors_json.json')
scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres, 
                fiducial = scale.factors$fiducial_diameter_fullres, 
		hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef)

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                datatype = "spatial", coordinates = spatial.locs, scale.factors = scale.factors)

CellChatDB <- CellChatDB.human # Use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                               distance.use = TRUE, interaction.length = 200, scale.distance = 0.01)
cellchat <- computeCommunProbPathway(cellchat)

saveRDS(cellchat,"cellchat_computeCommunProb_2.rds")

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

pdf("aggregated_cell-cell_communication_network.pdf", height=6, width=6 )
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

write.table(cellchat@net$count, file="Intaraction_Count.txt", quote=F, col.names=NA, sep="\t")
write.table(cellchat@net$weight, file="Intaraction_weight.txt", quote=F, col.names=NA, sep="\t")
write.table(subsetCommunication(cellchat), file="Intaraction_LR_network.txt", quote=F, col.names=NA, sep="\t")
write.table(subsetCommunication(cellchat, slot.name = "netP"), file="Intaraction_pathway.txt", quote=F, col.names=NA, sep="\t")
pw.df <- subsetCommunication(cellchat, slot.name = "netP")

pw.df2 <- pw.df %>% mutate(Pair = paste( !!!rlang::syms(c("source","target")), sep="-")) %>%
	select(c("Pair", "pathway_name", "prob")) %>%
	pivot_wider(names_from = pathway_name, values_from = prob, values_fill = list(value = 0)) %>% as.data.frame()
pw.df2[is.na(pw.df2)] <- 0
write.table(pw.df2, file="Intaraction_pathway_Table.txt", quote=F, col.names=NA, sep="\t")


pdf("each_cluser_network_weight.pdf", height=15, width=15)
mat <- cellchat@net$weight
par(mfrow = c(4,3), mar = c(1, 2, 1, 1), oma = c(0, 0, 0, 0))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i], margin = 0)
}
dev.off()

pathways.show <- cellchat@netP$pathways

pdf("Pathway_network_circle.pdf", height=7, width=7)
for (i in pathways.show) {
  print(i)
  par(mfrow = c(1,1), xpd=TRUE)
  netVisual_aggregate(cellchat, signaling = i, layout = "circle")
}
dev.off()

pdf("Pathway_network_chord.pdf", height=7, width=7)
for (i in pathways.show) {
  print(i)
  par(mfrow = c(1,1), xpd=TRUE)
  netVisual_aggregate(cellchat, signaling = i, layout = "chord")
}
dev.off()

pdf("Pathway_network_Visium.pdf", height=7, width=7)
for (i in pathways.show) {
  print(i)
  par(mfrow=c(1,1))
  p <- netVisual_aggregate(cellchat, signaling = i, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)
  print(p)
}
dev.off()

pdf("Pathway_network_heatmap.pdf", height=7, width=7)
for (i in pathways.show) {
  print(i)
  par(mfrow = c(1,1), xpd=TRUE)
  p <- netVisual_heatmap(cellchat, signaling = c(i), color.heatmap = "Reds")
  print(p)
}
dev.off()

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
pdf("Pathway_network_Centrality.pdf", height=10, width=10)
for (i in pathways.show[1:82]) {
	par(mfrow=c(1,1))
	netAnalysis_signalingRole_network(cellchat, signaling = i, width = 8, height = 2.5, font.size = 10)
}
dev.off()

pdf("Pathway_network_Visium2.pdf", height=7, width=7)
for (i in pathways.show[1:82]) {
  print(i)
  par(mfrow=c(1,1))
  p <- netVisual_aggregate(cellchat, signaling = i, layout = "spatial", edge.width.max = 2, alpha.image = 0.2, vertex.weight = "outgoing", vertex.size.max = 3, vertex.label.cex = 3.5)
  print(p)
}
dev.off()


pdf("bubble.pdf", height=20, width=15)
for (i in pathways.show) {
	par(mfrow = c(1,1) , xpd=TRUE )
	netVisual_bubble(cellchat, sources.use = levels(cellchat@idents),  targets.use = levels(cellchat@idents), remove.isolate = FALSE, font.size = 5)
}
dev.off()

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways

for (i in 1:length(pathways.show.all)) {
  gg <- netAnalysis_contribution( cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0("LR/",pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 6, height = 12, units = 'in', dpi = 300)
}

