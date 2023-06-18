library(reticulate)
path_to_python ="/path/to/python"
use_python(path_to_python)
py_config()

library(Giotto)
library(patchwork)
library(ggplot2)
library(dplyr)

set.seed(1234)
source("./SUB_function.R")

SRanger <-  commandArgs(trailingOnly=TRUE)[1]  # Space Ranger result
Dir     <-  ommandArgs(trailingOnly=TRUE)[2]   # output Dir name 
Pname   <-  commandArgs(trailingOnly=TRUE)[3]  # Sample name

xmax    <- as.numeric(commandArgs(trailingOnly=TRUE)[4])
xmin    <- asas.numeric(commandArgs(trailingOnly=TRUE)[5])
ymax    <- asas.numeric(commandArgs(trailingOnly=TRUE)[6])
ymin    <- asas.numeric(commandArgs(trailingOnly=TRUE)[7])


temp_dir = Dir
instrs = createGiottoInstructions(save_dir = temp_dir, save_plot = TRUE, show_plot = FALSE, python_path = path_to_python)
results_folder <- paste0(SRanger)

GiottoOBJ = createGiottoVisiumObject(visium_dir = results_folder, expr_data = 'raw',
                                         png_name = 'tissue_lowres_image.png',
                                         gene_column_index = 2, instructions = instrs)

GiottoOBJ = updateGiottoImage(GiottoOBJ, image_name = 'image',
                                  xmax_adj = xmax, xmin_adj = xmin,
                                  ymax_adj = ymax, ymin_adj = ymin)

spatPlot(gobject = GiottoOBJ, cell_color = 'in_tissue', show_image = T, point_alpha = 0.5,
         save_param = list(save_name = 'spatplot_image'))


## Subset on spots that were covered by the tissue
metadata = pDataDT(GiottoOBJ)
in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
GiottoOBJ = subsetGiotto(GiottoOBJ, cell_ids = in_tissue_barcodes)

## Spot filter
GiottoOBJ <- filterGiotto(gobject = GiottoOBJ,
                              expression_threshold = 1,
                              gene_det_in_min_cells = 1,
                              min_det_genes_per_cell = 10,
                              expression_values = c('raw'),
                              verbose = T)


## Normalization
GiottoOBJ <- normalizeGiotto(gobject = GiottoOBJ, scalefactor = 10000, verbose = T)


## Add gene & cell statistics
GiottoOBJ <- addStatistics(gobject = GiottoOBJ)
gene_metadata = fDataDT(GiottoOBJ)

# Read Marker data
cellMaker <- read.table("./Cancer_TypeMarker.txt", sep="\t", header= T)
cellMaker[!cellMaker$Gene %in% gene_metadata$gene_ID,]
cellMaker <- cellMaker[cellMaker$Gene %in% gene_metadata$gene_ID,]

sign_lsit <- list()
for(i in unique(cellMaker$Celltype) ) { sign_lsit <- c( sign_lsit , list(cellMaker[cellMaker$Celltype == i ,1 ] ))    }

# PAGE
signature_matrix = makeSignMatrixPAGE(sign_names = unique(cellMaker$Celltype), sign_list =sign_lsit )
GiottoOBJ = runPAGEEnrich(gobject = GiottoOBJ, sign_matrix = signature_matrix, min_overlap_genes =2)

File <- paste0(Dir, "/", Pname, "_PAGEscore_Type_Marker.txt")
write.table(as.data.frame(GiottoOBJ@spatial_enrichment$PAGE), file=File, sep="\t", col.names=NA, row.names=T, quote=F)

typeCol <- read.table("./Cancer_Type_Col.txt", sep="\t", header= T)
col_names <- typeCol$col
names(col_names) <- typeCol$CellType

Annot <- as.data.frame(GiottoOBJ@spatial_enrichment$PAGE)
Annot$Annotation <- apply(Annot, 1, function(x){names(x)[which.max(x)]})
GiottoOBJ@cell_metadata$CellType <- Annot$Annotation
GiottoOBJ@cell_metadata$sample <- Pname

# plot scores with custom function
p1 <- MyspatPlot2D_single(gobject = GiottoOBJ, cell_color_gradient=c("white","blue"), spat_enr_names = 'PAGE',
     gradient_limits = c(-4,5), gradient_midpoint = 0, point_border_col = "gray",
     return_plot = T, cell_color  = cell_types[1], color_as_factor = F , point_size = 3)
p2 <- MyspatPlot2D_single(gobject = GiottoOBJ ,  cell_color_gradient=c("white","red"), spat_enr_names = 'PAGE',
     gradient_limits = c(-4,5), gradient_midpoint = 0 , point_border_col ="gray",
     return_plot =T, cell_color = cell_types[2], color_as_factor = F, point_size = 3)
p3 <- MyspatPlot2D_single(gobject = GiottoOBJ,  cell_color_gradient = c("white","#005731"), spat_enr_names = 'PAGE',
     gradient_limits = c(-4,5), gradient_midpoint = 0 , point_border_col ="gray",
     return_plot = T, cell_color = cell_types[3], color_as_factor = F, point_size = 3)
p <- p1 + p2 + p3

filePath <- paste0(Dir, "/Mode_Page_result_", Pname, ".pdf")
ggsave(file = filePath, plot = p, dpi = 100, width = 15, height = 5)
