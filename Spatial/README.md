# Spatial transcriptome analysis
## Analysis of Visium data  
- [Visium_Seurat.R](./Visium_Seurat.R): basic analysis (dimentional reduction, clustering, etc.) by Seurat (R 4.2.1; Seurat 4.3.0).  
- [PAGE.R](./PAGE.R): PAGE analysis by Giotto for scoring Visium spots using the signature genes.  
  - [SUB_function.R](./SUB_function.R): coloring PAGE plots.  
  - [Cancer_TypeMarker.txt](./Cancer_TypeMarker.txt): the list of signature genes (Supplementary Table S6; Haga, Sakamoto, Kajiya et al.).  
  - [Cancer_Type_Col.txt](./Cancer_Type_Col.txt): the color list for each stage  
- [CellChat.R](./CellChat.R): ligand-receptor interaction analysis by CellChat.

## Analysis of Xenium data
- [Xenium_Seurat.R](./Xenium_Seurat.R): basic analysis (dimentional reduction, clustering, etc.) by Seurat (R 4.2.1; Seurat 4.3.0).

## Decovolution analysis of Visium data
- [public_scRNA_ana.R](./public_scRNA_ana.R): custom annotation of scRNA-seq data (datasets provided from Zhu et al. 2022 _Exp Mol Med_ PMID: 36434043).
- [Deconv_byScRNA.R](./Deconv_byScRNA.R): deconvolution of Visium data with scRNA-seq reference datasets by RCTD.
- [Deconv_byXenium.R](./Deconv_byXenium.R): deconvolution of Visium data with Xenium reference data by RCTD.
  - [VIS_topic.R](./VIS_topic.R): visualization of results of Visium deconvolution analysis.
