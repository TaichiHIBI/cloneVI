library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
set.seed(1234)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)


#PT1_pos
features <- readr::read_tsv('mtscATAC-seq/mtscATACpaper_large_data_files/source/cellranger_output/CLL_PT1_CD19pos_peak_bc_matrix/peaks.bed', col_names = F) %>% tidyr::unite(feature, c(X1,X2), sep=':')
features <- features %>% tidyr::unite(feature, c(feature, X3), sep='-')

barcodes <- readr::read_tsv('mtscATAC-seq/mtscATACpaper_large_data_files/source/cellranger_output/CLL_PT1_CD19pos_peak_bc_matrix/barcodes.tsv', col_names = F) %>% tidyr::unite(barcodes)

mtx <- Matrix::readMM('mtscATAC-seq/mtscATACpaper_large_data_files/source/cellranger_output/CLL_PT1_CD19pos_peak_bc_matrix/matrix.mtx') %>%
  magrittr::set_rownames(features$feature) %>%
  magrittr::set_colnames(barcodes$barcodes)

del_PT1pos <- readr::read_tsv('mtscATAC-seq/row_del/CLL_PT1_CD19pos_del.tsv')
del_PT1pos <- del_PT1pos$index +1

mtx <- mtx[,del_PT1pos]

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevelsStyle(annotations) <- 'UCSC'

chrom_assay <- CreateChromatinAssay(
  counts = mtx,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = 'mtscATAC-seq/mtscATACpaper_large_data_files/source/cellranger_output/CLL_PT1_CD19pos_v12-mtMask_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
  )

PT1_pos <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
  )

PT1_pos

PT1_pos[['peaks']]

Annotation(PT1_pos) <- annotations

PT1_pos <- FindTopFeatures(PT1_pos, min.cutoff = 'q75')

PT1_pos <- RunTFIDF(PT1_pos)
PT1_pos <- RunSVD(PT1_pos)
DepthCor(PT1_pos)
PT1_pos <- RunUMAP(object = PT1_pos, reduction = 'lsi', dims = 2:30)
PT1_pos <- FindNeighbors(object = PT1_pos, reduction = 'lsi', dims = 2:30)
PT1_pos <- FindClusters(object = PT1_pos, verbose = FALSE, algorithm = 3)
DimPlot(object = PT1_pos, label = TRUE) + NoLegend()


gene.activities <- GeneActivity(PT1_pos)

PT1_pos[['RNA']] <- CreateAssayObject(counts = gene.activities)

PT1_pos <- NormalizeData(
  object=PT1_pos,
  assay='RNA',
  normalization.method='LogNormalize',
  scale.factor=median(PT1_pos$nCount_RNA)
  )

DefaultAssay(pbmc) <- 'peaks'



SaveH5Seurat(PT1_pos, filename='mtscATAC-seq/h5_files/PT1_pos_pp.h5seurat')
Convert('mtscATAC-seq/h5_files/PT1_pos_pp.h5seurat', dest='h5ad')


#PT2_pos
features <- readr::read_tsv('mtscATAC-seq/mtscATACpaper_large_data_files/source/cellranger_output/CLL_PT2_CD19pos_peak_bc_matrix/peaks.bed', col_names = F) %>% tidyr::unite(feature, c(X1,X2), sep=':')
features <- features %>% tidyr::unite(feature, c(feature, X3), sep='-')

barcodes <- readr::read_tsv('mtscATAC-seq/mtscATACpaper_large_data_files/source/cellranger_output/CLL_PT2_CD19pos_peak_bc_matrix/barcodes.tsv', col_names = F) %>% tidyr::unite(barcodes)

mtx <- Matrix::readMM('mtscATAC-seq/mtscATACpaper_large_data_files/source/cellranger_output/CLL_PT2_CD19pos_peak_bc_matrix/matrix.mtx') %>%
  magrittr::set_rownames(features$feature) %>%
  magrittr::set_colnames(barcodes$barcodes)

del_PT2pos <- readr::read_tsv('mtscATAC-seq/row_del/CLL_PT2_CD19pos_del.tsv')
del_PT2pos <- del_PT2pos$index +1

mtx <- mtx[,del_PT2pos]

chrom_assay <- CreateChromatinAssay(
  counts = mtx,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = 'mtscATAC-seq/mtscATACpaper_large_data_files/source/cellranger_output/CLL_PT2_CD19pos_v12-mtMask_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

PT2_pos <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)

Annotation(PT2_pos) <- annotations

PT2_pos <- RunTFIDF(PT2_pos)
PT2_pos <- FindTopFeatures(PT2_pos, min.cutoff = 'q0')
PT2_pos <- RunSVD(PT2_pos)
DepthCor(PT2_pos)
PT2_pos <- RunUMAP(object = PT2_pos, reduction = 'lsi', dims = 2:30)
PT2_pos <- FindNeighbors(object = PT2_pos, reduction = 'lsi', dims = 2:30)
PT2_pos <- FindClusters(object = PT2_pos, verbose = FALSE, algorithm = 3)
DimPlot(object = PT2_pos, label = TRUE) + NoLegend()

PT2_pos <- FindTopFeatures(PT2_pos, min.cutoff = 'q75')

SaveH5Seurat(PT2_pos, filename='mtscATAC-seq/h5_files/PT2_pos_pp.h5seurat')
Convert('mtscATAC-seq/h5_files/PT2_pos_pp.h5seurat', dest='h5ad')
