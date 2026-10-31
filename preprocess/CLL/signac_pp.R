# Signac processing of CLL patient 1, CD19+ (CLL_PT1_CD19pos): the cellranger-atac peak-by-cell counts of the cells
# kept by get_filtering_num.R, peaks observed in at least 10 cells and cells with at least 200 detected peaks, then the
# top 25% of peaks by count (FindTopFeatures, q75) for TF-IDF, LSI, UMAP and clustering (the LSI clusters of
# Fig. 4a), and the gene activity.  Writes mtscATAC-seq/h5_files/PT1_pos_pp.h5seurat and its h5ad (counts of all
# kept peaks in raw.X).
#
# The study ran these lines for both patients from a script whose saved copies differ in the FindTopFeatures cutoff;
# the output used for the study has the 26,884 variable peaks of q75, as here.  Only patient 1 is used.
# Run from the repository root: Rscript preprocess/CLL/signac_pp.R
library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
set.seed(1234)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)

CR <- 'mtscATAC-seq/mtscATACpaper_large_data_files/source/cellranger_output'
OUT <- 'mtscATAC-seq/h5_files'
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

features <- readr::read_tsv(file.path(CR, 'CLL_PT1_CD19pos_peak_bc_matrix/peaks.bed'), col_names = F) %>% tidyr::unite(feature, c(X1,X2), sep=':')
features <- features %>% tidyr::unite(feature, c(feature, X3), sep='-')

barcodes <- readr::read_tsv(file.path(CR, 'CLL_PT1_CD19pos_peak_bc_matrix/barcodes.tsv'), col_names = F) %>% tidyr::unite(barcodes)

mtx <- Matrix::readMM(file.path(CR, 'CLL_PT1_CD19pos_peak_bc_matrix/matrix.mtx')) %>%
  magrittr::set_rownames(features$feature) %>%
  magrittr::set_colnames(barcodes$barcodes)

# cells kept by get_filtering_num.R (0-based positions)
del_PT1pos <- readr::read_tsv('preprocess/CLL/row_del/CLL_PT1_CD19pos_del.tsv')
del_PT1pos <- del_PT1pos$index +1

mtx <- mtx[,del_PT1pos]

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevelsStyle(annotations) <- 'UCSC'

chrom_assay <- CreateChromatinAssay(
  counts = mtx,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = file.path(CR, 'CLL_PT1_CD19pos_v12-mtMask_fragments.tsv.gz'),
  min.cells = 10,
  min.features = 200
  )

PT1_pos <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
  )

PT1_pos

Annotation(PT1_pos) <- annotations

PT1_pos <- FindTopFeatures(PT1_pos, min.cutoff = 'q75')

PT1_pos <- RunTFIDF(PT1_pos)
PT1_pos <- RunSVD(PT1_pos)
PT1_pos <- RunUMAP(object = PT1_pos, reduction = 'lsi', dims = 2:30)
PT1_pos <- FindNeighbors(object = PT1_pos, reduction = 'lsi', dims = 2:30)
PT1_pos <- FindClusters(object = PT1_pos, verbose = FALSE, algorithm = 3)

gene.activities <- GeneActivity(PT1_pos)

PT1_pos[['RNA']] <- CreateAssayObject(counts = gene.activities)

PT1_pos <- NormalizeData(
  object=PT1_pos,
  assay='RNA',
  normalization.method='LogNormalize',
  scale.factor=median(PT1_pos$nCount_RNA)
  )

DefaultAssay(PT1_pos) <- 'peaks'

SaveH5Seurat(PT1_pos, filename = file.path(OUT, 'PT1_pos_pp.h5seurat'), overwrite = TRUE)
Convert(file.path(OUT, 'PT1_pos_pp.h5seurat'), dest = 'h5ad', overwrite = TRUE)
