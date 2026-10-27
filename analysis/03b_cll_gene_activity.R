# CLL mtscATAC-seq (2/4): gene activity from the accessibility.
#
# Gene activity of each gene in each cell = the number of fragments in its gene body and 2 kb upstream (GeneActivity of
# Signac, protein-coding genes of Ensembl v75). Genes with fragments in more than 50 cells are kept.
#
# Input:  analysis/results/cll/r_input/ (peak counts, written by 03a_cll.ipynb) and the fragment file of the source study
#         (data/ATAC/flagment/CLL_PT1_CD19pos_v12-mtMask_fragments.tsv.gz)
# Output: analysis/results/cll/r_output/gene_activity/ (genes x cells, Matrix Market), read by 03c_cll_gene_activity.ipynb
# Environment: conda env Signac (R 4.1.1, Signac 1.4.0, Seurat 4.0.5)
# Run from the repository root: Rscript analysis/03b_cll_gene_activity.R

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(EnsDb.Hsapiens.v75)
  library(GenomeInfoDb)
  library(Matrix)
})
set.seed(1234)

R_IN <- 'analysis/results/cll/r_input'
OUT <- 'analysis/results/cll/r_output/gene_activity'
FRAGMENTS <- 'data/ATAC/flagment/CLL_PT1_CD19pos_v12-mtMask_fragments.tsv.gz'
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# peak counts of the cells (peaks x cells)
features <- read.table(file.path(R_IN, 'features.tsv'))
barcodes <- read.table(file.path(R_IN, 'barcodes.tsv'))
mtx <- Matrix::readMM(file.path(R_IN, 'matrix.mtx'))
rownames(mtx) <- rownames(features)
colnames(mtx) <- barcodes$barcode
cat('peak counts:', nrow(mtx), 'peaks x', ncol(mtx), 'cells\n')

# chromatin assay with the gene annotation (hg19, UCSC chromosome names)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevelsStyle(annotations) <- 'UCSC'
chrom_assay <- CreateChromatinAssay(counts = mtx, sep = c(':', '-'), genome = 'hg19', fragments = FRAGMENTS,
                                    min.cells = 10, min.features = 200)
cll <- CreateSeuratObject(counts = chrom_assay, assay = 'peaks')
Annotation(cll) <- annotations
cat('after the filters of CreateChromatinAssay (peaks in >= 10 cells, cells with >= 200 peaks):',
    nrow(cll), 'peaks x', ncol(cll), 'cells\n')

# gene activity
ga <- GeneActivity(cll)
cat('gene activity:', nrow(ga), 'genes x', ncol(ga), 'cells\n')
ga <- ga[rowSums(ga > 0) > 50, ]
cat('genes with fragments in more than 50 cells:', nrow(ga), '\n')

Matrix::writeMM(ga, file.path(OUT, 'matrix.mtx'))
write.table(data.frame(gene = rownames(ga)), file.path(OUT, 'genes.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(data.frame(barcode = colnames(ga)), file.path(OUT, 'barcodes.tsv'), sep = '\t', quote = FALSE,
            row.names = FALSE)
cat('wrote', OUT, '\n')
