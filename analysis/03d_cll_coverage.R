# CLL mtscATAC-seq (4/4): accessibility at the genes of each clone in cells with high and low clonal probability.
#
# For every clone k, the cells are split into the "high" group (the cells with the largest pi that together account for
# 50% of the clone's probability, 03a_cll.ipynb) and the "low" group (the other cells).  For the three first-ranked genes
# of the clone (03c_cll_gene_activity.ipynb), the accessibility of the two groups is drawn with CoveragePlot (40 kb
# upstream, 20 kb downstream).  The adjusted p-value in the title is from the Wilcoxon rank-sum test of gene activity
# between the two groups (FindMarkers, Bonferroni over all genes); it is descriptive only, because the gene was selected
# by the ranking and the groups come from the same data (Methods).  The test uses the gene-activity counts as they are
# (the ACTIVITY assay is not normalized, as in the original analysis).
#   Fig. 4d   PAX9, Clone 4, the gene body only and without the p-value (coverage_fig4d_PAX9_clone4.pdf)
#   Fig. S4d  the gene of rank 2 of Clone 4: HIST1H2AE (coverage_clone4_2_HIST1H2AE.pdf)
#
# No region is highlighted (earlier versions passed the gene name as region.highlight, which Signac 1.4 ignores with a
# warning, so the plots are unchanged).
#
# Input:  analysis/results/cll/r_input/ (03a), r_output/gene_activity/ (03b), clone_top_genes.csv (03c), fragment file
# Output: analysis/results/cll/r_output/coverage/ (plots, FindMarkers tables, coverage_genes_pvalues.csv)
# Environment: conda env Signac (R 4.1.1, Signac 1.4.0, Seurat 4.0.5)
# Run from the repository root: Rscript analysis/03d_cll_coverage.R

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(EnsDb.Hsapiens.v75)
  library(GenomeInfoDb)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})
set.seed(1234)

R_IN <- 'analysis/results/cll/r_input'
GA <- 'analysis/results/cll/r_output/gene_activity'
OUT <- 'analysis/results/cll/r_output/coverage'
FRAGMENTS <- 'data/ATAC/flagment/CLL_PT1_CD19pos_v12-mtMask_fragments.tsv.gz'
N_TOP <- 3                 # genes per clone
W <- 6; H <- 3             # inches; CoveragePlot fixes the size of the gene labels, so the plots are drawn larger than
                           # print size and scaled down by about 0.55 when the figure is assembled (text about 6 pt)
AXIS <- theme(axis.title = element_text(size = 8),   # track titles, so that the long y title fits the coverage tracks
              plot.margin = margin(5.5, 20, 5.5, 5.5))   # room at the right for the half of the last position label
                                                       # that hangs beyond the axis (cut at the page edge otherwise)
W4D <- 6.5; H4D <- 2.6     # Fig. 4d: the aspect of its place in Fig. 4 (two columns wide, 2.5 : 1); not lower than
                           # 2.6 in, where the y title no longer fits the coverage tracks
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# the cells and their peak counts, as in 03b
features <- read.table(file.path(R_IN, 'features.tsv'))
barcodes <- read.table(file.path(R_IN, 'barcodes.tsv'))
mtx <- Matrix::readMM(file.path(R_IN, 'matrix.mtx'))
rownames(mtx) <- rownames(features)
colnames(mtx) <- barcodes$barcode
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevelsStyle(annotations) <- 'UCSC'
chrom_assay <- CreateChromatinAssay(counts = mtx, sep = c(':', '-'), genome = 'hg19', fragments = FRAGMENTS,
                                    min.cells = 10, min.features = 200)
cll <- CreateSeuratObject(counts = chrom_assay, assay = 'peaks')
Annotation(cll) <- annotations
DefaultAssay(cll) <- 'peaks'

# gene activity from 03b
ga <- Matrix::readMM(file.path(GA, 'matrix.mtx'))
rownames(ga) <- read.table(file.path(GA, 'genes.tsv'), header = TRUE)$gene
colnames(ga) <- read.table(file.path(GA, 'barcodes.tsv'), header = TRUE)$barcode
stopifnot(identical(colnames(ga), colnames(cll)))
cll[['ACTIVITY']] <- CreateAssayObject(counts = as(ga, 'dgCMatrix'))

# high / low groups (column clone_{k+1} = Clone k) and the ranked genes of each clone
flags <- read.table(file.path(R_IN, 'clone_flags.tsv'))
stopifnot(nrow(flags) == ncol(cll))
top <- read.csv('analysis/results/cll/clone_top_genes.csv', check.names = FALSE, colClasses = 'character')

set_groups <- function(k) {
  grp <- factor(ifelse(flags[[k + 1]] == 1, 'high', 'low'), levels = c('low', 'high'))
  names(grp) <- colnames(cll)
  Idents(cll) <- grp
  cll
}

# ONLY_FIG4D=1 draws Fig. 4d alone (the other outputs of a full run are kept)
ONLY_FIG4D <- Sys.getenv('ONLY_FIG4D') == '1'
pvals <- list()
for (k in if (ONLY_FIG4D) integer(0) else 0:(ncol(flags) - 1)) {
  cll <- set_groups(k)
  deg <- FindMarkers(cll, ident.1 = 'high', min.pct = 0.1, logfc.threshold = 0.25, assay = 'ACTIVITY')
  write.csv(deg, file.path(OUT, sprintf('findmarkers_clone%d.csv', k)))
  genes <- top[[sprintf('Clone_%d', k)]]
  genes <- genes[!is.na(genes) & genes != '']
  cat(sprintf('Clone %d: %d high cells, %d genes ranked\n', k, sum(flags[[k + 1]]), length(genes)))
  for (r in seq_len(min(N_TOP, length(genes)))) {
    g <- genes[r]
    p <- if (g %in% rownames(deg)) deg[g, 'p_val_adj'] else NA   # NA: not tested (below min.pct or logfc.threshold)
    cat(sprintf('  rank %d: %s  adjusted p = %s\n', r, g, format(p, digits = 15)))
    pvals[[length(pvals) + 1]] <- data.frame(clone = k, rank = r, gene = g, p_val_adj = p)
    title <- sprintf('Clone %d, %s (adjusted p = %s)', k, g, if (is.na(p)) 'not tested' else format(p, digits = 3))
    plt <- CoveragePlot(cll, region = g, extend.upstream = 40000, extend.downstream = 20000) +
      plot_annotation(title = title)
    ggsave(file.path(OUT, sprintf('coverage_clone%d_%d_%s.pdf', k, r, g)), plt & AXIS, width = W, height = H)
  }
}
if (!ONLY_FIG4D) write.csv(do.call(rbind, pvals), file.path(OUT, 'coverage_genes_pvalues.csv'), row.names = FALSE)

# Fig. 4d: PAX9, Clone 4, the gene body only (22,148 bp; extended by 2 bp at the 3' end because Signac 1.4 samples 10%
# of the positions and current dplyr rejects a non-integer sample size), without the p-value
cll <- set_groups(4)
plt <- CoveragePlot(cll, region = 'PAX9', extend.upstream = 0, extend.downstream = 2)
ggsave(file.path(OUT, 'coverage_fig4d_PAX9_clone4.pdf'), plt & AXIS, width = W4D, height = H4D)
cat('wrote', OUT, '\n')
