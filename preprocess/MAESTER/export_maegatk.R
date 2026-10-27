# Writes the maegatk allele counts of the MAESTER sample BPDCN712 (GEO GSM5534706) as .mtx files (cells in columns),
# with the cells (barcodes.tsv) and the positions (genes.tsv), in MAESTER/rds_output/, for build_model_inputs.ipynb.
# The GEO file, once decompressed (gunzip; the .rds inside is compressed again, as R saves it), is byte-identical to
# the maegatk object used for the study (TenX_BPDCN712_All_mr3_maegatk.rds of the source study).  These are the export lines of the study (output.R) restricted to the
# allele counts and the coverage; the base qualities are not used.
# Run from the repository root: Rscript preprocess/MAESTER/export_maegatk.R
library(SummarizedExperiment)
library(Matrix)

rna_mgv <- readRDS('MAESTER/GSM5534706_BPDCN712_10x_MAESTER_maegatk.rds')     # the GEO file after gunzip
OUT <- 'MAESTER/rds_output'
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

for (n in c(paste0(rep(c('A', 'C', 'G', 'T'), each = 2), '_counts_', c('fw', 'rev')), 'coverage')) {
  writeMM(assay(rna_mgv, n), file.path(OUT, paste0(n, '.mtx')))
}
write(colnames(rna_mgv), file.path(OUT, 'barcodes.tsv'))
write.table(rowRanges(rna_mgv), file.path(OUT, 'genes.tsv'), sep = '\t')
cat(ncol(rna_mgv), 'cells,', nrow(rna_mgv), 'positions\n')
