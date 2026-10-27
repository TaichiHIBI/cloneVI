# Writes the inputs of build_model_inputs.ipynb that come from R objects of the source study:
# - the mgatk allele counts of CLL_PT1_CD19pos as .mtx files (cells in columns), with the cells (barcodes.tsv) and
#   the positions (genes.tsv), in mtscATAC-seq/rds_output/PT1_pos_mg/;
# - the mtDNA clones of patient 1 defined by the source study (figure_CLL/output/PT1_clone_definition.rds of its
#   reproducibility code) as mtscATAC-seq/rds_output/clone_definition/PT1.tsv.
# The export used for the study was not kept; these lines follow the MAESTER export (preprocess/MAESTER/export_maegatk.R)
# and reproduce the files read by the study's notebooks (build_model_inputs.ipynb compares the results).
# Run from the repository root: Rscript preprocess/CLL/export_mgatk.R
library(SummarizedExperiment)
library(Matrix)

L <- 'mtscATAC-seq/mtscATACpaper_large_data_files/source/mgatk_output'
OUT <- 'mtscATAC-seq/rds_output'
dir.create(file.path(OUT, 'PT1_pos_mg'), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, 'clone_definition'), recursive = TRUE, showWarnings = FALSE)

mg <- readRDS(file.path(L, 'CLL_PT1_CD19pos_v12-mtMask_mgatk.rds'))
# the allele counts on each strand and the coverage (the base qualities are not used)
for (n in c(paste0(rep(c('A', 'C', 'G', 'T'), each = 2), '_counts_', c('fw', 'rev')), 'coverage')) {
  writeMM(assay(mg, n), file.path(OUT, 'PT1_pos_mg', paste0(n, '.mtx')))
}
write(colnames(mg), file.path(OUT, 'PT1_pos_mg', 'barcodes.tsv'))
write.table(rowRanges(mg), file.path(OUT, 'PT1_pos_mg', 'genes.tsv'), sep = '\t')

clones <- readRDS('mtscATAC-seq/mtscATACpaper_reproducibility/figure_CLL/output/PT1_clone_definition.rds')
write.table(clones, file.path(OUT, 'clone_definition', 'PT1.tsv'), sep = '\t')
cat(ncol(mg), 'cells,', nrow(mg), 'positions;', nrow(clones), 'cells with a clone\n')
