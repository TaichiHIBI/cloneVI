# Ovarian tumor cohort (4/4): the mutation-based motif test of 04d_ovarian_motif.ipynb against Seurat.
#
# Methods: the weighted Wilcoxon test with weights of only 0 and 1 gives the p-values of Seurat's FindMarkers.  This runs
# FindMarkers on the motif scores and mutant calls of the cells tested for the clone of 12879T>C (Clone 11 of this fit, number read from r_input/clone.txt; Signac vignette settings:
# assay "chromvar", mean.fxn = rowMeans, fc.name = "avg_diff"; min.pct 0.1 and no fold-change filter, as in the paper)
# and compares it with the Python implementation: the same tested and significant motifs, and the same p-values,
# effect sizes and detection fractions.
#
# Input:  analysis/results/ovarian/motif/r_input/ (z.csv: motifs x cells, groups.csv) and clone<k>_motifs.csv (04d; k = 11 for this fit)
# Output: analysis/results/ovarian/motif/clone<k>_seurat_findmarkers.csv, clone<k>_seurat_check.csv
# Environment: conda env Signac (Seurat 4.0.5; limma installed, used by FindMarkers for the Wilcoxon test)
# Run from the repository root: Rscript analysis/04e_ovarian_motif_seurat_check.R

suppressMessages({ library(Seurat); library(Matrix) })
options(Seurat.limma.wilcox.msg = FALSE)
DIR <- 'analysis/results/ovarian/motif'
K <- readLines(file.path(DIR, 'r_input', 'clone.txt'))       # the clone of 12879T>C, written by 04d

z <- as.matrix(read.csv(file.path(DIR, 'r_input', 'z.csv'), row.names = 1, check.names = FALSE))      # motifs x cells
groups <- read.csv(file.path(DIR, 'r_input', 'groups.csv'), row.names = 1, check.names = FALSE)
stopifnot(identical(rownames(groups), colnames(z)))
obj <- CreateSeuratObject(counts = CreateAssayObject(data = z), assay = 'chromvar')
Idents(obj) <- factor(groups[colnames(obj), 'genotype'], levels = c('pos', 'neg'))
cat('cells:', ncol(obj), '(mutant', sum(Idents(obj) == 'pos'), '), motifs in the assay:', nrow(obj), '\n')

res <- FindMarkers(obj, ident.1 = 'pos', ident.2 = 'neg', assay = 'chromvar', mean.fxn = rowMeans, fc.name = 'avg_diff',
                   logfc.threshold = 0, verbose = FALSE)
write.csv(res, file.path(DIR, sprintf('clone%s_seurat_findmarkers.csv', K)))

py <- read.csv(file.path(DIR, sprintf('clone%s_motifs.csv', K)), check.names = FALSE)
rownames(py) <- py$motif_id
tested_py <- py$motif_id[py$genotype__tested == 'True']          # Python writes booleans as True / False
sig_py <- py$motif_id[py$genotype__significant == 'True']
m <- rownames(res)
fin <- res$p_val > 0
check <- data.frame(
  tested_motifs_seurat = nrow(res), tested_motifs_python = length(tested_py),
  same_tested_motifs = setequal(m, tested_py),
  significant_seurat = sum(res$p_val_adj < 0.05), significant_python = length(sig_py),
  same_significant_motifs = setequal(m[res$p_val_adj < 0.05], sig_py),
  max_abs_diff_log10_p = max(abs(log10(res$p_val[fin]) - py[m[fin], 'genotype__log10p'])),
  max_abs_diff_avg_diff = max(abs(res$avg_diff - py[m, 'genotype__avg_diff'])),
  max_abs_diff_pct = max(abs(c(res$pct.1 - py[m, 'genotype__pct_1'], res$pct.2 - py[m, 'genotype__pct_2']))))
check$same <- check$same_tested_motifs & check$same_significant_motifs & check$max_abs_diff_log10_p < 1e-8 &
  check$max_abs_diff_avg_diff < 1e-10 & check$max_abs_diff_pct == 0
write.csv(check, file.path(DIR, sprintf('clone%s_seurat_check.csv', K)), row.names = FALSE)
print(t(check))
