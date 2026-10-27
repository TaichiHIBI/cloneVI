#!/usr/bin/env Rscript
# mgatk RDS (RangedSummarizedExperiment) -> MTX群 + genes.tsv + barcodes.tsv
# rds2h5ad.ipynb が期待する形式に合わせる。
suppressMessages({ library(SummarizedExperiment); library(Matrix) })
args <- commandArgs(trailingOnly = TRUE)
rds <- args[1]; out <- args[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

se <- readRDS(rds)
cat("assayNames:", paste(assayNames(se), collapse=", "), "\n")
cat("dim (positions x cells):", paste(dim(se), collapse=" x "), "\n")

# 各 assay を positions x cells の MTX で出力（rds2h5ad は読み込み後に .T する）
for (a in assayNames(se)) {
  m <- as(assay(se, a), "dgCMatrix")
  Matrix::writeMM(m, file = file.path(out, paste0(a, ".mtx")))
}

# 行(=mtDNA position)メタ: seqnames,start,end,width,strand,refAllele
gr <- rowRanges(se)
genes <- as.data.frame(gr)
if (!is.null(rowData(se)$refAllele)) genes$refAllele <- rowData(se)$refAllele
write.table(genes, file = file.path(out, "genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 列(=cell barcode)
writeLines(colnames(se), con = file.path(out, "barcodes.tsv"))
cat("exported to", out, "\n")
