#!/usr/bin/env Rscript
# 1サンプルの FeatureMatrix（統一peakで再カウント）を独立Rプロセスで実行し RDS 保存。
# build_atac_cohort_signac.R から並列起動される想定（fork不使用→Rsamtools/tabix安全）。
# 使い方: Rscript fm_one_sample.R <COHORT> <GSM>
suppressMessages({ library(Signac); library(GenomicRanges); library(Matrix); library(readr) })
a <- commandArgs(trailingOnly = TRUE); COHORT <- a[1]; GSM <- a[2]
std_chr <- paste0("chr", c(1:22, "X", "Y"))
PEAK_MIN_W <- 20; PEAK_MAX_W <- 10000
MAP <- "../gsm_donor_map.csv"; INCLUDE_PBMC <- TRUE
CACHEDIR <- paste0("cache_atac_", COHORT); dir.create(CACHEDIR, showWarnings = FALSE)
outf <- file.path(CACHEDIR, paste0(GSM, ".rds"))
if (file.exists(outf) && !is.null(readRDS(outf))) { cat("cached:", GSM, "\n"); quit(save="no") }

map <- readr::read_csv(MAP, show_col_types = FALSE)
map$cohort <- ifelse(grepl("^SU-L", map$donor), "Lung", "Ovarian")
m <- map[map$cohort == COHORT, ]; if (!INCLUDE_PBMC) m <- m[m$tissue != "PBMC", ]
samples <- m$gsm
cr_of <- function(s) file.path(paste0("pack_", s), "cellranger")

# 統一peak（全サンプルのpeaksをunion→reduce→幅フィルタ; 決定的なので各プロセスで同一）
read_peaks <- function(s) {
  bed <- readr::read_tsv(file.path(cr_of(s), "peaks.bed"), col_names = FALSE,
                         show_col_types = FALSE, comment = "#")
  bed <- bed[bed$X1 %in% std_chr, ]
  GRanges(bed$X1, IRanges(bed$X2 + 1, bed$X3))
}
combined <- GenomicRanges::reduce(do.call(c, lapply(samples, read_peaks)))
w <- width(combined); combined <- combined[w >= PEAK_MIN_W & w <= PEAK_MAX_W]
combined <- GenomeInfoDb::keepSeqlevels(combined, std_chr, pruning.mode = "coarse")

cr <- cr_of(GSM)
bcs <- readr::read_tsv(file.path(cr, "filtered_peak_bc_matrix", "barcodes.tsv"),
                       col_names = FALSE, show_col_types = FALSE)$X1
cells <- setNames(bcs, paste0(GSM, "#", bcs))
frag <- CreateFragmentObject(file.path(cr, "fragments.tsv.gz"), cells = cells,
                             validate.fragments = FALSE, verbose = FALSE)
t0 <- Sys.time()
# cells はオブジェクト側セル名(GSM#bc=names)を渡す。名前付きベクトル全体だと無次元になる。
mtx <- FeatureMatrix(fragments = frag, features = combined, cells = Cells(frag), verbose = FALSE)
stopifnot(!is.null(dim(mtx)), ncol(mtx) == length(bcs))
saveRDS(mtx, outf)
cat(sprintf("done %s: %d peaks x %d cells (%.1f min) -> %s\n", GSM, nrow(mtx), ncol(mtx),
            as.numeric(difftime(Sys.time(), t0, units = "mins")), outf))
