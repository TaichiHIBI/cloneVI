#!/usr/bin/env Rscript
# ドナー単位の mtSNV 検出 + Signac FindClonotypes（**1プロセス=1ドナー**）。
# 集約単位を per-donor にする以外は Signac mito vignette
#   (https://stuartlab.org/signac/articles/mito) を忠実に踏襲する。
#   IdentifyVariants → high.conf を subset → AlleleFreq → FindClonotypes。
# クローンは「同一ドナーの細胞系譜」内にしか存在しないので、検出・変異選択・
# クローノタイプ識別をドナー単位（そのドナーの全GSMをプール）で行う。
#
# 大ドナー(35k–47k細胞)では 132,552行×細胞数 の密化で 38–50GB になるため、
# **ドナーごとに別Rプロセス**で実行し（run_clonotypes_perdonor.sh）、プール行列は
# 逐次 cbind + rm/gc、obj 生成後に疎行列を解放して密コピーを1個に抑える。
#
# 使い方: conda run -n Signac Rscript clonotypes_perdonor_signac.R <COHORT> <DONOR>
# 入力 : cohort_<COHORT>_atac_cells.csv (cell=GSM#bc, gsm, donor, tissue), pack_<GSM>/mgatk/*
# 出力(per-donor):
#   cohort_<COHORT>_<DONOR>_variant_stats.tsv   IdentifyVariants 統計(+donor,high_conf)
#   cohort_<COHORT>_<DONOR>_clonotype.csv       cell(GSM#bc), donor, clonotype, mtDNA_depth
#   figs_clonotypes_<COHORT>/<DONOR>_variantplot.png / _clonotype_heatmap.png
suppressMessages({
  library(Seurat); library(Signac); library(Matrix); library(readr); library(ggplot2)
  library(data.table)
})
set.seed(1234)
CACHEDIR <- "cache_clono"   # per-sample(cohort-subset) 8ブロック行列の RDS キャッシュ

args   <- commandArgs(trailingOnly = TRUE)
COHORT <- if (length(args) >= 1) args[1] else "Ovarian"
DONOR  <- if (length(args) >= 2) args[2] else stop("DONOR を指定してください")
CELLS  <- sprintf("cohort_%s_atac_cells.csv", COHORT)
OUT_VS <- sprintf("cohort_%s_%s_variant_stats.tsv", COHORT, DONOR)
OUT_CL <- sprintf("cohort_%s_%s_clonotype.csv", COHORT, DONOR)
FIGDIR <- sprintf("figs_clonotypes_%s", COHORT)
dir.create(FIGDIR, showWarnings = FALSE)

# ===== CONFIG = Signac mito vignette の high.conf / FindClonotypes 既定 =====
HC <- list(min_n_cells_conf_detected = 5, min_strand_correlation = 0.65, min_vmr = 0.01)
CL <- list(metric = "cosine", resolution = 1, k = 10, algorithm = 3)
# ==========================================================================

pick    <- function(...) { for (p in c(...)) if (file.exists(p)) return(p); stop("not found: ", ...) }
# mtx を data.table::fread で高速読込（readMM と完全一致・約5–10倍速）
read_mm <- function(p) {
  dims <- scan(p, what=integer(), skip=1, nlines=1, quiet=TRUE)   # rows cols nnz
  dt <- data.table::fread(p, skip=2, header=FALSE, col.names=c("i","j","x"),
                          colClasses=c("integer","integer","double"), showProgress=FALSE)
  Matrix::sparseMatrix(i=dt$i, j=dt$j, x=dt$x, dims=c(dims[1], dims[2]))
}

# 1 サンプル(GSM)の 8ブロック counts（行=A,C,T,G fwd→A,C,T,G rev × position, 列=GSM#bc）。
# keep_cells で読込直後に列を絞りメモリ節約。cohort-subset 済みを RDS キャッシュ（再実行が一瞬）。
build_sample_counts <- function(gsm, keep_cells, cohort) {
  cf <- file.path(CACHEDIR, sprintf("%s_%s.rds", cohort, gsm))
  if (file.exists(cf)) { cat(sprintf("  %s: cache hit\n", gsm)); return(readRDS(cf)) }
  MG <- file.path(paste0("pack_", gsm), "mgatk")
  bc <- readr::read_tsv(pick(file.path(MG,"barcodes.tsv")), col_names=FALSE, show_col_types=FALSE)$X1
  cn <- paste0(gsm, "#", bc)
  keepidx <- which(cn %in% keep_cells)
  strand_block <- function(strand) {
    do.call(rbind, lapply(c("A","C","T","G"), function(b) {
      m <- read_mm(pick(file.path(MG, sprintf("%s_counts_%s.mtx", b,
                        ifelse(strand=="fwd","fw","rev")))))
      m <- m[, keepidx, drop = FALSE]
      rownames(m) <- paste(b, seq_len(nrow(m)), strand, sep = "-")
      m
    }))
  }
  counts <- rbind(strand_block("fwd"), strand_block("rev"))
  colnames(counts) <- cn[keepidx]
  dir.create(CACHEDIR, showWarnings = FALSE)
  saveRDS(counts, cf)
  counts
}

refallele_of <- function(gsm) {
  MG <- file.path(paste0("pack_", gsm), "mgatk")
  g  <- readr::read_tsv(pick(file.path(MG,"genes.tsv")), show_col_types=FALSE)
  g  <- g[order(g$start), ]
  data.frame(pos = as.integer(g$start), ref = toupper(as.character(g$refAllele)),
             stringsAsFactors = FALSE)
}

# --- コホート細胞リスト（このドナーのみ） ---
cells_df <- readr::read_csv(CELLS, show_col_types = FALSE)
cells_df$cell  <- as.character(cells_df$cell)
cells_df$gsm   <- as.character(cells_df$gsm)
cells_df$donor <- as.character(cells_df$donor)
sub  <- cells_df[cells_df$donor == DONOR, ]
gsms <- sort(unique(sub$gsm))
if (nrow(sub) == 0) stop("donor に該当細胞なし: ", DONOR)
cat(sprintf("COHORT=%s DONOR=%s : %d samples, %d cohort cells\n",
            COHORT, DONOR, length(gsms), nrow(sub)))

refallele <- refallele_of(gsms[1])
maxpos    <- nrow(refallele)

# --- プール行列を逐次 cbind（毎回 rm + gc で密コピーを1個に抑える） ---
counts <- NULL
for (g in gsms) {
  keep <- sub$cell[sub$gsm == g]
  cnt  <- build_sample_counts(g, keep, COHORT)
  cat(sprintf("  %s: %d cells\n", g, ncol(cnt)))
  counts <- if (is.null(counts)) cnt else cbind(counts, cnt)
  rm(cnt); gc(verbose = FALSE)
}
stopifnot(nrow(counts) == 8 * maxpos)
cat(sprintf("  pooled: %d cells, %.1f GB sparse\n",
            ncol(counts), as.numeric(object.size(counts)) / 1e9))

depth <- Matrix::colSums(counts) / maxpos
obj <- CreateSeuratObject(counts = CreateAssayObject(counts = counts, min.cells = 0,
                          min.features = 0), assay = "mito")
obj$mtDNA_depth <- depth
obj$donor <- DONOR
rm(counts); gc(verbose = FALSE)   # ★ 疎行列を解放してから IdentifyVariants の密化に備える

# --- 変異同定（vignette）---
variants <- IdentifyVariants(obj, refallele = refallele)
variants$donor <- DONOR
hc <- with(variants,
           n_cells_conf_detected >= HC$min_n_cells_conf_detected &
           strand_correlation    >= HC$min_strand_correlation &
           vmr                   >  HC$min_vmr)
hc[is.na(hc)] <- FALSE
variants$high_conf <- hc
write_tsv(variants, OUT_VS)
high.conf <- variants[hc, , drop = FALSE]
cat(sprintf("  high.conf 変異: %d / %d\n", nrow(high.conf), nrow(variants)))

# --- 図1: VariantPlot（必須・vignette）---
tryCatch({
  p <- VariantPlot(variants = variants,
                   min.cells = HC$min_n_cells_conf_detected,
                   concordance.threshold = HC$min_strand_correlation,
                   vmr.threshold = HC$min_vmr) + ggtitle(sprintf("%s (%s)", DONOR, COHORT))
  ggsave(file.path(FIGDIR, sprintf("%s_variantplot.png", DONOR)), p, width = 6, height = 5, dpi = 150)
}, error = function(e) message("  VariantPlot skip: ", conditionMessage(e)))

if (nrow(high.conf) < 2) {
  cat("  → high.conf 変異が不足。AlleleFreq/FindClonotypes はスキップ。\n")
  write_csv(data.frame(cell = colnames(obj), donor = DONOR,
                       clonotype = NA_character_, mtDNA_depth = obj$mtDNA_depth,
                       row.names = NULL), OUT_CL)
  quit(save = "no", status = 0)
}

# --- アレル頻度 → クローノタイプ識別（vignette 既定値）---
obj <- AlleleFreq(obj, variants = high.conf$variant, assay = "mito", new.assay.name = "alleles")
DefaultAssay(obj) <- "alleles"
obj <- FindClonotypes(obj, assay = "alleles", metric = CL$metric,
                      resolution = CL$resolution, k = CL$k, algorithm = CL$algorithm)
ct <- paste0(DONOR, "_", as.character(Idents(obj)))   # 例: SU-O-005_3
cat(sprintf("  FindClonotypes: %d クローノタイプ\n", length(unique(ct))))
print(table(clonotype = Idents(obj)))
write_csv(data.frame(cell = colnames(obj), donor = DONOR,
                     clonotype = ct, mtDNA_depth = obj$mtDNA_depth,
                     row.names = NULL), OUT_CL)

# --- 図2: クローノタイプ別アレル頻度ヒートマップ（vignette）---
# 全細胞を列にすると大ドナーで激重(かつOOM危険)なので、表示は最大 HM_CELLS 細胞に
# ダウンサンプル(データには影響なし・図のみ)。クラスタごとの構造は保たれる。
HM_CELLS <- 4000
tryCatch({
  set.seed(1)
  cells_h <- colnames(obj)
  if (length(cells_h) > HM_CELLS) cells_h <- sample(cells_h, HM_CELLS)
  h <- DoHeatmap(obj, cells = cells_h, features = rownames(obj),
                 slot = "data", disp.max = 0.1) + scale_fill_viridis_c() +
       ggtitle(sprintf("%s allele freq by clonotype (%d/%d cells shown)",
                       DONOR, length(cells_h), ncol(obj)))
  ggsave(file.path(FIGDIR, sprintf("%s_clonotype_heatmap.png", DONOR)), h,
         width = 8, height = 6, dpi = 150, limitsize = FALSE)
}, error = function(e) message("  heatmap skip: ", conditionMessage(e)))

cat(sprintf("DONE %s/%s -> %s , %s\n", COHORT, DONOR, OUT_VS, OUT_CL))
