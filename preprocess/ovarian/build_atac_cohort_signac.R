#!/usr/bin/env Rscript
# コホート単位 ATAC h5ad ビルダー（Lung / Ovarian）— cloneVI 入力用。
# per-GSM ではなく、コホート内の全GSMを統合して特徴量選択(q75)・LSI まで1回で行う。
#   統一peak set → 各GSMを FeatureMatrix で再カウント → merge(barcode=GSM#bc)
#   → QC → TFIDF → FindTopFeatures(q75) → LSI/UMAP/クラスタ
# 患者(donor)=batch は obs に記録（バッチ補正は cloneVI 側、または harmony 導入時に任意）。
# クローンコールは別途ドナー単位（mtSNV側）で行う想定。
# 使い方: conda run -n Signac Rscript build_atac_cohort_signac.R <COHORT>   COHORT=Ovarian|Lung
suppressMessages({
  library(Seurat); library(Signac); library(SeuratDisk)
  library(GenomeInfoDb); library(GenomicRanges); library(Matrix)
  library(readr); library(dplyr)
})
set.seed(1234)

args   <- commandArgs(trailingOnly = TRUE)
COHORT <- if (length(args) >= 1) args[1] else "Ovarian"

# ===== CONFIG =====================================================
INCLUDE_PBMC <- TRUE     # FALSE で tumor/normal/met のみ（論文の細胞状態アトラス流）
# GeneActivity はメモリ大。大規模コホート(Lung)は別プロセス build_geneactivity_cohort.R に分離。
# 環境変数で制御: DO_GENEACTIVITY=FALSE で本体ではスキップ（既定TRUE=Ovarian再現）。
DO_GENEACTIVITY <- toupper(Sys.getenv("DO_GENEACTIVITY", "TRUE")) == "TRUE"
USE_HARMONY  <- TRUE     # harmony 導入済みなら donor で batch 補正（無ければ自動skip）
MAP    <- "../gsm_donor_map.csv"
QC <- list(min_nCount = 1000, max_nCount = 50000, min_pct_reads_in_peaks = 40,
           max_blacklist_ratio = 0.05, min_TSS = 3, max_nucleosome = 4)
FEAT_MIN_CUTOFF <- "q75"; DIMS <- 2:50; CLUSTER_RES <- 0.7
PEAK_MIN_W <- 20; PEAK_MAX_W <- 10000   # 統一peak の幅フィルタ(極端な結合peak除去)
std_chr <- paste0("chr", c(1:22, "X", "Y"))
# =================================================================

map <- readr::read_csv(MAP, show_col_types = FALSE)
map$cohort <- ifelse(grepl("^SU-L", map$donor), "Lung", "Ovarian")
m <- map[map$cohort == COHORT, ]
if (!INCLUDE_PBMC) m <- m[m$tissue != "PBMC", ]
samples <- m$gsm
cr_of <- function(s) file.path(paste0("pack_", s), "cellranger")
stopifnot(all(file.exists(file.path(cr_of(samples), "fragments.tsv.gz"))))
cat(sprintf("COHORT=%s  samples=%d  donors=%d  (PBMC=%s)\n", COHORT, length(samples),
            length(unique(m$donor)), INCLUDE_PBMC))

# --- 1. 統一 peak set: 全GSMのpeaksを結合→reduce→幅フィルタ（標準染色体のみ） ---
read_peaks <- function(s) {
  bed <- readr::read_tsv(file.path(cr_of(s), "peaks.bed"), col_names = FALSE,
                         show_col_types = FALSE, comment = "#")
  bed <- bed[bed$X1 %in% std_chr, ]
  GenomicRanges::GRanges(bed$X1, IRanges::IRanges(bed$X2 + 1, bed$X3))  # bed 0-based → 1-based
}
combined <- GenomicRanges::reduce(do.call(c, lapply(samples, read_peaks)))
w <- GenomicRanges::width(combined)
combined <- combined[w >= PEAK_MIN_W & w <= PEAK_MAX_W]
combined <- GenomeInfoDb::keepSeqlevels(combined, std_chr, pruning.mode = "coarse")
cat(sprintf("統一peak: %d (幅%d-%dbp, 標準染色体)\n", length(combined), PEAK_MIN_W, PEAK_MAX_W))

# --- 2. 各GSMを統一peakで FeatureMatrix 再カウント（barcode=GSM#bc）→ merge ---
#   FeatureMatrix は単一スレッドで重い(統一peak大)ので mclapply でサンプル並列。
#   per-sample に RDS キャッシュし、進捗は逐次フラッシュ(cache_atac_<COHORT>/progress.log)。
library(parallel)
CACHEDIR <- paste0("cache_atac_", COHORT); dir.create(CACHEDIR, showWarnings = FALSE)
PROG <- file.path(CACHEDIR, "progress.log"); cat("", file = PROG)
frag_list <- list(); cells_list <- list(); sc_list <- list()
for (s in samples) {        # Fragment/cells/singlecell は軽いので親で作る
  cr <- cr_of(s)
  bcs <- readr::read_tsv(file.path(cr, "filtered_peak_bc_matrix", "barcodes.tsv"),
                         col_names = FALSE, show_col_types = FALSE)$X1
  cells_list[[s]] <- setNames(bcs, paste0(s, "#", bcs))   # name=object側, value=file内barcode
  frag_list[[s]] <- CreateFragmentObject(path = file.path(cr, "fragments.tsv.gz"),
                               cells = cells_list[[s]], validate.fragments = FALSE, verbose = FALSE)
  sc <- as.data.frame(readr::read_csv(file.path(cr, "singlecell.csv"), show_col_types = FALSE))
  rownames(sc) <- paste0(s, "#", sc$barcode); sc_list[[s]] <- sc
}
# per-sample FeatureMatrix キャッシュを読む。未生成のものは直列で計算(fork不使用=安全)。
#   並列化したい場合は事前に `fm_one_sample.R <COHORT> <GSM>` を別プロセスで流す(下記driver)。
valid_cache <- function(s) { f <- file.path(CACHEDIR, paste0(s, ".rds"))
  file.exists(f) && !is.null(tryCatch(readRDS(f), error = function(e) NULL)) }
for (s in samples) if (!valid_cache(s)) {
  cat(sprintf("[%s] FeatureMatrix %s (cache無し→計算)...\n", format(Sys.time(),"%H:%M:%S"), s))
  mtx <- FeatureMatrix(fragments = frag_list[[s]], features = combined,
                       cells = Cells(frag_list[[s]]), verbose = FALSE)  # 名前付きvec全体は不可
  saveRDS(mtx, file.path(CACHEDIR, paste0(s, ".rds")))
}
counts_list <- lapply(samples, function(s) readRDS(file.path(CACHEDIR, paste0(s, ".rds"))))
merged <- do.call(cbind, counts_list)   # 行=統一peak(同順), 列=全GSMの GSM#bc
rm(counts_list); gc()                   # 27個の元行列(~10GB)を解放(大規模コホートのOOM回避)
feat_str <- GRangesToString(combined, sep = c("-", "-"))   # FeatureMatrixは行名を付けないので明示
stopifnot(nrow(merged) == length(feat_str))
rownames(merged) <- feat_str
cat(sprintf("merged: %d peaks x %d cells\n", nrow(merged), ncol(merged)))

# --- 3. ChromatinAssay（ranges=combined を明示し行名パースを回避。
#         genome= は GenomeInfoDb 1.30.1 バグ回避で渡さない）---
chrom_assay <- CreateChromatinAssay(counts = merged, ranges = combined, sep = c("-", "-"),
                                    fragments = frag_list, min.cells = 10, min.features = 200)
obj <- CreateSeuratObject(counts = chrom_assay, assay = "peaks")
rm(merged, chrom_assay); gc()           # merged複製(~10GB)を解放

# --- 4. メタデータ（GSM#bc から donor/tissue/sorting を付与）---
cellmeta <- data.frame(cell = colnames(obj),
                       gsm = sub("#.*$", "", colnames(obj)), stringsAsFactors = FALSE)
cellmeta <- dplyr::left_join(cellmeta, m[, c("gsm","donor","tissue","sorting")], by = "gsm")
rownames(cellmeta) <- cellmeta$cell
obj <- AddMetaData(obj, cellmeta[colnames(obj), c("gsm","donor","tissue","sorting")])
obj$cohort <- COHORT
# cellranger singlecell.csv の per-cell 指標
sc_all <- do.call(rbind, unname(sc_list))   # unname: 名前付きlistだとrownameにGSM接頭辞が二重付与され照合失敗する
sc_all <- sc_all[colnames(obj), , drop = FALSE]
obj <- AddMetaData(obj, sc_all[, c("passed_filters","peak_region_fragments")])  # 必要2列のみ(版間で列差あり)
obj$pct_reads_in_peaks <- obj$peak_region_fragments / obj$passed_filters * 100
stopifnot(mean(is.na(obj$pct_reads_in_peaks)) < 0.5)   # 結合失敗(全NA)を早期検知
data("blacklist_hg38_unified", package = "Signac")
obj$blacklist_ratio <- FractionCountsInRegion(obj, assay = "peaks", regions = blacklist_hg38_unified)

# --- 5. annotation（標準染色体のみ手動リネームでバグ回避）→ Nucleosome/TSS ---
ok_anno <- tryCatch({
  library(EnsDb.Hsapiens.v86)
  ann <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  std <- c(as.character(1:22), "X", "Y", "MT")
  ann <- ann[as.character(seqnames(ann)) %in% std]
  ann <- keepSeqlevels(ann, std, pruning.mode = "coarse")
  ann <- renameSeqlevels(ann, setNames(c(paste0("chr", c(1:22,"X","Y")), "chrM"), std))
  Annotation(obj) <- ann; TRUE
}, error = function(e){ message("annotation skip: ", conditionMessage(e)); FALSE })
obj <- NucleosomeSignal(obj)
ok_tss <- tryCatch({ obj <- TSSEnrichment(obj, fast = TRUE); TRUE },
                   error = function(e){ message("TSS skip: ", conditionMessage(e)); FALSE })

# --- 6. QC 分布表示 → セルフィルタ ---
qc_metrics <- c("nCount_peaks","pct_reads_in_peaks","blacklist_ratio","nucleosome_signal",
                if (ok_tss) "TSS.enrichment")
cat("QC 分位点[0/25/50/75/100%] (フィルタ前):\n")
for (mx in qc_metrics)
  cat(sprintf("  %-20s %s\n", mx, paste(signif(quantile(obj[[mx]][,1], c(0,.25,.5,.75,1),
      na.rm=TRUE), 3), collapse="  ")))
n0 <- ncol(obj)
qc_pass <- obj$nCount_peaks > QC$min_nCount & obj$nCount_peaks < QC$max_nCount &
           obj$pct_reads_in_peaks > QC$min_pct_reads_in_peaks &
           obj$blacklist_ratio < QC$max_blacklist_ratio &
           obj$nucleosome_signal < QC$max_nucleosome
if (ok_tss) qc_pass <- qc_pass & obj$TSS.enrichment > QC$min_TSS
qc_pass[is.na(qc_pass)] <- FALSE
cat(sprintf("QC: %d → %d cells 通過 (%.1f%%)\n", n0, sum(qc_pass), 100*sum(qc_pass)/n0))
obj <- subset(obj, cells = colnames(obj)[qc_pass])
cat("donor別 通過セル数:\n"); print(table(obj$donor, obj$tissue))

# --- 7. 正規化 → 特徴量選択(q75, コホート全体で) → LSI → (任意 harmony) → UMAP/クラスタ ---
PF <- paste0("progress_atac_", COHORT, ".log")
pf <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), sprintf(...)),
                          file = PF, append = TRUE)
  cat(sprintf("  mem: %.1f GB\n", sum(gc()[,2])/1024), file = PF, append = TRUE) }
pf("QC後 %d cells. RunTFIDF...", ncol(obj))
obj <- RunTFIDF(obj); pf("TFIDF done")
obj <- FindTopFeatures(obj, min.cutoff = FEAT_MIN_CUTOFF)
cat(sprintf("FindTopFeatures(%s): VariableFeatures = %d / %d peaks\n",
            FEAT_MIN_CUTOFF, length(VariableFeatures(obj)), nrow(obj)))
pf("FindTopFeatures done (%d VF). RunSVD...", length(VariableFeatures(obj)))
obj <- RunSVD(obj); pf("RunSVD done")
emb <- "lsi"
if (USE_HARMONY && requireNamespace("harmony", quietly = TRUE)) {
  obj <- harmony::RunHarmony(obj, group.by.vars = "donor", reduction = "lsi",
                             assay.use = "peaks", project.dim = FALSE, reduction.save = "harmony")
  emb <- "harmony"; cat("harmony: donorでbatch補正 → reduction='harmony'\n")
} else cat("harmony: 未使用（donorはobsに記録, 補正はcloneVI側）\n")
pf("RunUMAP...")
obj <- RunUMAP(obj, reduction = emb, dims = DIMS); pf("UMAP done. FindNeighbors...")
obj <- FindNeighbors(obj, reduction = emb, dims = DIMS); pf("FindNeighbors done. FindClusters...")
obj <- FindClusters(obj, resolution = CLUSTER_RES, verbose = FALSE, algorithm = 3); pf("FindClusters done")
obj$lsi_clusters <- obj$seurat_clusters

has_rna <- FALSE
if (DO_GENEACTIVITY && ok_anno) tryCatch({
  ga <- GeneActivity(obj); obj[["RNA"]] <- CreateAssayObject(counts = ga)
  obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize",
                       scale.factor = median(obj$nCount_RNA)); has_rna <- TRUE
}, error = function(e) message("GeneActivity skip: ", conditionMessage(e)))

DefaultAssay(obj) <- "peaks"
OUT <- paste0("cohort_", COHORT)
pf("SaveH5Seurat...")
SaveH5Seurat(obj, filename = paste0(OUT, "_atac.h5seurat"), overwrite = TRUE); pf("SaveH5Seurat done. Convert...")
Convert(paste0(OUT, "_atac.h5seurat"), dest = "h5ad", overwrite = TRUE); pf("Convert done")
cat("DONE ->", paste0(OUT, "_atac.h5ad"), "\n")
write.csv(data.frame(cell = colnames(obj), gsm = obj$gsm, donor = obj$donor,
          tissue = obj$tissue), paste0(OUT, "_atac_cells.csv"), row.names = FALSE)
print(obj)

if (has_rna) tryCatch({
  rna <- CreateSeuratObject(counts = GetAssayData(obj, assay="RNA", slot="counts"),
                            assay = "RNA", meta.data = obj[[]])
  rna <- SetAssayData(rna, assay="RNA", slot="data",
                      new.data = GetAssayData(obj, assay="RNA", slot="data"))
  if ("umap" %in% Reductions(obj)) { red <- obj[["umap"]]; red@assay.used <- "RNA"; rna[["umap"]] <- red }
  SaveH5Seurat(rna, filename = paste0(OUT, "_geneactivity.h5seurat"), overwrite = TRUE)
  Convert(paste0(OUT, "_geneactivity.h5seurat"), dest = "h5ad", overwrite = TRUE)
  cat("DONE ->", paste0(OUT, "_geneactivity.h5ad"), "\n")
}, error = function(e) message("GeneActivity h5ad skip: ", conditionMessage(e)))
