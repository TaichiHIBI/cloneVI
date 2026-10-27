# Cells of CLL patient 1, CD19+ (CLL_PT1_CD19pos) kept for the model: called by cellranger-atac, at least 1,000
# fragments, at least 50% of them in the fixed-width peaks of the source study, and a mitochondrial coverage of at
# least 20.  Writes preprocess/CLL/row_del/CLL_PT1_CD19pos_del.tsv: one row per kept cell, `index` = its 0-based
# position among the cells of the sample (the order of the cellranger barcodes, the mgatk cells and the fragment
# table, all sorted alike).  Despite the name of the file, these are the cells that are kept.
#
# The counting functions are those of the source study's reproducibility code; the version used for the study ran
# over all six CLL libraries from the uncompressed fragment files, of which only CLL_PT1_CD19pos is used.
# Run from the repository root (needs about 30 GB of memory): Rscript preprocess/CLL/get_filtering_num.R
library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(dplyr)

L <- "mtscATAC-seq/mtscATACpaper_large_data_files/source"
R <- "mtscATAC-seq/mtscATACpaper_reproducibility/figure_CLL/data"
peaks <- diffloop::bedToGRanges(file.path(R, "CLL_heme.fixedwidthpeaks.bed"))


# function to get counts
getCountsFromFrags <- function(frag_gz_file,
                               peaks_gr,
                               barcodes){

  # Make GRanges of fragments that are solid for the cells that we care about
  frags_valid <- data.table::fread(paste0(frag_gz_file)) %>%
    data.frame() %>% filter(V4 %in% barcodes) %>%  # filter for barcodes in our search set
    GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE)

  # Get a denominator, per cell
  denom <- table(GenomicRanges::mcols(frags_valid)$V4)
  barcodes_found <- names(denom)

  # Get the overlaps with peaks
  ovPEAK <- GenomicRanges::findOverlaps(peaks_gr, frags_valid)

  # Establish a numeric index for the barcodes for sparse matrix purposes
  id <- factor(as.character(GenomicRanges::mcols(frags_valid)$V4), levels = barcodes_found)

  # Make sparse matrix with counts with peaks by  unique barcode
  countdf <- data.frame(peaks = S4Vectors::queryHits(ovPEAK),
                        sample = as.numeric(id)[S4Vectors::subjectHits(ovPEAK)]) %>%
    dplyr::group_by(peaks,sample) %>% dplyr::summarise(count = n()) %>% data.matrix()

  m <- Matrix::sparseMatrix(i = c(countdf[,1], length(peaks_gr)),
                            j = c(countdf[,2], length(barcodes_found)),
                            x = c(countdf[,3],0))
  colnames(m) <- barcodes_found

  # Make a polished colData
  colData <- data.frame(
    sample = barcodes_found,
    depth = as.numeric(denom),
    FRIP = Matrix::colSums(m)/as.numeric(denom)
  )
  # Make sure that the SE can be correctly constructed
  stopifnot(all(colData$sample == colnames(m)))

  # Make summarized Experiment
  SE <- SummarizedExperiment::SummarizedExperiment(
    rowRanges = peaks_gr,
    assays = list(counts = m),
    colData = colData
  )
  return(SE)
}

# Count the fragments of the cells called by cellranger-atac, add the mitochondrial coverage of mgatk, and write the
# cells that pass the filters
importExperiment <- function(exp, peaks, frip_threshold = 0.5){

  qcdf <- fread(file.path(R, "singlecell_sumstats", paste0(exp, "_v12-mtMask_singlecell.csv.gz")), header = TRUE, sep = ",") %>%
    data.frame() %>% filter(cell_id != "None")
  bc <- as.character(qcdf$barcode)
  SE <- getCountsFromFrags(file.path(L, "cellranger_output", paste0(exp, "_v12-mtMask_fragments.tsv.gz")), peaks, bc)

  df <- data.frame(colData(SE))
  index <- seq(0,nrow(df)-1,by=1)
  df <- cbind(df, index)

  # Import mtDNA DF
  cov_mtDNA <- data.frame(colData(readRDS(file.path(L, "mgatk_output", paste0(exp, "_v12-mtMask_mgatk.rds")))))
  vec <- as.numeric(cov_mtDNA$depth); names(vec) <- as.character(cov_mtDNA$sample)
  df$mtDNAcoverage <- vec[as.character(df$sample)] %>% unname()

  del_row <- df[log10(df$depth) >= 3 & df$FRIP >= frip_threshold & df$mtDNAcoverage >= 20 ,]
  dir.create("preprocess/CLL/row_del", showWarnings = FALSE)
  write.table(del_row, file = paste0("preprocess/CLL/row_del/", exp, "_del.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  return(del_row)
}

kept <- importExperiment("CLL_PT1_CD19pos", peaks, 0.5)
cat(nrow(kept), "cells kept\n")
