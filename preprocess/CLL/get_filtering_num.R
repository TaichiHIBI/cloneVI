library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
#library(ggrastr)
library(dplyr)
#library(BuenColors)

peaks <- diffloop::bedToGRanges("mtscATAC-seq/mtscATACpaper_reproducibility/figure_CLL/data/CLL_heme.fixedwidthpeaks.bed")


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

# Import experiment to make an RDS of the SE and a basic qc plot
importExperiment <- function(exp, peaks, frip_threshold = 0.5){
  
  qcdf <- fread(paste0("mtscATAC-seq/mtscATACpaper_reproducibility/figure_CLL/data/singlecell_sumstats/", exp,"_v12-mtMask_singlecell.csv.gz"), header = TRUE, sep = ",") %>% 
    data.frame() %>% filter(cell_id != "None")
  bc <- as.character(qcdf$barcode)
  SE <- getCountsFromFrags(paste0("mtscATAC-seq/mtscATACpaper_large_data_files/source/cellranger_output/", exp, "_v12-mtMask_fragments.tsv"), peaks, bc)
  
  df <- data.frame(colData(SE))
  index <- seq(0,nrow(df)-1,by=1)
  df <- cbind(df, index)
  
  if(FALSE){
    df$density <- get_density(log10(df$depth + 1), df$FRIP*100)
    p1 <- ggplot(df %>% arrange(density), aes(x = log10(depth + 1),  y = FRIP*100, color = density)) +
      geom_point_rast(size = 0.1) + scale_color_gradientn(colors = jdb_palette("brewer_celsius")) +
      pretty_plot() + L_border() + labs(x = "log10 # Fragments", y = "FRIP") +
      theme(legend.position = "bottom") + ggtitle(exp) +
      geom_vline(xintercept = 3.5, linetype = 2) +
      geom_hline(yintercept = 40, linetype = 2) +
      scale_y_continuous(limits = c(0,100))
  }
  # Import mtDNA DF
  cov_mtDNA <- data.frame(colData(readRDS(paste0("mtscATAC-seq/mtscATACpaper_large_data_files/source/mgatk_output/", exp, "_v12-mtMask_mgatk.rds"))))
  vec <- as.numeric(cov_mtDNA$depth); names(vec) <- as.character(cov_mtDNA$sample)
  colData(SE)$mtDNAcoverage <- vec[as.character(df$sample)] %>% unname()
  df$mtDNAcoverage <- vec[as.character(df$sample)] %>% unname()
  
  del_row <- df[log10(df$depth) >= 3 & df$FRIP >= frip_threshold & df$mtDNAcoverage >= 20 ,] # modify
  #SE2 <- SE[, log10(df$depth) >= 3 & df$FRIP >= frip_threshold & df$mtDNAcoverage >= 20]
  #saveRDS(SE2, file = paste0("mtscATAC-seq/mtscATACpaper_large_data_files/intermediate/", exp, ".rds"))
  if(!dir.exists("preprocess/CLL/row_del")){
    dir.create("preprocess/CLL/row_del")
    }
  write.table(del_row, file = paste0("preprocess/CLL/row_del/", exp, "_del.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  return(del_row)
}

# Import and process all samples

# gunzip flagmentfile!!

importExperiment("CLL_PT1_bulk1", peaks, 0.5)
importExperiment("CLL_PT1_bulk2", peaks, 0.5)
importExperiment("CLL_PT1_CD19pos", peaks, 0.5)
importExperiment("CLL_PT1_CD19neg", peaks, 0.5)
importExperiment("CLL_PT2_CD19pos", peaks, 0.5)
importExperiment("CLL_PT2_CD19neg", peaks, 0.5)

# Manually process public data that doesn't really have mitochondrial DNA content
#bc_5k <- fread("mtscATAC-seq/mtscATACpaper_large_data_files/source/other/public_atac_data_10x/atac_pbmc_5k_nextgem_barcodes.tsv", header = FALSE)[[1]]
#frag_file <- "mtscATAC-seq/mtscATACpaper_large_data_files/source/other/public_atac_data_10x/atac_pbmc_5k_nextgem_fragments.tsv"
#SE <- getCountsFromFrags(frag_file, peaks,bc_5k)
#SE2 <- SE[, colData(SE)$FRIP >= 0.5 & colData(SE)$depth >= 1000]
#dim(SE2)
#saveRDS(SE2, file = paste0("mtscATAC-seq/mtscATACpaper_large_data_files/intermediate/public_atac_pbmc_5k_ng.rds"))

