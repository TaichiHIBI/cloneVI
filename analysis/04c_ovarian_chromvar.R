# Ovarian tumor cohort (3/4, part 1): motif activity of every cell (chromVAR), for the motif tests of 04d_ovarian_motif.ipynb.
#
# The Signac standard workflow on the whole cohort (Methods): the 746 vertebrate CORE motifs of JASPAR2020 are matched to
# the peaks (AddMotifs) and the deviation z score of each motif in each cell is computed (RunChromVAR; peaks with no
# counts are dropped, GC bias added, background peaks chosen with the chromVAR defaults).
#
# Input:  cohort_Ovarian_atac.h5seurat of build_atac_cohort_signac.R (preprocessing; only the peak counts are loaded)
# Output: analysis/results/ovarian/chromvar/ chromvar_z.h5 (dataset "z": cells x motifs; "motif_id", "cell"),
#         motifs.csv (motif_id, tf, family, species, n_motif_peaks, n_nonfinite_cells), peaks_used.txt, session_info.txt
# Environment: conda env Signac (R 4.1.1, Seurat 4.0.5, Signac 1.4.0, chromVAR 1.16.0, JASPAR2020, BSgenome hg38)
# Run from the repository root: Rscript analysis/04c_ovarian_chromvar.R   (about one hour with 12 cores)

suppressMessages({
  library(Seurat); library(Signac); library(SeuratDisk); library(Matrix)
  library(JASPAR2020); library(TFBSTools); library(BSgenome.Hsapiens.UCSC.hg38)
  library(BiocParallel); library(hdf5r)
})
register(MulticoreParam(12))
set.seed(1234)

# Signac object of the preprocessing (cohort_Ovarian_atac.h5seurat, deposited in data/mitotrek/); CLONEVI_OVARIAN_PRE sets another directory
PRE <- Sys.getenv('CLONEVI_OVARIAN_PRE', 'data/mitotrek')
OUT <- 'analysis/results/ovarian/chromvar'
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time()
log <- function(...) cat(sprintf('[%s min] ', format(difftime(Sys.time(), t0, units = 'mins'), digits = 3)), ..., '\n')

obj <- LoadH5Seurat(file.path(PRE, 'cohort_Ovarian_atac.h5seurat'), assays = list(peaks = 'counts'), reductions = FALSE,
                    graphs = FALSE, neighbors = FALSE, misc = FALSE, tools = FALSE, verbose = FALSE)
DefaultAssay(obj) <- 'peaks'
log('loaded:', ncol(obj), 'cells,', nrow(obj), 'peaks')

# the saved peak ranges carry genome NA, which AddMotifs cannot merge with BSgenome hg38; the peaks were called on hg38
gr <- obj[['peaks']]@ranges
GenomeInfoDb::genome(gr) <- 'hg38'
obj[['peaks']]@ranges <- gr

pfm <- getMatrixSet(x = JASPAR2020, opts = list(collection = 'CORE', tax_group = 'vertebrates', all_versions = FALSE))
obj <- AddMotifs(object = obj, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
log('AddMotifs:', length(pfm), 'motifs')
obj <- RunChromVAR(object = obj, genome = BSgenome.Hsapiens.UCSC.hg38)
log('RunChromVAR done')

z <- GetAssayData(obj, assay = 'chromvar', slot = 'data')
ids <- rownames(z)
mm <- GetMotifData(obj, assay = 'peaks', slot = 'data')
kept_peaks <- Matrix::rowSums(GetAssayData(obj, assay = 'peaks', slot = 'counts')) > 0     # RunChromVAR drops zero-count peaks
motifs <- data.frame(motif_id = ids,
                     tf = vapply(ids, function(i) name(pfm[[i]]), ''),
                     family = vapply(ids, function(i) paste(tags(pfm[[i]])$family, collapse = ';'), ''),
                     species = vapply(ids, function(i) paste(tags(pfm[[i]])$species, collapse = ';'), ''),
                     n_motif_peaks = as.integer(Matrix::colSums(mm[kept_peaks, ids, drop = FALSE])),
                     n_nonfinite_cells = as.integer(rowSums(!is.finite(as.matrix(z)))))
write.csv(motifs, file.path(OUT, 'motifs.csv'), row.names = FALSE)
writeLines(rownames(obj)[kept_peaks], file.path(OUT, 'peaks_used.txt'))
h5 <- H5File$new(file.path(OUT, 'chromvar_z.h5'), mode = 'w')
h5[['z']] <- as.matrix(z)
h5[['motif_id']] <- ids
h5[['cell']] <- colnames(z)
h5$close_all()
writeLines(capture.output(sessionInfo()), file.path(OUT, 'session_info.txt'))
log('saved', nrow(z), 'motifs x', ncol(z), 'cells; motifs with >= 50 peaks:', sum(motifs$n_motif_peaks >= 50),
    '; motifs with non-finite z:', sum(motifs$n_nonfinite_cells > 0))
