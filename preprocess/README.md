# From the public raw data to the model inputs

The scripts expect the raw data in folders at the repository root (`MAESTER/`, `mtscATAC-seq/`, `GSE302113/`) and are run
from the repository root, the notebooks from their own folder. The files they produce are deposited on Zenodo (`data/`),
so the results do not require rerunning them. Each dataset ends with a notebook `build_model_inputs.ipynb` that writes
the model inputs to a `rebuilt/` folder next to the deposited ones and compares the two; the notebooks are kept with
the outputs of that comparison. The ovarian scripts are the ones used for the study, with only their paths changed and
run on a cluster. For MAESTER and CLL the study's code, written for one machine and for more samples than were used,
was gathered into the scripts below for the one sample used, and the steps whose code was not kept were rewritten;
they reproduce the deposited model inputs from the public data, as the notebooks show.

## MAESTER (GEO GSE182685; Miller et al. 2022)
Place the GEO files `GSM5534705_BPDCN712_10x_scRNAseq_Expression.txt.gz`, `GSM5534705_BPDCN712_10x_scRNAseq_MetaData.txt.gz`
and `GSM5534706_BPDCN712_10x_MAESTER_maegatk.rds.gz` in `MAESTER/`, and decompress the last one (`gunzip`; the `.rds`
inside is the maegatk object of the source study used for the study, byte for byte).
1. `MAESTER/export_maegatk.R` — the maegatk allele counts and coverage as `.mtx` files in `MAESTER/rds_output/`
   (identical to the files read for the study).
2. `MAESTER/build_model_inputs.ipynb` — `rna_mgv_counts.h5ad` from these files and `rna_exp.h5ad` from the expression and
   metadata of GEO, written to `data/h5ad_files_RNA/rebuilt/`. `rna_mgv_counts.h5ad` is identical to the deposited
   file. For the study `rna_exp.h5ad` was converted from the source study's Seurat object (`BPDCN712_Seurat_with_TCR.rds`),
   distributed on a web site that is no longer online; the GEO files give the same counts and cell metadata, but not the
   columns of the source study's T/NK, TCR and mtDNA-clone analyses (`MT_clone`, `MT_clone_num` and five others).
   `retrain/prep_maester.py` reads `MT_clone_num` to color a plot drawn after the fit, so it runs on the deposited file.
The 26 variants used are those analysed in the original study (`retrain/prep_maester.py`).

## CLL mtscATAC-seq (Lareau et al. 2021)
Place in `mtscATAC-seq/`:
- from https://osf.io/bupge/, `mtscATACpaper_large_data_files/source/cellranger_output/` with
  `CLL_PT1_CD19pos_peak_bc_matrix/` and `CLL_PT1_CD19pos_v12-mtMask_fragments.tsv.gz` (make its index with
  `tabix -p bed`), and `mtscATACpaper_large_data_files/source/mgatk_output/CLL_PT1_CD19pos_v12-mtMask_mgatk.rds`;
- the source study's reproducibility repository (https://github.com/caleblareau/mtscATACpaper_reproducibility) as
  `mtscATACpaper_reproducibility/` (its `figure_CLL/data/` and `figure_CLL/output/` are read).
1. `CLL/get_filtering_num.R` — the cells kept: called by cellranger-atac, at least 1,000 fragments, at least 50% of them
   in the fixed-width peaks of the source study, mitochondrial coverage of at least 20 (`CLL/row_del/`, included; the
   script writes the same file).
2. `CLL/signac_pp.R` — Signac: the peak-by-cell counts of these cells, peaks observed in at least 10 cells and cells with
   at least 200 detected peaks, then TF-IDF, LSI, UMAP and clustering on the top 25% of peaks, and the gene activity →
   `mtscATAC-seq/h5_files/PT1_pos_pp.h5seurat`/`.h5ad`. The model input is the raw counts of all kept peaks.
3. `CLL/export_mgatk.R` — the mgatk allele counts as `.mtx` files (`mtscATAC-seq/rds_output/PT1_pos_mg/`; the export used
   for the study was not kept) and the clones that the source study defined (`PT1_clone_definition.rds`) as
   `mtscATAC-seq/rds_output/clone_definition/PT1.tsv`.
4. `CLL/build_model_inputs.ipynb` — `PT1_pos.h5ad` and `PT1_pos_mg.h5ad`, written to `data/ATAC/rebuilt/`.
   `PT1_pos_mg.h5ad` is identical to the deposited file. `PT1_pos.h5ad` has the same counts, peaks, cells, clones and
   LSI clusters, and the same LSI up to the sign of its components; the UMAP of Signac, which is not used, differs. Its
   `raw.var` lists the peaks of patient 1, whereas the deposited file lists there those of patient 2, by an error of
   the study's notebook (its counts are those of patient 1, and the model reads the peaks from `var`).

## Ovarian tumor mtscATAC-seq cohort (GEO GSE302113)
1. Download from SRA, alignment and peak calling (cellranger-atac 2.1.0, GRCh38), mitochondrial allele counts
   (mgatk 0.7.0: `mgatk tenx -g hg38` on the cells of the filtered peak-barcode matrix, in the biocontainer
   `quay.io/biocontainers/mgatk:0.7.0--pyhdfd78af_2` pulled by `ovarian/verify_biocontainer.sh`), export to `.mtx`
   (`03b_export_mgatk_mtx.R`) and packing of the files used below. GSM9096478 was run first with the single-library
   scripts `ovarian/01_sra_download.sh` (its SRA run downloaded with `par_get.sh`), `02_cellranger_atac.sh`,
   `03_mgatk_apptainer.sh` and `04_pack_results.sh`; the other ten libraries with the job arrays in `ovarian/array/`
   (`01_sra_array.sh` → `02_cellranger_array.sh` → `03_mgatk_array.sh` → `04_pack_array.sh`, one task per line of
   `samples_ovarian_todo.tsv`, which lists the SRA run and the library; `submit_lung.sh` submits the four arrays for a
   given sheet). The scripts are as they were on the cluster after the runs. Written for a cluster with Univa Grid
   Engine; `<user>`, the bind paths and the scheduler headers must be adapted.
2. `ovarian/build_atac_cohort_signac.R` (with `fm_one_sample.R`; `gsm_donor_map.csv` maps libraries to donors) —
   cohort peak set, recount of every library, cell QC, top 25% of peaks, LSI, gene activity →
   `cohort_Ovarian_atac.h5seurat`/`.h5ad` (the `.h5seurat` is deposited as `data/mitotrek/cohort_Ovarian_atac.h5seurat`
   for the chromVAR step of `analysis/`).
3. `ovarian/run_clonotypes_perdonor.sh`, `clonotypes_perdonor_signac.R` — per-donor variant calling following the Signac
   mitochondrial vignette → `cohort_Ovarian_perdonor_variant_stats.tsv`, `cohort_Ovarian_clonotypes.csv`.
4. `ovarian/build_model_inputs.ipynb` — the model inputs from the outputs of steps 1–3 (in `GSE302113/test_h5ad/`):
   the top-25% peaks with raw counts (`cohort_Ovarian_atac_top25.h5ad`), the depth and mutant counts of the 64 variants
   called in at least one donor (`cohort_Ovarian_mtsnv_perdonor_frac0.5pct.h5ad`, from the mgatk output of every
   library), and the cell types of the source study (`tableS2_ovarian_tumor_labels.csv`, from its Supplementary Table S2;
   place `NIHMS2184724-supplement-Table_S2.xlsx` in `data/mitotrek/paper_supp/`). The files are written to
   `data/mitotrek/rebuilt/` and compared with the deposited ones (all identical); the notebook is kept with its outputs.

## Simulation
`simulation/simulation_datasets.ipynb` generates the simulated data from the Scanpy PBMC3k data. The clone of each cell
and the read counts were drawn without a fixed seed, so the deposited data (`data/simulation/back_0.5_dropout/`) are the
reference; the notebook writes a new draw to `data/simulation/regenerated/`. `analysis/01_simulation.ipynb` checks that
the deposited data follow the rules of this notebook.
