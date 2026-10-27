# From the public raw data to the model inputs

The scripts are the ones used for the study; only their paths were changed. They expect the raw data in folders at the
repository root (`MAESTER/`, `mtscATAC-seq/`, `GSE302113/`) and were run on different machines (the ovarian raw data on a
cluster with a job scheduler). The files they produce are deposited on Zenodo (`data/`), so the results do not require
rerunning them. Steps whose code was not preserved are marked; for these the deposited files are the reference.

## MAESTER (GEO GSE182685; Miller et al. 2022)
Place `BPDCN712_Seurat_with_TCR.rds` and `TenX_BPDCN712_All_mr3_maegatk.rds` from GEO in `MAESTER/clonalhemato/`.
1. `MAESTER/output.R` — writes the maegatk allele counts as `.mtx` files (`MAESTER/rds_output/`) and converts the Seurat
   object to `rna_exp.h5seurat`/`.h5ad` with SeuratDisk. The deposited `data/h5ad_files_RNA/rna_exp.h5ad` is the reference:
   its `MT_clone_num` column holds numbers, whereas this script would store it as true/false, so the file was made with a
   slightly different version of these lines that was not kept.
2. `MAESTER/MAESTER_RNAmgv2h5ad.ipynb` (cells 1–20) — `.mtx` files → `rna_mgv_counts.h5ad` (the output is identical to
   the deposited file).
The 26 variants used are those analysed in the original study (`retrain/prep_maester.py`).

## CLL mtscATAC-seq (Lareau et al. 2021; https://osf.io/bupge/)
Place `mtscATACpaper_large_data_files/` (source data, cellranger and mgatk outputs) and the study's reproducibility
repository `mtscATACpaper_reproducibility/` in `mtscATAC-seq/`.
1. `CLL/get_filtering_num.R` — lists the cells to remove (`CLL/row_del/`, included).
2. `CLL/signac_pp.R` — Signac: peak-by-cell counts, cell and peak filters, top 25% of peaks → `PT1_pos_pp`.
3. `CLL/rds2h5ad.ipynb`, `CLL/rds_setup_h5ad.ipynb`, `CLL/filtering_call.ipynb` — mgatk allele counts → `PT1_pos_mg.h5ad`.
   The step that exported the mgatk `.rds` object to the `.mtx` files read here was not preserved, and which of these
   notebooks wrote the deposited file is not certain.
4. `CLL/clone_definition.ipynb` — `PT1_pos_pp` and the clone table of the source study → `PT1_pos.h5ad`.

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
