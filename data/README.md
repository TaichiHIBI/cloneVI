# data/

The files of the Zenodo record (https://doi.org/10.5281/zenodo.23687068) go here. At the repository root:

```bash
tar -xzf cloneVI_data.tar.gz                              # everything below marked "Zenodo"
mv cohort_Ovarian_atac.h5seurat data/mitotrek/            # a separate file of the record, needed only by 04c
sha256sum -c SHA256SUMS.txt                               # every file of the record, after unpacking
```

The files marked "public" are not in the record; download them from the source.

| Path | Content | Source | Made by | Used by |
|---|---|---|---|---|
| `simulation/filtered_gene_bc_matrices/hg19/` | PBMC3k raw counts (10x Genomics) | public (https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz, unpacked in `simulation/`) | | training, `analysis/01_simulation.ipynb` |
| `simulation/back_0.5_dropout/{adata.h5ad, adata_mg.h5ad, cluster_probs.json}` | simulated data of Fig. 2 (the clone draws were not seeded, so the deposited data are the reference) | Zenodo | `preprocess/simulation/simulation_datasets.ipynb` | training, `analysis/01_simulation.ipynb` |
| `h5ad_files_RNA/rna_exp.h5ad` | MAESTER expression | Zenodo | `preprocess/MAESTER/` | training, `analysis/02_maester.ipynb` |
| `h5ad_files_RNA/rna_mgv_counts.h5ad` | MAESTER mtDNA allele counts | Zenodo | `preprocess/MAESTER/` | training |
| `ATAC/PT1_pos.h5ad` | CLL accessibility (raw counts of the 107,497 peaks observed in at least 10 cells; LSI clusters) | Zenodo | `preprocess/CLL/` | training, `analysis/03a_cll.ipynb`, `analysis/03e_cll_figure4.ipynb` |
| `ATAC/PT1_pos_mg.h5ad` | CLL mtDNA allele counts | Zenodo | `preprocess/CLL/` | training |
| `ATAC/flagment/CLL_PT1_CD19pos_v12-mtMask_fragments.tsv.gz(.tbi)` | fragments of the source study | public (https://osf.io/bupge/, `source/cellranger_output/`; make the index with `tabix -p bed`) | | `analysis/03b_cll_gene_activity.R`, `analysis/03d_cll_coverage.R` |
| `mitotrek/cohort_Ovarian_atac_top25.h5ad` | ovarian cohort accessibility (top 25% of peaks, raw counts) | Zenodo | `preprocess/ovarian/build_model_inputs.ipynb` | training |
| `mitotrek/cohort_Ovarian_mtsnv_perdonor_frac0.5pct.h5ad` | ovarian cohort mtDNA allele counts, 64 variants | Zenodo | `preprocess/ovarian/build_model_inputs.ipynb` | training |
| `mitotrek/cohort_Ovarian_atac.h5seurat` | ovarian cohort accessibility, all peaks (Signac object) | Zenodo (separate file) | `preprocess/ovarian/build_atac_cohort_signac.R` | `analysis/04c_ovarian_chromvar.R` |
| `mitotrek/chromvar/{chromvar_z.h5, motifs.csv, peaks_used.txt, session_info.txt}` | motif activity of every cell, as computed for the paper | Zenodo | `analysis/04c_ovarian_chromvar.R` | `analysis/04d_ovarian_motif.ipynb` (see below) |
| `mitotrek/paper_supp/tableS2_ovarian_tumor_labels.csv` | cell-type labels of the source study, matched to the cells | Zenodo | `preprocess/ovarian/build_model_inputs.ipynb` from the source study's Supplementary Table S2 | `analysis/04a_ovarian.ipynb`, `analysis/04d_ovarian_motif.ipynb`, `analysis/04f_ovarian_figure5.ipynb` |

`analysis/04d_ovarian_motif.ipynb` reads the motif activity from `analysis/results/ovarian/chromvar/`, where
`analysis/04c_ovarian_chromvar.R` writes it (about one hour). To run the motif tests without R, copy the deposited one
there instead: `mkdir -p analysis/results/ovarian && cp -r data/mitotrek/chromvar analysis/results/ovarian/`.

The CLL gene activity and the inputs of the CLL R scripts are computed by `analysis/` (`analysis/results/cll/`).
The notebooks `preprocess/{MAESTER,CLL,ovarian}/build_model_inputs.ipynb` write the model inputs to
`h5ad_files_RNA/rebuilt/`, `ATAC/rebuilt/` and `mitotrek/rebuilt/` and compare them with the deposited ones.
