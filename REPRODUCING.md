# Reproducing the paper

How to reproduce the results of the manuscript from the deposited data and fits, and how those were made.

## Data
The processed input data, the outputs of the fits used in the paper and the intermediate files that cannot be
regenerated are deposited on Zenodo (https://doi.org/10.5281/zenodo.23687068). Unpack `cloneVI_data.tar.gz` and
`cloneVI_output.tar.gz` at the repository root, which fills `data/` and `output/`, and put `cohort_Ovarian_atac.h5seurat`
in `data/mitotrek/` (layouts, the two files to download from their sources and the checksums in
[data/README.md](data/README.md) and [output/README.md](output/README.md)). The raw data are public:
PBMC3k (10x Genomics, via Scanpy), MAESTER (GEO GSE182685), CLL mtscATAC-seq (https://osf.io/bupge/) and the ovarian
tumor mtscATAC-seq cohort (GEO GSE302113).

## Results: `analysis/`
The notebooks and R scripts in `analysis/` compute the figure panels, tables and numbers of the Results and Methods from
the deposited data and fits, except the schematics drawn by hand (Fig. 1, Supplementary Fig. S1 and Fig. 2a) and the
training settings of Supplementary Table 2 (in `retrain/`).
Each notebook says at each step what is computed and why, shows the values
it computes, and ends with a table that compares the numbers quoted in the paper with the computed ones. The notebooks
are kept with their outputs. They are run from `analysis/` (each starts with `%cd ../src`); the R scripts are run from
the repository root (`Rscript analysis/<script>.R`). Everything is written to `analysis/results/<dataset>/`: one file per
figure panel, each figure assembled from its panels (`figure<N>.pdf`, `supfig<N>.pdf`, and `.png`), and the tables behind
them as CSV files. Panels drawn for one clone are written for every clone. Assembling Fig. 4 and Supplementary Fig. S4
needs `pdftoppm` (poppler) to place the coverage plots drawn by R.

| Item | Code (in order) | Output in `analysis/results/` |
|---|---|---|
| Fig. 1, Supplementary Fig. S1 | schematics (drawn by hand) | |
| Fig. 2, Supplementary Fig. S2, Supplementary Table 3 | `01_simulation.ipynb` (the schematic of Fig. 2a, drawn by hand, is `analysis/figure_parts/fig2a_schematic.png`) | `simulation/` |
| Fig. 3, Supplementary Fig. S3 | `02_maester.ipynb` | `maester/` |
| Fig. 4, Supplementary Fig. S4 | `03a_cll.ipynb` (Fig. 4a–c, S4a–c) → `03b_cll_gene_activity.R` (gene activity) → `03c_cll_gene_activity.ipynb` (genes of each clone, Fig. 4e) → `03d_cll_coverage.R` (coverage plots, Fig. 4d and S4d) → `03e_cll_figure4.ipynb` (the two figures) | `cll/`, `cll/r_output/` |
| Fig. 5a–e, Supplementary Figs. S5–S7 and S8d, Supplementary Table 1 | `04a_ovarian.ipynb` | `ovarian/` |
| Supplementary Fig. S9 | `04b_ovarian_within_donor.ipynb` | `ovarian/` |
| Fig. 5f–h, Supplementary Fig. S8a–c | `04c_ovarian_chromvar.R` (chromVAR) → `04d_ovarian_motif.ipynb` (motif tests, for Clone 11 as in the paper and for every clone) → `04e_ovarian_motif_seurat_check.R` (the comparison with Seurat stated in Methods) | `ovarian/chromvar/`, `ovarian/motif/` |
| Fig. 5 and Supplementary Fig. S8 as figures | `04f_ovarian_figure5.ipynb` (after the four steps above) | `ovarian/` |
| Supplementary Table 2 and the other training settings | `retrain/` (and `used_parameters.txt` of each fit in `output/`) | |

`03a_cll.ipynb` must be run before the CLL R scripts, `03e_cll_figure4.ipynb` after them, `04c_ovarian_chromvar.R`
before `04d_ovarian_motif.ipynb` (or the deposited motif activity copied in its place, see
[data/README.md](data/README.md)), and `04f_ovarian_figure5.ipynb` last.
The long steps take about 4 minutes (`03b`), 75 minutes (`03d`), 15 minutes (`04b`) and one hour (`04c`, 12 cores); the
others take about a minute or less.

## Training
The fits of the paper are trained by the scripts in `retrain/` ([retrain/README.md](retrain/README.md)) with the code
in `src/`. The results above are computed from the deposited fits.

## Preprocessing
From the public raw data to the model inputs: [preprocess/README.md](preprocess/README.md). The model inputs of
MAESTER, CLL and the ovarian cohort are rebuilt and compared with the deposited ones by the notebooks
`preprocess/{MAESTER,CLL,ovarian}/build_model_inputs.ipynb`.
