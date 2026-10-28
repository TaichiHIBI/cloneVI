# CloneVI

CloneVI is a deep generative model for single-cell multi-omics data in which mitochondrial DNA (mtDNA) mutations are
measured together with gene expression (scRNA-seq, e.g. MAESTER) or chromatin accessibility (e.g. mtscATAC-seq) in the
same cells. It infers a latent cell state from the expression or accessibility counts and then estimates, from that cell
state, a clonal probability for each cell and a clone-specific allele frequency for each clone by reconstructing the
observed mtDNA mutation counts. The clones it forms are mutation-defined populations associated with cell states.

This repository contains the model (`src/`) and the code behind the results of the accompanying manuscript
(`analysis/`, `retrain/`, `preprocess/`); see [Reproducing the paper](#reproducing-the-paper).

## Installation

```bash
git clone https://github.com/TaichiHIBI/cloneVI.git
conda create -n cloneVI python=3.11
conda activate cloneVI
pip install -r environment/requirements.txt      # PyTorch: choose the build for your CUDA version
```

The versions used for the paper are listed in `environment/` (Python 3.11, PyTorch 2.2.1 with CUDA 12.1, scanpy 1.9.4,
anndata 0.9.2). The analyses of chromatin accessibility (coverage plots, gene activity, chromVAR) use R 4.1.1 with
Signac 1.4.0 and Seurat 4.0.5 (`environment/R_packages.txt`).

## Usage

```python
import sys
sys.path.insert(0, 'cloneVI/src')
import utils, workflow

loss_hparams = {'kld_temp': utils.DecreasingTemp(dec_temp_steps=300, init_temp=100),
                'binom_temp': utils.DecreasingTemp(dec_temp_steps=300, init_temp=10000),
                'sp_temp': utils.DecreasingTemp(dec_temp_steps=300, init_temp=100),
                'lmd_pi': 1e-2, 'lmd_phi': 1}
adata, adata_mg = workflow.estimate_probability(adata, adata_mg, mutations, save_dir, loss_hparams,
                                                model_params={'use_phi_common': True})
```

### Required parameters
* `adata` : AnnData of the cell-state input; `X` holds raw counts (genes or peaks) of the cells.
* `adata_mg` : AnnData of the mtDNA mutations for the same cells in the same order; `X` holds the read depth and
  `layers['all_mut']` the counts of the mutant allele at each position.
* `mutations` : positions (columns of `adata_mg`) to use.
* `save_dir` : output directory; results are written to `save_dir/experiment_name/YYYYMM/DD/HHMMSS/`.
* `loss_hparams` : temperatures and regularization weights of the loss.
    * `kld_temp`, `binom_temp`, `sp_temp` : temperatures of the KL divergence, the mutation likelihood and the
      regularization term, e.g. `utils.DecreasingTemp(dec_temp_steps=300, init_temp=100)`; each decreases to 1.
    * `lmd_pi`, `lmd_phi` : weights of the group-lasso regularization on the clonal probabilities and the allele
      frequencies. The regularization stage is run only when one of them is non-zero.

### Optional parameters
* `model_params` : encoder and decoder.
    * `xz_dim` : dimension of the latent cell state (default `10`, increased automatically for large inputs).
    * `enc_z_h_dim`, `dec_z_h_dim` : hidden units (default `50`, increased automatically for large inputs).
    * `num_enc_z_layers`, `num_dec_z_layers` : hidden layers (default `2`).
    * `k` : initial number of clones (default `20`).
    * `use_phi_common` : decompose the allele frequency of each clone into a component shared by all clones and a
      clone-specific component (default `False`; the CLL and ovarian fits of the paper used `True`).
    * `freeze_z_after_x` : keep the encoder fixed after the first stage (default `True`).
* `epoch_num` : maximum epochs per stage (default `20000`); each stage stops early when its validation loss stops improving
  (`stage_patience_x`, `stage_patience_pi`).
* `lr` : learning rate (default `4e-3`), with `ReduceLROnPlateau` when `lr_scheduler=True` (`lr_patience`, `lr_factor`).
* `warmup_epochs` : learning-rate warm-up of the first stage (default `3000`).
* `val_ratio`, `test_ratio` : fractions of cells held out for validation and test (defaults `0.05`, `0.1`).
* `batch_size`, `num_workers` : mini-batches (defaults `300`, `2`).
* `n_neighbors`, `min_dist` : UMAP of the latent cell state (defaults `30`, `0.5`).
* `celltype`, `genotype` : columns of `adata.obs` used to colour the diagnostic plots.
* `train_mode` : `train_all` (all stages, default), `train_x` (cell state only) or `train_pi` (clonal probabilities from a
  fixed cell state).
* `sample_info` : column of `adata.obs` with the sample of origin (e.g. donor), given to the encoder and decoders.
* `qz` : a given q(z|x) of the cells (posterior means and standard deviations, e.g. from `workflow.qz_from_checkpoint`);
  the later stages then start from it without training the first stage (`qz_source` is recorded in
  `used_parameters.txt`).
* `use_tensorboard`, `experiment_name`, `seed` (default `42`; the train/validation/test split always uses 42).

### Output
`cloneVI_result.h5ad` (all cells) and `cloneVI_test_result.h5ad` (test cells): `obsm['z']` (latent cell state, posterior
mean), `obsm['pi_sorted']` and `obsm['pi_deleted']` (clonal probabilities of all clones and of the retained clones),
`varm['phi_sorted']` and `varm['phi_deleted']` (allele frequencies), `obsm['X_cloneVI_umap']`. At the end of each stage
the weights with the lowest validation loss over the stage are restored, so the results and the saved model
(`models/<loss>.pth`) hold the same weights. `used_parameters.txt` records the settings, whether every stage ended before
`epoch_num` (`completed`) and the best and last epoch of each stage; `split_idx.npz` holds the train/validation/test
split. With the decomposition of the allele frequency, `varm['phi_0']`, `varm['phi_k']`
and `uns['alpha']` hold its components in the clone order of `phi_sorted` (`phi_sorted = alpha * phi_0 + (1 - alpha) *
phi_k`); `varm['phi_specific']` holds the clone-specific component sorted by the probabilities of the test cells, whose
clone order can differ.

## Reproducing the paper

### Data
The processed input data, the outputs of the fits used in the paper and the intermediate files that cannot be
regenerated are deposited on Zenodo (https://doi.org/10.5281/zenodo.23687068). Unpack `cloneVI_data.tar.gz` and
`cloneVI_output.tar.gz` at the repository root, which fills `data/` and `output/`, and put `cohort_Ovarian_atac.h5seurat`
in `data/mitotrek/` (layouts, the two files to download from their sources and the checksums in
[data/README.md](data/README.md) and [output/README.md](output/README.md)). The raw data are public:
PBMC3k (10x Genomics, via Scanpy), MAESTER (GEO GSE182685), CLL mtscATAC-seq (https://osf.io/bupge/) and the ovarian
tumor mtscATAC-seq cohort (GEO GSE302113).

### Results: `analysis/`
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

### Training
The fits of the paper are trained by the scripts in `retrain/` ([retrain/README.md](retrain/README.md)) with the code
in `src/`. The results above are computed from the deposited fits.

### Preprocessing
From the public raw data to the model inputs: [preprocess/README.md](preprocess/README.md). The model inputs of
MAESTER, CLL and the ovarian cohort are rebuilt and compared with the deposited ones by the notebooks
`preprocess/{MAESTER,CLL,ovarian}/build_model_inputs.ipynb`.

## License
MIT (see `LICENSE`).
