# CloneVI

CloneVI is a deep generative model for single-cell multi-omics data in which mitochondrial DNA (mtDNA) mutations are
measured together with gene expression (scRNA-seq, e.g. MAESTER) or chromatin accessibility (e.g. mtscATAC-seq) in the
same cells. It infers a latent cell state from the expression or accessibility counts and then estimates, from that cell
state, a clonal probability for each cell and a clone-specific allele frequency for each clone by reconstructing the
observed mtDNA mutation counts. The clones it forms are mutation-defined populations associated with cell states.

This repository contains the model (`src/`) and the code behind the results of the accompanying manuscript
(`analysis/`, `retrain/`, `preprocess/`).

* To run CloneVI on your own data, start from the template [template/run_cloneVI.py](template/run_cloneVI.py).
* The inputs, the parameters and the outputs are described in [Usage](#usage).
* To reproduce the results of the paper, see [REPRODUCING.md](REPRODUCING.md).

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

To run CloneVI on your own data, start from [template/run_cloneVI.py](template/run_cloneVI.py): with a comment at
each step, it builds the two inputs from your files, checks them and trains the model.

### Required parameters
* `adata` : AnnData of the cell-state input; `X` holds raw counts (genes or peaks) of the cells. Every column enters
  the model, so the features are chosen beforehand (e.g. the 2,000 most variable genes).
* `adata_mg` : AnnData of the mtDNA mutations for the same cells in the same order (the order is not checked); `X`
  holds the read depth and `layers['all_mut']` the counts of the mutant allele at each position. Every column enters
  the model, so it holds only the positions to use.
* `mutations` : the names of the columns of `adata_mg`, in the same order.
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
* `epoch_num` : maximum number of epochs, counted over all stages (default `20000`); each stage stops early when its
  validation loss stops improving (`stage_patience_x`, `stage_patience_pi`).
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
Written to `save_dir/experiment_name/YYYYMM/DD/HHMMSS/`.
* `cloneVI_result.h5ad` : `adata_mg` of all cells with the results.
    * `obsm['z']` : latent cell state (posterior mean); `obsm['X_cloneVI_umap']` : its UMAP.
    * `obsm['pi_sorted']` : clonal probabilities of all clones (cells x clones).
    * `obsm['pi_deleted']` : clonal probabilities of the retained clones (mean probability above 0.01, or a cell with a
      probability of at least 0.2).
    * `varm['phi_sorted']`, `varm['phi_deleted']` : allele frequencies of all clones and of the retained clones
      (variants x clones).
    * `layers['reconst_VAF']` : VAF expected for each cell and variant (`pi_sorted @ phi_sorted`).
    * With the decomposition of the allele frequency: `varm['phi_0']`, `varm['phi_k']` and `uns['alpha']` hold its
      components in the clone order of `phi_sorted` (`phi_sorted = alpha * phi_0 + (1 - alpha) * phi_k`);
      `varm['phi_specific']` holds the clone-specific component sorted by the probabilities of the test cells, whose
      clone order can differ.
* `cloneVI_test_result.h5ad` : the same for the test cells (without `phi_specific`), with the clones sorted by their
  probabilities in these cells.
* `models/<loss>.pth` : the weights. At the end of each stage the weights with the lowest validation loss over the stage
  are restored, so the results and the saved model hold the same weights.
* `used_parameters.txt` : the settings, whether every stage ended before `epoch_num` (`completed`) and the best and last
  epoch of each stage.
* `split_idx.npz` : the train/validation/test split.
* `plots/` : diagnostic plots; `logs/` : TensorBoard logs (with `use_tensorboard=True`).

## Reproducing the paper

The data, the code and the order of the steps that reproduce the results of the manuscript are described in
[REPRODUCING.md](REPRODUCING.md): the deposited inputs and fits (Zenodo, https://doi.org/10.5281/zenodo.23687068),
the notebooks and R scripts behind each figure and table (`analysis/`), the training of the fits (`retrain/`) and
the preprocessing from the public raw data (`preprocess/`).

## License
MIT (see `LICENSE`).
