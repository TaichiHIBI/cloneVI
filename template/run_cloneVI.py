"""Template for running CloneVI on your own data.

Copy this file, replace every TODO and run it (python run_cloneVI.py). It builds the two inputs of the model from your
files, checks them and trains the model, which writes the clonal probabilities and the allele frequencies of the
clones to cloneVI_result.h5ad.

CloneVI needs, for the same cells:
  1. the counts of the cell-state modality, genes (scRNA-seq) or peaks (scATAC-seq): raw counts, cells x features;
  2. for each mtDNA variant used, the read depth at its position and the reads with the mutant allele: cells x variants
     (e.g. from mgatk or maegatk).
Any file format will do and the names of the files do not matter. The scripts of the paper's fits are worked examples:
retrain/prep_maester.py with retrain/train_maester.py (scRNA-seq), retrain/train_cll.py (mtscATAC-seq) and
retrain/train_ovarian.py (cells of several donors).
"""
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

# the src/ folder of this repository; change it if you copy this file elsewhere
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
import utils, workflow

SEED = 42
SAVE_DIR = 'TODO_output'       # results go to SAVE_DIR/EXPERIMENT/YYYYMM/DD/HHMMSS/
EXPERIMENT = 'my_data'
utils.set_seed(SEED)

# ---- 1. Cell-state counts ---------------------------------------------------------------------------------------------
# adata: cells x features with the raw counts in adata.X (sparse or dense). The model normalizes the counts itself, so
# do not give it normalized or log-transformed values.
adata = sc.read_h5ad('TODO_counts.h5ad')     # or sc.read_10x_h5('filtered_feature_bc_matrix.h5'), sc.read_10x_mtx(...)
# if the raw counts are not in adata.X, put them there, e.g.
# adata = adata.raw.to_adata()
# adata.X = adata.layers['counts']
adata.var_names_make_unique()
# filter the cells here with your usual quality control

# Every column of adata.X enters the model, so choose the features here. More features and more cells make the training
# slower: the paper's fits, run several at a time on one GPU, took about 30 minutes for 9,346 cells x 2,000 genes
# (MAESTER) and about 5 hours for 5,610 cells x 107,497 peaks (CLL).
MODALITY = 'rna'               # 'rna' or 'atac'
if MODALITY == 'rna':
    # as for the MAESTER fit: cells with at least 200 genes, genes seen in at least 3 cells, the 2,000 most variable genes
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    lognorm = adata.copy()                   # normalized only to choose the genes; the model gets the raw counts
    sc.pp.normalize_total(lognorm, target_sum=1e4)
    sc.pp.log1p(lognorm)
    sc.pp.highly_variable_genes(lognorm, n_top_genes=2000)
    adata = adata[:, lognorm.var['highly_variable']].copy()
else:
    # as for the CLL fit: the peaks observed in at least 10 cells (the ovarian fit used the top 25% of the peaks)
    sc.pp.filter_genes(adata, min_cells=10)

X = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
assert np.all(X >= 0) and np.allclose(X, np.round(X)), 'adata.X must hold raw counts (non-negative integers)'

# ---- 2. mtDNA variant counts ------------------------------------------------------------------------------------------
# Two tables, cells x variants, with the cell barcodes as the index and the variants as the columns:
#   depth:  reads covering the position of the variant
#   mutant: reads with the mutant allele at that position (at most depth)
# preprocess/ shows how the paper's tables were made from the outputs of maegatk (MAESTER/) and mgatk (CLL/, ovarian/).
depth = pd.read_csv('TODO_depth.csv', index_col=0)
mutant = pd.read_csv('TODO_mutant.csv', index_col=0)

# The variants to model. Every column of adata_mg enters the model, so keep only these, and choose them beforehand with
# a variant caller: the MAESTER and CLL fits used the variants of the source studies, the ovarian fit the high-confidence
# variants of Signac's IdentifyVariants called per donor (preprocess/ovarian/clonotypes_perdonor_signac.R).
MUTATIONS = ['TODO_variant_1', 'TODO_variant_2']
depth, mutant = depth[MUTATIONS], mutant[MUTATIONS]

# ---- 3. The same cells in the same order ------------------------------------------------------------------------------
# The model does not match the cells by barcode: row i of adata and row i of adata_mg must be the same cell.
cells = adata.obs_names.intersection(depth.index.astype(str))
print(f'{len(cells)} cells in both inputs ({adata.n_obs} with counts, {len(depth)} with variant counts)')
adata = adata[cells].copy()
depth.index, mutant.index = depth.index.astype(str), mutant.index.astype(str)
adata_mg = ad.AnnData(X=depth.loc[cells].to_numpy(dtype=np.float32),
                      obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=[str(m) for m in MUTATIONS]))
adata_mg.layers['all_mut'] = mutant.loc[cells].to_numpy(dtype=np.float32)

assert (adata.obs_names == adata_mg.obs_names).all()
assert np.all(adata_mg.layers['all_mut'] >= 0) and np.all(adata_mg.layers['all_mut'] <= adata_mg.X), \
    'the mutant reads must lie between 0 and the depth'

# ---- 4. Training ------------------------------------------------------------------------------------------------------
loss_hparams = {
    # temperatures of the KL divergence, the mutation likelihood and the regularization, each decreasing to 1
    'kld_temp': utils.DecreasingTemp(dec_temp_steps=300, init_temp=100),
    'binom_temp': utils.DecreasingTemp(dec_temp_steps=300, init_temp=10000),
    'sp_temp': utils.DecreasingTemp(dec_temp_steps=300, init_temp=100),
    # weights of the group-lasso regularization on the clonal probabilities and the allele frequencies: the CLL and
    # ovarian fits used 1e-2 and 1; with both 0 the regularization stage is skipped (the MAESTER fit)
    'lmd_pi': 1e-2,
    'lmd_phi': 1,
}
model_params = {
    'k': 20,                   # initial number of clones
    # decompose the allele frequency of each clone into a part shared by all clones and a clone-specific part: True
    # where variants shared across clones (background variants) are expected (the CLL and ovarian fits), False for
    # variants known to mark clones (the MAESTER fit)
    'use_phi_common': True,
}
adata, adata_mg = workflow.estimate_probability(
    adata, adata_mg, list(adata_mg.var_names), SAVE_DIR, loss_hparams, model_params=model_params,
    experiment_name=EXPERIMENT,
    celltype=None,             # a column of adata.obs to colour the diagnostic plots (e.g. cell types), if any
    genotype=None,             # likewise (e.g. clones called by another method)
    sample_info=None,          # a column of adata.obs with the sample (e.g. donor) when the cells come from several
    seed=SEED)
# Other settings (number of epochs, learning rate, warm-up, batch size, ...): see Optional parameters in README.md.

# ---- 5. Results -------------------------------------------------------------------------------------------------------
# The run folder holds cloneVI_result.h5ad (all cells; the adata_mg returned above), cloneVI_test_result.h5ad (test
# cells), the settings, the model and diagnostic plots (Output in README.md). In cloneVI_result.h5ad:
#   obsm['pi_deleted']   cells x clones: clonal probability of the clones kept (mean probability above 0.01 or a cell
#                        with a probability of at least 0.2); obsm['pi_sorted'] has all clones
#   varm['phi_deleted']  variants x clones: allele frequency of the clones kept; varm['phi_sorted'] has all clones
#   obsm['z'], obsm['X_cloneVI_umap']   latent cell state and its UMAP
# Read it later with sc.read_h5ad('<run folder>/cloneVI_result.h5ad').
print(f"{adata_mg.obsm['pi_deleted'].shape[1]} clones kept")
