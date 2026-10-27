"""Simulation fit of Fig. 2.

Settings: K = 10, learning rate 1e-3, 3,000 warm-up epochs, at most 20,000 epochs, initial temperature 100 for the KL
and mutation terms, no regularization, seed 42, without the phi decomposition (the simulation has no background
variant; see retrain/README.md).

Run from the repository root: python retrain/train_simulation.py
Writes output/v2/simulation/<date>/<time>/.
"""
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(ROOT, 'src'))
import numpy as np
import scanpy as sc
import utils, workflow

D = f'{ROOT}/data/simulation'
adata = sc.read_h5ad(f'{D}/back_0.5_dropout/adata.h5ad')
adata_mg = sc.read_h5ad(f'{D}/back_0.5_dropout/adata_mg.h5ad')
raw = sc.read_10x_mtx(f'{D}/filtered_gene_bc_matrices/hg19', var_names='gene_symbols', cache=False)
raw.var_names_make_unique()
adata.X = raw[adata.obs_names, adata.var_names].X          # raw counts as the model input
adata.obs['clone'] = adata_mg.obs['clone']
adata.var_names_make_unique()

seed = 42
utils.set_seed(seed)
loss_hparams = {'kld_temp': utils.DecreasingTemp(dec_temp_steps=300, init_temp=100),
                'binom_temp': utils.DecreasingTemp(dec_temp_steps=300, init_temp=100),
                'sp_temp': utils.DecreasingTemp(dec_temp_steps=300, init_temp=100),
                'lmd_pi': 0, 'lmd_phi': 0}
workflow.estimate_probability(adata, adata_mg, np.arange(adata_mg.n_vars).tolist(), f'{ROOT}/output',
                              loss_hparams, model_params={'num_enc_z_layers': 2, 'num_dec_z_layers': 2, 'k': 10,
                                                          'use_phi_common': False},
                              epoch_num=20000, lr=1e-3, celltype='louvain', genotype='clone', use_tensorboard=True,
                              experiment_name='v2/simulation', seed=seed)
