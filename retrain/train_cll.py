"""CLL fit of Fig. 4.

Settings: learning rate 4e-3, 3,000 warm-up epochs, at most 20,000 epochs, batch 300, initial temperature 100 (KL) and
10,000 (mutation term), regularization lmd_pi = 0.01 and lmd_phi = 1, phi decomposition on, encoder fixed after
stage 1, seed 42.

Run from the repository root: python retrain/train_cll.py
Writes output/v2/cll/<date>/<time>/.
"""
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(ROOT, 'src'))
import scanpy as sc
import utils, workflow

seed = 42
utils.set_seed(seed)
adata = sc.read_h5ad(f'{ROOT}/data/ATAC/PT1_pos.h5ad')
adata.X = adata.raw.X   # h5seurat2h5ad で、ATACのカウントがXではなくraw.Xに保存されるため、それで置換
adata_mg = sc.read_h5ad(f'{ROOT}/data/ATAC/PT1_pos_mg.h5ad')
mutations = ['3179', '12213', '4716', '16519', '12067', '5809', '2702', '11711', '14858', '1260', '16390', '1872',
             '5140', '2587', '9950', '16052', '5835', '5274']
adata_mg = adata_mg[:, mutations]
adata.obs.rename(columns={"cluster_id": "genotype", "seurat_clusters": "lsi_clusters"}, inplace=True)
adata.obs['genotype'] = adata.obs['genotype'].astype('category')
adata.obs['lsi_clusters'] = adata.obs['lsi_clusters'].astype('category')
adata.var_names_make_unique()

loss_hparams = {'kld_temp': utils.DecreasingTemp(dec_temp_steps=300, init_temp=100),
                'binom_temp': utils.DecreasingTemp(dec_temp_steps=300, init_temp=10000),
                'sp_temp': utils.DecreasingTemp(dec_temp_steps=300, init_temp=100),
                'lmd_pi': 1e-2, 'lmd_phi': 1}
workflow.estimate_probability(adata, adata_mg, mutations, f'{ROOT}/output', loss_hparams,
                              model_params={'use_phi_common': True, 'freeze_z_after_x': True},
                              epoch_num=20000, lr=4e-3, lr_scheduler=True, celltype='lsi_clusters', genotype='genotype',
                              train_mode='train_all', use_tensorboard=True, experiment_name='v2/cll', seed=seed)
