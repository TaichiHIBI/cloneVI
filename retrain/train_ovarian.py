"""Ovarian cohort fits of Fig. 5.

bkg: with the phi decomposition (the fit of the paper).  nobkg: without it, trained from scratch; the fit without the
decomposition compared in analysis/04a_ovarian.ipynb is made by retrain/train_ovarian_nobkg.py instead, from the first
stage of the bkg fit.
Settings: learning rate 4e-3, 150 warm-up epochs, at most 10,000 epochs, batch 1,024, donor as covariate, stage
patience 20, initial temperature 100 (KL) and 10,000 (mutation term), regularization lmd_pi = 0.01 and lmd_phi = 1 with
its temperature decreasing from 100 to 1 over 30 epochs (DecreasingTemp_adapt; the training loop does not pass it a
loss, so it is not adaptive and behaves as DecreasingTemp), seed 42.

Run from the repository root: python retrain/train_ovarian.py [bkg|nobkg]      (default bkg = phi decomposition on)
Writes output/v2/ovarian_bkg/ or output/v2/ovarian_nobkg/<date>/<time>/.
"""
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(ROOT, 'src'))
import scanpy as sc
import utils, workflow

variant = sys.argv[1] if len(sys.argv) > 1 else 'bkg'
assert variant in ('bkg', 'nobkg')
seed = 42
utils.set_seed(seed)
adata = sc.read_h5ad(f'{ROOT}/data/mitotrek/cohort_Ovarian_atac_top25.h5ad')
adata_mg = sc.read_h5ad(f'{ROOT}/data/mitotrek/cohort_Ovarian_mtsnv_perdonor_frac0.5pct.h5ad')
mutations = adata_mg.var_names
adata.obs['lsi_clusters'] = adata.obs['lsi_clusters'].astype('category')
adata.var_names_make_unique()

loss_hparams = {'kld_temp': utils.DecreasingTemp(dec_temp_steps=300, init_temp=100),
                'binom_temp': utils.DecreasingTemp(dec_temp_steps=300, init_temp=10000),
                'sp_temp': utils.DecreasingTemp_adapt(dec_temp_steps=30, init_temp=100, monitor_ratio=True),
                'lmd_pi': 1e-2, 'lmd_phi': 1}
workflow.estimate_probability(adata, adata_mg, mutations, f'{ROOT}/output', loss_hparams,
                              model_params={'use_phi_common': variant == 'bkg', 'freeze_z_after_x': True},
                              epoch_num=10000, lr=4e-3, lr_scheduler=True, batch_size=1024, num_workers=8,
                              celltype='lsi_clusters', genotype=None, use_tensorboard=True,
                              experiment_name=f'v2/ovarian_{variant}', warmup_epochs=150, sample_info='donor',
                              stage_patience_x=20, seed=seed)
