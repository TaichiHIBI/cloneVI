"""Ovarian cohort fit without the phi decomposition (compared in analysis/04a_ovarian.ipynb), from the stage 1 of the
fit with the decomposition.

Stage 1 (the cell state) does not involve phi, so the fit without the decomposition starts from the q(z|x) of the
finished fit with it (retrain/train_ovarian.py bkg) and trains stages 2-3 only.  The two fits then share the cell
state, and their comparison isolates the decomposition.  Other settings as retrain/train_ovarian.py; seed 42.

Run from the repository root: python retrain/train_ovarian_nobkg.py [bkg run directory]
(default: the latest run under output/v2/ovarian_bkg).  Writes output/v2/ovarian_nobkg/<date>/<time>/.
"""
import glob
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(ROOT, 'src'))
import scanpy as sc
import utils, workflow

bkg_dir = sys.argv[1] if len(sys.argv) > 1 else sorted(glob.glob(f'{ROOT}/output/v2/ovarian_bkg/*/*/*'))[-1]
assert 'completed: True' in open(f'{bkg_dir}/used_parameters.txt').read(), f'{bkg_dir} has not finished'
pth = glob.glob(f'{bkg_dir}/models/*.pth')
assert len(pth) == 1, pth

seed = 42
utils.set_seed(seed)
adata = sc.read_h5ad(f'{ROOT}/data/mitotrek/cohort_Ovarian_atac_top25.h5ad')
adata_mg = sc.read_h5ad(f'{ROOT}/data/mitotrek/cohort_Ovarian_mtsnv_perdonor_frac0.5pct.h5ad')
mutations = adata_mg.var_names
adata.obs['lsi_clusters'] = adata.obs['lsi_clusters'].astype('category')
adata.var_names_make_unique()

# q(z|x) from the encoder of the fit with the decomposition (its model_params, so that the weights load)
qz = workflow.qz_from_checkpoint(adata, pth[0], mutations, model_params={'use_phi_common': True, 'freeze_z_after_x': True},
                                 sample_info='donor')

loss_hparams = {'kld_temp': utils.DecreasingTemp(dec_temp_steps=300, init_temp=100),
                'binom_temp': utils.DecreasingTemp(dec_temp_steps=300, init_temp=10000),
                'sp_temp': utils.DecreasingTemp_adapt(dec_temp_steps=30, init_temp=100, monitor_ratio=True),
                'lmd_pi': 1e-2, 'lmd_phi': 1}
workflow.estimate_probability(adata, adata_mg, mutations, f'{ROOT}/output', loss_hparams,
                              model_params={'use_phi_common': False, 'freeze_z_after_x': True},
                              epoch_num=10000, lr=4e-3, lr_scheduler=True, batch_size=1024, num_workers=8,
                              celltype='lsi_clusters', genotype=None, use_tensorboard=True,
                              experiment_name='v2/ovarian_nobkg', warmup_epochs=150, sample_info='donor',
                              stage_patience_x=20, seed=seed, qz=qz, qz_source=pth[0])
