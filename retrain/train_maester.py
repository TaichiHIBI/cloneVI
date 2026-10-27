"""MAESTER fit of Fig. 3.

Settings: seed 2024, learning rate 4e-3, 1,000 warm-up epochs, at most 20,000 epochs, initial temperature 100 (KL) and
10,000 (mutation term), no regularization, without the phi decomposition (the 26 variants are clone markers of the
original study, no background variant is expected; see retrain/README.md).

Run from the repository root: python retrain/train_maester.py
Writes output/v2/maester/<date>/<time>/.
"""
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
os.environ.setdefault('CLONEVI_ROOT', ROOT)              # prep_maester reads data from here
sys.path.insert(0, os.path.join(ROOT, 'src'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils, workflow, prep_maester

seed = 2024
utils.set_seed(seed)
adata, adata_mg = prep_maester.load()
loss_hparams = {'kld_temp': utils.DecreasingTemp(dec_temp_steps=300, init_temp=100),
                'binom_temp': utils.DecreasingTemp(dec_temp_steps=300, init_temp=10000),
                'sp_temp': utils.DecreasingTemp(dec_temp_steps=300, init_temp=100),
                'lmd_pi': 0, 'lmd_phi': 0}
workflow.estimate_probability(adata, adata_mg, prep_maester.MUTATIONS, f'{ROOT}/output', loss_hparams,
                              model_params={'use_phi_common': False},
                              epoch_num=20000, lr=4e-3, lr_scheduler=True, celltype='CellType', genotype='MT_clone',
                              train_mode='train_all', use_tensorboard=True,
                              experiment_name='v2/maester', warmup_epochs=1000, seed=seed)
