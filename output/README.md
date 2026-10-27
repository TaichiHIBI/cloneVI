# output/

The fits of the Zenodo record (https://doi.org/10.5281/zenodo.23687068) go here: unpack `cloneVI_output.tar.gz` at the
repository root (`tar -xzf cloneVI_output.tar.gz`). Each fit is a folder `<experiment>/<YYYYMM>/<DD>/<HHMMSS>/` with
`cloneVI_result.h5ad`, `cloneVI_test_result.h5ad`, `used_parameters.txt`, `split_idx.npz`, `models/` (the weights of
the results), `plots/` and `logs/` (TensorBoard).

| Folder | Fit | Trained by |
|---|---|---|
| `v2/simulation/202610/11/080706` | simulation (Fig. 2, Supplementary Fig. S2) | `retrain/train_simulation.py` |
| `v2/maester/202610/11/085606` | MAESTER, seed 2024 (Fig. 3, Supplementary Fig. S3) | `retrain/train_maester.py` |
| `v2/cll/202610/11/080706` | CLL (Fig. 4, Supplementary Fig. S4) | `retrain/train_cll.py` |
| `v2/ovarian_bkg/202610/11/080711` | ovarian cohort (Fig. 5, Supplementary Figs. S5–S9, Supplementary Table 1) | `retrain/train_ovarian.py bkg` |
| `v2/ovarian_nobkg/202610/11/174848` | ovarian cohort without the decomposition of the allele frequency, from the first stage of the fit above (`qz.npz`; compared in `analysis/04a_ovarian.ipynb`) | `retrain/train_ovarian_nobkg.py` |
