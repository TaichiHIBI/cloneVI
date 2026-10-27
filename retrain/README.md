# Training the fits of the paper

Run from the repository root with the data in `data/` (see `data/README.md`); results are written to
`output/v2/<name>/<YYYYMM>/<DD>/<HHMMSS>/`, next to the deposited fits. Training on a GPU is not bit-identical between
machines or library versions; the deposited fits (`output/README.md`) are the ones analysed in the paper.

| Fit (paper) | Script | Deposited fit |
|---|---|---|
| Simulation (Fig. 2) | `train_simulation.py` | `output/v2/simulation/202610/11/080706` |
| MAESTER (Fig. 3) | `train_maester.py` (input prepared by `prep_maester.py`) | `output/v2/maester/202610/11/085606` |
| CLL (Fig. 4) | `train_cll.py` | `output/v2/cll/202610/11/080706` |
| Ovarian cohort (Fig. 5) | `train_ovarian.py bkg` | `output/v2/ovarian_bkg/202610/11/080711` |
| Ovarian cohort, no decomposition of the allele frequency (`analysis/04a_ovarian.ipynb`) | `train_ovarian_nobkg.py` (the later stages from the first stage of the fit above) | `output/v2/ovarian_nobkg/202610/11/174848` |

The settings of each fit are written in its script and recorded in its `used_parameters.txt`.

## Decomposition of the allele frequency

The decomposition is used where variants shared across clones (background variants) are expected: the CLL and
ovarian fits. The simulation (no background by design) and MAESTER (26 clone-marker variants of the original study)
are fitted without it.

## Seeds

The model seeds are simulation 42, MAESTER 2024, CLL 42 and ovarian 42.
