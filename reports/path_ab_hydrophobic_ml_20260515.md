# Path A+B: Hydrophobic microenvironment features + ML classifier

**Time**: 2026-05-15 18:27-18:33

## Path A: Feature extraction

Extracted 26-dim local hydrophobic microenvironment features for 173 K sites across 55 AlphaFold structures:

- Sequence window features: hydrophobic/polar/aromatic/charge density (±10 residues)
- 3D spatial neighbor features: residue type distribution within 8A CA distance
- Extended motif features: psiKxDE signature + charge balance within ±5 range
- Local exposure estimation from neighbor density

## Path B: ML classification

8 feature sets × 3 classifiers × 3-fold cross-validation (24 combinations total):

| Feature set | Best classifier | Balanced accuracy | AUC |
|-------------|----------------|------------------|-----|
| motif_only | RF | **75.6%** | 0.764 |
| motif+window | LR | 69.0% | **0.810** |
| motif+spatial | GB | 72.0% | 0.779 |
| all (20 feats) | LR | 71.5% | **0.818** |
| spatial_only | LR | 71.0% | 0.768 |

### Feature importance (RF, all features)

1. motif_extended — 13%
2. motif_score — 10%
3. n_polar_window — 7%
4. spatial_hydro_density — 6%
5. neg_nearby_5 — 6%

### Comparison to original

| Method | Accuracy |
|--------|---------|
| Original motif-based prediction | ~53% |
| Best ML (RF, motif_only) | 75.6% |
| Best ML (LR, all features) | 71.5% (AUC=0.818) |

## Conclusion

Hydrophobic microenvironment features are **informative but not complementary** to motif features. The bottleneck in Ub vs SUMO discrimination lies outside our computational feature space — molecular docking or MD would be needed for further progress.
