# Hydrophobic effect analysis report

**Time**: 2026-05-15 17:56-18:02

## Approach

Three-step SASA-based hydrophobic analysis to test whether hydrophobic surface properties differentiate Ub from SUMO modification sites.

## Step 1: Single structure SASA (CA neighborhood method)

| Metric | Ub (1UBQ) | SUMO-1 (1A5R) | Diff |
|--------|-----------|---------------|------|
| Total SASA | 4403 A² | 6487 A² | +47% |
| Hydrophobic SASA | 1480 A² | 1880 A² | +27% |
| Hydrophobic fraction | 33.6% | 29.0% | -14% |
| Hydrophobic/residue | 19.47 A² | 18.26 A² | -6% |

## Step 2: Bulk baseline statistics

- **Ub baseline** (16 structures): hydrophobic/residue = 17.52 ± 2.36 A²
- **SUMO** (8 structures, all isoforms): hydrophobic/residue = 17.21 ± 1.94 A²
- **No SUMO structure had z-score > 2** — no statistical significance

## Step 3: Hydrophobic patch analysis

- Ub main patch: 4 residues (LEU×3 + GLY), 331 A², no aromatics
- SUMO-1 main patch: 3 residues (GLY×2 + VAL), 234 A², PHE(100 A²) in secondary patch
- SIM-binding score: Ub=70, SUMO-1=53

## Conclusion

The hydrophobic hypothesis is **falsified** — no significant difference in hydrophobic surface area, density, or patch continuity between Ub and SUMO. The only observable difference is aromatic patch position specificity, not total amount.

## Data files

- `ubsumo_em_simulator/data/results/hydrophobic_sasa_baselines.json`
- `ubsumo_em_simulator/data/results/hydrophobic_patches.json`
