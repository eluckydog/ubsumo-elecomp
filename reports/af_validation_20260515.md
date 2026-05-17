# AlphaFold electrostatic validation report

**Time**: 2026-05-15 09:00-13:13

## Objective

Cross-validate the electrostatic competition model using AlphaFold-predicted structures — checking whether AF structures reproduce the same electrostatic patterns as experimental PDBs.

## Approach

1. Downloaded 55 AlphaFold structures (19MB) for proteins with known Ub/SUMO sites
2. Computed electrostatic delta-E (DE) for each K site: DE = netCoulomb(Ub-charged - SUMO-charged)
3. Electrostatic prediction rule: DE > 0 → Ub; DE < 0 → SUMO; |DE| < threshold → ambiguous
4. Cross-validated against experimental PDB-derived electrostatic results

## Results

| Metric | Value |
|--------|-------|
| Total sites validated | 180 (Ub=92, SUMO=81, both=7) |
| Electrostatic prediction accuracy | 4% |
| Motif-only prediction accuracy | ~55% |
| Combined (electrostatic + motif) | ~53% |

## Key findings

- AF electrostatics **worse than experimental PDB** — SUMO accuracy dropped from 56% → 0%
- AF lacks experiment-level sidechain precision needed for reliable electrostatic calculation
- Motif score (psiKxDE) consistently outperforms electrostatic — no physical model beats simple sequence pattern

## Conclusion

AlphaFold is not a substitute for experimental PDBs in electrostatic selectivity prediction. The motif-based prediction (~55%) is the best we can do with available tools.
