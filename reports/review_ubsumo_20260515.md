# Review: Ub-SUMO Competition Summary

**Time**: 2026-05-15 13:02-13:15

## Approach

Following 12403 electrostatic evaluations (Ub site accuracy: 3%, SUMO site accuracy: 56%), we converged on an asymmetric competition framework: Ub as the basal modification, SUMO as the special case.

## Key findings

- **Electrostatics is not the discriminator** — Delta-E distributions overlap completely
- **Ub**: ~600 E3 ligases → universal coverage of nearly all K sites (basal channel)
- **SUMO**: UBC9 only → only modifies sites with canonical consensus motifs (special case channel)
- "Ub-SUMO competition" is better understood as SUMO intervention into the Ub pathway under specific conditions

## Written output

`ubsumo_review.docx` — Chinese review document containing:
- Abstract with core asymmetric competition conclusion
- Structural basis (beta-grasp fold, charge differences, enzymatic asymmetry)
- Data & methods: 2693 PDBs (2.2GB), DBREF parsing, electrostatic model
- Results: coverage statistics, electrostatic accuracy, asymmetry discovery
- Discussion: redefining competition, model limitations
- Closing remarks

## Data status
- `data/pdb_all/`: 2693 PDB files, 2.2GB
- `data/results/electrostatic_results.json`: full analysis results
