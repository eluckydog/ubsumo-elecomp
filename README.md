# ubsumo-elecomp

**Ub/SUMO modification site electrostatic competition analysis · We failed.**

This repo documents our attempt to computationally predict ubiquitin (Ub) vs SUMO competition at lysine (K) sites. The bottom line: **it didn't work.**

All three physical hypotheses were falsified:

- **Electrostatic model** — Experimental PDB + AlphaFold electrostatic calculations (84 validation sites). Accuracy: 4%. Delta-E values completely overlapped.
- **Hydrophobic effect** — SASA analysis (16 Ub + 8 SUMO structures). Hydrophobic per-residue SASA: Ub=17.52±2.36, SUMO=17.21±1.94. No significant difference.
- **ML classifier** — 26-dim hydrophobic microenvironment features × 4 models × 24 combinations. Best balanced accuracy: 71.5% — no better than simple motif scoring.

## What we did

Through intensive work:
- Built UniProt→PDB mapping for 54 proteins (509KB, 2778 PDB entries)
- Downloaded 2693 PDB files (2.2GB)
- Extracted 55 AlphaFold predicted structures (19MB)
- Extracted 191 K sites across 90 proteins from DeepPCT dataset
- Three complete validation rounds: electrostatic → motif analysis → hydrophobic effect
- Built 26-dim local hydrophobic microenvironment features
- Ran 24 ML model combinations (3-fold cross-validation)

## Repository structure

```
ubsumo-elecomp/
├── code/            # Python source code
│   ├── core/        — EM force kernel, competition model
│   ├── analysis/    — Electrostatics, binding, pipelines
│   ├── io/          — PDB reader, trajectory
│   ├── tools/       — Batch validation, RCSB API
│   └── validation/  — Benchmarks (PCNA etc.)
├── data/
│   ├── datasets/    — Clean datasets & method references
│   │   └── ubsumo_191k_sites.csv  — 180 K sites (92 Ub, 81 SUMO, 7 both)
│   ├── pdb_sample/  — 9 curated PDBs (3.8MB)
│   │   — Ub: 1UBQ; SUMO-1: 1A5R; SUMO-2/3: 1WM3/2BFW
│   │   — Validation targets: p53, PCNA, DNMT1, ER-alpha, PML
│   ├── alphafold/   — AF2 structure listing (PDB files gitignored)
│   └── results/     — All intermediate JSON results
├── reports/         — Reports & review
└── docs/            — TBD
```

## Why we stopped

What comes next — molecular docking, molecular dynamics, FreeSASA — requires specialized software and computing power beyond what a human collaborator with a laptop and an AI assistant can provide.

All code and intermediate results are here. Take them, reproduce, or push further.

The human collaborator hopes Ub-SUMO secrets will be uncovered sooner, and Alzheimer's cured sooner.

---

MIT
