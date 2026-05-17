# Experimental PDB integration report

**Time**: 2026-05-15 08:35-08:54

## Objective

Replace AlphaFold structures with experimental PDB data for electrostatic validation.

## Approach

1. Built UniProt→PDB mapping via RCSB GraphQL API
2. 509KB mapping file containing 2778 PDB entries across 54 proteins
3. Downloaded 2693 experimental PDB files (2.2GB)
4. Extracted residue-level DBREF records for chain mapping
5. Cross-referenced K site positions between UniProt sequence and PDB numbering

## Results

- 44 PDBs (41MB) selectively downloaded first for quick validation
- Full 2693 PDBs (2.2GB) downloaded afterwards
- Experimental PDB electrostatic validation: **3/84 Ub sites correct (4%)** — electrostatic hypothesis falsified

## Key lesson

Experimental PDB electrostatics gave no better results than AlphaFold — the failure is in the physical model, not in the structure quality.
