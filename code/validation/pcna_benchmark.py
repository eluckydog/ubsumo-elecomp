"""
PCNA benchmark: validate Ub/SUMO competition predictions against experimental data.

PCNA (Proliferating Cell Nuclear Antigen) is a well-characterized substrate
for both Ub and SUMO modification.

Key site: K164 — known competition site (both Ub and SUMO observed).

NOTE: K127 is GLY in PCNA (1VYM), not lysine. Earlier reports of "K127 Ub-only"
likely refer to sequence-based predictions (DeepPCT), not structural positions.
This benchmark focuses on the structurally verified K164 site.
"""

import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Optional

# PCNA site data based on structural analysis of PDB 1VYM
@dataclass
class PCNASite:
    """Known PCNA modification site."""
    site_id: int
    position: np.ndarray  # CA position from PDB
    chain: str
    charge: float  # Lysine charge
    solvent_accessibility: float  # 0-1
    experimental: str  # "competition" or "ambiguous"
    expected_preference: str  # "SUMO" (under no DNA damage conditions)


def load_pcna_sites_from_pdb(pdb_path: str = "data/pdb/pcna_1VYM.pdb") -> List[PCNASite]:
    """Load K164 from PCNA PDB structure.
    
    Parameters
    ----------
    pdb_path : str
        Path to PCNA PDB file
        
    Returns
    -------
    list of PCNASite
        K164 positions (3 instances for trimer)
    """
    try:
        from ubsumo_em_simulator.io.pdb_reader import read_pdb_ca
        residues, positions, res_ids = read_pdb_ca(pdb_path)
    except ImportError:
        return get_default_pcna_sites()
    
    sites = []
    chain_labels = ["A", "B", "C"]
    k164_count = 0
    
    for rid, pos in zip(res_ids, positions):
        if rid.strip() == "164":
            chain = chain_labels[k164_count % 3] if k164_count < 3 else "?"
            sites.append(PCNASite(
                site_id=164,
                position=pos,
                chain=chain,
                charge=1.0,
                solvent_accessibility=0.75,
                experimental="competition",
                expected_preference="SUMO"
            ))
            k164_count += 1
    
    if len(sites) == 0:
        return get_default_pcna_sites()
    
    return sites


def get_default_pcna_sites() -> List[PCNASite]:
    """Fallback PCNA sites when PDB is not available."""
    return [
        PCNASite(
            site_id=164,
            position=np.array([19.1, 36.9, 20.1]),
            chain="A",
            charge=1.0,
            solvent_accessibility=0.75,
            experimental="competition",
            expected_preference="SUMO"
        ),
    ]


def build_ub_from_pdb(pdb_path: str = "data/pdb/ubiquitin_1UBQ.pdb") -> Optional[Dict]:
    """Build Ub modifier from real PDB structure."""
    try:
        from ubsumo_em_simulator.io.pdb_reader import read_pdb_ca, assign_charges
        residues, positions, _ = read_pdb_ca(pdb_path)
        charges = assign_charges(residues)
        center = np.mean(positions, axis=0)
        
        return {
            "name": "Ub",
            "active_site_position": center,
            "charges": charges.tolist(),
            "positions": [p.tolist() for p in positions],
            "n_residues": len(residues),
            "net_charge": float(sum(charges)),
        }
    except (ImportError, FileNotFoundError):
        return None


def build_sumo_from_pdb(pdb_path: str = "data/pdb/sumo1_1A5R.pdb") -> Optional[Dict]:
    """Build SUMO-1 modifier from real PDB structure."""
    try:
        from ubsumo_em_simulator.io.pdb_reader import read_pdb_ca, assign_charges
        residues, positions, _ = read_pdb_ca(pdb_path)
        charges = assign_charges(residues)
        center = np.mean(positions, axis=0)
        
        return {
            "name": "SUMO",
            "active_site_position": center,
            "charges": charges.tolist(),
            "positions": [p.tolist() for p in positions],
            "n_residues": len(residues),
            "net_charge": float(sum(charges)),
        }
    except (ImportError, FileNotFoundError):
        return None


def build_ub_modifier() -> Dict:
    """Build Ub modifier: try real PDB first, fallback to approximate."""
    real = build_ub_from_pdb()
    if real is not None:
        return real
    
    # Fallback: approximate charges
    charges = [
        1.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,
        1.0,  0.0, -1.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,
        0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  1.0,  0.0,
        0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,
        1.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0,  0.0,
       -1.0, -1.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
        0.0,  0.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
        0.0,  1.0,  0.0,  1.0,  0.0,  0.0
    ]
    np.random.seed(42)
    positions = np.random.normal(0, 5.0, (76, 3))
    return {
        "name": "Ub (approx)",
        "active_site_position": np.array([0.0, 0.0, 0.0]),
        "charges": charges,
        "positions": positions,
        "n_residues": 76,
        "net_charge": float(sum(charges)),
    }


def build_sumo_modifier() -> Dict:
    """Build SUMO modifier: try real PDB first, fallback to approximate."""
    real = build_sumo_from_pdb()
    if real is not None:
        return real
    
    # Fallback: approximate charges
    charges = [
        0.0,  0.0,  0.0, -1.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,
        1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,
       -1.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,
        1.0,  0.0,  0.0,  0.0, -1.0,  0.0,  1.0,  0.0,  0.0,  0.0,
        0.0,  1.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,
       -1.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
        0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,
        0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
        0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
        0.0
    ]
    np.random.seed(24)
    positions = np.random.normal(0, 5.0, (91, 3))
    return {
        "name": "SUMO (approx)",
        "active_site_position": np.array([1.0, 1.0, 1.0]),
        "charges": charges,
        "positions": positions,
        "n_residues": 91,
        "net_charge": float(sum(charges)),
    }


def benchmark_simulation(simulator_result: List, 
                        sites: List[PCNASite] = None,
                        tolerance_kj: float = 5.0
                        ) -> Dict:
    """Benchmark simulation results against experimental PCNA data.
    
    Parameters
    ----------
    simulator_result : list of CompetitionResult
        Results from CompetitionSimulator
    sites : list of PCNASite, optional
        PCNA sites for comparison. If None, uses default K164.
    tolerance_kj : float
        Energy tolerance for "match" (kJ/mol)
        
    Returns
    -------
    dict
        Benchmark metrics
    """
    if sites is None:
        sites = load_pcna_sites_from_pdb()
    
    correct = 0
    total = min(len(sites), len(simulator_result))
    details = []
    
    for i in range(total):
        site = sites[i]
        result = simulator_result[i]
        
        # K164 is a competition site where SUMO is slightly favored
        # (no DNA damage conditions)
        if site.expected_preference == "SUMO":
            match = result.delta_e > 0  # SUMO preferred
        elif site.expected_preference == "Ub":
            match = result.delta_e < 0  # Ub preferred
        else:
            match = abs(result.delta_e) < tolerance_kj  # ambiguous
        
        if match:
            correct += 1
        
        details.append({
            "site": site.site_id,
            "chain": site.chain,
            "experimental": site.experimental,
            "expected": site.expected_preference,
            "predicted": result.predicted_modifier,
            "delta_e": result.delta_e,
            "match": match,
        })
    
    return {
        "accuracy": correct / total if total > 0 else 0.0,
        "correct": correct,
        "total": total,
        "details": details,
    }