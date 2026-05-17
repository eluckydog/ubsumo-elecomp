"""
PDB file reader and charge assignment for protein structures.

Supports reading PDB format files and assigning charges to amino acid residues
based on standard pKa values at physiological pH.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from collections import OrderedDict


# Standard amino acid charges at pH 7.4
AA_CHARGES: Dict[str, float] = {
    # Positively charged
    "LYS": 1.0,   # Lysine (pKa ~10.5)
    "ARG": 1.0,   # Arginine (pKa ~12.5)
    "HIS": 0.5,   # Histidine (pKa ~6.0, partially charged at 7.4)
    # Negatively charged
    "ASP": -1.0,  # Aspartate (pKa ~3.9)
    "GLU": -1.0,  # Glutamate (pKa ~4.3)
    # N-terminal / C-terminal (charges at termini)
    # Neutral
    "ALA": 0.0, "ASN": 0.0, "CYS": 0.0, "GLN": 0.0,
    "GLY": 0.0, "ILE": 0.0, "LEU": 0.0, "MET": 0.0,
    "PHE": 0.0, "PRO": 0.0, "SER": 0.0, "THR": 0.0,
    "TRP": 0.0, "TYR": 0.0, "VAL": 0.0,
}

# Three-letter to one-letter mapping
AA_3TO1: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def read_pdb(pdb_path: str) -> Tuple[List[str], List[np.ndarray], List[str]]:
    """Read all atom coordinates from a PDB file.
    
    Parameters
    ----------
    pdb_path : str
        Path to PDB file
        
    Returns
    -------
    tuple
        (residue_names, positions, atom_names)
    """
    residue_names = []
    positions = []
    atom_names = []
    
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Parse residue name (columns 18-20)
                res_name = line[17:20].strip()
                
                # Parse coordinates (columns 31-54)
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue
                
                # Parse atom name
                atom_name = line[12:16].strip()
                
                residue_names.append(res_name)
                positions.append(np.array([x, y, z]))
                atom_names.append(atom_name)
    
    return residue_names, positions, atom_names


def read_pdb_ca(pdb_path: str) -> Tuple[List[str], List[np.ndarray], List[str]]:
    """Read CA (alpha-carbon) atoms only from a PDB file.
    
    One CA atom per residue = one entry per residue. Handles:
    - NMR ensembles: stops at first ENDMDL record
    - Multi-chain: all chains included
    
    Parameters
    ----------
    pdb_path : str
        Path to PDB file
        
    Returns
    -------
    tuple
        (residue_names, positions, residue_ids)
    """
    residue_names = []
    positions = []
    residue_ids = []
    
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ENDMDL"):
                break
            if line.startswith("ATOM") and " CA " in line[12:16]:
                res_name = line[17:20].strip()
                res_id = line[22:26].strip()
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue
                
                residue_names.append(res_name)
                positions.append(np.array([x, y, z]))
                residue_ids.append(res_id)
    
    return residue_names, positions, residue_ids


def assign_charges(residue_names: List[str]) -> np.ndarray:
    """Assign charges to residues based on amino acid type.
    
    Parameters
    ----------
    residue_names : list of str
        Three-letter amino acid names
        
    Returns
    -------
    np.ndarray
        Charges in elementary units (e)
    """
    charges = []
    for res in residue_names:
        charges.append(AA_CHARGES.get(res, 0.0))
    return np.array(charges)


def extract_lysine_sites(residue_names: List[str], positions: List[np.ndarray],
                         charges: np.ndarray) -> List[Dict]:
    """Extract lysine (K) sites from parsed PDB data.
    
    Parameters
    ----------
    residue_names : list of str
        Three-letter residue names
    positions : list of np.ndarray
        Residue positions (CA atoms)
    charges : np.ndarray
        Residue charges
        
    Returns
    -------
    list of dict
        Lysine site info: {residue_id, position, charge, ...}
    """
    sites = []
    residue_id = 0
    
    for i, res in enumerate(residue_names):
        if res == "LYS":
            sites.append({
                "residue_id": residue_id,
                "position": positions[i],
                "charge": charges[i],
            })
        residue_id += 1
    
    return sites


def build_modifier_from_pdb(pdb_path: str, name: str) -> Dict:
    """Build a ModifierProtein from a PDB file.
    
    Parameters
    ----------
    pdb_path : str
        Path to PDB file
    name : str
        Modifier name ("Ub" or "SUMO")
        
    Returns
    -------
    dict
        Modifier protein info with charges and positions
    """
    residues, positions, atoms = read_pdb(pdb_path)
    charges = assign_charges(residues)
    
    # Use center of mass as active site (approximation)
    center = np.mean(positions, axis=0)
    
    return {
        "name": name,
        "active_site_position": center,
        "charges": charges.tolist(),
        "positions": [p.tolist() for p in positions],
    }