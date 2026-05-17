"""
Ub/SUMO competition simulator based on electrostatic energy comparison.

The core logic: Ub and SUMO compete for the same lysine site.
The one with lower electrostatic interaction energy wins.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
from .em_force import EMForce


@dataclass
class ProteinSite:
    """A lysine site on a target protein."""
    residue_id: int
    position: np.ndarray  # 3D coordinates in Angstroms
    charge: float
    solvent_accessibility: float  # 0-1
    known_modification: Optional[str] = None  # "Ub", "SUMO", or None


@dataclass
class ModifierProtein:
    """A modifying protein (Ubiquitin or SUMO) with charge distribution."""
    name: str  # "Ub" or "SUMO"
    active_site_position: np.ndarray  # active site residue position
    charges: List[float]  # per-residue charges
    positions: List[np.ndarray]  # per-residue positions


@dataclass
class CompetitionResult:
    """Result of Ub/SUMO competition at a single site."""
    site_id: int
    e_ub: float  # kJ/mol
    e_sumo: float  # kJ/mol
    delta_e: float  # E_Ub - E_SUMO
    predicted_modifier: str  # "Ub", "SUMO", or "ambiguous"


class CompetitionSimulator:
    """Simulate Ub/SUMO competition on protein substrates.
    
    Parameters
    ----------
    epsilon_r : float
        Solvent dielectric constant (default: 80.0)
    """
    
    def __init__(self, epsilon_r: float = 80.0, cutoff: float = 1.2):
        self.epsilon_r = epsilon_r
        self.cutoff = cutoff
        self.em_force = EMForce(epsilon_r=epsilon_r, cutoff=cutoff)
        
    def evaluate_site(self, site: ProteinSite, ub: ModifierProtein, 
                     sumo: ModifierProtein) -> CompetitionResult:
        """Evaluate competition at a single lysine site.
        
        Parameters
        ----------
        site : ProteinSite
            Target lysine site
        ub : ModifierProtein
            Ubiquitin modifier
        sumo : ModifierProtein
            SUMO modifier
            
        Returns
        -------
        CompetitionResult
            Predicted competition outcome
        """
        # Calculate interaction energy between site and Ub
        e_ub = self._interaction_energy(site, ub)
        
        # Calculate interaction energy between site and SUMO
        e_sumo = self._interaction_energy(site, sumo)
        
        delta_e = e_ub - e_sumo
        
        # Determine winner
        if delta_e < -1.0:  # Significantly favors Ub
            predicted = "Ub"
        elif delta_e > 1.0:  # Significantly favors SUMO
            predicted = "SUMO"
        else:
            predicted = "ambiguous"
        
        return CompetitionResult(
            site_id=site.residue_id,
            e_ub=e_ub,
            e_sumo=e_sumo,
            delta_e=delta_e,
            predicted_modifier=predicted
        )
    
    def _interaction_energy(self, site: ProteinSite, modifier: ModifierProtein) -> float:
        """Calculate electrostatic interaction energy between site and modifier.
        
        Parameters
        ----------
        site : ProteinSite
            Target lysine site
        modifier : ModifierProtein
            Modifying protein
            
        Returns
        -------
        float
            Interaction energy in kJ/mol
        """
        energy = 0.0
        
        # Site charge interacts with all modifier charges
        for mod_charge, mod_pos in zip(modifier.charges, modifier.positions):
            r = np.linalg.norm(site.position - mod_pos) / 10.0  # A -> nm
            if r > 0.01:
                e = EMForce.COULOMB_CONST * site.charge * mod_charge / (self.epsilon_r * r)
                energy += e
        
        # Solvent accessibility correction
        energy *= site.solvent_accessibility
        
        return energy
    
    def evaluate_sites(self, sites: List[ProteinSite], ub: ModifierProtein,
                      sumo: ModifierProtein) -> List[CompetitionResult]:
        """Evaluate competition at multiple sites.
        
        Parameters
        ----------
        sites : list of ProteinSite
            Target lysine sites
        ub : ModifierProtein
            Ubiquitin modifier
        sumo : ModifierProtein
            SUMO modifier
            
        Returns
        -------
        list of CompetitionResult
            Predicted competition outcomes for all sites
        """
        return [self.evaluate_site(site, ub, sumo) for site in sites]
    
    def score_to_probability(self, results: List[CompetitionResult], 
                            temperature: float = 2.5) -> Dict[str, List[float]]:
        """Convert energy differences to probabilities using softmax.
        
        Parameters
        ----------
        results : list of CompetitionResult
            Competition results
        temperature : float
            Softmax temperature (lower = more decisive)
            
        Returns
        -------
        dict
            Probabilities for Ub and SUMO at each site
        """
        p_ub = []
        p_sumo = []
        
        for r in results:
            # Softmax: P(Ub) = exp(-E_Ub/T) / (exp(-E_Ub/T) + exp(-E_SUMO/T))
            exp_ub = np.exp(-r.e_ub / temperature)
            exp_sumo = np.exp(-r.e_sumo / temperature)
            total = exp_ub + exp_sumo
            
            p_ub.append(exp_ub / total)
            p_sumo.append(exp_sumo / total)
        
        return {"Ub": p_ub, "SUMO": p_sumo}