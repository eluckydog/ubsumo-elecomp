"""
Binding free energy estimation for Ub/SUMO-target interactions.
"""

import numpy as np
from typing import Dict, Optional


class BindingEnergyEstimator:
    """Estimate binding free energy using simple electrostatic model.
    
    Parameters
    ----------
    epsilon_r : float
        Dielectric constant
    """
    
    def __init__(self, epsilon_r: float = 80.0):
        self.epsilon_r = epsilon_r
        self.K = 138.935456  # kJ*nm/(mol*e^2)
    
    def estimate_delta_g(self, interaction_energy: float, 
                        desolvation_penalty: float = 10.0,
                        entropy_loss: float = 15.0
                        ) -> Dict[str, float]:
        """Estimate binding free energy.
        
        dG = E_electrostatic + dG_desolvation - T*dS
        
        Parameters
        ----------
        interaction_energy : float
            Electrostatic interaction energy (kJ/mol)
        desolvation_penalty : float
            Desolvation free energy cost (kJ/mol)
        entropy_loss : float
            Entropy loss upon binding (T*dS, kJ/mol)
            
        Returns
        -------
        dict
            Energy components and total
        """
        dG_total = interaction_energy + desolvation_penalty + entropy_loss
        
        return {
            "E_electrostatic": interaction_energy,
            "dG_desolvation": desolvation_penalty,
            "TdS": entropy_loss,
            "dG_total": dG_total,
            "binding_affinity": np.exp(-dG_total / 2.479)  # Kd in M, kT=2.479 at 298K
        }
    
    def compare_binding(self, e_ub: float, e_sumo: float
                       ) -> Dict[str, float]:
        """Compare Ub and SUMO binding energies.
        
        Parameters
        ----------
        e_ub : float
            Ub interaction energy (kJ/mol)
        e_sumo : float
            SUMO interaction energy (kJ/mol)
            
        Returns
        -------
        dict
            Binding comparison metrics
        """
        dG_ub = e_ub + 25.0  # Solvation + entropy penalty
        dG_sumo = e_sumo + 25.0
        
        ddG = dG_ub - dG_sumo
        
        # Selectivity
        Kd_ratio = np.exp(ddG / 2.479)
        
        return {
            "dG_ub": dG_ub,
            "dG_sumo": dG_sumo,
            "ddG": ddG,
            "selectivity": Kd_ratio,  # >1 favors SUMO, <1 favors Ub
            "preferred": "SUMO" if ddG > 0 else "Ub"
        }