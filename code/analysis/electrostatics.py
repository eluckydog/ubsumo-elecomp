"""
Electrostatic analysis for protein modifier interactions.
"""

import numpy as np
from typing import List, Dict, Tuple


class ElectrostaticAnalyzer:
    """Analyze electrostatic properties of protein-modifier interactions.
    
    Parameters
    ----------
    epsilon_r : float
        Dielectric constant
    """
    
    def __init__(self, epsilon_r: float = 80.0):
        self.epsilon_r = epsilon_r
        # Coulomb constant
        self.K = 138.935456  # kJ*nm/(mol*e^2)
    
    def pairwise_energy_matrix(self, positions_a: np.ndarray, charges_a: np.ndarray,
                               positions_b: np.ndarray, charges_b: np.ndarray) -> np.ndarray:
        """Compute pairwise electrostatic energy matrix between two groups.
        
        Parameters
        ----------
        positions_a, positions_b : np.ndarray
            Residue positions in Angstroms
        charges_a, charges_b : np.ndarray
            Residue charges in e
            
        Returns
        -------
        np.ndarray
            Energy matrix E[i,j] in kJ/mol
        """
        n_a, n_b = len(positions_a), len(positions_b)
        E = np.zeros((n_a, n_b))
        
        for i in range(n_a):
            for j in range(n_b):
                # Distance in nm
                r = np.linalg.norm(positions_a[i] - positions_b[j]) / 10.0
                if r > 0.01:
                    E[i, j] = self.K * charges_a[i] * charges_b[j] / (self.epsilon_r * r)
        
        return E
    
    def electrostatic_potential(self, positions: np.ndarray, charges: np.ndarray,
                               grid_points: np.ndarray) -> np.ndarray:
        """Compute electrostatic potential at given grid points.
        
        Parameters
        ----------
        positions : np.ndarray
            Source charge positions in Angstroms
        charges : np.ndarray
            Source charges in e
        grid_points : np.ndarray
            Grid points where potential is evaluated, in Angstroms
            
        Returns
        -------
        np.ndarray
            Potential at each grid point in kJ/(mol*e)
        """
        potential = np.zeros(len(grid_points))
        
        for i, gp in enumerate(grid_points):
            for pos, charge in zip(positions, charges):
                r = np.linalg.norm(gp - pos) / 10.0  # A -> nm
                if r > 0.01:
                    potential[i] += self.K * charge / (self.epsilon_r * r)
        
        return potential
    
    def charge_complementarity(self, site_charge: float, modifier_charges: np.ndarray,
                              modifier_positions: np.ndarray, site_position: np.ndarray
                              ) -> Dict[str, float]:
        """Compute charge complementarity between a site and a modifier.
        
        Parameters
        ----------
        site_charge : float
            Target site charge
        modifier_charges : np.ndarray
            Modifier residue charges
        modifier_positions : np.ndarray
            Modifier residue positions
        site_position : np.ndarray
            Target site position
            
        Returns
        -------
        dict
            Complementarity metrics
        """
        distances = np.linalg.norm(modifier_positions - site_position, axis=1) / 10.0
        contacts = np.sum(distances < 0.5)  # Residues within 5 Angstrom
        
        # Local charge around site
        local_charge = 0.0
        for charge, dist in zip(modifier_charges, distances):
            if dist < 0.8:
                local_charge += charge / (dist + 0.1)
        
        return {
            "contacts": contacts,
            "local_charge": local_charge,
            "complementarity": -site_charge * local_charge  # Higher = more complementary
        }