"""
Electromagnetic interaction force for UbSUMO competition simulation.

Provides both OpenMM-based and pure NumPy implementations.
"""

import numpy as np
from typing import Optional

# OpenMM is optional - only needed for MD engine integration
try:
    from openmm import CustomNonbondedForce, System
    from openmm.unit import nanometer, elementary_charge, kilojoule_per_mole
    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False
    CustomNonbondedForce = None
    System = None


class EMForce:
    """Custom Coulomb force for protein electrostatic interactions.
    
    Implements: E = q1*q2 * coulomb / (epsilon_r * r)
    
    Parameters
    ----------
    epsilon_r : float
        Relative dielectric constant (default: 80.0 for water)
    cutoff : float
        Cutoff distance in nm (default: 1.2)
    """
    
    COULOMB_CONST = 138.935456  # kJ*nm/(mol*e^2)
    
    def __init__(self, epsilon_r: float = 80.0, cutoff: float = 1.2):
        self.epsilon_r = epsilon_r
        self.cutoff = cutoff
        self.force = None
    
    def create_force(self, system, charges: np.ndarray):
        """Create CustomNonbondedForce for electrostatic interactions.
        
        Requires OpenMM. Falls back to pure NumPy if not available.
        
        Parameters
        ----------
        system : openmm.System or None
            OpenMM System object
        charges : np.ndarray
            Array of charges in elementary units (e)
            
        Returns
        -------
        CustomNonbondedForce or None
            Configured force object, or None if OpenMM not available
        """
        if not HAS_OPENMM:
            raise ImportError(
                "OpenMM is required for CustomNonbondedForce. "
                "Install with: pip install openmm"
            )
        
        expression = "coulomb*charge1*charge2/(eps_r*r)"
        
        force = CustomNonbondedForce(expression)
        force.addPerParticleParameter("charge")
        force.addGlobalParameter("coulomb", self.COULOMB_CONST)
        force.addGlobalParameter("eps_r", self.epsilon_r)
        
        force.setCutoffDistance(self.cutoff * nanometer)
        force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
        
        for charge in charges:
            force.addParticle([charge])
        
        self.force = force
        system.addForce(force)
        
        return force
    
    @staticmethod
    def calculate_energy(positions: np.ndarray, charges: np.ndarray, 
                        epsilon_r: float = 80.0) -> float:
        """Calculate total electrostatic energy (pure NumPy).
        
        Parameters
        ----------
        positions : np.ndarray
            Particle positions in Angstroms
        charges : np.ndarray
            Particle charges in elementary units
        epsilon_r : float
            Relative dielectric constant
            
        Returns
        -------
        float
            Total energy in kJ/mol
        """
        n = len(positions)
        energy = 0.0
        
        for i in range(n):
            for j in range(i + 1, n):
                r = np.linalg.norm(positions[i] - positions[j]) / 10.0
                if r > 0.01:
                    e = EMForce.COULOMB_CONST * charges[i] * charges[j] / (epsilon_r * r)
                    energy += e
        
        return energy
    
    @staticmethod
    def pairwise_energy(positions_a: np.ndarray, charges_a: np.ndarray,
                       positions_b: np.ndarray, charges_b: np.ndarray,
                       epsilon_r: float = 80.0) -> np.ndarray:
        """Calculate pairwise electrostatic energy between two groups.
        
        Parameters
        ----------
        positions_a, positions_b : np.ndarray
            Residue positions in Angstroms
        charges_a, charges_b : np.ndarray
            Residue charges in e
        epsilon_r : float
            Dielectric constant
            
        Returns
        -------
        np.ndarray
            Energy matrix E[i,j] in kJ/mol
        """
        n_a, n_b = len(positions_a), len(positions_b)
        E = np.zeros((n_a, n_b))
        
        for i in range(n_a):
            for j in range(n_b):
                r = np.linalg.norm(positions_a[i] - positions_b[j]) / 10.0
                if r > 0.01:
                    E[i, j] = EMForce.COULOMB_CONST * charges_a[i] * charges_b[j] / (epsilon_r * r)
        
        return E


class DistanceDielectricEMForce(EMForce):
    """EM force with distance-dependent dielectric.
    
    Uses: epsilon_r = r to account for dielectric saturation.
    """
    
    def create_force(self, system, charges: np.ndarray):
        if not HAS_OPENMM:
            raise ImportError("OpenMM required. Install with: pip install openmm")
        
        expression = "coulomb*charge1*charge2/(r*r)"
        
        force = CustomNonbondedForce(expression)
        force.addPerParticleParameter("charge")
        force.addGlobalParameter("coulomb", self.COULOMB_CONST)
        
        force.setCutoffDistance(self.cutoff * nanometer)
        force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
        
        for charge in charges:
            force.addParticle([charge])
        
        self.force = force
        system.addForce(force)
        
        return force