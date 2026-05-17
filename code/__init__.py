"""
ubsumo-em-simulator: Ubiquitin/SUMO Electromagnetic Competition Simulator

A notebook-compatible Python package for simulating Ub/SUMO competition
on protein substrates based on electrostatic interaction physics.
"""

__version__ = "0.1.0"
__author__ = "math-science agent"

# Core modules work without OpenMM
from .core.competition import CompetitionSimulator, ProteinSite, ModifierProtein, CompetitionResult
from .core.em_force import EMForce, HAS_OPENMM

# I/O modules
from .io.pdb_reader import read_pdb, read_pdb_ca, assign_charges, extract_lysine_sites
from .io.trajectory import TrajectoryWriter, SummaryWriter

# Analysis modules
from .analysis.electrostatics import ElectrostaticAnalyzer
from .analysis.binding import BindingEnergyEstimator

# Validation
from .validation.pcna_benchmark import benchmark_simulation, load_pcna_sites_from_pdb

__all__ = [
    "CompetitionSimulator",
    "ProteinSite", 
    "ModifierProtein",
    "CompetitionResult",
    "EMForce",
    "HAS_OPENMM",
    "read_pdb",
    "assign_charges",
    "extract_lysine_sites",
    "TrajectoryWriter",
    "SummaryWriter",
    "ElectrostaticAnalyzer",
    "BindingEnergyEstimator",
    "PCNA_SITES",
    "benchmark_simulation",
    "load_pcna_sites_from_pdb",
]