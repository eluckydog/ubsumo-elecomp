"""
Trajectory writer for UbSUMO competition simulation output.

Supports HDF5 format for efficient storage of simulation results.
h5py is optional - falls back to text output if not available.
"""

import numpy as np
from typing import List, Optional
from datetime import datetime

try:
    import h5py
    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False


class TrajectoryWriter:
    """Write simulation trajectories to HDF5 format.
    
    Requires h5py. Falls back gracefully if not available.
    
    Parameters
    ----------
    filename : str
        Output filename
    mode : str
        File mode ("w" for write, "a" for append)
    """
    
    def __init__(self, filename: str, mode: str = "w"):
        if not HAS_H5PY:
            raise ImportError("h5py is required for TrajectoryWriter. Install with: pip install h5py")
        self.filename = filename
        self.mode = mode
        self._file: Optional["h5py.File"] = None
    
    def open(self):
        """Open the HDF5 file."""
        self._file = h5py.File(self.filename, self.mode)
        self._file.attrs["created"] = datetime.now().isoformat()
        self._file.attrs["software"] = "ubsumo-em-simulator"
        return self
    
    def close(self):
        """Close the HDF5 file."""
        if self._file is not None:
            self._file.close()
            self._file = None
    
    def __enter__(self):
        return self.open()
    
    def __exit__(self, *args):
        self.close()
    
    def write_competition_results(self, results: List, group_name: str = "competition"):
        """Write competition results to HDF5."""
        if self._file is None:
            raise RuntimeError("File not open. Use `open()` or context manager.")
        
        group = self._file.create_group(group_name)
        
        n = len(results)
        site_ids = np.zeros(n, dtype=int)
        e_ub = np.zeros(n)
        e_sumo = np.zeros(n)
        delta_e = np.zeros(n)
        predictions = []
        
        for i, r in enumerate(results):
            site_ids[i] = r.site_id
            e_ub[i] = r.e_ub
            e_sumo[i] = r.e_sumo
            delta_e[i] = r.delta_e
            predictions.append(r.predicted_modifier)
        
        group.create_dataset("site_id", data=site_ids)
        group.create_dataset("e_ub", data=e_ub)
        group.create_dataset("e_sumo", data=e_sumo)
        group.create_dataset("delta_e", data=delta_e)
        group.create_dataset("prediction", data=np.array(predictions, dtype="S10"))
    
    def write_probabilities(self, probs: dict, group_name: str = "probabilities"):
        """Write softmax probabilities to HDF5."""
        if self._file is None:
            raise RuntimeError("File not open.")
        
        group = self._file.create_group(group_name)
        group.create_dataset("p_ub", data=np.array(probs["Ub"]))
        group.create_dataset("p_sumo", data=np.array(probs["SUMO"]))


class SummaryWriter:
    """Write text summary of competition results. No dependencies needed."""
    
    @staticmethod
    def write_summary(results: List, filepath: str):
        """Write human-readable summary.
        
        Parameters
        ----------
        results : list of CompetitionResult
            Competition results
        filepath : str
            Output file path
        """
        with open(filepath, "w", encoding="utf-8") as f:
            f.write("UbSUMO Competition Results Summary\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"{'Site':>6} {'E_Ub':>10} {'E_SUMO':>10} {'Delta':>10} {'Prediction':>12}\n")
            f.write("-" * 50 + "\n")
            
            for r in results:
                f.write(f"{r.site_id:>6} {r.e_ub:>10.2f} {r.e_sumo:>10.2f} "
                       f"{r.delta_e:>10.2f} {r.predicted_modifier:>12}\n")