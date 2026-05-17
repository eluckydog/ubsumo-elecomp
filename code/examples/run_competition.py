"""
Run UbSUMO EM competition simulation on PCNA benchmark.

This example demonstrates the full pipeline:
1. Read/modify protein structures
2. Run electrostatic competition simulation
3. Validate against experimental data
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
from ubsumo_em_simulator.core.competition import (
    CompetitionSimulator, ProteinSite, ModifierProtein, CompetitionResult
)
from ubsumo_em_simulator.validation.pcna_benchmark import (
    PCNA_SITES, benchmark_simulation, build_ub_modifier, build_sumo_modifier
)
from ubsumo_em_simulator.io.trajectory import TrajectoryWriter, SummaryWriter
from ubsumo_em_simulator.analysis.binding import BindingEnergyEstimator


def main():
    print("=" * 60)
    print("UbSUMO EM Competition Simulator - PCNA Benchmark")
    print("=" * 60)
    
    # 1. Build modifiers
    print("\n[1/4] Building Ub and SUMO modifiers...")
    ub_data = build_ub_modifier()
    sumo_data = build_sumo_modifier()
    
    ub = ModifierProtein(
        name=ub_data["name"],
        active_site_position=ub_data["active_site_position"],
        charges=ub_data["charges"],
        positions=ub_data["positions"]
    )
    
    sumo = ModifierProtein(
        name=sumo_data["name"],
        active_site_position=sumo_data["active_site_position"],
        charges=sumo_data["charges"],
        positions=sumo_data["positions"]
    )
    
    print(f"  Ub: {len(ub.charges)} residues, net charge {sum(ub.charges):+.1f}")
    print(f"  SUMO: {len(sumo.charges)} residues, net charge {sum(sumo.charges):+.1f}")
    
    # 2. Prepare PCNA sites
    print("\n[2/4] Preparing PCNA sites...")
    sites = []
    for s in PCNA_SITES:
        sites.append(ProteinSite(
            residue_id=s.site_id,
            position=s.position,
            charge=s.charge,
            solvent_accessibility=s.solvent_accessibility,
            known_modification=s.experimental
        ))
    print(f"  {len(sites)} sites: {[s.residue_id for s in sites]}")
    
    # 3. Run competition simulation
    print("\n[3/4] Running electrostatic competition simulation...")
    sim = CompetitionSimulator(epsilon_r=80.0)
    results = sim.evaluate_sites(sites, ub, sumo)
    
    # Print results
    print(f"\n{'Site':>6} {'E_Ub':>10} {'E_SUMO':>10} {'DeltaE':>10} {'Winner':>10}")
    print("-" * 50)
    for r in results:
        winner = "SUMO" if r.delta_e > 0 else "Ub"
        print(f"{r.site_id:>6} {r.e_ub:>10.2f} {r.e_sumo:>10.2f} "
              f"{r.delta_e:>10.2f} {winner:>10}")
    
    # Softmax probabilities
    probs = sim.score_to_probability(results, temperature=2.5)
    print(f"\n  Site  P(Ub)    P(SUMO)")
    for i, r in enumerate(results):
        print(f"  {r.site_id:>4}  {probs['Ub'][i]:.3f}    {probs['SUMO'][i]:.3f}")
    
    # 4. Benchmark against experimental data
    print("\n[4/4] Benchmarking against experimental data...")
    benchmark = benchmark_simulation(results)
    
    print(f"\n  Accuracy: {benchmark['accuracy']:.0%} ({benchmark['correct']}/{benchmark['total']})")
    print(f"\n  {'Site':>6} {'Experimental':>14} {'Predicted':>10} {'DeltaE':>10} {'Match':>6}")
    print("-" * 55)
    for d in benchmark["details"]:
        print(f"  {d['site']:>6} {d['experimental']:>14} {d['predicted']:>10} "
              f"{d['delta_e']:>10.2f} {'YES' if d['match'] else 'NO':>6}")
    
    # Binding energy analysis
    print("\n[Bonus] Binding energy comparison:")
    binder = BindingEnergyEstimator()
    for r in results:
        comparison = binder.compare_binding(r.e_ub, r.e_sumo)
        print(f"\n  Site {r.site_id}:")
        print(f"    dG_Ub = {comparison['dG_ub']:.1f} kJ/mol")
        print(f"    dG_SUMO = {comparison['dG_sumo']:.1f} kJ/mol")
        print(f"    ddG = {comparison['ddG']:.1f} kJ/mol")
        print(f"    Selectivity = {comparison['selectivity']:.2f}")
        print(f"    Preferred: {comparison['preferred']}")
    
    # Save results
    print("\nSaving results...")
    
    # Text summary
    SummaryWriter.write_summary(results, "pcna_competition_summary.txt")
    print("  -> pcna_competition_summary.txt")
    
    # HDF5 trajectory (if h5py available)
    try:
        with TrajectoryWriter("pcna_competition.h5") as writer:
            writer.write_competition_results(results)
            writer.write_probabilities(probs)
        print("  -> pcna_competition.h5")
    except ImportError:
        print("  (h5py not available, skipping HDF5 output)")
    
    print("\nDone!")
    return results, benchmark


if __name__ == "__main__":
    main()