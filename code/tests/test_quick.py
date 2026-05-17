"""Quick end-to-end test for UbSUMO EM simulator."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'ubsumo_em_simulator'))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from ubsumo_em_simulator.core.competition import CompetitionSimulator, ProteinSite, ModifierProtein
from ubsumo_em_simulator.validation.pcna_benchmark import (
    PCNA_SITES, benchmark_simulation, build_ub_modifier, build_sumo_modifier
)

# Build modifiers
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

print(f"Ub: {len(ub.charges)} residues, net charge {sum(ub.charges):+.1f}")
print(f"SUMO: {len(sumo.charges)} residues, net charge {sum(sumo.charges):+.1f}")

# Prepare sites
sites = [ProteinSite(s.site_id, s.position, s.charge, s.solvent_accessibility, s.experimental)
         for s in PCNA_SITES]

# Run simulation
sim = CompetitionSimulator(epsilon_r=80.0)
results = sim.evaluate_sites(sites, ub, sumo)
probs = sim.score_to_probability(results)

print("\nResults:")
print(f"{'Site':>6} {'E_Ub':>10} {'E_SUMO':>10} {'DeltaE':>10} {'P(Ub)':>8} {'P(SUMO)':>8} {'Winner':>10}")
print("-" * 68)

for i, r in enumerate(results):
    p_ub = probs["Ub"][i]
    p_sumo = probs["SUMO"][i]
    print(f"K{r.site_id:>4} {r.e_ub:>10.2f} {r.e_sumo:>10.2f} {r.delta_e:>10.2f} {p_ub:>8.3f} {p_sumo:>8.3f} {r.predicted_modifier:>10}")

# Benchmark
bench = benchmark_simulation(results)
print(f"\nBenchmark Accuracy: {bench['accuracy']:.0%} ({bench['correct']}/{bench['total']})")

for d in bench["details"]:
    status = "PASS" if d["match"] else "FAIL"
    print(f"  K{d['site']:>4}: exp={d['experimental']:>10}, pred={d['predicted']:>10}, deltaE={d['delta_e']:.2f}, {status}")

print("\nAll tests complete!")