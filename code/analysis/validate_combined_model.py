"""Validation: combined electrostatic + enzyme specificity model on 4 benchmark proteins.

Compares predictions against DeepPCT experimental labels.
"""
import sys, os, csv, re, math
import numpy as np
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from ubsumo_em_simulator.core.competition import CompetitionSimulator, ProteinSite, ModifierProtein
from ubsumo_em_simulator.io.pdb_reader import assign_charges

# ===== PDB Reader (inline to avoid import issues) =====
def read_pdb_ca(pdb_path):
    residues, positions, residue_ids = [], [], []
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ENDMDL"): break
            if line.startswith("ATOM") and " CA " in line[12:16]:
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                except ValueError: continue
                residues.append(line[17:20].strip())
                positions.append(np.array([x, y, z]))
                residue_ids.append(line[22:26].strip())
    return residues, positions, residue_ids

AA3TO1 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
          "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
          "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

# ===== Enzyme Specificity Filter =====
class EnzymeFilter:
    SUMO_MOTIF = re.compile(r"[ILVFM]K.[DE]")
    
    def score(self, seq, k_pos):
        s = 0.0
        # SUMO consensus
        for i in range(len(seq)-3):
            if i+1 == k_pos and self.SUMO_MOTIF.match(seq[i:i+4]):
                s += 2.0
                break
        # Local charge
        neg = sum(1 for i,aa in enumerate(seq) if aa in "DE" and abs(i-k_pos)<=5)
        pos = sum(1 for i,aa in enumerate(seq) if aa in "KRH" and abs(i-k_pos)<=5 and i!=k_pos)
        s += 0.3 * (neg - pos)
        # Hydrophobic context (favors Ub)
        hydro = sum(1 for i,aa in enumerate(seq) if aa in "ILVFM" and abs(i-k_pos)<=5)
        s -= 0.2 * hydro
        return s
    
    def combined_prob(self, seq, k_pos, de_electrostatic, T=2.5):
        raw = de_electrostatic / T + self.score(seq, k_pos)
        p = 1.0 / (1.0 + math.exp(-raw))
        return np.clip(p, 0.01, 0.99)

# ===== Load DeepPCT labels =====
PTM_NORM = {"ubiquitination":"Ub","Ubiquitination":"Ub","sumoylation":"SUMO",
            "Sumoylation":"SUMO","SUMOylation":"SUMO"}
def parse_site(s): return (s[0], int(s[1:]))

def load_labels(csv_path):
    labels = {}
    with open(csv_path, "r", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            for site_col, ptm_col in [("Site1","PTM1"), ("Site2","PTM2")]:
                ptm = PTM_NORM.get(row[ptm_col].strip(), "")
                aa, pos = parse_site(row[site_col].strip())
                if ptm in ("Ub","SUMO") and aa == "K":
                    key = (row["UniProt ID"].strip(), pos)
                    if key not in labels:
                        labels[key] = {"Ub": 0, "SUMO": 0}
                    labels[key][ptm] += 1
    return labels

def load_all_labels():
    all_l = {}
    for csv_path in [
        "archives/subagent_artifacts/temp_deeppct/datasets/training-set.csv",
        "archives/subagent_artifacts/temp_deeppct/datasets/indepentent-test-set.csv",
    ]:
        for key, counts in load_labels(csv_path).items():
            if key not in all_l: all_l[key] = {"Ub": 0, "SUMO": 0}
            all_l[key]["Ub"] += counts["Ub"]
            all_l[key]["SUMO"] += counts["SUMO"]
    return all_l

# ===== Build modifiers =====
def build_modifier(pdb_path, name):
    residues, positions, _ = read_pdb_ca(pdb_path)
    charges = assign_charges(residues)
    return ModifierProtein(name, np.mean(positions, axis=0),
                          list(charges), [p.tolist() for p in positions])

# ===== Main validation =====
print("=" * 65)
print("COMBINED MODEL VALIDATION: Electrostatic + Enzyme Specificity")
print("=" * 65)

ub = build_modifier("data/pdb/ubiquitin_1UBQ.pdb", "Ub")
sumo = build_modifier("data/pdb/sumo1_1A5R.pdb", "SUMO")
sim = CompetitionSimulator(epsilon_r=80.0)
ef = EnzymeFilter()
labels = load_all_labels()
sim_total = 0
sim_match = 0
combined_total = 0
combined_match = 0
targets = [
    ("p53", "P04637", "data/pdb/p53_1TUP.pdb"),
    ("ER-alpha", "P03372", "data/pdb/ER-alpha_1A52.pdb"),
    ("DNMT1", "P26358", "data/pdb/DNMT1_3PTA.pdb"),
    ("PML", "P29590", "data/pdb/PML_1BOR.pdb"),
]

for name, uniprot, pdb_path in targets:
    residues, positions, res_ids = read_pdb_ca(pdb_path)
    lys_sites = [(rid.strip(), res, pos, i) for i, (rid, res, pos) in 
                 enumerate(zip(res_ids, residues, positions)) if res == "LYS"]
    # Build full sequence
    seq_parts = []
    seq_ids = []
    for rn, rid in zip(residues, res_ids):
        if rid.strip() not in [x[0] for x in seq_ids]:
            seq_parts.append(AA3TO1.get(rn, "X"))
            seq_ids.append((rid.strip(), len(seq_parts)-1))
    full_seq = "".join(seq_parts)
    
    print("\n--- %s (%s) | %d residues, %d lysines ---" % (name, uniprot, len(residues), len(lys_sites)))
    
    # Count ground truth labels for this protein
    ub_sites = sum(1 for (uid, pos), c in labels.items() if uid == uniprot and c["Ub"] > c["SUMO"])
    sumo_sites = sum(1 for (uid, pos), c in labels.items() if uid == uniprot and c["SUMO"] > c["Ub"])
    both_sites = sum(1 for (uid, pos), c in labels.items() if uid == uniprot and c["Ub"] > 0 and c["SUMO"] > 0)
    total_labeled = ub_sites + sumo_sites
    print("DeepPCT labels: %d Ub-favoring, %d SUMO-favoring, %d both (total=%d)" % (
        ub_sites, sumo_sites, both_sites, total_labeled))
    
    # Test each lysine with E-only and E+Seq
    results = []
    for site_id, res_name, pos, idx in lys_sites:
        # Need Lys at K position = site_id in PCNA numbering
        # Actually for these proteins, site_id is the residue number from PDB
        protein_site = ProteinSite(int(site_id) if site_id.isdigit() else 0, pos, 1.0, 0.75)
        r = sim.evaluate_site(protein_site, ub, sumo)
        
        # Get sequence window around this K
        k_pos_in_seq = None
        for sid, si in seq_ids:
            if sid == site_id:
                k_pos_in_seq = si
                break
        
        seq_score = 0.0
        combined_p = 0.5
        if k_pos_in_seq is not None:
            start = max(0, k_pos_in_seq - 7)
            end = min(len(full_seq), k_pos_in_seq + 8)
            window = full_seq[start:end]
            k_in_win = k_pos_in_seq - start
            seq_score = ef.score(window, k_in_win)
            combined_p = ef.combined_prob(window, k_in_win, r.delta_e)
        
        # Get ground truth
        gt_key = (uniprot, int(site_id)) if site_id.isdigit() else None
        gt_label = "N/A"
        if gt_key and gt_key in labels:
            counts = labels[gt_key]
            if counts["Ub"] > counts["SUMO"]: gt_label = "Ub"
            elif counts["SUMO"] > counts["Ub"]: gt_label = "SUMO"
            elif counts["Ub"] > 0 and counts["SUMO"] > 0: gt_label = "both"
            else: gt_label = "N/A"
        
        results.append({
            "site": site_id,
            "delta_e": r.delta_e,
            "seq_score": seq_score,
            "combined_p_sumo": combined_p,
            "pred_e": r.predicted_modifier,
            "pred_combined": "SUMO" if combined_p > 0.55 else ("Ub" if combined_p < 0.45 else "ambiguous"),
            "gt": gt_label,
        })
    
    # Display results (top and bottom by combined P)
    results.sort(key=lambda x: x["combined_p_sumo"], reverse=True)
    
    # Count matches
    e_matches = 0
    c_matches = 0
    compared = 0
    for r in results:
        if r["gt"] in ("Ub", "SUMO"):
            if r["pred_e"] == r["gt"]:
                e_matches += 1
            if r["pred_combined"] == r["gt"]:
                c_matches += 1
            compared += 1
    
    print("%d lysines with labels: E-only acc=%.1f%% (%d/%d) | Combined acc=%.1f%% (%d/%d)" % (
        compared, 100*e_matches/compared if compared else 0, e_matches, compared,
        100*c_matches/compared if compared else 0, c_matches, compared))
    
    sim_total += compared
    sim_match += e_matches
    combined_total += compared
    combined_match += c_matches
    
    # Show top/bottom of each
    print("  Top-5 SUMO-predicted:")
    for r in results[:5]:
        print("    K%s: E=%+.2f Seq=%+.1f P(SUMO)=%.3f -> %s vs GT=%s" % (
            r["site"], r["delta_e"], r["seq_score"], r["combined_p_sumo"],
            r["pred_combined"], r["gt"]))
    print("  Bottom-5 (Ub-predicted):")
    for r in results[-5:]:
        print("    K%s: E=%+.2f Seq=%+.1f P(SUMO)=%.3f -> %s vs GT=%s" % (
            r["site"], r["delta_e"], r["seq_score"], r["combined_p_sumo"],
            r["pred_combined"], r["gt"]))

print("\n" + "=" * 65)
print("OVERALL ACCURACY")
print("=" * 65)
print("Electrostatic only: %.1f%% (%d/%d)" % (100*sim_match/sim_total if sim_total else 0, sim_match, sim_total))
print("E + Enzyme filter:  %.1f%% (%d/%d)" % (100*combined_match/combined_total if combined_total else 0, combined_match, combined_total))
print("Total compared sites: %d" % combined_total)
print("\nDone!")