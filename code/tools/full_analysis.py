"""Full analysis: electrostatic + DBREF matching across ALL local PDBs."""
import sys, os, csv, json
import numpy as np
from collections import defaultdict

WORKSPACE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, WORKSPACE)
from ubsumo_em_simulator.core.competition import CompetitionSimulator, ProteinSite, ModifierProtein
from ubsumo_em_simulator.io.pdb_reader import assign_charges

PDB_DIR = os.path.join(WORKSPACE, "data", "pdb_all")
RESULTS_DIR = os.path.join(WORKSPACE, "data", "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

def read_ca(pdb_path):
    """Parse CA atom positions from PDB file."""
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("MODEL") and not line.startswith("MODEL        1"):
                break
            if line.startswith("ATOM") and line[12:16].strip() in ("CA", "P"):
                try:
                    atoms.append((
                        line[17:20].strip(), line[21].strip(),
                        int(line[22:26].strip()) if line[22:26].strip().lstrip("-").isdigit() else 0,
                        line[26].strip(),
                        float(line[30:38]), float(line[38:46]), float(line[46:54])
                    ))
                except: pass
    return atoms

def parse_dbref(pdb_path):
    """Parse DBREF records from PDB file."""
    dbrefs = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("DBREF"):
                parts = line.split()
                if len(parts) >= 10:
                    try:
                        dbrefs.append({
                            "chain": parts[2],
                            "pdb_start": int(parts[3]), "pdb_end": int(parts[4]),
                            "uniprot": parts[6],
                            "unp_start": int(parts[8]), "unp_end": int(parts[9])
                        })
                    except: pass
    return dbrefs

def load_labels():
    PTM_NORM = {"ubiquitination":"Ub","sumoylation":"SUMO",
                "Ubiquitination":"Ub","Sumoylation":"SUMO","SUMOylation":"SUMO"}
    labels = defaultdict(lambda: {"Ub":0,"SUMO":0})
    for csv_path in ["training-set.csv", "indepentent-test-set.csv"]:
        fp = os.path.join(WORKSPACE, "archives", "subagent_artifacts", "temp_deeppct", "datasets", csv_path)
        with open(fp, "r", encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                for sc, pc in [("Site1","PTM1"),("Site2","PTM2")]:
                    ptm = PTM_NORM.get(row[pc].strip(), "")
                    site = row[sc].strip()
                    if ptm in ("Ub","SUMO") and site.startswith("K"):
                        labels[(row["UniProt ID"].strip(), int(site[1:]))][ptm] += 1
    return labels

# Load
mapping = json.load(open(os.path.join(WORKSPACE, "data", "uniprot_pdb_mapping.json")))
labels = load_labels()

# Build modifiers
def build_mod(pdb):
    atoms = read_ca(pdb)
    pos = [np.array([a[4],a[5],a[6]]) for a in atoms]
    ch = assign_charges([a[0] for a in atoms])
    return ModifierProtein(os.path.basename(pdb).replace(".pdb",""), np.mean(pos,axis=0), list(ch), pos)

ub = build_mod(os.path.join(PDB_DIR, "ubiquitin_1ubq.pdb"))
sumo = build_mod(os.path.join(PDB_DIR, "sumo1_1a5r.pdb"))
sim = CompetitionSimulator(80.0)

# For each labeled site, find ALL PDB positions and run electrostatic analysis
results = []

for (uid, k), gt in sorted(labels.items()):
    gt_label = "both" if gt["Ub"]>0 and gt["SUMO"]>0 else ("Ub" if gt["Ub"]>0 else "SUMO")
    
    # Find all PDB entries for this protein
    pdb_entries = mapping.get(uid, [])
    
    site_results = []  # (pdb_id, chain, unp_pos, de, pred)
    
    for entry in pdb_entries:
        pid = entry["pdb_id"].lower()
        pdb_path = os.path.join(PDB_DIR, "%s.pdb" % pid)
        if not os.path.exists(pdb_path):
            pdb_path = os.path.join(PDB_DIR, "%s.cif" % pid)
            if not os.path.exists(pdb_path):
                continue
        
        dbrefs = parse_dbref(pdb_path)
        atoms = read_ca(pdb_path)
        
        for dr in dbrefs:
            if dr["uniprot"] != uid:
                continue
            if not (dr["unp_start"] <= k <= dr["unp_end"]):
                continue
            
            offset = dr["unp_start"] - dr["pdb_start"]
            pdb_res = k - offset
            
            # Find this residue in PDB CA atoms
            for rname, chain, rnum, icode, x, y, z in atoms:
                if chain != dr["chain"]:
                    continue
                if rnum == pdb_res and rname == "LYS":
                    pos = np.array([x, y, z])
                    site = ProteinSite(k, pos, 1.0, 0.75)
                    r = sim.evaluate_site(site, ub, sumo)
                    site_results.append((pid, dr["chain"], pdb_res, r.delta_e, r.predicted_modifier))
                    break
    
    if not site_results:
        results.append({"uid": uid, "k": k, "gt": gt_label, "n_pdbs": 0, "de_list": [], "pred_list": []})
    else:
        de_vals = [sr[3] for sr in site_results]
        pred_vals = [sr[4] for sr in site_results]
        n_ub_pred = pred_vals.count("Ub")
        n_sumo_pred = pred_vals.count("SUMO")
        n_amb_pred = pred_vals.count("ambiguous")
        majority_pred = "Ub" if n_ub_pred > n_sumo_pred else ("SUMO" if n_sumo_pred > n_ub_pred else "ambiguous")
        
        results.append({
            "uid": uid, "k": k, "gt": gt_label,
            "n_pdbs": len(site_results),
            "de_mean": np.mean(de_vals),
            "de_std": np.std(de_vals),
            "de_min": min(de_vals),
            "de_max": max(de_vals),
            "majority_pred": majority_pred,
            "ub_pred": n_ub_pred, "sumo_pred": n_sumo_pred, "amb_pred": n_amb_pred,
            "sites": site_results
        })

# Save raw results
json.dump(results, open(os.path.join(RESULTS_DIR, "electrostatic_results.json"), "w"), indent=1)

# Summary
print("=" * 70)
print("ELECTROSTATIC ANALYSIS: %d labeled K sites across 54 proteins" % len(results))
print("=" * 70)

n_with_pdb = sum(1 for r in results if r["n_pdbs"] > 0)
n_without = sum(1 for r in results if r["n_pdbs"] == 0)
n_ub_gt = sum(1 for r in results if r["gt"] == "Ub")
n_sumo_gt = sum(1 for r in results if r["gt"] == "SUMO")
n_both_gt = sum(1 for r in results if r["gt"] == "both")

print("\nTotal: %d sites (Ub=%d, SUMO=%d, both=%d)" % (len(results), n_ub_gt, n_sumo_gt, n_both_gt))
print("Evidenced from PDB: %d sites (%.0f%%)" % (n_with_pdb, 100*n_with_pdb/len(results)))
print("No PDB evidence: %d sites" % n_without)

# Accuracy: majority vote vs ground truth
print("\n--- Majority Prediction vs Ground Truth ---")
for gt_type in ["Ub", "SUMO", "both"]:
    subset = [r for r in results if r["gt"] == gt_type and r["n_pdbs"] > 0]
    if not subset:
        continue
    correct = sum(1 for r in subset if r["majority_pred"] == gt_type)
    pred_ub = sum(1 for r in subset if r["majority_pred"] == "Ub")
    pred_sumo = sum(1 for r in subset if r["majority_pred"] == "SUMO")
    pred_amb = sum(1 for r in subset if r["majority_pred"] == "ambiguous")
    print("  GT=%s: %d sites, correct=%d/%.0f%%, pred Ub=%d, SUMO=%d, amb=%d" % (
        gt_type, len(subset), correct, 100*correct/len(subset), pred_ub, pred_sumo, pred_amb))

# Separatrix analysis
print("\n--- DE Distribution by Ground Truth ---")
for gt_type in ["Ub", "SUMO", "both"]:
    de_vals = [r["de_mean"] for r in results if r["gt"] == gt_type and r["n_pdbs"] > 0]
    if de_vals:
        print("  GT=%s: mean DE=%+.3f, std=%.3f, min=%+.3f, max=%+.3f" % (
            gt_type, np.mean(de_vals), np.std(de_vals), min(de_vals), max(de_vals)))

# Detailed output
print("\n--- Sites with strong signals (|DE| > 1.0) ---")
for r in results:
    if r["n_pdbs"] > 0 and abs(r["de_mean"]) >= 1.0:
        match = "[MATCH]" if r["majority_pred"] == r["gt"] else "[MISMATCH]"
        print("  %s %s K%d: n=%d DE=%.2f(+/-%.2f) majority=%s GT=%s %s" % (
            match, r["uid"], r["k"], r["n_pdbs"], r["de_mean"], r["de_std"],
            r["majority_pred"], r["gt"], "" if r["majority_pred"]==r["gt"] else ""))

# Sites found
print("\n--- All categorized ---")
for r in results:
    if r["n_pdbs"] > 0:
        match = "+" if r["majority_pred"] == r["gt"] else " "
        print("  %s %s K%d: n=%d DE=%.2f majority=%s GT=%s" % (
            match, r["uid"], r["k"], r["n_pdbs"], r["de_mean"], r["majority_pred"], r["gt"]))
    else:
        print("    %s K%d: NO PDB GT=%s" % (r["uid"], r["k"], r["gt"]))
