"""Validation: DBREF-filtered, electrostatic + enzyme specificity vs DeepPCT labels."""
import sys, os, re, math, csv, json
import numpy as np
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from ubsumo_em_simulator.core.competition import CompetitionSimulator, ProteinSite, ModifierProtein
from ubsumo_em_simulator.io.pdb_reader import assign_charges

AA3TO1 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
          "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
          "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

def read_pdb_ca(pdb_path):
    residues, positions, residue_ids = [], [], []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ENDMDL"): break
            if line.startswith("ATOM") and " CA " in line[12:16]:
                try: x=float(line[30:38]); y=float(line[38:46]); z=float(line[46:54])
                except: continue
                residues.append(line[17:20].strip())
                positions.append(np.array([x,y,z]))
                residue_ids.append(line[22:26].strip())
    return residues, positions, residue_ids

def parse_dbref(pdb_path):
    """DBREF 1TUP A 94 312 UNP P04637 P53_HUMAN 94 312"""
    ranges = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("DBREF") and "UNP" in line:
                parts = line.split()
                if len(parts) >= 10:
                    ranges.append({
                        "chain": parts[2], "pdb_start": int(parts[3]), "pdb_end": int(parts[4]),
                        "uniprot": parts[6], "uniprot_start": int(parts[8]), "uniprot_end": int(parts[9]),
                        "offset": int(parts[8]) - int(parts[3])
                    })
    return ranges

class EnzymeFilter:
    SUMO_MOTIF = re.compile(r"[ILVFM]K.[DE]")
    def score(self, seq, k_pos):
        s = 0.0
        for i in range(len(seq)-3):
            if i+1==k_pos and self.SUMO_MOTIF.match(seq[i:i+4]): s += 2.0; break
        neg = sum(1 for i,a in enumerate(seq) if a in "DE" and abs(i-k_pos)<=5)
        pos = sum(1 for i,a in enumerate(seq) if a in "KRH" and abs(i-k_pos)<=5 and i!=k_pos)
        s += 0.3*(neg-pos)
        s -= 0.2*sum(1 for i,a in enumerate(seq) if a in "ILVFM" and abs(i-k_pos)<=5)
        return s
    def prob(self, seq, k_pos, de, T=2.5):
        raw = de/T + self.score(seq, k_pos)
        return np.clip(1.0/(1.0+math.exp(-raw)), 0.01, 0.99)

def build_mod(pdb, name):
    r, pos, _ = read_pdb_ca(pdb)
    ch = assign_charges(r)
    center = np.mean(pos, axis=0)
    return ModifierProtein(name, center, list(ch), pos)

def load_deeppct_labels(csv_paths):
    PTM_NORM = {"ubiquitination":"Ub","Ubiquitination":"Ub","sumoylation":"SUMO",
                "Sumoylation":"SUMO","SUMOylation":"SUMO"}
    labels = defaultdict(lambda: {"Ub":0,"SUMO":0})
    for path in csv_paths:
        with open(path, "r", encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                for sc, pc in [("Site1","PTM1"),("Site2","PTM2")]:
                    ptm = PTM_NORM.get(row[pc].strip(), "")
                    site = row[sc].strip()
                    if ptm in ("Ub","SUMO") and site.startswith("K") and len(site) > 1:
                        key = (row["UniProt ID"].strip(), int(site[1:]))
                        labels[key][ptm] += 1
    return labels

def get_pdb_lyso_sequence(pdb_path):
    """Build unique residue 1-letter sequence from PDB."""
    res, _, ids = read_pdb_ca(pdb_path)
    seq_aa = []
    id_to_idx = {}
    seen = set()
    for rn, rid in zip(res, ids):
        if rid.strip() not in seen:
            seen.add(rid.strip())
            id_to_idx[rid.strip()] = len(seq_aa)
            seq_aa.append(AA3TO1.get(rn, "X"))
    return "".join(seq_aa), id_to_idx

# ===== MAIN =====
print("=" * 65)
print("CORRECTED VALIDATION: DBREF-filtered sites")
print("=" * 65)

labels = load_deeppct_labels([
    "archives/subagent_artifacts/temp_deeppct/datasets/training-set.csv",
    "archives/subagent_artifacts/temp_deeppct/datasets/indepentent-test-set.csv",
])

targets = [
    ("p53",       "P04637", "data/pdb/p53_1TUP.pdb"),
    ("ER-alpha",  "P03372", "data/pdb/ER-alpha_1A52.pdb"),
    ("DNMT1",     "P26358", "data/pdb/DNMT1_3PTA.pdb"),
]

ub = build_mod("data/pdb/ubiquitin_1UBQ.pdb", "Ub")
sumo = build_mod("data/pdb/sumo1_1A5R.pdb", "SUMO")
sim = CompetitionSimulator(80.0)
ef = EnzymeFilter()

all_results = []
sim_total=sim_match=comb_total=comb_match=0

for name, uniprot, pdb_path in targets:
    res, pos, ids = read_pdb_ca(pdb_path)
    dbrefs = [d for d in parse_dbref(pdb_path) if d["uniprot"] == uniprot]
    if not dbrefs:
        print("\n--- %s: No DBREF found for %s" % (name, uniprot))
        continue
    
    seq_aa, rid_to_idx = get_pdb_lyso_sequence(pdb_path)
    
    # Find all unique K positions in PDB that map to UniProt and have labels
    seen_up = set()
    labeled_sites = []
    for rid_s, aa, p in zip(ids, res, pos):
        if aa != "LYS": continue
        rid = int(rid_s.strip()) if rid_s.strip().isdigit() else 0
        
        up = None
        for dr in dbrefs:
            if dr["pdb_start"] <= rid <= dr["pdb_end"]:
                up = rid + dr["offset"]
                break
        if up is None or up in seen_up: continue
        seen_up.add(up)
        
        key = (uniprot, up)
        if key not in labels: continue
        gt = labels[key]
        
        # Run simulation
        site = ProteinSite(rid, p, 1.0, 0.75)
        r = sim.evaluate_site(site, ub, sumo)
        
        # Sequence
        seq_idx = rid_to_idx.get(rid_s.strip())
        ss = 0.0; cp = 0.5
        if seq_idx is not None:
            start = max(0, seq_idx-7); end = min(len(seq_aa), seq_idx+8)
            window = seq_aa[start:end]; k_win = seq_idx - start
            ss = ef.score(window, k_win)
            cp = ef.prob(window, k_win, r.delta_e)
        
        # Determine ground truth
        if gt["Ub"] > 0 and gt["SUMO"] > 0: gt_label = "both"
        elif gt["Ub"] > 0: gt_label = "Ub"
        elif gt["SUMO"] > 0: gt_label = "SUMO"
        else: gt_label = "none"
        
        pred_e = r.predicted_modifier
        pred_comb = "SUMO" if cp > 0.55 else ("Ub" if cp < 0.45 else "amb")
        
        e_match = (pred_e == gt_label and gt_label in ("Ub","SUMO"))
        c_match = (pred_comb == gt_label and gt_label in ("Ub","SUMO"))
        if gt_label in ("Ub","SUMO"):
            sim_total += 1; comb_total += 1
            if e_match: sim_match += 1
            if c_match: comb_match += 1
        
        labeled_sites.append((up, r.delta_e, ss, cp, pred_e, pred_comb, gt_label, e_match, c_match, window if seq_idx is not None else "N/A"))
    
    if not labeled_sites:
        print("\n--- %s (%s): 0 in-PDB labeled K sites" % (name, uniprot))
        continue
    
    print("\n--- %s (%s) | PDB %s | %d in-PDB labeled K sites ---" % (
        name, uniprot, os.path.basename(pdb_path), len(labeled_sites)))
    
    e_ok = sum(1 for x in labeled_sites if x[7])
    c_ok = sum(1 for x in labeled_sites if x[8])
    print("  E-only: %d/%d (%.0f%%) | Combined: %d/%d (%.0f%%)" % (
        e_ok, len(labeled_sites), 100*e_ok/max(1,len(labeled_sites)),
        c_ok, len(labeled_sites), 100*c_ok/max(1,len(labeled_sites))))
    
    for up, de, ss, cp, pe, pc, gt, em, cm, win in sorted(labeled_sites):
        flag = "[OK]" if cm else ("[NO]" if gt in ("Ub","SUMO") else "[--]")
        print("  K%d: DE=%+.2f Seq=%+.1f P(S)=%.3f E->%s E+S->%s GT=%s %s | %s" % (
            up, de, ss, cp, pe, pc, gt, flag, win))

print("\n" + "=" * 65)
print("FINAL ACCURACY (DBREF-filtered)")
print("=" * 65)
print("Electrostatic only:  %.0f%% (%d/%d)" % (
    100*sim_match/sim_total if sim_total else 0, sim_match, sim_total))
print("E + Enzyme filter:   %.0f%% (%d/%d)" % (
    100*comb_match/comb_total if comb_total else 0, comb_match, comb_total))
print("Total in-PDB labeled sites: %d" % comb_total)
print()

# Summary table
print("Summary:")
for name, uniprot, pdb_path in targets:
    dbrefs = [d for d in parse_dbref(pdb_path) if d["uniprot"] == uniprot]
    if dbrefs:
        print("  %s: UNP %d-%d (%d sites)" % (
            name, dbrefs[0]["uniprot_start"], dbrefs[0]["uniprot_end"],
            sum(1 for l in all_results if l.startswith(name))))