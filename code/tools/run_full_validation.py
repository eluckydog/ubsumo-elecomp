"""Full-Map Validation: DBREF-based, electrostatic + enzyme filter across all 54 proteins."""
import sys, os, re, math, csv, json
import numpy as np
from collections import defaultdict

WORKSPACE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, WORKSPACE)

from ubsumo_em_simulator.core.competition import CompetitionSimulator, ProteinSite, ModifierProtein
from ubsumo_em_simulator.io.pdb_reader import assign_charges

DATA_DIR = os.path.join(WORKSPACE, "data")
PDB_DIR = os.path.join(DATA_DIR, "pdb")
MAPPING_PATH = os.path.join(DATA_DIR, "uniprot_pdb_mapping.json")

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
    ranges = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("DBREF") and "UNP" in line:
                parts = line.split()
                if len(parts) >= 10:
                    ranges.append({
                        "pdb_chain": parts[2], "pdb_start": int(parts[3]), "pdb_end": int(parts[4]),
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
    return ModifierProtein(name, np.mean(pos, axis=0), list(ch), pos)

def load_labels():
    PTM_NORM = {"ubiquitination":"Ub","Ubiquitination":"Ub",
                "sumoylation":"SUMO","Sumoylation":"SUMO","SUMOylation":"SUMO"}
    labels = defaultdict(lambda: {"Ub":0,"SUMO":0})
    for csv_path in ["training-set.csv", "indepentent-test-set.csv"]:
        fp = os.path.join(WORKSPACE, "archives", "subagent_artifacts", "temp_deeppct", "datasets", csv_path)
        with open(fp, "r", encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                for sc, pc in [("Site1","PTM1"),("Site2","PTM2")]:
                    ptm = PTM_NORM.get(row[pc].strip(), "")
                    site = row[sc].strip()
                    if ptm in ("Ub","SUMO") and site.startswith("K") and len(site) > 1:
                        key = (row["UniProt ID"].strip(), int(site[1:]))
                        labels[key][ptm] += 1
    return labels

# Load
print("=" * 65)
print("FULL VALIDATION: %d proteins, local PDB cache" % 54)
print("=" * 65)

labels = load_labels()
mapping = json.load(open(MAPPING_PATH))

ub = build_mod(os.path.join(PDB_DIR, "ubiquitin_1UBQ.pdb"), "Ub")
sumo = build_mod(os.path.join(PDB_DIR, "sumo1_1A5R.pdb"), "SUMO")
sim = CompetitionSimulator(80.0)
ef = EnzymeFilter()

all_results = []
sim_total=sim_match=comb_total=comb_match=0

for uid in sorted(labels.keys()):
    # Get all (pdb, chain, range) that cover K sites
    ks_all = set(k for (u,k),v in labels.items() if u == uid)
    
    pdb_cover_maps = []  # (pdb_file, unp_start, unp_end, covered_Ks)
    for entry in mapping.get(uid, []):
        pid = entry["pdb_id"]
        for ch_info in entry.get("chains", []):
            s, e = ch_info.get("unp_start"), ch_info.get("unp_end")
            if s and e:
                covered = {k for k in ks_all if s <= k <= e}
                if covered:
                    pdb_path = os.path.join(PDB_DIR, "%s.pdb" % pid.lower())
                    if os.path.exists(pdb_path):
                        pdb_cover_maps.append((pdb_path, pid, ch_info["chains"], s, e, covered))
    
    if not pdb_cover_maps:
        continue
    
    # Pick best PDB per site (greedy, can use multiple PDBs)
    all_covered = set()
    used_pdbs = []
    pdb_cover_maps.sort(key=lambda x: -len(x[5]))
    
    for pdb_path, pid, chains, s, e, covered in pdb_cover_maps:
        # Read PDB, find which K sites we can actually evaluate
        try:
            res, pos, ids = read_pdb_ca(pdb_path)
            dbrefs = parse_dbref(pdb_path)
            
            # Find uniprot range from DBREF
            uniprot_dbrefs = [d for d in dbrefs if d["uniprot"] == uid]
            if not uniprot_dbrefs:
                continue
            
            for k in covered:
                if k in all_covered:
                    continue
                
                # Find PDB residue ID for this UniProt position
                pdb_residue = None
                for dr in uniprot_dbrefs:
                    # PDB position = k - offset
                    pdb_pos = k - dr["offset"]
                    if dr["pdb_start"] <= pdb_pos <= dr["pdb_end"]:
                        pdb_residue = pdb_pos
                        break
                
                if pdb_residue is None:
                    continue
                
                # Find this residue in PDB CA atoms
                for rid_s, aa, p in zip(ids, res, pos):
                    if rid_s.strip() == str(pdb_residue) or rid_s.strip().lstrip("ABCDEFGH") == str(pdb_residue):
                        # Actually residue IDs in PDB have insertion codes
                        rid_num = ''.join(c for c in rid_s.strip() if c.isdigit() or c == '-')
                        try:
                            if rid_num and int(rid_num) == pdb_residue:
                                if aa == "LYS":
                                    # Found it!
                                    site = ProteinSite(pdb_residue, p, 1.0, 0.75)
                                    r = sim.evaluate_site(site, ub, sumo)
                                    
                                    gt = labels[(uid, k)]
                                    if gt["Ub"] > 0 and gt["SUMO"] > 0: gt_label = "both"
                                    elif gt["Ub"] > 0: gt_label = "Ub"
                                    elif gt["SUMO"] > 0: gt_label = "SUMO"
                                    else: gt_label = "none"
                                    
                                    all_covered.add(k)
                                    used_pdbs.append((uid, k, pdb_path, pid, pdb_residue, r.delta_e, r.predicted_modifier, gt_label))
                        except ValueError:
                            pass
                    if k in all_covered:
                        break
        
        except Exception as e:
            pass
    
    if used_pdbs:
        print("  %s: evaluated %d/%d K sites via PDB" % (
            uid, len(all_covered), len(ks_all)))

# Print results
print("\nSites evaluated from PDB: %d" % len(all_results))
# ... (rest of analysis coming)
