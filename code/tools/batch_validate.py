"""Batch validate: for all proteins with local PDBs, extract K site positions and run simulation."""
import sys, os, re, math, csv, json
import numpy as np
from collections import defaultdict

WORKSPACE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, WORKSPACE)

from ubsumo_em_simulator.core.competition import CompetitionSimulator, ProteinSite, ModifierProtein
from ubsumo_em_simulator.io.pdb_reader import assign_charges

PDB_DIR = os.path.join(WORKSPACE, "data", "pdb")
AA3TO1 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
          "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
          "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
A1TO3 = {v:k for k,v in AA3TO1.items()}

def read_all_atoms(pdb_path):
    """Read CA atoms + residue info from PDB."""
    atoms = []  # (resname, chain, resnum, icode, x, y, z)
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("MODEL") and not line.startswith("MODEL        1"):
                break  # only first model
            if line.startswith("ENDMDL") or line.startswith("END"):
                continue
            if line.startswith("ATOM") and line[12:16].strip() in ("CA", "P"):
                try:
                    x=float(line[30:38]); y=float(line[38:46]); z=float(line[46:54])
                except: continue
                resname = line[17:20].strip()
                chain = line[21].strip()
                resnum = line[22:26].strip()
                icode = line[26].strip()
                atoms.append((resname, chain, resnum, icode, x, y, z))
    return atoms

def parse_dbrefs(pdb_path):
    """Parse DBREF + SEQRES to understand PDB->UniProt mapping."""
    db_wanted = {}
    seqres_chain = {}
    
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("DBREF"):
                parts = line.split()
                if len(parts) >= 10:
                    chain = parts[2]
                    pdb_start, pdb_end = int(parts[3]), int(parts[4])
                    uniprot = parts[6]
                    unp_start, unp_end = int(parts[8]), int(parts[9])
                    db_wanted[chain] = {
                        "pdb_start": pdb_start, "pdb_end": pdb_end,
                        "uniprot": uniprot, "unp_start": unp_start, "unp_end": unp_end,
                        "offset": unp_start - pdb_start
                    }
    return db_wanted

def build_modifier(name, pdb_path):
    """Build Ub or SUMO modifier from PDB."""
    atoms = read_all_atoms(pdb_path)
    resnames = [a[0] for a in atoms]
    positions = [np.array([a[4], a[5], a[6]]) for a in atoms]
    charges = assign_charges(resnames)
    center = np.mean(positions, axis=0)
    return ModifierProtein(name, center, list(charges), positions)

def laod_labels():
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

labels = laod_labels()
mapping = json.load(open(os.path.join(WORKSPACE, "data", "uniprot_pdb_mapping.json")))

ub = build_modifier("Ub", os.path.join(PDB_DIR, "ubiquitin_1UBQ.pdb"))
sumo = build_modifier("SUMO", os.path.join(PDB_DIR, "sumo1_1A5R.pdb"))
sim = CompetitionSimulator(80.0)

print("=" * 60)
print("FULL PDB-VALIDATION")
print("=" * 60)

total_all = 0
total_eval = 0
match_e = 0
match_e_sumo = 0
match_e_ub = 0
found = 0
not_found = 0

for uid in sorted(set(k[0] for k in labels)):
    # Collect all labeled K sites for this protein
    ks_all = {}
    for (u,k), v in labels.items():
        if u == uid:
            ks_all[k] = v
    
    # Find PDBs with range info from mapping
    covering_pdbs = []
    for entry in mapping.get(uid, []):
        pid = entry["pdb_id"].lower()
        pdb_path = os.path.join(PDB_DIR, "%s.pdb" % pid)
        if not os.path.exists(pdb_path):
            continue
        for ch_info in entry.get("chains", []):
            s, e = ch_info.get("unp_start"), ch_info.get("unp_end")
            if s and e:
                ks_covered = {k for k in ks_all if s <= k <= e}
                if ks_covered:
                    covering_pdbs.append((pdb_path, pid, ch_info["chains"], s, e, ks_covered))
    
    if not covering_pdbs:
        continue
    
    # For each site, try to find it in a PDB
    found_list = []
    missing_list = []
    
    for k in sorted(ks_all.keys()):
        k_found = False
        
        for pdb_path, pid, chain_info, s, e, ks_covered in sorted(covering_pdbs, key=lambda x: -len(x[5])):
            if k not in ks_covered:
                continue
            
            dbrefs = parse_dbrefs(pdb_path)
            if not dbrefs:
                continue
            
            # For each DBREF chain, check if this PDB covers our K
            for chain_label, dbref in dbrefs.items():
                if dbref["uniprot"] != uid:
                    continue
                if not (dbref["unp_start"] <= k <= dbref["unp_end"]):
                    continue
                
                # PDB residue number
                pdb_res = k - dbref["offset"]
                
                # Read CA atoms and find K at this position
                atoms = read_all_atoms(pdb_path)
                for resname, chain, resnum, icode, x, y, z in atoms:
                    if chain != chain_label:
                        continue
                    try:
                        rn = int(resnum)
                    except:
                        continue
                    if rn == pdb_res and resname == "LYS":
                        k_found = True
                        pos = np.array([x, y, z])
                        site = ProteinSite(k, pos, 1.0, 0.75)
                        
                        try:
                            r = sim.evaluate_site(site, ub, sumo)
                            
                            gt = ks_all[k]
                            if gt["Ub"] > 0 and gt["SUMO"] > 0: gt_label = "both"
                            elif gt["Ub"] > 0: gt_label = "Ub"
                            elif gt["SUMO"] > 0: gt_label = "SUMO"
                            else: gt_label = "none"
                            
                            found_list.append((uid, k, r.delta_e, r.predicted_modifier, gt_label, pid, chain_label, pdb_res))
                            total_eval += 1
                            
                            if gt_label == "Ub" and r.predicted_modifier == "Ub":
                                match_e_ub += 1
                            elif gt_label == "SUMO" and r.predicted_modifier == "SUMO":
                                match_e_sumo += 1
                            if (gt_label == "Ub" and r.predicted_modifier == "Ub") or \
                               (gt_label == "SUMO" and r.predicted_modifier == "SUMO"):
                                match_e += 1
                            
                        except Exception as ex:
                            pass
                        break
                if k_found:
                    break
            if k_found:
                break
        
        if not k_found:
            missing_list.append((uid, k))
    
    total_all += len(ks_all)
    found += len(found_list)
    not_found += len(missing_list)
    
    if found_list:
        n_ub = sum(1 for f in found_list if f[4] == "Ub")
        n_sumo = sum(1 for f in found_list if f[4] == "SUMO")
        n_both = sum(1 for f in found_list if f[4] == "both")
        m = sum(1 for f in found_list if f[3] == f[4] and f[4] in ("Ub","SUMO"))
        d = sum(1 for f in found_list if f[4] in ("Ub","SUMO"))
        print("  %s: %d/%d K from PDB, E-only=%d/%d=%.0f%% (Ub=%d SUMO=%d both=%d)" % (
            uid, len(found_list), len(ks_all), m, d, 100*m/d if d else 0,
            n_ub, n_sumo, n_both))
        
        # Show each site
        for uk, k, de, pred, gt, pid, ch, pr in found_list:
            m_flag = "[OK]" if (pred == gt and gt in ("Ub","SUMO")) else "[NO]"
            print("    K%d: DE=%+.2f pred=%s GT=%s %s (PDB %s %s@%d)" % (
                k, de, pred, gt, m_flag, pid, ch, pr))
    else:
        print("  %s: 0/%d K found in PDB" % (uid, len(ks_all)))

print()
print("=" * 60)
print("SUMMARY")
print("=" * 60)
print("Total labeled K sites: %d" % total_all)
print("Found in PDB: %d (%.0f%%)" % (found, 100*found/total_all))
print("Not found in PDB: %d" % not_found)
print()
d_total = match_e + (total_eval - match_e - match_e)
print("E-only accuracy: %.0f%% (%d/%d clear-label sites)" % (
    100*match_e/total_eval if total_eval else 0, match_e, total_eval))
print("  Ub correct: %d/%d" % (match_e_ub, total_eval))
print("  SUMO correct: %d/%d" % (match_e_sumo, total_eval))
