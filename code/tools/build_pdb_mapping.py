"""Batch download all PDB→UniProt mapping data to local JSON."""
import urllib.request, json, time, os, csv
from collections import defaultdict

WORKSPACE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(WORKSPACE, "data")
MAPPING_PATH = os.path.join(DATA_DIR, "uniprot_pdb_mapping.json")

def load_deeppct_proteins():
    PTM_NORM = {"ubiquitination":"Ub","Ubiquitination":"Ub",
                "sumoylation":"SUMO","Sumoylation":"SUMO","SUMOylation":"SUMO"}
    labels = defaultdict(lambda: {"Ub":set(), "SUMO":set()})
    csv_dir = os.path.join(WORKSPACE, "archives", "subagent_artifacts", "temp_deeppct", "datasets")
    for csv_name in ["training-set.csv", "indepentent-test-set.csv"]:
        fp = os.path.join(csv_dir, csv_name)
        with open(fp, "r", encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                for sc, pc in [("Site1","PTM1"),("Site2","PTM2")]:
                    ptm = PTM_NORM.get(row[pc].strip(), "")
                    site = row[sc].strip()
                    if ptm in ("Ub","SUMO") and site.startswith("K") and len(site) > 1:
                        key = row["UniProt ID"].strip()
                        labels[key][ptm].add(int(site[1:]))
    return labels, sorted(labels.keys())

# Step 1: Download UniProt→PDB mapping via UniProt REST API
# /uniprotkb/{accession}.json returns all cross-references including PDB
labels, all_ids = load_deeppct_proteins()
print("Loading PDB cross-references for %d UniProt proteins..." % len(all_ids))

mapping = {}  # {uniprot: [{pdb_id, resolution, chain, unp_start, unp_end}, ...]}

URL_TPL = "https://rest.uniprot.org/uniprotkb/%s.json"

for i, uid in enumerate(all_ids):
    print("[%d/%d] %s ... " % (i+1, len(all_ids), uid), end="", flush=True)
    url = URL_TPL % uid
    
    try:
        req = urllib.request.urlopen(url, timeout=15)
        data = json.loads(req.read())
        
        # Extract PDB cross-references with UniProt alignment info
        xrefs = data.get("uniProtKBCrossReferences", [])
        pdb_entries = [x for x in xrefs if x.get("database") == "PDB"]
        
        entry_list = []
        for entry in pdb_entries:
            props = {p["key"]: p["value"] for p in entry.get("properties", [])}
            pdb_id = entry["id"]
            
            # Extract chain info from properties
            # Typical PDB properties: "Chains", "Resolution", "Method"
            chains = props.get("Chains", "")
            resolution = props.get("Resolution", "")
            method = props.get("Method", "")
            
            # Parse chains: format "A/B/C=94-312"
            chain_info = []
            for part in chains.split(","):
                part = part.strip()
                if "=" in part:
                    chain_ids, range_str = part.split("=", 1)
                    if "-" in range_str:
                        try:
                            s, e = range_str.split("-")
                            chain_info.append({
                                "chains": chain_ids.strip(),
                                "unp_start": int(s) if s else None,
                                "unp_end": int(e) if e else None
                            })
                        except:
                            pass
            
            entry_list.append({
                "pdb_id": pdb_id,
                "resolution": resolution,
                "method": method,
                "chains": chain_info
            })
        
        mapping[uid] = entry_list
        print("%d PDBs" % len(entry_list))
        
    except Exception as e:
        print("FAILED: %s" % str(e)[:40])
        mapping[uid] = []
    
    if (i+1) % 10 == 0:
        with open(MAPPING_PATH + ".partial", "w") as f:
            json.dump(mapping, f)
    
    time.sleep(0.15)

# Save final
os.makedirs(DATA_DIR, exist_ok=True)
with open(MAPPING_PATH, "w") as f:
    json.dump(mapping, f, indent=1)

total_pdb = sum(len(v) for v in mapping.values())
with_pdb = sum(1 for v in mapping.values() if v)
print("\nDone! %s" % MAPPING_PATH)
print("Proteins with PDB: %d/%d" % (with_pdb, len(all_ids)))
print("Total PDB entries: %d" % total_pdb)

# Print coverage summary
for uid in all_ids:
    k_sites = labels[uid]["Ub"] | labels[uid]["SUMO"]
    if not k_sites:
        continue
    min_k = min(k_sites)
    max_k = max(k_sites)
    
    # Count PDBs that cover any of our K sites
    covering = []
    for entry in mapping.get(uid, []):
        for chain in entry.get("chains", []):
            s = chain.get("unp_start")
            e = chain.get("unp_end")
            if s and e:
                covered = [k for k in k_sites if s <= k <= e]
                if covered:
                    covering.append((entry["pdb_id"], entry["resolution"], entry["method"],
                                     chain["chains"], s, e, len(covered), sorted(covered)))
    
    if covering:
        covering.sort(key=lambda x: -x[6])  # by coverage count
        best = covering[0]
        print("  %s %s: PDBs=%d, best=%s (res=%s, %s, chains=%s, unp=%d-%d, covers %d/%d K=%s)" % (
            "✓" if best[6] >= len(k_sites) else " ",
            uid, len(covering), best[0], best[1], best[2], best[3],
            best[4], best[5], best[6], len(k_sites), best[7]))
    else:
        print("  %s %s: NO covering PDB (%d K sites: %d-%d)" % ("✗", uid, len(k_sites), min_k, max_k))
