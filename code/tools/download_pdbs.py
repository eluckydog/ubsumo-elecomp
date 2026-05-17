"""Download covering PDBs to local cache."""
import urllib.request, json, csv, os, sys, time
from collections import defaultdict

WORKSPACE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(WORKSPACE, "data")
PDB_DIR = os.path.join(DATA_DIR, "pdb")
MAPPING_PATH = os.path.join(DATA_DIR, "uniprot_pdb_mapping.json")

os.makedirs(PDB_DIR, exist_ok=True)

# Load mapping
m = json.load(open(MAPPING_PATH))

# Load labels
PTM_NORM = {'ubiquitination':'Ub','Ubiquitination':'Ub','sumoylation':'SUMO','Sumoylation':'SUMO','SUMOylation':'SUMO'}
labels = defaultdict(lambda: {'Ub':set(), 'SUMO':set()})
for csv_path in ['training-set.csv', 'indepentent-test-set.csv']:
    fp = os.path.join(WORKSPACE, 'archives', 'subagent_artifacts', 'temp_deeppct', 'datasets', csv_path)
    with open(fp, 'r', encoding='utf-8-sig') as f:
        for row in csv.DictReader(f):
            for sc,pc in [('Site1','PTM1'),('Site2','PTM2')]:
                ptm = PTM_NORM.get(row[pc].strip(), '')
                site = row[sc].strip()
                if ptm in ('Ub','SUMO') and site.startswith('K'):
                    labels[row['UniProt ID'].strip()][ptm].add(int(site[1:]))

# Build download list
download_list = []
for uid in sorted(labels.keys()):
    ks = labels[uid]['Ub'] | labels[uid]['SUMO']
    covering = []
    for entry in m.get(uid, []):
        for ch_info in entry.get('chains', []):
            s, e = ch_info.get('unp_start'), ch_info.get('unp_end')
            if s and e:
                c = {k for k in ks if s <= k <= e}
                if c:
                    covering.append((entry['pdb_id'], ch_info['chains'], s, e, sorted(c)))
    
    # Dedup by PDB
    dedup = {}
    for pid, ch, s, e, c in covering:
        if pid not in dedup or len(c) > len(dedup[pid][4]):
            dedup[pid] = (pid, ch, s, e, c)
    
    # Download best PDBs that incrementally cover new sites
    cum_covered = set()
    for pid, ch, s, e, c in sorted(dedup.values(), key=lambda x: -len(x[4])):
        new = [k for k in c if k not in cum_covered]
        if new:
            download_list.append((pid, uid, ch, s, e, new))
            cum_covered.update(c)
            if len(cum_covered) >= len(ks):
                break

print("Download plan: %d PDBs to download\n" % len(download_list))

# Download
success = 0
fail = 0
existing = 0

for i, (pid, uid, ch, s, e, covered_k) in enumerate(download_list):
    pdb_path = os.path.join(PDB_DIR, "%s.pdb" % pid.lower())
    
    if os.path.exists(pdb_path) and os.path.getsize(pdb_path) > 1000:
        existing += 1
        print("[%d/%d] %s.pdb (cached, %.0f KB)" % (i+1, len(download_list), pid, os.path.getsize(pdb_path)/1024))
        continue
    
    # Download from RCSB PDB
    # Use the PDB REST API: https://files.rcsb.org/download/{PDBID}.pdb
    url = "https://files.rcsb.org/download/%s.pdb" % pid.upper()
    
    try:
        req = urllib.request.urlopen(url, timeout=30)
        content = req.read()
        
        with open(pdb_path, "wb") as f:
            f.write(content)
        
        size_kb = len(content) / 1024
        print("[%d/%d] %s.pdb downloaded (%d KB)" % (i+1, len(download_list), pid, size_kb))
        success += 1
        
    except Exception as e:
        print("[%d/%d] %s.pdb FAILED: %s" % (i+1, len(download_list), pid, str(e)[:40]))
        fail += 1
    
    # Rate limit
    time.sleep(0.2)

print("\nDone!")
print("Downloaded: %d | Existing: %d | Failed: %d" % (success, existing, fail))
print("Total PDBs in cache: %d" % len([f for f in os.listdir(PDB_DIR) if f.endswith('.pdb')]))
print("PDB cache size: %.1f MB" % (sum(os.path.getsize(os.path.join(PDB_DIR, f)) for f in os.listdir(PDB_DIR) if f.endswith('.pdb')) / 1024 / 1024))
