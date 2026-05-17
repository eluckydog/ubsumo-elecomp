"""Batch download all PDBs in parallel batches, 300 per batch."""
import urllib.request, json, os, time, sys
import concurrent.futures

WORKSPACE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PDB_DIR = os.path.join(WORKSPACE, "data", "pdb_all")
MAPPING_PATH = os.path.join(WORKSPACE, "data", "uniprot_pdb_mapping.json")
os.makedirs(PDB_DIR, exist_ok=True)

m = json.load(open(MAPPING_PATH))
all_pdbs = sorted(set(e["pdb_id"].lower() for entries in m.values() for e in entries))
print("Total PDB IDs: %d" % len(all_pdbs))

existing = set(f.replace(".pdb","").replace(".cif","") for f in os.listdir(PDB_DIR))
to_dl = [p for p in all_pdbs if p not in existing]
print("Already cached: %d" % len(existing))
print("To download: %d" % len(to_dl))

if not to_dl:
    print("All done!")
    sys.exit(0)

def download_one(pid):
    path = os.path.join(PDB_DIR, "%s.pdb" % pid)
    try:
        url = "https://files.rcsb.org/download/%s.pdb" % pid.upper()
        req = urllib.request.urlopen(url, timeout=30)
        content = req.read()
        with open(path, "wb") as f:
            f.write(content)
        return (pid, len(content), None)
    except Exception as e1:
        try:
            url = "https://files.rcsb.org/download/%s.cif" % pid.upper()
            req = urllib.request.urlopen(url, timeout=30)
            content = req.read()
            with open(os.path.join(PDB_DIR, "%s.cif" % pid), "wb") as f:
                f.write(content)
            return (pid, len(content), None)
        except:
            return (pid, 0, str(e1)[:40])

BATCH = 300
total_bytes = 0
success = 0
fail = 0

for batch_start in range(0, len(to_dl), BATCH):
    batch = to_dl[batch_start:batch_start+BATCH]
    print("Batch %d/%d (%d PDBs)..." % (batch_start//BATCH+1, (len(to_dl)+BATCH-1)//BATCH, len(batch)), flush=True)
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        results = list(ex.map(download_one, batch))
    
    for pid, size, err in results:
        if err is None:
            success += 1
            total_bytes += size
        else:
            fail += 1
    
    print("  Progress: %d/%d downloaded, %.1f MB, %d failed" % (
        success, len(to_dl), total_bytes/1024/1024, fail), flush=True)
    time.sleep(1)  # pause between batches

total = total_bytes/1024/1024
print("\nDone! Downloaded: %d/%d - Failed: %d - Total: %.1f MB" % (
    success, len(to_dl), fail, total))

# Verification
pdb_files = [f for f in os.listdir(PDB_DIR) if f.endswith(".pdb")]
cif_files = [f for f in os.listdir(PDB_DIR) if f.endswith(".cif")]
final_size = sum(os.path.getsize(os.path.join(PDB_DIR,f)) for f in os.listdir(PDB_DIR))
print("PDB: %d, CIF: %d, Total: %d files, %.1f MB" % (
    len(pdb_files), len(cif_files), len(pdb_files)+len(cif_files), final_size/1024/1024))
