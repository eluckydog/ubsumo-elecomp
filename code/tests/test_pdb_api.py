"""Query PDBe API for UniProt-to-PDB mappings."""
import urllib.request, json, sys, time

def get_pdb_for_uniprot(uniprot_id):
    """Get PDB entries that cover this UniProt protein."""
    url = f'https://www.ebi.ac.uk/pdbe/graph-api/uniprot/{uniprot_id}'
    try:
        req = urllib.request.urlopen(url, timeout=15)
        data = json.loads(req.read())
        
        pdbs = data.get('PDB', {})
        results = []
        for pdb_id, info in pdbs.items():
            mappings = info.get('mappings', [])
            for m in mappings:
                # Get the start/end of this chain's UniProt coverage
                unp_start = m.get('unp_start', '')
                unp_end = m.get('unp_end', '')
                results.append({
                    'pdb_id': pdb_id.lower(),
                    'chain': m.get('chain_id', ''),
                    'unp_start': int(unp_start) if unp_start else 0,
                    'unp_end': int(unp_end) if unp_end else 0,
                })
        return results
    except Exception as e:
        print(f'  Error for {uniprot_id}: {e}', file=sys.stderr)
        return []

# Test with p53
print('Testing PDBe graph API...')
p53_pdbs = get_pdb_for_uniprot('P04637')
print(f'p53: {len(p53_pdbs)} PDB chain mappings')

# Sort by coverage width
p53_pdbs.sort(key=lambda x: -(x['unp_end'] - x['unp_start']))

# Show top results
seen_ids = set()
count = 0
for p in p53_pdbs:
    pid = p['pdb_id']
    if pid not in seen_ids:
        seen_ids.add(pid)
        count += 1
        print(f'  {pid} chain {p["chain"]}: UNP {p["unp_start"]}-{p["unp_end"]} ({p["unp_end"]-p["unp_start"]} aa)')
        if count >= 15:
            break

# Now test a few more and estimate total compute
test_ids = ['P04637', 'P01106', 'Q04206', 'P37840', 'P03372']
results = {}
for uid in test_ids:
    pdbs = get_pdb_for_uniprot(uid)
    print(f'\n{uid}: {len(set(p['pdb_id'] for p in pdbs))} unique PDBs')
    time.sleep(0.3)
