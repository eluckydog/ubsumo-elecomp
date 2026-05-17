"""Test RCSB GraphQL API for UniProt-to-PDB mapping."""
import urllib.request, json, sys

url = 'https://data.rcsb.org/graphql'

query_str = """
{
    uniprot(accession: "P04637") {
        rcsb_uniprot_container_identifiers {
            uniprot_id
        }
        rcsb_uniprot_related_protein_described_pdb_ids
    }
}
"""

try:
    payload = json.dumps({'query': query_str}).encode()
    headers = {'Content-Type': 'application/json'}
    req = urllib.request.Request(url, data=payload, headers=headers)
    resp = urllib.request.urlopen(req, timeout=15)
    data = json.loads(resp.read())
    
    if 'errors' in data:
        print('Errors:', data['errors'])
    else:
        result = data.get('data', {}).get('uniprot', {})
        print('Result keys:', list(result.keys()))
        pdbs = result.get('rcsb_uniprot_related_protein_described_pdb_ids', [])
        print(f'p53 PDB entries: {len(pdbs)}')
        print('Sample:', pdbs[:5])
        
        # Also get residue range coverage info
        # We need a more detailed query to get chain coverage
        print()
        print('Now querying for chain coverage...')
        
except Exception as e:
    print('GraphQL error:', e)
    import traceback; traceback.print_exc()

# Second query: get detailed coverage info for each PDB
query2 = """
{
    uniprot(accession: "P04637") {
        rcsb_uniprot_related_protein_described_pdb_ids
        rcsb_entry_info_derived_pdb_id
    }
}
"""
