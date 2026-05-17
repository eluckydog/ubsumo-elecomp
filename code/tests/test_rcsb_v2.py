"""Get UniProt-to-PDB mapping via RCSB GraphQL."""
import urllib.request, json, time

url = 'https://data.rcsb.org/graphql'

# Query with correct arg name
query = """
{
    uniprot(uniprot_id: "P04637") {
        rcsb_id
        rcsb_uniprot_container_identifiers {
            uniprot_id
        }
        rcsb_uniprot_alignments {
            ... on RcsbUniprotAlignments {
                pdb_chain_alignments {
                    align_end
                    align_start
                    pdb_model_num
                    pdb_chain_id
                    pdb_id
                }
                query_begin
                query_end
                sequence_identity
            }
        }
    }
}
"""

payload = json.dumps({'query': query}).encode()
headers = {'Content-Type': 'application/json'}
req = urllib.request.Request(url, data=payload, headers=headers)
resp = urllib.request.urlopen(req, timeout=15)
data = json.loads(resp.read())

if 'errors' in data:
    print('Errors:', data['errors'])
else:
    result = data.get('data', {}).get('uniprot', {})
    pdb_id = result.get('rcsb_id')
    print(f'UniProt ID: {pdb_id}')
    
    alignments = result.get('rcsb_uniprot_alignments', [])
    print(f'Alignment groups: {len(alignments)}')
    
    all_pdb_links = []
    for align_group in alignments:
        if not isinstance(align_group, dict): continue
        query_begin = align_group.get('query_begin', '?')
        query_end = align_group.get('query_end', '?')
        seq_id = align_group.get('sequence_identity', '?')
        
        pdb_aligns = align_group.get('pdb_chain_alignments', [])
        for pa in pdb_aligns:
            all_pdb_links.append({
                'pdb_id': pa.get('pdb_id'),
                'chain': pa.get('pdb_chain_id'),
                'align_start': pa.get('align_start'),
                'align_end': pa.get('align_end'),
                'unp_start': query_begin,
                'unp_end': query_end,
                'seq_id': seq_id
            })
    
    print(f'Total PDB chain alignments: {len(all_pdb_links)}')
    
    # Dedup by PDB ID and show coverage
    by_pdb = {}
    for link in all_pdb_links:
        pid = link['pdb_id']
        if pid not in by_pdb:
            by_pdb[pid] = []
        by_pdb[pid].append(link)
    
    print(f'Unique PDB entries: {len(by_pdb)}')
    print()
    
    # Show top 10 by coverage
    sorted_pdbs = sorted(by_pdb.items(), key=lambda x: -max(
        a['align_end'] - a['align_start'] for a in x[1]
    ))
    for pid, links in sorted_pdbs[:10]:
        max_range = max(l['align_end'] - l['align_start'] for l in links)
        chains = [l['chain'] for l in links]
        print(f'  {pid} (chains={",".join(chains[:3])}, range={max_range}aa, seq_id={links[0]["seq_id"]})')
