"""Get all PDB entries for a UniProt protein via full_text search + entity alignment."""
import urllib.request, json, time

RCSB_URL = "https://data.rcsb.org/graphql"
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
HEADERS = {"Content-Type": "application/json"}

def search_pdb_by_uniprot(uniprot_id, max_rows=500):
    """Full text search for PDB entries containing this UniProt ID."""
    payload = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {"value": uniprot_id}
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": max_rows}
        }
    }
    
    req = urllib.request.Request(
        SEARCH_URL,
        data=json.dumps(payload).encode(),
        headers=HEADERS
    )
    resp = urllib.request.urlopen(req, timeout=30)
    result = json.loads(resp.read())
    
    # Parse results: identifiers are like "1TUP.1", "1TUP.2", "1TUP.3"
    entries = {}
    for item in result.get("result_set", []):
        identifier = item.get("identifier", "")
        parts = identifier.split(".")
        if len(parts) >= 2:
            entry_id = parts[0].lower()
            entity_id = parts[1]
            if entry_id not in entries:
                entries[entry_id] = []
            entries[entry_id].append(entity_id)
    
    return entries

def get_entity_info(entry_id, entity_ids):
    """Get detailed alignment info for entities in a PDB entry."""
    # Batch query: get all polymer entities at once
    # Use entry-level query with polymer_entities
    q = """
    {
        entry(entry_id: "%s") {
            polymer_entities {
                entity_poly {
                    pdbx_strand_id
                    pdbx_seq_one_letter_code_can
                }
                uniprots {
                    rcsb_id
                }
                rcsb_polymer_entity_align {
                    aligned_regions {
                        entity_beg_seq_id
                        ref_beg_seq_id
                        length
                    }
                    reference_database_accession
                }
                rcsb_id
            }
        }
    }
    """ % entry_id.upper()
    
    payload = json.dumps({"query": q})
    req = urllib.request.Request(RCSB_URL, data=payload.encode(), headers=HEADERS)
    resp = urllib.request.urlopen(req, timeout=15)
    data = json.loads(resp.read())
    
    if "errors" in data:
        print("  GraphQL error for %s: %s" % (entry_id, data["errors"][0]["message"][:60]))
        return []
    
    entities = data.get("data", {}).get("entry", {}).get("polymer_entities", [])
    results = []
    
    for ent in entities:
        rcsb_id = ent.get("rcsb_id", "")
        ent_poly = ent.get("entity_poly") or {}
        chains = ent_poly.get("pdbx_strand_id", "")
        
        aligns = ent.get("rcsb_polymer_entity_align", [])
        if not aligns:
            continue
        
        # Dedup entities we didn't ask for but are aligned to same UniProt
        # Filter to only entities with our UniProt
        for a in aligns:
            acc = a.get("reference_database_accession", "")
            if acc not in entity_ids:  # Check if this is the protein we want
                # Actually entity_ids are internal PDB entity IDs (1, 2, 3...)
                pass
            
            regions = a.get("aligned_regions", [])
            for r in regions:
                results.append({
                    "pdb_id": entry_id,
                    "entity_id": rcsb_id,
                    "chain": chains,
                    "unp_start": r.get("ref_beg_seq_id", 0),
                    "pdb_start": r.get("entity_beg_seq_id", 0),
                    "length": r.get("length", 0),
                    "accession": acc
                })
    
    return results

# Test with p53
print("Testing with P04637 (p53)...")
entries = search_pdb_by_uniprot("P04637", max_rows=20)
print("Found %d unique PDB entries" % len(entries))
for eid, eids in list(entries.items())[:5]:
    print("  %s: entities=%s" % (eid, eids))

print()
print("Getting alignment info for 1TUP...")
results = get_entity_info("1tup", ["1", "2", "3"])
for r in results:
    print("  %s chain=%s unp=%d len=%d (acc=%s)" % (
        r["pdb_id"], r["chain"], r["unp_start"], r["length"], r["accession"]))

print("\nGetting for 7BWN...")
results2 = get_entity_info("7bwn", [])
for r in results2:
    print("  %s chain=%s unp=%d len=%d (acc=%s)" % (
        r["pdb_id"], r["chain"], r["unp_start"], r["length"], r["accession"]))
