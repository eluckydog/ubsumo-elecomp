"""Find entity IDs for chains."""
import urllib.request, json

url = "https://data.rcsb.org/graphql"
headers = {"Content-Type": "application/json"}

# Query entry for non-polymer and polymer entities
query = """
{
    entry(entry_id: "1TUP") {
        polymer_entities {
            entity_poly {
                rcsb_id
                pdbx_strand_id
                pdbx_seq_one_letter_code_can
            }
            rcsb_entity_polymer_type
            rcsb_polymer_entity_align {
                aligned_regions {
                    entity_beg_seq_id
                    ref_beg_seq_id
                    length
                }
                reference_database_accession
                reference_database_name
            }
        }
    }
}
"""
payload = json.dumps({"query": query})
req = urllib.request.Request(url, data=payload.encode(), headers=headers)
resp = urllib.request.urlopen(req, timeout=15)
data = json.loads(resp.read())

if "errors" in data:
    print("ERRORS:")
    for e in data["errors"]:
        print("  ", e["message"])
else:
    entities = data.get("data", {}).get("entry", {}).get("polymer_entities", [])
    print("Polymer entities in 1TUP:")
    for i, ent in enumerate(entities):
        poly = ent.get("entity_poly", {})
        entity_id = poly.get("rcsb_id", "?")
        chains = poly.get("pdbx_strand_id", "")
        seq = poly.get("pdbx_seq_one_letter_code_can", "")[:30]
        
        align = ent.get("rcsb_polymer_entity_align", [])
        db_match = ""
        if align and isinstance(align, list) and align:
            a = align[0]
            db_match = "%s (accession: %s)" % (a.get("reference_database_name",""), a.get("reference_database_accession",""))
            regions = a.get("aligned_regions", [])
            if regions:
                r = regions[0]
                db_match += " | entity_beg=%d ref_beg=%d len=%d" % (r["entity_beg_seq_id"], r["ref_beg_seq_id"], r["length"])
        
        print("  Entity %s: chains=%s | seq=%s... | %s" % (entity_id, chains, seq, db_match))
