"""Fix entity_id type and try RcsbPolymerEntityAlign."""
import urllib.request, json, pprint

url = "https://data.rcsb.org/graphql"
headers = {"Content-Type": "application/json"}

# entity_id is a String, not Int
query = """
{
    polymer_entity(entity_id: "1", entry_id: "1TUP") {
        rcsb_polymer_entity_align {
            aligned_regions {
                entity_beg_seq_id
                ref_beg_seq_id
                length
            }
            reference_database_accession
            reference_database_name
        }
        entity_poly {
            pdbx_seq_one_letter_code
            pdbx_seq_one_letter_code_can
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
    result = data.get("data", {}).get("polymer_entity", {})
    print("Result keys:", list(result.keys()))
    
    align = result.get("rcsb_polymer_entity_align")
    if align:
        print("\nrcsb_polymer_entity_align:")
        if isinstance(align, list):
            print("  type: list, length:", len(align))
            if align:
                a = align[0]
                print("  keys:", list(a.keys()))
                print("  db:", a.get("reference_database_name"), a.get("reference_database_accession"))
                regions = a.get("aligned_regions", [])
                print("  aligned_regions:", len(regions))
                for r in regions[:3]:
                    print("    entity_beg=%d, ref_beg=%d, len=%d" % (
                        r["entity_beg_seq_id"], r["ref_beg_seq_id"], r["length"]))
        elif isinstance(align, dict):
            print("  type: dict, keys:", list(align.keys()))
        else:
            print("  type:", type(align).__name__, "val:", align)
    
    seq = result.get("entity_poly", {}).get("pdbx_seq_one_letter_code", "")
    if seq:
        print("Sequence: len=%d, first 50: %s" % (len(seq), seq[:50]))
