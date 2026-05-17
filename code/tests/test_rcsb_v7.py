"""Query polymer_entities direct."""
import urllib.request, json

url = "https://data.rcsb.org/graphql"
headers = {"Content-Type": "application/json"}

query = """
{
    entry(entry_id: "1TUP") {
        polymer_entities {
            entity_poly {
                pdbx_strand_id
                pdbx_seq_one_letter_code_can
                pdbx_seq_one_letter_code
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
                reference_database_name
            }
            rcsb_id
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
    print("Entities in 1TUP: %d\n" % len(entities))
    
    for i, ent in enumerate(entities):
        rcsb_id = ent.get("rcsb_id", "?")
        poly = ent.get("entity_poly", {})
        chains = poly.get("pdbx_strand_id", "")
        seq = (poly.get("pdbx_seq_one_letter_code_can") or poly.get("pdbx_seq_one_letter_code") or "")[:40]
        align = ent.get("rcsb_polymer_entity_align", [])
        
        print("[%d] rcsb_id=%s chains=%s seq=%s..." % (i, rcsb_id, chains, seq))
        
        uniprots_list = ent.get("uniprots") or []
        uniprots = [u["rcsb_id"] for u in uniprots_list if u.get("rcsb_id")]
        if uniprots:
            print("    UniProt: %s" % ", ".join(uniprots))
        
        if align:
            for a in align:
                db = "%s: %s" % (a.get("reference_database_name",""), a.get("reference_database_accession",""))
                for r in a.get("aligned_regions", []):
                    print("    Align: %s | pdb_beg=%d uniprot_beg=%d len=%d" % (
                        db, r["entity_beg_seq_id"], r["ref_beg_seq_id"], r["length"]))
        print()
