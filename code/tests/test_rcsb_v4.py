"""Try various query paths."""
import urllib.request, json, pprint

url = "https://data.rcsb.org/graphql"
headers = {"Content-Type": "application/json"}

# Try querying polymer_entity alignment directly
query = """
{
    polymer_entity(entity_id: 1, entry_id: "1TUP") {
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
"""
payload = json.dumps({"query": query})
req = urllib.request.Request(url, data=payload.encode(), headers=headers)
resp = urllib.request.urlopen(req, timeout=15)
data = json.loads(resp.read())

if "errors" in data:
    print("Query 1 (entity align):")
    for e in data["errors"]:
        print("  ", e["message"])
else:
    result = data.get("data", {}).get("polymer_entity", {})
    print("Query 1 result keys:", list(result.keys()))
    align = result.get("rcsb_polymer_entity_align")
    if align:
        print("  1st align:", json.dumps(align, indent=2)[:400])
    print()

# Try different approach: use entry + entity composite alignments
query2 = """
{
    polymer_entity(entity_id: 1, entry_id: "1TUP") {
        entity_poly {
            pdbx_seq_one_letter_code
        }
    }
}
"""
payload2 = json.dumps({"query": query2})
req2 = urllib.request.Request(url, data=payload2.encode(), headers=headers)
resp2 = urllib.request.urlopen(req2, timeout=15)
data2 = json.loads(resp2.read())

if "errors" in data2:
    print("Query 2 (simple entity):")
    for e in data2["errors"]:
        print("  ", e["message"])
else:
    result = data2["data"]["polymer_entity"]
    print("Query 2 result keys:", list(result.keys()))
    pp = result.get("entity_poly")
    if pp:
        print("  sequence available, len=%d" % len(pp.get("pdbx_seq_one_letter_code", "")))
print()

# Try: what does the entry-level data for 1TUP say?
query3 = """
{
    entry(entry_id: "1TUP") {
        rcsb_entry_info {
            structure_determination_methodology
            resolution_combined
            molecular_weight
        }
        struct_keywords {
            pdbx_keywords
        }
    }
}
"""
payload3 = json.dumps({"query": query3})
req3 = urllib.request.Request(url, data=payload3.encode(), headers=headers)
resp3 = urllib.request.urlopen(req3, timeout=15)
data3 = json.loads(resp3.read())

if "errors" in data3:
    print("Query 3 (entry info):")
    for e in data3["errors"]:
        print("  ", e["message"])
else:
    result = data3["data"]["entry"]
    print("Query 3 - 1TUP resolution: %.2f A" % result.get("rcsb_entry_info", {}).get("resolution_combined", [0])[0])
