"""Debug RCSB GraphQL response for P04637."""
import urllib.request, json, pprint

url = "https://data.rcsb.org/graphql"
headers = {"Content-Type": "application/json"}

# Simplified query - avoid aliasing issue
query = """
{
    uniprot(uniprot_id: "P04637") {
        rcsb_uniprot_alignments {
            core_entity_alignments {
                aligned_regions {
                    query_begin
                    target_begin
                    length
                }
                core_entity_identifiers {
                    entry_id
                    entity_id
                }
                scores {
                    query_coverage
                    target_coverage
                }
            }
        }
    }
}
"""

payload = json.dumps({"query": query})
print("Query:")
print(payload[:500])
print()

req = urllib.request.Request(url, data=payload.encode(), headers=headers)
resp = urllib.request.urlopen(req, timeout=15)
data = json.loads(resp.read())

if "errors" in data:
    print("ERRORS:")
    for e in data["errors"]:
        print(" ", e["message"])
else:
    result = data.get("data", {}).get("uniprot", {})
    print("Keys in result:", list(result.keys()))
    
    alignments = result.get("rcsb_uniprot_alignments")
    if alignments is None:
        print("rcsb_uniprot_alignments is None!")
        print("Full result structure:")
        pprint.pprint(result, depth=3)
    elif isinstance(alignments, list):
        print("Type: list, length:", len(alignments))
        if alignments:
            a0 = alignments[0]
            print("First element type:", type(a0).__name__)
            if isinstance(a0, dict):
                print("First element keys:", list(a0.keys()))
                core = a0.get("core_entity_alignments")
                print("core_entity_alignments type:", type(core).__name__)
                if isinstance(core, list) and core:
                    print("Sample core entity:", json.dumps(core[0], indent=2)[:500])
    else:
        print("Type:", type(alignments).__name__)
        print("Value:", str(alignments)[:300])
