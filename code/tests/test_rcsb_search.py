"""Debug RCSB search API query format."""
import urllib.request, json

# RCSB search API v2 accepts structured queries
url = "https://search.rcsb.org/rcsbsearch/v2/query"

# Try the simplest possible search
payloads = [
    # Option 1: full text search for UniProt
    {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {
                "value": "P04637"
            }
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": 5}}
    },
    # Option 2: text search on specific attribute
    {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_align.reference_database_accession",
                "operator": "exact_match",
                "value": "P04637"
            }
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": 5}}
    },
    # Option 3: text search with different attribute
    {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entry_info.related_uniprot_accession",
                "operator": "exact_match",
                "value": "P04637"
            }
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": 5}}
    },
]

for i, payload in enumerate(payloads):
    print("Payload %d:" % (i+1))
    print(json.dumps(payload, indent=2)[:250])
    print()
    
    try:
        req = urllib.request.Request(
            url, data=json.dumps(payload).encode(),
            headers={"Content-Type": "application/json"}
        )
        resp = urllib.request.urlopen(req, timeout=15)
        result = json.loads(resp.read())
        
        identifiers = result.get("result_set", [])
        print("  Found: %d entries" % len(identifiers))
        for item in identifiers[:3]:
            print("  - %s (score=%.2f)" % (item.get("identifier",""), item.get("score",0)))
        print()
    except urllib.error.HTTPError as e:
        print("  HTTP %d: %s" % (e.code, e.read()[:200].decode()))
        print()
    except Exception as e:
        print("  Error: %s" % e)
        print()
