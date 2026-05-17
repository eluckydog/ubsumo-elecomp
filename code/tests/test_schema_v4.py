"""Find alignment-related types."""
import urllib.request, json

url = "https://data.rcsb.org/graphql"
headers = {"Content-Type": "application/json"}

# Get all type names
intro = {"query": "{ __schema { types { name kind } } }"}
req = urllib.request.Request(url, data=json.dumps(intro).encode(), headers=headers)
resp = urllib.request.urlopen(req, timeout=15)
data = json.loads(resp.read())

all_types = data["data"]["__schema"]["types"]

# Filter alignment-related types
align_types = [t for t in all_types if "align" in t["name"].lower()]
pdb_types = [t for t in all_types if "pdb" in t["name"].lower() and "x" in t["name"].lower()]

print("Alignment-related types ({})".format(len(align_types)))
for t in sorted(align_types, key=lambda x: x["name"]):
    print("  {}: {}".format(t["kind"], t["name"]))

print("\nPDBx alignment types ({})".format(len(pdb_types)))
for t in sorted(pdb_types, key=lambda x: x["name"]):
    print("  {}: {}".format(t["kind"], t["name"]))
