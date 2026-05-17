"""Get RcsbUniprotAlignmentsCoreEntityAlignments schema."""
import urllib.request, json

url = "https://data.rcsb.org/graphql"
headers = {"Content-Type": "application/json"}

types_to_check = [
    "RcsbUniprotAlignmentsCoreEntityAlignments",
    "CoreEntityAlignmentsAlignedRegions",
    "CoreEntityAlignmentsCoreEntityIdentifiers",
    "CoreEntityAlignmentsScores",
    "RcsbPolymerEntityAlign",
    "RcsbPolymerEntityAlignAlignedRegions",
]

for tname in types_to_check:
    q = '{{ __type(name: "{}") {{ name fields {{ name type {{ name kind }} }} }} }}'.format(tname)
    req = urllib.request.Request(url, data=json.dumps({"query": q}).encode(), headers=headers)
    resp = urllib.request.urlopen(req, timeout=10)
    data = json.loads(resp.read())
    
    t = data.get("data", {}).get("__type")
    if t:
        print("{}:".format(t["name"]))
        for f in t.get("fields", []):
            ft = f["type"]
            print("  {}: {} ({})".format(f["name"], ft.get("name"), ft.get("kind")))
        print()
