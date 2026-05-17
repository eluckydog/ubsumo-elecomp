"""Explore CorePolymerEntity schema."""
import urllib.request, json

url = "https://data.rcsb.org/graphql"
headers = {"Content-Type": "application/json"}

q = """
{
    __type(name: "CorePolymerEntity") {
        fields {
            name
            type { name kind }
        }
    }
}
"""
payload = json.dumps({"query": q})
req = urllib.request.Request(url, data=payload.encode(), headers=headers)
resp = urllib.request.urlopen(req, timeout=15)
data = json.loads(resp.read())

print("CorePolymerEntity fields:")
for f in data["data"]["__type"]["fields"]:
    ft = f["type"]
    print("  {}: {} ({})".format(f["name"], ft.get("name","?"), ft.get("kind","?")))
