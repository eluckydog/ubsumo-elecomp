"""Discover correct RCSB GraphQL schema."""
import urllib.request, json

url = 'https://data.rcsb.org/graphql'

# Introspection query
query = """
{
    __type(name: "CoreUniprot") {
        name
        fields {
            name
            type {
                name
                kind
            }
        }
    }
}
"""

payload = json.dumps({'query': query}).encode()
headers = {'Content-Type': 'application/json'}
req = urllib.request.Request(url, data=payload, headers=headers)
resp = urllib.request.urlopen(req, timeout=15)
data = json.loads(resp.read())

print('CoreUniprot type:')
for field in data.get('data', {}).get('__type', {}).get('fields', []):
    ftype = field.get('type', {})
    print(f'  {field["name"]}: {ftype.get("name")} ({ftype.get("kind")})')

print()

# Also check what arguments the uniprot query takes
intro = """
{
    __schema {
        queryType {
            fields {
                name
                args {
                    name
                    type {
                        name
                        kind
                    }
                }
            }
        }
    }
}
"""
payload2 = json.dumps({'query': intro}).encode()
req2 = urllib.request.Request(url, data=payload2, headers=headers)
resp2 = urllib.request.urlopen(req2, timeout=15)
data2 = json.loads(resp2.read())

for field in data2.get('data', {}).get('__schema', {}).get('queryType', {}).get('fields', []):
    if field['name'] == 'uniprot':
        print('uniprot query args:')
        for arg in field.get('args', []):
            atype = arg.get('type', {})
            print(f'  {arg["name"]}: {atype.get("name")} ({atype.get("kind")})')
        break
