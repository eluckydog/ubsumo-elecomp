"""Explore deeper schema."""
import urllib.request, json

url = 'https://data.rcsb.org/graphql'

query = """
{
    __type(name: "RcsbCoreEntityAlignments") {
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

if 'errors' in data:
    print('No RcsbCoreEntityAlignments type, checking Rcsb...')
    # Try with plural
    alt_types = ['RcsbCoreAlignment', 'CoreEntityAlignments', 'RcsbEntityAlignments', 'RcsbAlignment',
                 'RcsbPdbxUniprotAlignment', 'PdbxUniprotAlignment']
    for t in alt_types:
        q = "{{ __type(name: \"{}\") {{ name fields {{ name type {{ name kind }} }} }} }}".format(t)
        p = json.dumps({'query': q}).encode()
        try:
            r2 = urllib.request.Request(url, data=p, headers=headers)
            r2_resp = urllib.request.urlopen(r2, timeout=10)
            d2 = json.loads(r2_resp.read())
            if d2.get('data', {}).get('__type'):
                print(f'{t} exists!')
                for f in d2['data']['__type']['fields']:
                    print(f'  {f["name"]}: {f["type"]["name"]}')
        except:
            pass
else:
    print('Fields:')
    for f in data['data']['__type']['fields']:
        print(f'  {f["name"]}: {f["type"]["name"]}')
