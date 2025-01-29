import requests, sys
 
server = "https://rest.ensembl.org"
ext = "/vep/human/hgvs/6:g.152219031G>T?"
 
r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
decoded = r.json()
print(repr(decoded))