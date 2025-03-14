import requests, sys
 
server = "https://rest.ensembl.org"

# Nuclear genome variation
ext = "/vep/human/hgvs/NC_000006.12:g.152219031G>T?" # chromosome 6
 
r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
decoded = r.json()
#print(repr(decoded))

# Mitochondrial genome variation
ext = "/vep/human/id/rs2853493?"
ext = "/vep/human/hgvs/NC_012920.1:m.11467A>G?"

r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

if not r.ok:
  r.raise_for_status()
  sys.exit()

decoded = r.json()
print(repr(decoded))