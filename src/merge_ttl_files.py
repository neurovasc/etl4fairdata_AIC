
import os
import glob
from rdflib import Graph
from vcfaggregate2rdf_v2 import create_rdfgraph_namespace

count = 0
# The turtle files should not be loaded into memory
# but rather serialized to the final file
temporary_merge = '/home/bodrug-a/Devlopment/etl4fairdata_AIC/temp/temp_merge.ttl'
# Empty file if it already exists
if os.path.exists(temporary_merge):
    open(temporary_merge, 'w').close()
# Get all the turtle files in temp
ttlfiles = glob.glob('/home/bodrug-a/Devlopment/etl4fairdata_AIC/temp/*_*_intermediate.ttl')
print(ttlfiles)

# Merge all the turtle files into one, without loading them into memory all at once
with open(temporary_merge, 'w') as f:
    for ttl in ttlfiles:
        print(ttl)
        count += 1
        if count%100 == 0 or count == len(ttlfiles):
            print(f"Merging files: {count}/{len(ttlfiles)}")
        graph = Graph()
        graph.parse(ttl, format='turtle')
        f.write(graph.serialize(format='turtle'))
exit()
# Create graph for final serialization
megag = create_rdfgraph_namespace() # merged graph
megag.parse(temporary_merge, format='turtle')
megag.serialize('/home/bodrug-a/Devlopment/etl4fairdata_AIC/temp/mergedfile.ttl', format="turtle")
