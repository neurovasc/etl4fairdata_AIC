import os
import glob
from rdflib import Graph
import xxxrdfgraphinfo as rgi
import argparse

# Argument parser
parser = argparse.ArgumentParser(description="Merge TTL files with shared prefixes and deduplication.")
parser.add_argument('-d', '--temp-folder', type=str, required=True,
                    help='Temporary folder containing *_intermediate.ttl files to merge.')
args = parser.parse_args()

# Function
def create_rdfgraph_namespace():
    '''
    '''
    namespace = rgi.namespace
    g = Graph()
    for key, value in namespace.items():
        g.bind(key, value)
    return g

# Paths
temp_folder = args.temp_folder
temporary_merge = os.path.join(temp_folder, 'temp_merge.ttl')
final_merged_file = os.path.join(temp_folder, 'mergedfile.ttl')
# Empty the temporary merge file if it already exists
open(temporary_merge, 'w').close()

# Get all the intermediate turtle files
ttlfiles = glob.glob(os.path.join(temp_folder, '*_*_intermediate.ttl'))
print(f"Found {len(ttlfiles)} files to merge.")

# Step 1: Extract prefixes from the first file
prefix_lines = []
first_file = ttlfiles[0]
with open(first_file, 'r') as f:
    for line in f:
        if line.strip().startswith('@prefix'):
            prefix_lines.append(line)

print(f"Collected {len(prefix_lines)} prefixes from the first file.")

# Step 2: Merge
with open(temporary_merge, 'w') as outfile:
    # Write all prefixes first
    outfile.writelines(prefix_lines)
    outfile.write('\n')  # Add a newline after prefixes

    # Now process all files
    for idx, ttl in enumerate(ttlfiles, 1):
        with open(ttl, 'r') as infile:
            for line in infile:
                if not line.strip().startswith('@prefix'):
                    outfile.write(line)

        if idx % 1 == 0 or idx == len(ttlfiles):
            print(f"Merged {idx}/{len(ttlfiles)} files.")

print(f"Finished merging all files into {temporary_merge}.")

# Step 3: Deduplicate triples by parsing temp_merge.ttl into a clean graph
print("Deduplicating merged TTL file...")
megag = create_rdfgraph_namespace()  # Create fresh graph with namespaces
print("Parsing temporary merge file...")
megag.parse(temporary_merge, format='turtle')
print("Removing duplicates... (serializing to a new graph)")
megag.serialize(destination=final_merged_file, format='turtle')
print(f"Final merged TTL file written to: {final_merged_file}")
