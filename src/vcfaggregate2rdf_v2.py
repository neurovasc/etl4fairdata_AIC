# author: bodrug-a, chatgpt, github Copilot
# 2025-01-20
# This script takes in a VCF file and converts it to turtle format
# It (tried to) follows the semantic data model described here:
# https://gitlab.univ-nantes.fr/ican/variants-kg-schema/-/blob/main/linkML/ica-variants.yaml?ref_type=heads
# Unithread processing is (should be) limited to 100 variants
# Multithreading is used for processing more than 100 variants

import os
import io
import glob
import gzip
import time
import argparse
import logging
import pandas as pd
import concurrent.futures
import xxxutilitary as ut
import xxxrdfgraphinfo as rgi
from rdflib import Graph
from pathlib import Path

######################################################
# arguments
agp = argparse.ArgumentParser(description='Convert VCF to RDF-ttl')
agp.add_argument('-r', '--rdf', type=str, help='Path to the output RDF/ttl file')
agp.add_argument('-v', '--vcf', type=str, help='Path to the VCF.gz file') 
agp.add_argument('-c', '--chunksize', type=int, 
                 help='Number of variants to process in a single thread',
                 default=100)
agp.add_argument('-t', '--threads', type=int, 
                 help='Number of threads to use',
                 default=1)
agp.add_argument('-d', '--tempdir', type=str, 
                 help='Temporary folder for storing subgraphs',
                 default='temp/')
agp.add_argument('-T', '--testing', action='store_true',
                 help='Use only 150 variants for testing purposes')

args = agp.parse_args()

######################################################
# logging
loggy = logging.getLogger(__name__)

logging.basicConfig(
    level=logging.DEBUG,  
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',  
    datefmt='%Y-%m-%d %H:%M:%S', 
)
######################################################
# timer
class Timer:
    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.interval = self.end - self.start
        print(f"Code block executed in {self.interval:.4f} seconds")

######################################################
# Check validity of arguments
def vcfvalidity(file):
    ''' Ensures that the file exists, is a valid VCF, is zipped, and ends with .vcf.gz
    '''
    # exists
    if not os.path.exists(file):
        loggy.error(f"{file} VCF.gz file does not exist")
        return False
    # ends with vcf.gz
    if not args.vcf.endswith('.vcf.gz'):
        loggy.error(f"{file} is does not end with VCF.gz")
        return False
    # is a gz
    try:
        with gzip.open(file, 'r') as f:
            f.read(1)  # Attempt to read the first byte
    except (OSError, gzip.BadGzipFile):
        loggy.error(f"{file} does not seem to be a gz file")
        return False
    # is a vcf.gz
    try:
        loggy.info(f"Checking if {file} is a VCF with bcftools view")
        os.system(f"bcftools view {file} -h > /dev/null")
    except:
        loggy.error("{file} is not a VCF")
        return False
    return True
#
def process_variants(variants, chunk):
    '''
    '''
    '''
    Process the query line from the bcftools query
    For each bcftools query line a small ttl file is created
    With the namespace and the info of the one variant
    '''
    #loggy.debug(f"Processing variants in chunk {chunk} : initiate mini graph.")
    minig = create_rdfgraph_namespace()
    output = args.tempdir + '/' +str(os.path.basename(args.rdf))+'_'+str(chunk)+'_intermediate.ttl'

    df = pd.DataFrame()
    #loggy.debug(f"Processing variants in chunk {chunk} : build pandas df.")
    for line in variants:
        #print(line)
        dfl = pd.read_csv(io.StringIO(line), sep='\t', \
                        names=['chromosome', 'position', 'reference', 'alternate', 'info'])
        df = pd.concat([df, dfl])
    #
    #loggy.debug(f"Processing variants in chunk {chunk} : build mini graph.")     
    minig = rgi.build_rdfgraph(minig, df)
    #
    #loggy.debug(f"Processing variants in chunk {chunk} : writing to output files: {output}.")
    minig.serialize(destination=output, format="turtle")
    #
    #loggy.debug(f"Processing variants in chunk {chunk} : finished.")
    return minig
#
def create_rdfgraph_namespace():
    '''
    '''
    namespace = rgi.namespace
    g = Graph()
    for key, value in namespace.items():
        g.bind(key, value)
    return g
#
def check_turtle(file):
    '''
    Check if the file is a valid turtle file
    '''
    if not os.path.isfile(file):
        loggy.debug(f"File does not exist: {file}")
        return False
    try:
        g = Graph()
        g.parse(file, format='turtle') # this will load everything into memory
        return True
    except Exception as e:
        loggy.debug(f"Error: {e}")
        return False
#
def clean_temp_ttlfiles():
    ''' Deleting intermediate turtle files created during multithreading
    '''
    temp_folder = args.tempdir
    temp_files = glob.glob(os.path.join(temp_folder, '*.ttl'))
    for file in temp_files:
        os.remove(file)
#
def merge_turtle_files_clean(ttl_files, output_path, loggy=None):
    """
    Merge Turtle files as text, keeping only the first file's prefixes.
    Assumes all TTL files are valid and use consistent prefixes.
    """
    prefix_lines = []
    body_lines = []
    prefix_seen = False

    for count, filepath in enumerate(ttl_files, 1):
        if loggy and (count % 100 == 0 or count == len(ttl_files)):
            loggy.debug(f"Merging file {count}/{len(ttl_files)}: {filepath}")

        with open(filepath, 'r') as infile:
            lines = infile.readlines()

        local_prefix = []
        local_body = []

        for line in lines:
            if line.strip().startswith('@prefix'):
                if not prefix_seen:
                    local_prefix.append(line)
            else:
                local_body.append(line)

        if not prefix_seen:
            prefix_lines = local_prefix
            prefix_seen = True

        body_lines.extend(local_body)
    # Ensure there are no duplicates

    # Write the final output
    with open(output_path, 'w') as outfile:
        outfile.writelines(prefix_lines)
        outfile.write('\n')
        outfile.writelines(body_lines)

    if loggy:
        loggy.info(f"Merged {len(ttl_files)} TTL files into {output_path} with single prefix block.")

# Main
if __name__ == "__main__":
    # Check validity of the VCF file
    loggy.info("### Checking validity of file containing the variants: vcf.gz ###")
    okvcf = ut.color_truefalse(vcfvalidity(args.vcf))
    loggy.info(f"VCF.gz file is valid: {okvcf}")

    # Convert vcf.gz or bcf to bcftools query file
    queryfile = Path(f"{args.vcf}.query")
    if not queryfile.exists():
        loggy.info("### Converting VCF.gz to bcftools query file ###")
        start = time.time()
        os.system(f"bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%INFO\n' \
                {args.vcf} > {args.vcf}.query")
        end = time.time()
        loggy.info(f"New file available: {args.vcf}.query (execution time: {end - start:.2f} seconds)")
    else:
        loggy.info(f"Bcftools query file already available: {args.vcf}.query - NOT COMPUTED AGAIN")

    # Parsing the bcftools query file
    loggy.info("### Reading the bcftools query file ###")
    with open(f"{args.vcf}.query", 'r') as file:
        lines = file.readlines()
    # For testing purposes, shorten the number of lines
    loggy.info(f"There are {len(lines)} variants to process.")
    if args.testing:
        loggy.debug("### Processing the 150 variants, for testing purposes ###")
        lines = lines[:150]
    # Parallel processing of the variants
    chunksize = args.chunksize
    loggy.info("### Initiating multi-threading with futures ###")
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        loggy.debug(f"Submitting to concurrent executor with {args.threads} threads.")
        futures = [executor.submit(process_variants, lines[chunk:chunk+chunksize], 
                                   chunk) for chunk in range(0, len(lines), chunksize)]
        loggy.debug(f"Finished building futures. Waiting for results. Number of futures: {len(futures)}")
        results = [future.result() for future in concurrent.futures.as_completed(futures)]
    nb_chunks = len(results)
    loggy.info(f"### Finished processing all {nb_chunks} vcf chunks into subgraphs. ###")
    # Merging all the pocessed variants into one turtle file
    loggy.info(f"Merging all temp .ttl files. Temporary folder: {args.tempdir}")
    #
    ttlfiles = glob.glob(os.path.join(args.tempdir, f"{os.path.basename(args.rdf)}_*_intermediate.ttl"))
    merge_turtle_files_clean(ttlfiles, args.rdf, loggy=None)
    #
    # Ensure no duplicates
    loggy.info(f"Ensuring no duplicates in {args.rdf} (serialization). Be mindful of memory usage.")
    g = Graph()
    g.parse(args.rdf, format="ttl")
    g.serialize(destination=args.rdf, format="turtle")
    #
    loggy.info("### Finished merging temp .ttl files. ###")
    loggy.info(f"Checking if valid turtle file: {args.rdf}.")
    okttl = ut.color_truefalse(check_turtle(args.rdf))
    loggy.debug(f"turtle file is valid: {okttl}")

    # Clean up
    if args.threads:
        loggy.info("### Cleaning up temp files ###")
        clean_temp_ttlfiles()
        loggy.info("Done.")

    

