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
import argparse
import logging
import pandas as pd
import concurrent.futures
import xxxutilitary as ut
import xxxrdfgraphinfo as rgi
from rdflib import Graph

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
# Main
if __name__ == "__main__":
    # Check validity of the VCF file
    loggy.info("### Checking validity of file containing the variants: vcf.gz ###")
    okvcf = ut.color_truefalse(vcfvalidity(args.vcf))
    loggy.info(f"VCF.gz file is valid: {okvcf}")

    # Convert vcf.gz or bcf to bcftools query file
    loggy.info("### Converting VCF.gz to bcftools query file ###")
    os.system(f"bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%INFO\n' \
              {args.vcf} > {args.vcf}.query")
    loggy.info(f"New file available: {args.vcf}.query")

    # Parsing the bcftools query file
    loggy.info("### Reading the bcftools query file ###")
    with open(f"{args.vcf}.query", 'r') as file:
        lines = file.readlines()
    # For testing purposes, shorten the number of lines
    loggy.info(f"There are {len(lines)} variants to process.")
    lines = lines[:1000]
    # Parallel processing of the variants
    chunksize = args.chunksize
    loggy.info("### Initiating multi-threading with futures ###")
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        loggy.debug(f"Submitting to concurrent executor with {args.threads} threads.")
        futures = [executor.submit(process_variants, lines[chunk:chunk+chunksize], 
                                   chunk) for chunk in range(0, len(lines), chunksize)]
        loggy.debug(f"Finished building futures. Waiting for results.")
        results = [future.result() for future in concurrent.futures.as_completed(futures)]
    nb_chunks = len(results)
    loggy.info(f"### Finished processing all {nb_chunks} vcf chunks into subgraphs. ###")
    # Merging all the pocessed variants into one turtle file
    loggy.info(f"Merging all temp .ttl files.")
    count = 0
    # The turtle files should not be loaded into memory
    # but rather serialized to the final file
    temporary_merge = args.tempdir + '/temp_merge.ttl'
    # Empty file if it already exists
    if os.path.exists(temporary_merge):
        open(temporary_merge, 'w').close()
    # Get all the turtle files in temp
    ttlfiles = glob.glob(args.tempdir + '/' +str(os.path.basename(args.rdf))+'_*_intermediate.ttl')
    # Merge all the turtle files into one, without loading them into memory all at once
    with open(temporary_merge, 'w') as f:
        for ttl in ttlfiles:
            count += 1
            if count%100 == 0 or count == nb_chunks:
                loggy.debug(f"Merging files: {count}/{nb_chunks}")
            graph = Graph()
            graph.parse(ttl, format='turtle')
            f.write(graph.serialize(format='turtle'))
    # Create graph for final serialization
    megag = create_rdfgraph_namespace() # merged graph
    megag.parse(temporary_merge, format='turtle')
    megag.serialize(args.rdf, format="turtle")
    loggy.info("### Finished merging temp .ttl files. ###")
    loggy.info(f"Checking if valid turtle file: {args.rdf}.")
    okttl = ut.color_truefalse(check_turtle(args.rdf))
    loggy.debug(f"turtle file is valid: {okttl}")

    # Clean up
    if args.threads:
        loggy.info("### Cleaning up temp files ###")
        clean_temp_ttlfiles()
        loggy.info("Done.")

    

