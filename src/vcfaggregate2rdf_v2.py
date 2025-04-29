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
                 default=None)
agp.add_argument('-T', '--testing', action='store_true',
                 help='Use only 150 variants for testing purposes')
agp.add_argument('-L', '--testinglines', type=int,
                 nargs="+",
                 help='Use variants in the interval n:m for testing purposes')
agp.add_argument('-C', '--cleanup', action='store_true',
                 help='Cleanup temporary files after processing')

args = agp.parse_args()

######################################################
# temporary directory

if args.tempdir is None:
    # Make a temp folder based on timestamp and basename
    timestamp = time.strftime('%Y%m%d-%H%M%S')
    vcf_basename = os.path.splitext(os.path.basename(args.vcf))[0]
    temp_folder_name = f"temp/vcfaggregate2rdf_{timestamp}_{vcf_basename}"
    
    # Actually create the directory
    os.makedirs(temp_folder_name, exist_ok=True)

    args.tempdir = temp_folder_name

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
    Merge Turtle (.ttl) files into a single file.
    
    - Only keep @prefix declarations from the first file.
    - Write triples from all files.
    - Process files line-by-line to minimize memory usage.
    - Optionally logs progress.
    """
    prefix_written = False

    with open(output_path, 'w') as outfile:
        for count, filepath in enumerate(ttl_files, 1):
            if loggy and (count % 100 == 0 or count == len(ttl_files)):
                loggy.debug(f"Merging file {count}/{len(ttl_files)}: {filepath}")

            with open(filepath, 'r') as infile:
                for line in infile:
                    stripped_line = line.strip()
                    if stripped_line.startswith('@prefix'):
                        if not prefix_written:
                            outfile.write(line)
                    else:
                        outfile.write(line)

            if not prefix_written:
                # After finishing the first file, write a blank line and mark prefixes as done
                outfile.write('\n')
                prefix_written = True
    if loggy:
        loggy.info(f"Merged {len(ttl_files)} TTL files into {output_path} with a single prefix block.")


def chunk_lines(file_path, chunksize, start=0, end=None):
    """
    Generator yielding chunks of lines from a file without loading the entire file in memory.
    """
    with open(file_path, 'r') as f:
        chunk = []
        for idx, line in enumerate(f):
            if idx < start:
                continue
            if end is not None and idx >= end:
                break
            chunk.append(line)
            if len(chunk) >= chunksize:
                yield chunk, idx - len(chunk) + 1  # return chunk and its starting line number
                chunk = []
        if chunk:
            yield chunk, idx - len(chunk) + 1

# ---------------- Main Script ----------------

if __name__ == "__main__":
    # 1. Check validity of VCF file
    loggy.info("### Checking validity of file containing the variants: vcf.gz ###")
    okvcf = ut.color_truefalse(vcfvalidity(args.vcf))
    loggy.info(f"VCF.gz file is valid: {okvcf}")

    # 2. Convert vcf.gz or bcf to bcftools query file
    queryfile = Path(f"{args.vcf}.query")
    if not queryfile.exists():
        loggy.info("### Converting VCF.gz to bcftools query file ###")
        start = time.time()
        os.system(f"bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%INFO\n' {args.vcf} > {args.vcf}.query")
        end = time.time()
        loggy.info(f"New file available: {args.vcf}.query (execution time: {end - start:.2f} seconds)")
    else:
        loggy.info(f"Bcftools query file already available: {args.vcf}.query - NOT COMPUTED AGAIN")

    # 3. Count total number of lines (lightweight)
    loggy.info("### Counting the number of variants ###")
    total_lines = 0
    with open(f"{args.vcf}.query", 'r') as file:
        for _ in file:
            total_lines += 1
    loggy.info(f"There are {total_lines} variants to process.")

    # 4. Testing mode adjustments
    startat = 0
    endat = total_lines
    if args.testing:
        if args.testinglines is not None:
            startat = args.testinglines[0]
            endat = args.testinglines[1]
            loggy.debug(f"### Processing {startat}:{endat} variants, for testing purposes ###")
        else:
            loggy.debug("### Processing 200 variants, for testing purposes ###")
            endat = min(200, total_lines)

    # 5. Parallel processing of variants
    chunksize = args.chunksize
    max_pending_futures = args.threads  # Limit in-flight futures
    loggy.info("### Initiating multi-threading with futures ###")

    futures = set()
    chunk_generator = chunk_lines(f"{args.vcf}.query", chunksize, startat, endat)

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        loggy.debug(f"Throttled submission to concurrent executor with {args.threads} threads, max pending futures = {max_pending_futures}.")
        
        # Pre-fill futures queue
        for _ in range(max_pending_futures):
            try:
                chunk, chunk_start = next(chunk_generator)
                future = executor.submit(process_variants, chunk, chunk_start)
                futures.add(future)
                loggy.debug(f"Submitted new chunk starting at line {chunk_start}.")
            except StopIteration:
                break

        # Process completed futures and keep submitting new chunks
        while futures:
            done, futures = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)
            for fut in done:
                _ = fut.result()  # No accumulation of large results

                try:
                    chunk, chunk_start = next(chunk_generator)
                    future = executor.submit(process_variants, chunk, chunk_start)
                    futures.add(future)
                    loggy.debug(f"Submitted new chunk starting at line {chunk_start}.")
                except StopIteration:
                    pass

    loggy.info(f"### Finished processing all VCF chunks into subgraphs. ###")

    '''
    # 6. Merge temp turtle files
    loggy.info(f"Merging all temp .ttl files. Temporary folder: {args.tempdir}")
    ttlfiles = glob.glob(os.path.join(args.tempdir, f"{os.path.basename(args.rdf)}_*_intermediate.ttl"))
    loggy.debug(f"Number of temporary files to merge: {len(ttlfiles)}")
    merge_turtle_files_clean(ttlfiles, args.rdf, loggy=loggy)

    # 7. Optionally ensure no duplicates (careful, memory heavy)
    loggy.info(f"Ensuring no duplicates in {args.rdf} (serialization). Be mindful of memory usage.")
    g = Graph()
    g.parse(args.rdf, format="ttl")
    g.serialize(destination=args.rdf, format="turtle")

    loggy.info("### Finished merging temp .ttl files. ###")

    # 8. Check validity of final TTL file
    loggy.info(f"Checking if valid turtle file: {args.rdf}.")
    oktl = ut.color_truefalse(check_turtle(args.rdf))
    loggy.debug(f"Turtle file is valid: {oktl}")

    # 9. Cleanup
    if args.cleanup:
        loggy.info("### Cleaning up temp files ###")
        clean_temp_ttlfiles()
        loggy.info("Done.")
    '''
    

