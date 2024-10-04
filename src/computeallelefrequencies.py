# authors: bodrug-a, chatgpt, github Copilot
# 25/07/2024
# This script takes in phenotypic data for AIC patients
# and corresponding genotypes from the VCF file.
# The genotypes are in 0/1 format after a bcftools query. 
# The multiallelic sites were split into biallelic sites.
# The positions without alternative alleles were removed.
# Normally these positions without an alt are eliminated with previous
# filterings in the pipeline: bcftools view -c1:alt1 
# This script creates a new df for each variation, and concatenates it
# to the phenotypic data csv/df. Then it computes the allele frequencies
# according to different features.
# The output of this script is a VCF aggregate file with allele frequencies
# and no individuals genotypes.
# The output of this script is a VCF that retains all the annotations of the input VCF
# The input VCF should be annotated with VEP, snpEff, contain CADD scores
# The input VCF usually used for this script was annotated by Raphael Blanchet
# for example: QCed.VEP.AFctrls.GND.CADD.vcf.gz

# debug pyarrow: export LD_LIBRARY_PATH=/home/bodrug-a/miniconda3/envs/etl4fair/lib/:$LD_LIBRARY_PATH

# example of command line:
# python3 src/computeallelefrequencies.py \
# --genotype data-intermediate/aicdataset-querygenotype.tsv \
# --sequences data-intermediate/aicdataset-contigs.txt \
# --clinical data-intermediate/aicdataset-extraction_GAIA_ICAN_26-09-2023.reordered.csv \
# --outfile data-intermediate/test.vcf.gz --threads 8 --batchsize 10

import os
import subprocess
import argparse
import fireducks.pandas as pd
import xxxutilitary as ut
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import logging

# Logging   
# Logger 
logger = logging.getLogger(__name__)
logging.basicConfig(filename='computeallelefrequency.log', encoding='utf-8', level=logging.DEBUG)

# Argument parser
agp = argparse.ArgumentParser(description='Compute allele frequencies')
agp.add_argument('-g', '--genotypes', type=str, help='Path to the query genotype file',\
                required=True,\
                default='data-intermediate/aicdataset-querygenotype.tsv')
agp.add_argument('-s', '--sequences', type=str, help='Contig name',\
                required=True,\
                default='data-intermediate/aicdataset-contigs.txt')
agp.add_argument('-i', '--info', type=str, help='Path to the file holding info field of the OG input vcf',\
                required=False,\
                default='data-intermediate/aicdataset-infofields.txt')
agp.add_argument('-c', '--clinical', type=str, help='Path to the phenotype file csv',\
                required=True,\
                default='data-intermediate/aicdataset-extraction_GAIA_ICAN_26-09-2023.reordered.csv')
agp.add_argument('-o', '--outfile', type=str,\
                help='Path to the output file (should end with .vcf.gz as it will be a bgzipeed vcf file)',\
                required=True,\
                default='data-intermediate/aicdataset-aggregate-testing.vcf.gz')
agp.add_argument('-t', '--threads', type=int,\
                help='Number of threads to use (one thread per batch of lines)',\
                default=8)
agp.add_argument('-b', '--batchsize', type=int,\
                help='Number of lines (variants) to put into a batch',\
                default=5)

args = agp.parse_args()

# Functions
def build_phenotypedf_for_variant(clinical, genotypes):
    # Build a new df with the genotypes for the variant
    # The columns are the samples and the values are the genotypes
    genotypesarray = genotypes.split(',')
    genotypesarray.pop() # The last element of this array is empty because there is a comma at the end of
    # the line in the genotypes file obtained with bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT,\t]\n'
    # sanity check
    if len(clinical) != len(genotypesarray):
        raise ValueError("Length of clinical dataframe does not match length of genotypes array")
    clinical['genotypes'] = genotypesarray
    return clinical
#
def process_batchlines(lines):
    ''' Process a genotype line from the VCF
    One line per thread. For parallelization purposes.
    '''
    variant_lines = []
    for line in lines:
        chromosome, position, reference, alternative, genotypes, info = line.strip().split('\t')
        # Build a new df with the genotypes for the variant
        clinicalgenotyped = build_phenotypedf_for_variant(clinical, genotypes) # this is a pandas df
        ##################################################
        # Compute allele frequencies                     #
        # Compute allele counts                          #
        # Build corresponding VCF header lines           #
        ##################################################
        #
        # Getting frequencies and info fields
        
        frequencies, counts, infofields = ut.compute_frequencies(clinicalgenotyped)
        #logging.debug(f'frequencies: {frequencies}')
        #logging.debug(f'counts: {counts}')
        #logging.debug(f'infofields: {infofields}')
        # Setting up the info field for the variant
        info_field_for_variant = ''
        #
        for i in range(0, len(infofields)):
            af_field = infofields[i][0] # frequency label
            ac_field = infofields[i][1] # count label
            af = frequencies[i] # frequency value
            ac = counts[i] # count value 
            info_field_for_variant += f'{af_field}={af};'
            info_field_for_variant += f'{ac_field}={ac};'
        #
        info = ut.reshape_infofield(info)
        info_field_for_variant += info
        #
        variant_line = f'{chromosome}\t{position}\t.\t{reference}\t{alternative}\t.\t.\t{info_field_for_variant}\n'
        variant_lines.append(variant_line)
        #
        #logging.debug(variant_line)
        #
    return ''.join(variant_lines)
    
# Main
if __name__ == "__main__":

    # Load clinical data
    clinical = pd.read_csv(args.clinical, sep=',')

    # Uncompressed vcf 
    vcfout = args.outfile.replace('.gz', '')

    # Read genotype file line by line
    # Do not load it with pandas because the genotypes file is
    # potentially very large and it's better to process it line by line
    
    # Allele frequencies and allele counts should be specified in the vcf header
    # Choosing which frequencies to add to VCF:
    group = ['whole', 'male', 'female']

    # Write the variant positions
    with open(args.genotypes, 'r') as genotypes_file, open(vcfout, 'w') as outvcf_file:
        # write header in output file (vcf.gz) 
        outvcf_file.write(ut.write_headervcf(group, args.sequences, args.info))
        # skip header in input file (bcftools query)
        next(genotypes_file) 
        # process each line in a separate thread
        with ProcessPoolExecutor(max_workers=args.threads) as executor:
            futures = []
            batchcount = 0
            linecount = 0
            while True:
                # Read the next batch of X lines
                lines = [genotypes_file.readline() for _ in range(args.batchsize)]
                # Stop if we reached the end of the file
                #print(f"lines: {lines}")
                #print(f"lines[0]: {lines[0]}")
                # End of file was reached, there is nothing more p
                if any(line == '' for line in lines):
                    non_empty_lines = [line for line in lines if line != '']
                    if len(non_empty_lines) > 0:
                        print("\nEnd of file reached, there are a couple of lines to submit to processing.\n")
                        futures.append(executor.submit(process_batchlines, non_empty_lines))
                    else:
                        print("\nEnd of file reached.\n")
                    break

                # Submit each batch of lines for processing in parallel
                futures.append(executor.submit(process_batchlines, lines))
                batchcount += 1
                linecount += args.batchsize # in batchsize, we have the number of lines
                print(f"linecount: {linecount}\r", end='')
            # 
            print("\nChecking completed batched: \n")
            completed_lines = 0
            # Loop through completed futures and track progress
            
            for future in as_completed(futures):
                # The results variable is a string
                # It is a concatenation of all the vcf like
                # variant lines passed to
                # the function that processes the query lines in a batch
                results = future.result()
                # keep count
                completed_lines += args.batchsize
                # write to output file
                outvcf_file.write(results)


                # Print progress
                print(f"Progress: {completed_lines}/{linecount} lines completed\r", end='')
            
    # Compress the file using bgzip
    try:
        subprocess.run(['bgzip', '-f', vcfout], check=True)
        print(f"File compressed successfully to {args.outfile}")
    except subprocess.CalledProcessError as e:
        print(f"Error compressing file: {e}")
        # Clean up temporary file if compression fails
        if os.path.exists(args.outfile):
            os.remove(args.outfile)

