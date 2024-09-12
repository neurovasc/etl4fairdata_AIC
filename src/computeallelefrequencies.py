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

import os
import subprocess
import argparse
import fireducks.pandas as pd
import utilitary as ut
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

# Argument parser
agp = argparse.ArgumentParser(description='Compute allele frequencies')
agp.add_argument('-g', '--genotypes', type=str, help='Path to the query genotype file',\
                  required=True)
agp.add_argument('-s', '--sequences', type=str, help='Contig name',\
                  required=True)
agp.add_argument('-c', '--clinical', type=str, help='Path to the phenotype file csv',\
                  required=True)
agp.add_argument('-o', '--outfile', type=str,\
                  help='Path to the output file (should end with .vcf.gz as it will be a bgzipeed vcf file)',\
                  required=True)
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
        chromosome, position, reference, alternative, genotypes = line.strip().split('\t')
        # Build a new df with the genotypes for the variant
        clinicalgenotyped = build_phenotypedf_for_variant(clinical, genotypes) # this is a pandas df
        ##################################################
        # Compute allele frequencies                     #
        # Compute allele counts                          #
        # Build corresponding VCF header lines           #
        ##################################################
        #
        # Whole population
        AF, AC = ut.compute_AFAC_inpop(clinicalgenotyped) # allele frequency in population
        # Frequencies by sex: male, female, other
        AFm, AFf, ACm, ACf = ut.compute_AFAC_bysex(clinicalgenotyped)

        # Frequencies by case type: familial, sporadic, uncertain
        AFfam, AFspo, AFunc, ACfam, ACspo, ACunc = ut.compute_AFAC_bytype(clinicalgenotyped)
        # Frequencies by onset: early, late, other
        AFearly, AFlate, AFother, ACearly, AClate, ACother = ut.compute_AFAC_byonset(clinicalgenotyped)
        # Frequencies by number of stroked: multiple, single, duo, other
        # Frequencies by bmi: underweight, normal, overweight, obese1, obese2, other

        # overall data
        info_field_for_variant = f'AF_whole={AF};AC_whole={AC};'
        # by sex data
        info_field_for_variant += f'AF_female={AFf};AC_female={ACf};'
        info_field_for_variant += f'AF_male={AFm};AC_male={ACm};'
        # by type (sporaic, familail, uncertain)
        info_field_for_variant += f'AF_familial={AFfam};AC_familial={ACfam};'
        info_field_for_variant += f'AF_sporadic={AFspo};AC_sporadic={ACspo};'
        info_field_for_variant += f'AF_uncertain={AFunc};AC_uncertain={ACunc};'
        # by age of onset 
        info_field_for_variant += f'AF_earlyonset={AFearly};AC_earlyonset={ACearly};'
        info_field_for_variant += f'AF_lateonset={AFlate};AC_lateonset={AClate};'
        #info_field_for_variant += f'AF_otheronset={AFother};AC_otheronset={ACother};'
        #
        variant_line = f'{chromosome}\t{position}\t.\t{reference}\t{alternative}\t.\t.\t{info_field_for_variant}\n'
        variant_lines.append(variant_line)
        #
        #print(variant_line)
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
    info = ['whole', 'male', 'female', 'familialcase', 'uncertaincase', 
            'sporadiccase', 'earlyonset', 'lateonset']
    # Write the header
    #with open(vcfout, 'w') as f:
    #    f.write(ut.write_headervcf(info, args.sequences))
    # Write the variant positions
    with open(args.genotypes, 'r') as genotypes_file, open(vcfout, 'w') as outvcf_file:
        # write header in output file (vcf.gz) 
        outvcf_file.write(ut.write_headervcf(info, args.sequences))
        # skip header in input file (bcftools query)
        next(genotypes_file) 
        # process each line in a separate thread
        with ProcessPoolExecutor(max_workers=args.threads) as executor:
            futures = []
            batchcount = 0
            while True:
                # Read the next batch of 100 lines
                lines = [genotypes_file.readline() for _ in range(args.batchsize)]
                # Stop if we reached the end of the file
                if not any(line.strip() for line in lines):
                    break

                # Submit each batch of lines for processing in parallel
                futures.append(executor.submit(process_batchlines, lines))
                batchcount += args.batchsize
                print("linecount: ", batchcount, end='\r')
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
                print(f"Progress: {completed_lines}/{batchcount} lines completed\r")

                
            
    # Compress the file using bgzip
    try:
        subprocess.run(['bgzip', '-f', vcfout], check=True)
        print(f"File compressed successfully to {args.outfile}")
    except subprocess.CalledProcessError as e:
        print(f"Error compressing file: {e}")
        # Clean up temporary file if compression fails
        if os.path.exists(args.outfile):
            os.remove(args.outfile)

