# authors: bodrug-a, chatgpt, github Copilot
# 25/07/2024
# This script takes in phenotypic data for AIC patients
# and corresponding genotypes from the VCF file.
# The genotypes are in 0/1 format after a bcftools query. 
# The multiallelic sites were split into biallelic sites.
# The positions without alternative alleles were removed.
# This script creates a new df for each variation, and concatenates it
# to the phenotypic data csv/df. Then it computes the allele frequencies
# according to different features.
# The output of this script is a VCF aggregate file with allele frequencies
# and no individuals genotypes.

import os
import subprocess
import argparse
import pandas as pd
import utilitary as ut

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
    info = ['whole', 'male', 'female']
    # Write the header
    with open(vcfout, 'w') as f:
        f.write(ut.write_headervcf(info, args.sequences))
    # Write the variant positions
    with open(args.genotypes, 'r') as genotypes_file:
        next(genotypes_file) # Skip the header
        for line in genotypes_file:
            chromosome, position, reference, alternative, genotypes = line.strip().split('\t')
            # Build a new df with the genotypes for the variant
            clinicalgenotyped = build_phenotypedf_for_variant(clinical, genotypes)
            ##################################################
            # Compute allele frequencies                     #
            # Compute allele counts                          #
            # Build corresponding VCF header lines           #
            ##################################################
            #
            # Whole population
            AF, AC = ut.compute_AFAC_inpop(clinicalgenotyped) # allele frequency in population
            # Frequencies by sex: male, female, other
            # TODO: add frequencies, write function in utilitary.py
            # Frequencies by case type: familial certain, familial uncertain, sporadic, other
            # TODO: add frequencies, write function in utilitary.py
            # Frequencies by onset: early, late, other
            # Frequencies by number of stroked: multiple, single, duo, other
            # Frequencies by bmi: underweight, normal, overweight, obese1, obese2, other

            info_field_for_variant = f'AF_whole={AF};AC_whole={AC}'
            variant_line = f'{chromosome}\t{position}\t.\t{reference}\t{alternative}\t.\t.\t{info_field_for_variant}'
            with open(vcfout, 'a') as f:
                f.write(variant_line + "\n")
            
    # Compress the file using bgzip
    try:
        subprocess.run(['bgzip', '-f', vcfout], check=True)
        print(f"File compressed successfully to {args.outfile}")
    except subprocess.CalledProcessError as e:
        print(f"Error compressing file: {e}")
        # Clean up temporary file if compression fails
        if os.path.exists(args.outfile):
            os.remove(args.outfile)

