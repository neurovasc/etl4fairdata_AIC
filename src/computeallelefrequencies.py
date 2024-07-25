# authors: bodrug-a, chatgpt, github Copilot
# 25/07/2024
# This script takes in phenotypic data for AIC patients
# and corresponding genotypes from the VCF file.
# The genotypes are in 0/1 format after a bcftools query. 
# This script create a new df for each variation, and concatenates it
# to the phenotypic data csv/df. Then it computes the allele frequencies
# according to different features.

import argparse
import pandas as pd

# Argument parser
agp = argparse.ArgumentParser(description='Compute allele frequencies')
agp.add_argument('-g', '--genotypes', type=str, help='Path to the query genotype file', required=True)
agp.add_argument('-c', '--clinical', type=str, help='Path to the phenotype file csv', required=True)
agp.add_argument('-o', '--outfile', type=str, help='Path to the output file', default='stdout')

args = agp.parse_args()

# Functions
def build_phenotypedf_for_variant(clinical, chromosome, position, reference, alternative, genotypes):
    # Build a new df with the genotypes for the variant
    # The columns are the samples and the values are the genotypes
    genotypesarray = genotypes.split(',')
    genotypesarray.pop() # The last element of this array is empty because there is a comma at the end of
    # the genotypes string. It was created with bcftools query -f '[%GT,]' (see snakemake file)
    samples = clinical["N°ADN IRT 1"].values
    genotypesdf = pd.DataFrame(data=genotypesarray, index=samples, columns=['genotype'])
    genotypesdf['chromosome'] = chromosome
    genotypesdf['position'] = position
    genotypesdf['reference'] = reference
    genotypesdf['alternative'] = alternative
    return genotypesdf

# Load clinical data
clinical = pd.read_csv(args.clinical, sep=',')

# Read genotype file line by line
# Do not load it with pandas because the genotypes file is
# potentially very large and it's better to process it line by line
with open(args.genotypes, 'r') as genotypes_file:
    for line in genotypes_file:
        chromosome, position, reference, alternative, genotypes = line.strip().split('\t')


