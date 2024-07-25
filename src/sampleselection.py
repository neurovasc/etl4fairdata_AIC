# authors: bodrug-a, chatgpt, github Copilot
# 25/07/2024
# This script takes in phenotypic and clinical data
# from patients of ICA and relatives. It will look for
# patients that have a genomic sample as well.
# The output will be a list of biosamples that have both
# phenotypic and genomic data and that are ICA or related patients.

import argparse
import pandas as pd
import csv

# Argument parser
agp = argparse.ArgumentParser(description='Convert VCF to RDF-ttl')
agp.add_argument('-s', '--samples', type=str, help='Path to the list of samples in VCF', required=True )
agp.add_argument('-c', '--clinical', type=str, help='Path to the phenotype file csv', required=True)
agp.add_argument('-o', '--outsamples', type=str, help='Path to the output file (list of samples)', default='stdout')
agp.add_argument('-f', '--outfile', type=str, help='Path to the output file (filtered clinical csv)', default='stdout')

args = agp.parse_args()

# Load the samples
samples = pd.read_csv(args.samples, header=None)
clinical = pd.read_csv(args.clinical, sep=',')
# Filter out rows with no data or 'aucun' in "N°ADN IRT 1"
clinical = clinical[~clinical["N°ADN IRT 1"].isin(["", "AUCUN", "aucun", "NaN"])].dropna(subset=["N°ADN IRT 1"])
# Extract rows from clinical df where column "N°ADN IRT 1" or column "N°ADN IRT 2" contains a value present in the samples df
filtered_clinical = clinical[clinical["N°ADN IRT 1"].isin(samples[0]) | clinical["N°ADN IRT 2"].isin(samples[0])]
# Write filtered_clinical to output file
filtered_clinical.to_csv(args.outfile, sep=',', quoting=csv.QUOTE_ALL, index=False)
# Create a list of elements found in column "N°ADN IRT 1" and "N°ADN IRT 2"
aic_samples = []
# Merged samples: CD21425_CD23878 CD21678_CD21237 CD26637_CD22832
# These samples should be eliminated as there was an issue during sampleing (info source: Raphael)
merged_samples = ["CD21425", "CD23878", "CD21678", "CD21237", "CD26637", "CD22832"]
for index, row in filtered_clinical.iterrows():
    if row['N°ADN IRT 1'] not in merged_samples and row['N°ADN IRT 1'] not in aic_samples:
        aic_samples.append(row["N°ADN IRT 1"])

# Write aic_samples to outfile parameter path
with open(args.outsamples, 'w') as f: f.write('\n'.join(aic_samples))

