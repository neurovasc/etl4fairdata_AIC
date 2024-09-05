# # author bodrug-a, chatgpt, github Copilot
# 2024-08-06
#
import argparse
import json
import pandas as pd
import build_bff_metadata4beacon as bbm

# Arguments 
agp = argparse.ArgumentParser(description='Build BFF genomicVariantions file')
agp.add_argument('-g', '--genotypes', type=str,
                help='Path to the file containing genotypes (VCF)',
                required=True)
agp.add_argument('-o', '--output', type=str,
                help='Path to the output directory',
                default='data-deliverable/bff/')
args = agp.parse_args()
#
# TO DO
#
def write_genomicVariations_bff(vcf, output):
    # TO DO 
    #
    genomicVariations = []
    with open(output+'/genomicVariations.json', 'w') as f:
        json.dump(genomicVariations, f, indent=4)
    return 0    
#
# Main
if __name__ == "__main__":
    # Write bff
    write_genomicVariations_bff(args.genotypes, args.output)

