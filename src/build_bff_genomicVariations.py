# # author bodrug-a, chatgpt, github Copilot
# 2024-08-06
#
import argparse
import json
import pandas as pd
import build_bff_metadata4beacon as bbm

# https://github.com/EGA-archive/beacon2-ri-tools-v2
# sudo docker exec -it ri-tools python genomicVariations_vcf.py 
# with the vcf.gz to convert in files/vcf/files_to_read/
# > sudo docker exec ri-tools-mongo mongoexport --jsonArray --uri 
# "mongodb://root:example@127.0.0.1:27017/beacon?authSource=admin" 
# --collection genomicVariations > ~/bin/beacon2-ri-api/deploy/data/ICAN-exome-634/genomicVariations.json


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

