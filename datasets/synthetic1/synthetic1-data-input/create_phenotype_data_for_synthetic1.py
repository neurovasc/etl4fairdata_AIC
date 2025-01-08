import subprocess 
import argparse
import pandas as pd

# arguments
agp = argparse.ArgumentParser(description='Build BFF for synthetic1')
agp.add_argument('-v', '--vcf', type=str,
                help='Path to the VCF file',
                default="exon_subset_CCDG_14151_B01_GRM_WGS_2020-08-05_allchrs.filtered.shapeit2-duohmm-phased.vcf.gz")
agp.add_argument('-c', '--clinical', 
                help='Path to the clinical data',
                default="../test-data-input/fake-GAIA-extraction-clean.csv")
args = agp.parse_args()

cmd = subprocess.run(["bcftools", "query", "-l", args.vcf], capture_output=True, text=True, check=True)
samples = cmd.stdout.splitlines()
first_line = pd.read_csv(args.clinical, nrows=1)
df = pd.concat([first_line] * len(samples), ignore_index=True)
df['N°ADN IRT 1'] = samples
df['PATIENT'] = ['ivid-' + item for item in samples]

# save
df.to_csv('synthetic1_clinical-data.csv', index=False) 