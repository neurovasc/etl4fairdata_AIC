# author bodrug-a, chatgpt, github Copilot
# 2024-08-05
#
import argparse
import json
import build_bff_metadata4beacon as bbm 

# Arguments 
agp = argparse.ArgumentParser(description='Build BFF analyses file')
agp.add_argument('-p', '--phenotypes', type=str,
                help='Path to the file containing phenotypes (one individual per line)',
                required=True)
agp.add_argument('-o', '--output', type=str,
                help='Path to the output directory',
                default='data-deliverable/bff/')
args = agp.parse_args()