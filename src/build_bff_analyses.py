# author bodrug-a, chatgpt, github Copilot
# 2024-08-01

import argparse

# Arguments 
agp = argparse.ArgumentParser(description='Build BFF analyses file')
agp.add_argument('-s', '--samples', type=str, \
                 help='Path to the file containing samples', required=True)

args = parser.parse_args()

# Use args.sample_file in your code to access the file containing samples