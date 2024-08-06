# author bodrug-a, chatgpt, github Copilot
# 2024-08-05
#
import argparse
import json
import build_bff_metadata4beacon as bbm
from build_bff_cohorts import get_ids

# Arguments 
agp = argparse.ArgumentParser(description='Build BFF datasets file')
agp.add_argument('-p', '--phenotypes', type=str,
                help='Path to the file containing phenotypes (one individual per line)',
                required=True)
agp.add_argument('-o', '--output', type=str,
                help='Path to the output directory',
                default='data-deliverable/bff/')
args = agp.parse_args()

#
def write_datasets_bff(phenotypefile, output):
    '''
    '''
    datasetsinfo = bbm.datasets # dictionary of info about cohorts
    ids = get_ids(phenotypefile)

    dataset = bbm.Dataset(datasetsinfo['createDateTime'], 
                        datasetsinfo['description'], 
                        datasetsinfo['id'], 
                        datasetsinfo['name'], 
                        ids)
    datasets = [dataset.to_dict()] # there is only one dataset in my data
    with open(output+'/datasets.json', 'w') as f:
        json.dump(datasets, f, indent=4)
#
# Main
if __name__ == "__main__":

    # Write bff
    write_datasets_bff(args.phenotypes, args.output) 