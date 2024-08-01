# author bodrug-a, chatgpt, github Copilot
# 2024-08-01
# This script is coded in object oriented Python

import argparse
import json
import build_bff_metadata4beacon as bbm 

# Arguments 
agp = argparse.ArgumentParser(description='Build BFF analyses file')
agp.add_argument('-s', '--samples', type=str,
                help='Path to the file containing samples (one sample per line)',
                required=True)
agp.add_argument('-o', '--output', type=str,
                help='Path to the output directory',
                default='data-deliverable/bff/')
args = agp.parse_args()
#
def load_samples(samples):
    ''' Load samples from a file
    '''
    with open(samples, 'r') as f:
        samples = [line.strip() for line in f if line.strip()]
    return samples
#
# Main
if __name__ == "__main__":
    # Load samples
    slist = load_samples(args.samples)
    count = 0
    analyses = []
    for biosample in slist:
        count += 1
        analysis = bbm.Analysis(bbm.analyses['aligner'], 
                            bbm.analyses['analysisDate'],
                            bbm.analyses['pipelineName'],
                            bbm.analyses['pipelineRef'], 
                            bbm.analyses['variantCaller'],
                            'bsid-' + biosample, 
                            'anid-' + bbm.analyses['idbasename'] +'_'+str(count), 
                            'ivid-' + biosample,
                            'runId' + bbm.analyses['runbasename']+'_'+str(count)
                    )
        analyses.append(analysis.to_dict())
    with open(args.output+'/analyses.json', 'w') as f:
        json.dump(analyses, f, indent=4)        
    