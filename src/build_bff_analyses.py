# author bodrug-a, chatgpt, github Copilot
# 2024-08-01
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
#
def write_analysis_bff(slist, output):
    ''' Write the analyses json in a json file
    '''
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
                            'anid-' + bbm.analyses['idbasename'] +'-'+str(count), 
                            'ivid-' + biosample,
                            'runId-' + bbm.analyses['runbasename']+'-'+str(count)
                    )
        analyses.append(analysis.to_dict())
    with open(output+'/analyses.json', 'w') as f:
        json.dump(analyses, f, indent=4)
#
# Main
if __name__ == "__main__":
    # Load samples
    slist = bbm.load_samples(args.phenotypes)
    # Write bff
    write_analysis_bff(slist, args.output)        
    