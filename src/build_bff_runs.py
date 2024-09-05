# author bodrug-a, chatgpt, github Copilot
# 2024-08-01
#
import argparse
import json
import build_bff_metadata4beacon as bbm 

# Arguments 
agp = argparse.ArgumentParser(description='Build BFF runs file')
agp.add_argument('-p', '--phenotypes', type=str,
                help='Path to the file containing phenotypes (one individual per line)',
                required=True)
agp.add_argument('-o', '--output', type=str,
                help='Path to the output directory',
                default='data-deliverable/bff/')
args = agp.parse_args()
#
def get_runDate(biosample):
    ''' In phenotype file: column "date (gene)"
    '''
    return '2017-01-01'
def get_libraryStrategy(biosample):
    ''' In phenotype file: column "technique utilisee (gene)"
    '''
    # exome, Genetitan PMRA, Agilent SureSelect, unknown
    return 'unknown'
#
def write_runs_bff(samples, output):
    '''
    '''
    count = 0
    runs = []
    for biosample in samples:
        count += 1 
        runDate = get_runDate(biosample)
        libraryStrategy = get_libraryStrategy(biosample)
        run = bbm.Run(bbm.runs['libraryLayout'], 
                      bbm.runs['librarySelection'], 
                      bbm.runs['librarySource'], 
                      libraryStrategy, 
                      bbm.runs['platform'], 
                      bbm.runs['platformModel'],
                      biosample, 
                      'runId-' + bbm.analyses['runbasename']+'-'+str(count), 
                      'ivid-' + biosample, 
                      runDate)
        runs.append(run.to_dict())
    with open(output+'/runs.json', 'w') as f:
        json.dump(runs, f, indent=4)

#
# Main
if __name__ == "__main__":
    # Load samples
    slist = bbm.load_samples(args.phenotypes)
    # Write bff
    write_runs_bff(slist, args.output)      