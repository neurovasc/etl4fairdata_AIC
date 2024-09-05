# author bodrug-a, chatgpt, github Copilot
# 2024-08-01
#
import argparse
import json
import build_bff_metadata4beacon as bbm 

# Arguments 
agp = argparse.ArgumentParser(description='Build BFF cohorts file')
agp.add_argument('-p', '--phenotypes', type=str,
                help='Path to the file containing phenotypes (one individual per line)',
                required=True)
agp.add_argument('-o', '--output', type=str,
                help='Path to the output directory',
                default='data-deliverable/bff/')
args = agp.parse_args()
#
def get_cohortSize(phenotypefile):
    '''
    '''
    return 628
def get_collectionEvents(phenotypefile):
    '''
    '''
    return []
def get_ids(phenotypefile):
    '''
    '''
    
    biosampleids = []
    individualids = []
    slist = bbm.load_samples(phenotypefile)
    for id in slist:
        biosampleids.append(id)
        individualids.append('ivid-' + id)
    ids = {'biosampleIds' : biosampleids,
           'individualIds' : individualids}
    return ids
#
def write_cohorts_bff(phenotypefile, output):
    '''
    '''
    cohortsinfo = bbm.cohorts # dictionary of info about cohorts
    cohortSize = get_cohortSize(phenotypefile)
    collectionEvents = get_collectionEvents(phenotypefile)
    ids = get_ids(phenotypefile)

    cohort = bbm.Cohort(cohortsinfo['cohortDesign'], 
                         cohortsinfo['cohortId'], 
                         cohortsinfo['cohortType'], 
                         cohortsinfo['cohortName'], 
                         cohortsinfo['id'], 
                         cohortsinfo['inclusionCriteria'], 
                         cohortSize, 
                         collectionEvents, 
                         ids)
    cohorts = [cohort.to_dict()] # there is only one cohort in my data
    with open(output+'/cohorts.json', 'w') as f:
        json.dump(cohorts, f, indent=4)
#
# Main
if __name__ == "__main__":

    # Write bff
    write_cohorts_bff(args.phenotypes, args.output)   