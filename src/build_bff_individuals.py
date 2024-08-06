# author bodrug-a, chatgpt, github Copilot
# 2024-08-06
#
import argparse
import json
import pandas as pd
import build_bff_metadata4beacon as bbm

# Arguments 
agp = argparse.ArgumentParser(description='Build BFF individuals file')
agp.add_argument('-p', '--phenotypes', type=str,
                help='Path to the file containing phenotypes (one individual per line)',
                required=True)
agp.add_argument('-o', '--output', type=str,
                help='Path to the output directory',
                default='data-deliverable/bff/')
args = agp.parse_args()
#
def get_individualId(row):
    # TODO: Implement logic to get individual ID from row
    return 0

def get_sex(row):
    # TODO: Implement logic to get sex from row
    return 0

def get_measures(row):
    # TODO: Implement logic to get measures from row
    return 0

def get_enthnicity(row):
    # TODO: Implement logic to get ethnicity from row
    return 0

def get_geographicOrigin(row):
    # TODO: Implement logic to get geographic origin from row
    return 0

def get_pedigree(row):
    # TODO: Implement logic to get pedigree from row
    return 0

def get_exposures(row):
    # TODO: Implement logic to get exposures from row
    return 0

def get_diseases(row):
    # TODO: Implement logic to get diseases from row
    return 0

def get_interventionsOrProcedures(row):
    # TODO: Implement logic to get interventions or procedures from row
    return 0

def get_treatments(row):
    # TODO: Implement logic to get treatments from row
    return 0
#
def write_individuals_bff(phenotypefile, output):
    '''
    '''
    individuals = []
    individualsinfo = bbm.individuals # there is nothing here
    df = pd.read_csv(phenotypefile)
    # id, sex, measures, ethnicity, geographicOrigin, pedigree
    # exposures, diseases, interventionsOrProcedures, treatments
    for index, row in df.iterrows():
        id = get_individualId(row)
        sex = get_sex(row)
        measures = get_measures(row)
        ethnicity = get_enthnicity(row)
        geographicOrigin = get_geographicOrigin(row)
        pedigree = get_pedigree(row)
        exposures = get_exposures(row)
        diseases = get_diseases(row)
        interventionsOrProcedures = get_interventionsOrProcedures(row)
        treatments = get_treatments(row)
        individual = bbm.Individual(id, 
                                   sex, 
                                   measures, 
                                   ethnicity, 
                                   geographicOrigin, 
                                   pedigree,
                                   exposures,
                                   diseases, 
                                   interventionsOrProcedures, 
                                   treatments)
        individuals.append(individual.to_dict())
    #
    with open(output+'/individuals.json', 'w') as f:
        json.dump(individuals, f, indent=4)
    return 0
#
# Main
if __name__ == "__main__":
    # Write bff
    write_individuals_bff(args.phenotypes, args.output) 