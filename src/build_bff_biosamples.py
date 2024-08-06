# author bodrug-a, chatgpt, github Copilot
# 2024-08-06
#
import argparse
import json
import pandas as pd
import build_bff_metadata4beacon as bbm

# Arguments 
agp = argparse.ArgumentParser(description='Build BFF biosamples file')
agp.add_argument('-p', '--phenotypes', type=str,
                help='Path to the file containing phenotypes (one individual per line)',
                required=True)
agp.add_argument('-o', '--output', type=str,
                help='Path to the output directory',
                default='data-deliverable/bff/')
args = agp.parse_args()

# Info: for individuals and biosamples, the same base id is used from the phentoype file
# it is the DNA sample id (so technically the biosampleId) 
# It is therefore assumed that there is only one biosample per individual
# This should be the case after processing the phenotype file extracted from GAIA medical
# software through the Snakemake pipeline, as there is only one sample per individual 
# in the genotype file (VCF). Should this change, the ivid- indentifiers shoule be changed to
# correspond to the patient number. I did not do this as the patient number is potentially a
# very sentitive piece of data that I do not want to load into a Beacon. 
#
def get_biosampleStatus(row):
    '''
    '''
    biosampleStatus = {"id": "OBI:0000093","label": "patient"}
    #biosampleStatus = {"id": "???","label": "relative"}
    # NB relatives are not sequenced in the ICAN cohort for now
    return biosampleStatus
#
def get_obtentionProcedure(row):
    '''
    '''
    #
    if row['Type de prélèvement'] == 'veineux':
        obtentionProcedure = {"id": "OBI:2000014", "label": "venous blood collection"}
    elif row['Type de prélèvement'] == 'artériel':
        obtentionProcedure = {"id": "OBI:2000013", "label": "arterial blood collection"}
    else:
        obtentionProcedure = {"id": "OBI:0000655", "label": "blood specimen"}

    # for ICAN, all specimens are blood
    return obtentionProcedure
#
def get_sampleOriginType(row):
    '''
    '''
    #
    if row['Type de prélèvement'] == 'veineux':
        sampleOriginType = {"id": "UBERON:0013756", "label": "venous blood"}
    elif row['Type de prélèvement'] == 'artériel':
        sampleOriginType = {"id": "UBERON:0013755", "label": "arterial blood"}
    else:
        sampleOriginType = {"id": "UBERON:0000178", "label": "blood"}
    #
    return sampleOriginType
#
def get_collectionDate(row):
    '''
    '''
    #
    collectionDate = '2024-10-10'
    #
    return collectionDate
#
def get_collectionMoment(row):
    #
    collectionMoment = 'P40Y'
    #
    return collectionMoment
#
def write_biosamples_bff(phenotypefile, output):
    '''
    '''
    biosamples = []
    biosamplesinfo = bbm.biosamples
    df = pd.read_csv(phenotypefile)
    # info, biosampleStatus, obtentionProcedure, sampleOriginType, 
    # collectionDate, collectionMoment, id, individualId
    for index, row in df.iterrows():
        id = 'bsid-'+row['N°ADN IRT 1']
        individualId = 'ivid-'+row['N°ADN IRT 1']
        biosampleStatus = get_biosampleStatus(row)
        obtentionProcedure = get_obtentionProcedure(row)
        sampleOriginType = get_sampleOriginType(row)
        collectionDate = get_collectionDate(row)
        collectionMoment = get_collectionMoment(row)
        biosample = bbm.Biosample(biosamplesinfo['info'], 
                                  biosampleStatus, 
                                  obtentionProcedure, 
                                  sampleOriginType, 
                                  collectionDate, 
                                  collectionMoment, 
                                  id, 
                                  individualId)
        biosamples.append(biosample.to_dict())
    #
    with open(output+'/biosamples.json', 'w') as f:
        json.dump(biosamples, f, indent=4)
    return 0
#
# Main
if __name__ == "__main__":
    # Write bff
    write_biosamples_bff(args.phenotypes, args.output) 