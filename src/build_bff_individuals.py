# author bodrug-a, chatgpt, github Copilot
# 2024-08-06
#
import argparse
import json
import pandas as pd
import build_bff_metadata4beacon as bbm
import datetime
import logging

# Logger 
logger = logging.getLogger(__name__)
logging.basicConfig(filename='example.log', encoding='utf-8', level=logging.DEBUG)

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
    ''' Construct individualId as ivid- and the code of the ADN sample
    This code is also found in the variant calling format file.
    This function should be changed in the case of multiple samples per individual
    and add another id code (for example 'PATIENT' column in the GAIA extraction)
    '''
    id = ''
    try:
        id = 'ivid-'+row['N°ADN IRT 1']
        if str(row['N°ADN IRT 1']).lower == 'aucun':
            logger.debug("Issue with individual's id: 'aucun'")
    except TypeError:
        logger.debug("Issue with individual's id: empty field, 'aucun', non str value?")

    return id

def get_sex(row):
    '''
    '''
    malestr = ['M', 'm', 'Male', 'homme', 'Homme']
    femalestr = ['F', 'f', 'female', 'femme', 'Femme']
    sex = { "id": "", "label": ""}
    if row['sexe'] in malestr:
        sex['id'] = "NCIT:C20197"
        sex['label'] = "male"
    elif row['sexe'] in femalestr:
        sex['id'] = "NCIT:C16576"
        sex['label'] = "female"
    return sex


def get_measures(row):
    ''' BMI, Aneurysm size (1st), Number of aneurysms
    '''
    # inner functions
    def build_bmi(bmi, observationdate, dob):
        ''' Build dictionary structure for beacon BFF individuals bmi
        to welcome values for an individual and return the filled in structure.
        For individuals with no BMI information
        '''
        # bmi
        try:
            bmi = round(float(str(bmi).replace(',', '.')), 2)
            if bmi != bmi: # check if bmi is nan
                logger.debug("There is no BMI for individual (dob: ", dob, " )")
                return False
        except ValueError:
            bmi = float('nan')
            logger.debug("There is no BMI for individual (dob: ", dob, " )")
            return False
        # observation date
        try:
            observationdate = observationdate.split(" ")[0]
            observationdate = datetime.datetime.strptime(observationdate, "%d/%m/%Y").date()
            dob = datetime.datetime.strptime(dob, "%d/%m/%Y").date()
        except AttributeError:
            logger.debug("Date of observation for BMI is absent.")
            observationdate = ''
        # age at observation: ObservationMoment
        print(bmi, observationdate, dob)
        age = observationdate - dob
        years = age.days // 365
        months = (age.days % 365) // 30
        age_of_observation = f"P{years}Y{months}M"
        bmistructure = {
                "assayCode": {
                    "id": "LOINC:35925-4",
                    "label": "BMI"
                },
                "date": observationdate.strftime("%Y-%m-%d"),
                "measurementValue": {
                    "unit": {
                        "id": "NCIT:C49671",
                        "label": "Kilogram per Square Meter"
                    },
                    "value": bmi
                },
                "observationMoment": {
                    "age": {
                        "iso8601duration": age_of_observation
                    }
                }
            }
        return bmistructure

    bmi = build_bmi(row['imc'], row['valide le /conformite signature cst biocoll'], row['date de naissance'])
    measures = [ bmi ]
            
    return measures

def get_enthnicity(row):
    ''' Get ethnicity. In the ICAN cohort as registered in GAIA
    the ethnicity/origin is registed in 'Origine Ethnique' column
    and Précision autre origine ethnique
    '''
    ethnicitydict = bbm.ethnicitydict
    print(ethnicitydict)
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
        print(individual.to_dict())
    #
    with open(output+'/individuals.json', 'w') as f:
        json.dump(individuals, f, indent=4)
    return 0
#
# Main
if __name__ == "__main__":
    # Write bff
    write_individuals_bff(args.phenotypes, args.output) 