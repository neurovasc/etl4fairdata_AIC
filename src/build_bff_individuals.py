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
    id = 'ivid-'
    try:
        id = 'ivid-'+row['N°ADN IRT 1']
        if str(row['N°ADN IRT 1']).lower == 'aucun':
            logger.debug("Issue with individual's id: 'aucun'")
    except TypeError:
        logger.debug("Issue with individual's id: empty field, 'aucun', non str value?")

    return str(id)

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
    def build_bmi(row):
        ''' Build dictionary structure for beacon BFF individuals bmi
        to welcome values for an individual and return the filled in structure.
        For individuals with no BMI information, nothing is added in the measures (no
        empty structure, just no structure). But this does not happen when working with the
        rows/individuals that have a genotyping (so in data-intermediate of the snakemake
        pipeline). The try excepts are not useful in that case, and are just an additional
        safeguard.
        '''
        #
        bmi = row['imc']
        observationdate = row['valide le /conformite signature cst biocoll']
        dob = row['date de naissance']
        individualId = get_individualId(row)
        bmistructure = {}
        # bmi
        try:
            bmi = round(float(str(bmi).replace(',', '.')), 2)
            if bmi != bmi: # check if bmi is nan
                logger.debug("There is no BMI for individual "+individualId)
                return bmistructure
            bmistructure = {
                "assayCode": {
                    "id": "LOINC:35925-4",
                    "label": "BMI"
                },
                "measurementValue": {
                    "unit": {
                        "id": "NCIT:C49671",
                        "label": "Kilogram per Square Meter"
                    },
                    "value": bmi
                }
            }
        except ValueError:
            bmi = float('nan')
            logger.debug("There is no BMI for individual "+individualId)
            return bmistructure
        # dob
        try :
            dob = datetime.datetime.strptime(dob, "%d/%m/%Y").date()
            if dob != dob: # check if dob is nan
                logger.debug("There is no DOB for individual "+individualId)
        except ValueError:
            dob = float('nan')
            logger.debug("There is no DOB for individual "+individualId)
        except TypeError:
            logger.debug("The DOB value is 'nan' for individual "+individualId)
        # observation date
        try:
            observationdate = observationdate.split(" ")[0]
            observationdate = datetime.datetime.strptime(observationdate, "%d/%m/%Y").date()
            bmistructure["date"] = observationdate.strftime("%Y-%m-%d")
        except AttributeError:
            logger.debug("Date of observation for BMI is absent for individual "+individualId)
        # moment of observation and stucture build
        try:
            age = observationdate - dob
            years = age.days // 365
            months = (age.days % 365) // 30
            age_of_observation = f"P{years}Y{months}M"
            bmistructure["observationMoment"] = {"age": {"iso8601duration": age_of_observation}}
        except:
            logger.debug("Impossible to compute moment of observation for individual "+individualId)
            
        return bmistructure

    bmi = build_bmi(row)
    measures = []
    if bmi != {}:
        measures.append(bmi)
            
    return measures

def get_enthnicity(row):
    ''' Get ethnicity. In the ICAN cohort as registered in GAIA
    the ethnicity/origin is registed in 'Origine Ethnique' column
    and 'Précision autre origine ethnique'. Only Causasian and 
    Asians are processed, too complicated to choose ontology terms
    between race / ethnicity / geographic origin ontologies. 
    '''
    #
    ethnicity = {}
    #
    origin1 = str(row['Origine ethnique']).lower()
    origin2 = str(row['Précision autre origine ethnique']).lower()
    #
    ethnicitydict = bbm.ethnicitydict
    #
    for o in [origin1, origin2]:
        try:
            e = ethnicitydict[o]
            ethnicity = {"id": e[0], "label": e[1]}
            return ethnicity
        except KeyError:
            if o != 'nan':
                logger.debug("There is no ethnicity corresponding to "+o+" in metadata ethnicitydict.")
            return ethnicity

def get_geographicOrigin(row):
    ''' Geographic origin of individual
    Depends on the meaning: where do they live? where we they born? which hospital are they treated in?
    I put everything in France.
    '''
    return {"id" : "NCIT:C16592", "label" : "France"}

def get_pedigree(row):
    '''In the ICAN study, with the GAIA extraction, 
    there are no related individuals'''
    return {}

def get_exposures(row):
    ''' Exposures such as smoking, alcohol
    NB: In beacon exposures is an array
    https://github.com/ga4gh-beacon/beacon-v2/blob/main/models/json/beacon-v2-default-model/common/exposure.json
    No examples provided as of 09-2024
    Some studies put the number of smoked packets per year in the measures array
    DOI: 10.3233/SHTI240636 
    '''
    exposures = []
    def build_smokingstatus(row):
        ''' Jamais fumé, abstinence >3, abstinence <3, tabagisme actif
        '''
        cigsinapaq = 20 # nombre de cigarettes dans un paquet? idk
        consomation = row['consommation de tabac O/N']
        duration = row['duree en annees']
        smokingstatus = {}
        if consomation == 'Tabagisme actif':
            paqperday = float(str(row['nb de paquets/jour']).replace(',','.'))
            cigsperday = cigsinapaq * paqperday
            smokingstatus = {'exposureCode': {'id' : 'NCIT:C19796', 'label' : 'Smoking Status'},
                             'exposureDescription' : {'id' : 'NCIT:C67147', 'label' : 'Current Smoker'}, 
                                'unit' : {'label': 'cigarettes per day', 'id' : 'EFO:0006525'}, 
                                'value' : cigsperday, 
                                'duration' : duration}
        return smokingstatus
    
    smokingstatus = build_smokingstatus(row)
    if smokingstatus != {}:
        exposures.append(smokingstatus)
    return exposures

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