# author bodrug-a, chatgpt, github Copilot
# 2024-08-06
#
import argparse
import json
import pandas as pd
import build_bff_metadata4beacon as bbm
import datetime
import logging

#https://docs.genomebeacons.org/schemas-md/individuals_defaultSchema/

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

def get_biosampleId(row):
    return 0

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
    return []

def get_exposures(row):
    ''' Exposures such as smoking, alcohol
    NB: In beacon exposures is an array
    https://github.com/ga4gh-beacon/beacon-v2/blob/main/models/json/beacon-v2-default-model/common/exposure.json
    No examples provided as of 09-2024
    Some studies put the number of smoked packets per year in the measures array
    DOI: 10.3233/SHTI240636 
    '''
    exposures = []
    ##################
    # smoking status
    ##################
    def build_smokingstatus(row):
        ''' Jamais fumé, abstinence >3, abstinence <3, tabagisme actif
        '''
        consomation = row['consommation de tabac O/N']
        smokingstatus = handle_consumption(consomation, row)
        return smokingstatus
    def handle_consumption(consomation, row):
        ''' Building blocks of smoking status dictionnary according to consumption info in GAIA extraction
            Sevrage plus de 3 and et moins de 3 ans sont tous les deux --> Former Smoker
        '''
        consomation = str(consomation)
        cigsinapaq = 20 # nombre de cigarettes dans un paquet? idk
        duration = row['duree en annees']
        paqperday = float(str(row['nb de paquets/jour']).replace(',','.'))
        cigsperday = cigsinapaq * paqperday
        
        action = {
                'Tabagisme actif': {'exposureCode': {'id' : 'NCIT:C19796', 'label' : 'Smoking Status'},
                    'exposureDescription' : {'id' : 'NCIT:C67147', 'label' : 'Current Smoker'}, 
                    'unit' : {'label': 'cigarettes per day', 'id' : 'EFO:0006525'}, 
                    'value' : str(cigsperday), 
                    'duration' : str(duration)}, 

                'Jamais fumé' : {'exposureCode': {'id' : 'NCIT:C19796', 'label' : 'Smoking Status'},
                    'exposureDescription' : {'id' : 'NCIT:C65108', 'label' : 'Never Smoked'}},

                'nan' : {'exposureCode': {'id' : 'NCIT:C19796', 'label' : 'Smoking Status'},
                    'exposureDescription' : {'id' : 'NCIT:C163971', 'label' : 'Smoking Status Not Documented'}},

                'Sevrage > 3ans' : {'exposureCode': {'id' : 'NCIT:C19796', 'label' : 'Smoking Status'},
                    'exposureDescription' : {'id' : 'NCIT:C67148', 'label' : 'Former Smoker'}, 
                    'unit' : {'label': 'cigarettes per day', 'id' : 'EFO:0006525'}, 
                    'value' : str(cigsperday), 
                    'duration' : str(duration)},

                'Sevrage < 3 ans' : {'exposureCode': {'id' : 'NCIT:C19796', 'label' : 'Smoking Status'},
                    'exposureDescription' : {'id' : 'NCIT:C67148', 'label' : 'Former Smoker'}, 
                    'unit' : {'label': 'cigarettes per day', 'id' : 'EFO:0006525'}, 
                    'value' : str(cigsperday), 
                    'duration' : str(duration)}
                }
        return action[consomation]
    
    ##################
    # alcohol consumption
    ##################
    def build_alcoholconsumption(row):
        '''
        https://evsexplore.semantics.cancer.gov/evsexplore/hierarchy/ncit/C16273
        Difficult to find an ontology and described alcohol consumption per week...
        and draws a limit at 150g per week... 
        TODO
        '''
        alcoholconsumption = {}
        perweekalcohol = row['quantite alcool par semaine (g)']

        return alcoholconsumption
    #############
    # EXPOSURES
    #############
    smokingstatus = build_smokingstatus(row)
    alcoholconsumption = build_alcoholconsumption(row)
    if smokingstatus != {}:
        exposures.append(smokingstatus)
    if alcoholconsumption != {}:
        exposures.append(alcoholconsumption)
    return exposures

def get_diseases(row):
    ''' like {"diseases": [{"diseaseCode": {"id": "ICD10:D70", "label": "agranulocytosis"}}]
    '''
    diseases = []
    #
    def get_aneurysmtatus(row):
        ''' If patient is in the ICAN cohort, then it surely has a history of
        aneurysm, but not necessarily a stroke episode. However in the GAIA extraction
        there are individuals without aneurysms, but those were never sequenced in the 
        context of the ICAN cohort.
        '''
        individualid = get_individualId(row)
        # There is an aneurysm
        status = {'diseaseCode' : {'id' : 'NCIT:C27208', 'label' : 'Brain Aneurysm'}}
        presence = row['statut phenotypique']
        if presence == 'pas d\'anévrisme':
            status['phenotype'] = {'id' : 'NCIT:C41133', 'label' : 'Healed'}
        if presence == 'certain':
            status['phenotype'] = {'id' : 'NCIT:C25626', 'label' : 'Present'}
        if presence == 'incertain':
            status['phenotype'] = {'id' : 'NCIT:C47944', 'label' : 'Uncertain'}
        # Is it a sporadic case or familal history?
        sporadic = row['cas sporadique']
        familial1 = row['ATCD familial d\'AIC (1er degré)']
        familial2 = row['ATCD familial d\'AIC (2ème degré ou plus)']
        if familial1 == 'Oui certain' or familial2 == 'Oui certain':
            status['context'] = {'id' : 'SNOMEDCT:275104002', 'label' : 'Family history of stroke'}
        elif sporadic == 'Oui':
            status['context'] = {'id' : 'SNOMEDCT:75741005', 'label' : 'Sporadic'}
        # How was it discovered? -> This is also an info that can end up in Procedures, if there was
        # a screening
        '''
        discoverycircumstances = row['circonstances de decouverte']
        actions = {
            'Fortuite' : {'id' : 'SNOMED:', 'label' :  'Description'}, 
            'Dépistage familial' : {'id' : 'SNOMED:', 'label' :  'Description'},
            'Rupture AIC' : {'id' : 'SNOMED:', 'label' :  'Description'}, 
            'Compressif ou ischémique' :  {'id' : 'SNOMED:', 'label' :  'Description'}
                   }
        try:
            status['discovery'] =  actions[discoverycircumstances]
        except:
            logger.warning("No discovery context known for individual "+individualid)
        '''
        return status
    def get_strokestatus(row):
        individualid = get_individualId(row)
        status = {}
        # Was there a rupture? Often is discovery is 'Rupture d'AIC' but can also be something else
        # so we consider that there was a rupture if there is a rupture date
        rupturedate = row['date de rupture']
        try:
            rupturedate = datetime.datetime.strptime(rupturedate, "%d/%m/%Y").date()
            status = {'diseaseCode' : {'id' : 'HP:0001297', 'label' : 'Stroke'}, 
                      'date' : str(rupturedate)}
        except:
            logger.info('No rupture for individual '+individualid)
        #
        return status
    def get_htastatus(hta):
        status = {}
        if 'HTA' in hta:
            status = {'diseaseCode' : {'id' : 'HP:0000822', 'label' : 'Hypertension'}}
        if 'HTA gravidique' in hta :
            status = {'diseaseCode' : {'id' : 'HP:0008071', 'label' : 'Maternal hypertension'}}
        return status
    def get_diabetesstatus(diab):
        status = {}
        if 'DNID' in diab: # type 2 insuline resistant
            status = {'diseaseCode' : {'id' : 'HP:0005978', 'label' : 'Type II diabetes mellitus'}}
        if 'DID' in diab: # type 1 insulino dépendant
            status = {'diseaseCode' : {'id' : 'HP:0100651', 'label' : 'Type I diabetes mellitus'}}
        if 'gravidique' in diab:
            status = {'diseaseCode' : {'id' : 'HP:0009800', 'label' : 'Maternal diabetes'}}
        return status
    def get_dyslipidemiastatus(dislip):
        status = {}
        # TODO
        return status
    def get_asthmastatus(dislip):
        status = {}
        # TODO
        return status
    def get_allergystatus(dislip):
        status = {}
        # TODO
        return status
    def get_eczemastatus(dislip):
        status = {}
        # TODO
        return status
    #
    hta = str(row['hypertension arterielle'])
    diabetes = str(row['diabete'])
    htastatus = get_htastatus(hta)
    diabetesstatus = get_diabetesstatus(diabetes)
    aneurysmstatus = get_aneurysmtatus(row)
    strokestatus = get_strokestatus(row)
    #
    diseases.append(aneurysmstatus)
    if strokestatus != {}:
        diseases.append(strokestatus)
    if htastatus != {}:
        diseases.append(htastatus)
    if diabetesstatus != {}:
        diseases.append(diabetesstatus)
    #
    return diseases

def get_interventionsOrProcedures(row):
    # TODO: Implement logic to get interventions or procedures from row
    # Need to loop for approatiate ontology
    return []

def get_treatments(row):
    # TODO: Implement logic to get treatments from row
    # Need to look for appropriate ontology, no idea how to encode statin treatment for example
    return []
#
def write_individuals_bff(phenotypefile, output):
    '''
    '''
    individuals = []
    keeptrack = [] # some individuals have several lines (usually Exome/GWAS, different sequencing centers)
    # I keep track of processed individuals in this array 
    individualsinfo = bbm.individuals # there is nothing here
    df = pd.read_csv(phenotypefile)
    # id, sex, measures, ethnicity, geographicOrigin, pedigree
    # exposures, diseases, interventionsOrProcedures, treatments
    for index, row in df.iterrows():
        id = get_individualId(row)
        if id not in keeptrack:
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
            keeptrack.append(id)
            #print(individual.to_dict())
    #
    with open(output+'/individuals.json', 'w') as f:
        json.dump(individuals, f, indent=4)
    return 0
#
# Main
if __name__ == "__main__":
    # Write bff
    write_individuals_bff(args.phenotypes, args.output) 