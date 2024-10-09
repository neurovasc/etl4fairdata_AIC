# allele frequencies and allele counts in the genotype file
# according to the formula: AF = count(alt) / (count(ref|alt2) + count(alt))
# and taking into account phenotypic traits. 

# WARNING: For now, these functions do not fully handle missing genotypes
# for multiallelic sites. It needs to be reworked with a VCF having ./. genotypes
# 

# authors: bodrug-a, chatgpt, github Copilot
# 30/07/2024

# Get the frequency and the count of the alleles for
# the variant in the whole AIC population
# i.e. individuals with intracranial aneurysms, patients or relatives
# of the ICAN cohort

import sys
import logging
import fireducks.pandas as pd
import xxxvcfheaderinfo as vhi

# Logger 
logger = logging.getLogger(__name__)
logging.basicConfig(filename='utilitary.log', encoding='utf-8', level=logging.DEBUG)

#################################
# Write header for the VCF file #
#################################
def write_headervcf(info, contigfile, infofile):
    ''' Write the header for the VCF file
    Contig file contains the header portion of the input VCF with config information
    Info file contains the header portion of the input VCF with info fields, among them
    frequencies and functionnal annotations.
    Alter, possible confusion: info contains the info fields to add to the new vcf, according
    to new phenotypic traits, and infofile contains the info fields from the input vcf.
    '''
    # file format
    headervcf = vhi.fileformat # starts with the file format
    # contigs
    contig = open(contigfile, 'r').read()
    headervcf += contig
    headervcf += vhi.reference
    for i in info:
        try:
            headervcf += vhi.info_headerchunk[i]
        except KeyError as e:
            print(f"KeyError: {i} not found in dictionary for INFO field of vcf header.")
            print(f"KeyError: {e} not found in dictionary .")
            sys.exit(1)
    # annotations in the INFO field of the input VCF
    annotation = open(infofile, 'r').read().split("\n")
    tokeep = ['##INFO=<ID=CADD_RAW', '##INFO=<ID=PHRED', '##INFO=<ID=gnomad_AF', '##INFO=<ID=CSQ']
    #
    for a in annotation:
        if any(keep in a for keep in tokeep):
            headervcf += a+"\n"
    # final line in VCF header
    headervcf += vhi.columns # header ends with columns CHR POS ID REF ALT QUAL FILTER INFO
    return headervcf
#
#################################
# Reshape existing INFO field ###
#################################
def reshape_infofield(info):
    ''' Reshape the INFO field to be compatible with the VCF format
    '''
    # Keep only certain fields from the original INFO field of the input vcf file
    # Original annotated file: QCed.VEP.AFctrls.GND.CADD.vcf.gz by R Blanchet
    fields = info.split(';')
    csq = 'CSQ=nan'
    gnomad_AF = 'gnomad_AF=nan'
    gnomad_AF_AFR = 'gnomad_AF_AFR=nan'
    gnomad_AF_EAS = 'gnomad_AF_EAS=nan'
    gnomad_AF_NFE = 'gnomad_AF_NFE=nan'
    cadd_raw = 'CADD_RAW=nan'
    cadd_phred = 'PHRED=nan'

    for field in fields:
        if 'CSQ' in field:
            csq = field
        if 'gnomad_AF=' in field:
            gnomad_AF = field
        if 'gnomad_AF_AFR=' in field:
            gnomad_AF_AFR = field
        if 'gnomad_AF_EAS=' in field:
            gnomad_AF_EAS = field
        if 'gnomad_AF_NFE=' in field:
            gnomad_AF_NFE = field
        if 'CADD_RAW=' in field:
            cadd_raw = field
        if 'PHRED=' in field:
            cadd_phred = field
    reshaped = ";".join([gnomad_AF, gnomad_AF_AFR, gnomad_AF_EAS, gnomad_AF_NFE, cadd_raw, cadd_phred, csq])
    return reshaped
#
##################################
# Compute frequencies and counts #
##################################
def compute_frequencies(df, group):
    ''' df is a clinical dataframe with the genotypes for a given variant
    Here we compute frequencies and counts for one variant according to
    different clinical traits.
    female, male, familial case, sporadic case, early onset, late onset
    Exemple: what is the frequency of the alternative allele within the subset of
    the cohort that has a familial case of aneurysm? How many people does this represent?
    In group (array),there are all the types of frequencies and counts to compute.
    The output of this function is as follows:
    frequencies : array : floats : len N : definition - frequency of a variant ALLELE in a given group
    counts : array : str(integer/integegrs) : len N : definition - count of a variant ALLELE in a given group (/total)
    labels : array of tuples : (str, str) : len N : definition - labels for the frequencies and counts in the INFO field
    '''
    frequencies = ['nan']*len(group)
    counts = ['nan']*len(group)
    infofields = get_infofieldslabels(group)
    #
    for i, g in enumerate(group):
        if 'whole' == g: # whole population in the OC vcf file
            genotypes = ''.join(df['genotypes'].tolist())
            refc, altc, unkc = genotypes.count('0'), genotypes.count('1'), genotypes.count('.')
            frequencies[i] = calcfreq({'0': refc, '1': altc})
            counts[i] = f"{altc}/{refc+altc}"
        if "female" == g: # female population in the OC vcf file
            femalesex_list = ['F', 'female', 'femme']
            female_df = df[df['sexe'].isin(femalesex_list)]
            female_genotypes = ''.join(female_df['genotypes'].tolist())
            female_refc, female_altc, female_unkc = female_genotypes.count('0'), female_genotypes.count('1'), female_genotypes.count('.')
            frequencies[i] = calcfreq({'0': female_refc, '1': female_altc})
            counts[i] = f"{female_altc}/{female_refc+female_altc}"
        if "male" == g: # male population in the OC vcf file
            malesex_list = ['M', 'male', 'homme', 'H']
            male_df = df[df['sexe'].isin(malesex_list)]
            male_genotypes = ''.join(male_df['genotypes'].tolist())
            male_refc, male_altc, male_unkc = male_genotypes.count('0'), male_genotypes.count('1'), male_genotypes.count('.')
            frequencies[i] = calcfreq({'0': male_refc, '1': male_altc})
            counts[i] = f"{male_altc}/{male_refc+male_altc}"
        if "earlyonset" in g: # early onset population in the OC vcf file
            df = add_age_of_onset(df)
            earlyonset_df = df[df['age of onset']<35]
            earlyonset_genotypes = ''.join(earlyonset_df['genotypes'].tolist())
            earlyonset_refc, earlyonset_altc, earlyonset_unkc = earlyonset_genotypes.count('0'), earlyonset_genotypes.count('1'), earlyonset_genotypes.count('.')
            frequencies[i] = calcfreq({'0': earlyonset_refc, '1': earlyonset_altc})
            counts[i] = f"{earlyonset_altc}/{earlyonset_refc+earlyonset_altc}"
        if "obese" in g:
            df['imc'] = pd.to_numeric(df['imc'], errors='coerce')
            obese_df = df[(df['imc'] >= 30) & (df['imc'] < 100)]
            obese_genotypes = ''.join(obese_df['genotypes'].tolist())
            obese_refc, obese_altc, obese_unkc = obese_genotypes.count('0'), obese_genotypes.count('1'), obese_genotypes.count('.')
            frequencies[i] = calcfreq({'0': obese_refc, '1': obese_altc})
            counts[i] = f"{obese_altc}/{obese_refc+obese_altc}"
        if "overweight" in g:
            df['imc'] = pd.to_numeric(df['imc'], errors='coerce')
            overweight_df = df[(df['imc'] >= 25) & (df['imc'] < 30)]
            overweight_genotypes = ''.join(overweight_df['genotypes'].tolist())
            overweight_refc, overweight_altc, overweight_unkc = overweight_genotypes.count('0'), overweight_genotypes.count('1'), overweight_genotypes.count('.')
            frequencies[i] = calcfreq({'0': overweight_refc, '1': overweight_altc})
            counts[i] = f"{overweight_altc}/{overweight_refc+overweight_altc}"
        if "familialcase" in g:
            familial_df = df[(df['ATCD familial d\'AIC (1er degré)'] == 'Oui certain')]
            familial_genotypes = ''.join(familial_df['genotypes'].tolist())
            familial_refc, familial_altc, familial_unkc = familial_genotypes.count('0'), familial_genotypes.count('1'), familial_genotypes.count('.')
            frequencies[i] = calcfreq({'0': familial_refc, '1': familial_altc})
            counts[i] = f"{familial_altc}/{familial_refc+familial_altc}"
        if "sporadiccase" in g:
            sporadic_df = df[(df['cas sporadique'] == 'Oui')]
            sporadic_genotypes = ''.join(sporadic_df['genotypes'].tolist())
            sporadic_refc, sporadic_altc, sporadic_unkc = sporadic_genotypes.count('0'), sporadic_genotypes.count('1'), sporadic_genotypes.count('.')
            frequencies[i] = calcfreq({'0': sporadic_refc, '1': sporadic_altc})
            counts[i] = f"{sporadic_altc}/{sporadic_refc+sporadic_altc}"
        if "discoveryincidental" in g:
            discovery_df = df[(df['circonstances de decouverte'] == 'Fortuite')]
            discovery_genotypes = ''.join(discovery_df['genotypes'].tolist())
            discovery_refc, discovery_altc, discovery_unkc = discovery_genotypes.count('0'), discovery_genotypes.count('1'), discovery_genotypes.count('.')
            frequencies[i] = calcfreq({'0': discovery_refc, '1': discovery_altc})
            counts[i] = f"{discovery_altc}/{discovery_refc+discovery_altc}"
        if "discoveryfamilialscreening" in g:
            discovery_df = df[(df['circonstances de decouverte'] == 'Dépistage familial')]
            discovery_genotypes = ''.join(discovery_df['genotypes'].tolist())
            discovery_refc, discovery_altc, discovery_unkc = discovery_genotypes.count('0'), discovery_genotypes.count('1'), discovery_genotypes.count('.')
            frequencies[i] = calcfreq({'0': discovery_refc, '1': discovery_altc})
            counts[i] = f"{discovery_altc}/{discovery_refc+discovery_altc}"
        if "discoveryruptured" in g:
            discovery_df = df[(df['circonstances de decouverte'] == 'Rupture AIC')]
            discovery_genotypes = ''.join(discovery_df['genotypes'].tolist())
            discovery_refc, discovery_altc, discovery_unkc = discovery_genotypes.count('0'), discovery_genotypes.count('1'), discovery_genotypes.count('.')
            frequencies[i] = calcfreq({'0': discovery_refc, '1': discovery_altc})
            counts[i] = f"{discovery_altc}/{discovery_refc+discovery_altc}"
        if "discoveryischemic" in g:
            discovery_df = df[(df['circonstances de decouverte'] == 'Compressif ou ischémique')]
            discovery_genotypes = ''.join(discovery_df['genotypes'].tolist())
            discovery_refc, discovery_altc, discovery_unkc = discovery_genotypes.count('0'), discovery_genotypes.count('1'), discovery_genotypes.count('.')
            frequencies[i] = calcfreq({'0': discovery_refc, '1': discovery_altc})
            counts[i] = f"{discovery_altc}/{discovery_refc+discovery_altc}"
        if "ruptured" in g:
            ruptured_df = df[((df['AIC 1 Rompu'] == 'Oui')) | ((df['AIC 2 Rompu'] == 'Oui')) | ((df['AIC 3 Rompu'] == 'Oui'))
                             | ((df['AIC 4 Rompu'] == 'Oui')) | ((df['AIC 5 Rompu'] == 'Oui'))]
            ruptured_genotypes = ''.join(ruptured_df['genotypes'].tolist())
            ruptured_refc, ruptured_altc, ruptured_unkc = ruptured_genotypes.count('0'), ruptured_genotypes.count('1'), ruptured_genotypes.count('.')
            frequencies[i] = calcfreq({'0': ruptured_refc, '1': ruptured_altc})
            counts[i] = f"{ruptured_altc}/{ruptured_refc+ruptured_altc}"
        if "multipleica" in g:
            multipleica_df = df[(df['nb d anevrismes'] > 1)]
            multipleica_genotypes = ''.join(multipleica_df['genotypes'].tolist())
            multipleica_refc, multipleica_altc, multipleica_unkc = multipleica_genotypes.count('0'), multipleica_genotypes.count('1'), multipleica_genotypes.count('.')
            frequencies[i] = calcfreq({'0': multipleica_refc, '1': multipleica_altc})
            counts[i] = f"{multipleica_altc}/{multipleica_refc+multipleica_altc}"
        if "aht" in g:
            aht_df = df[(df['hypertension arterielle'] ==  'HTA traitée') | (df['hypertension arterielle'] == 'HTA non traitée')]
            aht_genotypes = ''.join(aht_df['genotypes'].tolist())
            aht_refc, aht_altc, aht_unkc = aht_genotypes.count('0'), aht_genotypes.count('1'), aht_genotypes.count('.')
            frequencies[i] = calcfreq({'0': aht_refc, '1': aht_altc})
            counts[i] = f"{aht_altc}/{aht_refc+aht_altc}"
        if "diabetes" in g:
            diabetes_df = df[(df['diabete'] == 'Oui')]
            diabetes_genotypes = ''.join(diabetes_df['genotypes'].tolist())
            diabetes_refc, diabetes_altc, diabetes_unkc = diabetes_genotypes.count('0'), diabetes_genotypes.count('1'), diabetes_genotypes.count('.')
            frequencies[i] = calcfreq({'0': diabetes_refc, '1': diabetes_altc})
            counts[i] = f"{diabetes_altc}/{diabetes_refc+diabetes_altc}"
        if "dyslipidemia" in g:
            dyslipidemia_df = df[(df['dyslipidemie'] == 'Oui')]
            dyslipidemia_genotypes = ''.join(dyslipidemia_df['genotypes'].tolist())
            dyslipidemia_refc, dyslipidemia_altc, dyslipidemia_unkc = dyslipidemia_genotypes.count('0'), dyslipidemia_genotypes.count('1'), dyslipidemia_genotypes.count('.')
            frequencies[i] = calcfreq({'0': dyslipidemia_refc, '1': dyslipidemia_altc})
            counts[i] = f"{dyslipidemia_altc}/{dyslipidemia_refc+dyslipidemia_altc}"
        #
    return frequencies, counts, infofields
#
def get_infofieldslabels(group):
    infofields = []
    for g in group:
        # the labels
        try:
            # here we split the vcf header info fields corresponding to the group
            # of interest, and extract only the labels for the frequencies and counts
            # to appear in the line of the variant in the VCF file
            labels = vhi.info_headerchunk[g].split('\n')
            l1 = labels[0].split('ID=')[1].split(',')[0]
            l2 = labels[1].split('ID=')[1].split(',')[0]
            infofields.append((l1, l2))   
        except KeyError as e:
            print(f"KeyError: {g} not found in dictionary for INFO field of vcf header.")
            print(f"KeyError: {e}")
            sys.exit(1)
    return infofields
#
def calcfreq(count_dict):
    ''' Calculate the allele frequency from a dictionary with counts of alleles
    This function assumes that the alleles are coded as 0 and 1
    And that there are no multialelic sites.
    Multiallelic sites are supposed to be split into biallelic sites with
    bcftools norm -m -any. This will create a new line for each allele
    and for genotypes that are like 1/2 it will replace them with 1/0 
    This will increase the amount of reference like alleles however it will not
    change the allele frequency of the allele that is coded as 1.
    Missing alleles coded as '.' or '*' are not taken into account for the
    allele frequency calculation.
    '''
    try:
        AF = count_dict['1'] / (count_dict['0'] + count_dict['1'])
    except ZeroDivisionError:
        AF = 0
    except KeyError:
        AF = 0
    return round(AF,6)
#
def add_age_of_onset(df):
    '''
    '''
    df['date de naissance'] = pd.to_datetime(df['date de naissance'], dayfirst=True)
    df['date du 1er diagnostic'] = pd.to_datetime(df['date du 1er diagnostic'], errors='coerce', dayfirst=True)
    df['age of onset'] = df.apply(lambda row: (row['date du 1er diagnostic'].year - row['date de naissance'].year 
                                           - ((row['date du 1er diagnostic'].month, row['date du 1er diagnostic'].day) 
                                              < (row['date de naissance'].month, row['date de naissance'].day))) 
                              if pd.notna(row['date du 1er diagnostic']) else pd.NA, axis=1)
    return df
def check_values_equal(val1, val2):
    '''
    '''
    if val1 != val2:
        print(f"Values are not equal: {val1} != {val2}")
        sys.exit(1)
    print("Values are equal.")
#
def compute_AFAC_bytype(df):
    ''' Compute allele frequency by type: familial, sporadic, uncertain

    '''
    # Familial cases are non sporadic and have a confirmed first OR second degree family
    # member that had a confirmed/certain aneurysm
    
    familial_df = df[((df['ATCD familial d\'AIC (1er degré)'] == 'Oui certain') | \
                     (df['ATCD familial d\'AIC (2ème degré ou plus)'] == 'Oui certain')) & \
                     (df['cas sporadique'] == 'Non')]
    # Sporadic cases are marked as such and no not have unconfirmed familal cases of first AND
    # second degree link.
    sporadic_df = df[(df['cas sporadique'] ==  'Oui') & \
                     (df['ATCD familial d\'AIC (1er degré)'] != 'Oui non confirmé') & \
                     (df['ATCD familial d\'AIC (2ème degré ou plus)'] != 'Oui non confirmé')]
    # Uncertain cases are cases where neither the familial or the sporadicity type was
    # confirmed. Could be a case marked sporadic but with an unconfirmed familial background zB
    # for cases where sporadic is set to 'Non' but nothing in ATCD fields
    indices_to_drop = pd.concat([familial_df, sporadic_df]).index
    uncertain_df = df.drop(indices_to_drop)
    #
    nbfam = len(familial_df) # ok, checked visually 11092024
    nbspo = len(sporadic_df) # ok, check visually 11092024
    nbunc = len(uncertain_df) # the remaining lines in df afer familial and sporadic extraction
    #print(uncertain_df[['cas sporadique', 'ATCD familial d\'AIC (1er degré)', 'ATCD familial d\'AIC (2ème degré ou plus)']].to_string())
    #logger.debug(check_values_equal(nbfam+nbspo+nbunc, len(df)))
    AFf_among, ACf = compute_AFAC_inpop(familial_df)
    AFs_among, ACs = compute_AFAC_inpop(sporadic_df)
    AFu_among, ACu = compute_AFAC_inpop(uncertain_df)
    AFf = get_round_af(ACf, df)
    AFs = get_round_af(ACs, df)
    AFu = get_round_af(ACu, df) 

    return AFf, AFs, AFu, ACf, ACs, ACu
#
def compute_AFAC_byonset(df):
    ''' Onseat can be early, late, or unknown?
    age of onset -> date of birth - L'age du premier diag: age
    Si moins de 35 ans, early onset.
    See: 
    '''
    df['date de naissance'] = pd.to_datetime(df['date de naissance'], dayfirst=True)
    df['date du 1er diagnostic'] = pd.to_datetime(df['date du 1er diagnostic'], errors='coerce', dayfirst=True)
    df['age of onset'] = df.apply(lambda row: (row['date du 1er diagnostic'].year - row['date de naissance'].year 
                                           - ((row['date du 1er diagnostic'].month, row['date du 1er diagnostic'].day) 
                                              < (row['date de naissance'].month, row['date de naissance'].day))) 
                              if pd.notna(row['date du 1er diagnostic']) else pd.NA, axis=1)
    earlyonsetthreshold = 35 # to be confirmed
    early_df = df[df['age of onset']<earlyonsetthreshold]
    late_df = df[df['age of onset']>=earlyonsetthreshold]
    indices_to_drop = pd.concat([early_df, late_df]).index
    unknown_df = df.drop(indices_to_drop)
    #
    AFe_among, ACe = compute_AFAC_inpop(early_df)
    AFl_among, ACl = compute_AFAC_inpop(late_df)
    AFo_among, ACo = compute_AFAC_inpop(unknown_df)
    #
    AFe = get_round_af(ACe, df)
    AFl = get_round_af(ACl, df)
    AFo = get_round_af(ACo, df)
    #
    return AFe, AFl, AFo, ACe, ACl, ACo


if __name__ == "__main__":
    print("This script is intended to be used as a module.")
