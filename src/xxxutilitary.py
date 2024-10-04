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

def write_headervcf(info, contigfile, infofile):
    ''' Write the header for the VCF file
    Contig file contains the header portion of the input VCF with config information
    Info file contains the header portion of the input VCF with info fields, among them
    frequencies and functionnal annotations.
    '''
    headervcf = vhi.fileformat # starts with the file format
    contig = open(contigfile, 'r').read()
    headervcf += contig
    headervcf += vhi.reference
    for i in info:
        try:
            headervcf += vhi.info_headerchunk[i]
        except KeyError as e:
            print(f"KeyError: {e} not found in dictionary.")
            sys.exit(1)
    headervcf += vhi.columns # header ends with columns CHR POS ID REF ALT QUAL FILTER INFO
    return headervcf
#
def write_headervcf_info(type):
    ''' Write the header for the VCF file: INFO field
    '''
    try:
        headervcf = vhi.info_headerchunk[type]
    except KeyError as e:
        print(f"KeyError: {e} not found in dictionary.")
        sys.exit(1) 
    return headervcf
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
    return AF
#
def compute_AFAC_inpop(df):
    ''' Input phenotype data with a column 'genotypes' that contains the genotypes
    for a given variant. The genotypes are in 0/1 or 0|1 format. Biallelic sites split.
    '''
    genotypes = df['genotypes'].tolist()
    count_dict = {}
    for genotype in genotypes:
        for char in genotype:
            if char != '/' and char != '|': 
                if char in count_dict:
                    count_dict[char] += 1
                else:
                    count_dict[char] = 1
    AF = round(calcfreq(count_dict), 8)
    try:
        AC = count_dict['1']
    except:
        #print(f"No alternative allele for selection.")
        AC = 0
    return AF, AC
#
def check_values_equal(val1, val2):
    '''
    '''
    if val1 != val2:
        print(f"Values are not equal: {val1} != {val2}")
        sys.exit(1)
    print("Values are equal.")
#
def get_round_af(count, df):
    '''
    '''
    AF = round(count / (len(df)*2), 8)
    return AF
#
def compute_AFAC_bysex(df):
    ''' Compute allele frequencies by sex. 
    The output is the allele frequency of the alternative allele occuring in females/males
    The AFf_among keeps track of the allele frequency of the alternative allele 
    among females i.e. answers the questions: how frequent is the allele in the female
    population. It's not quite what we want for the aggregate file, but we keep it in mind
    because it may be an useful metric down the line.
    '''
    femalesex_list = ['F', 'female', 'femme']
    malesex_list = ['M', 'male', 'homme', 'H']
    female_df = df[df['sexe'].isin(femalesex_list)]
    male_df = df[df['sexe'].isin(malesex_list)]
    AFf_among, ACf = compute_AFAC_inpop(female_df)
    AFm_among, ACm = compute_AFAC_inpop(male_df)
    AFf = get_round_af(ACf, df)
    AFm = get_round_af(ACm, df)
    #print('females: ',female_df['genotypes'].tolist())
    #print('males: ',male_df['genotypes'].tolist())
    #print('AFm, AFf, ACm, ACf',AFm, AFf, ACm, ACf)
    #print(len(df), len(female_df), len(male_df))
    
    return AFm, AFf, ACm, ACf
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

def reshape_infofield(info):
    ''' Reshape the INFO field to be compatible with the VCF format
    '''
    # Keep only certain fields from the original INFO field of the input vcf file
    # Original annotated file: QCed.VEP.AFctrls.GND.CADD.vcf.gz by R Blanchet
    return info

def compute_frequencies(df):
    ''' df is a clinical dataframe with the genotypes for a given variant
    Here we compute frequencies and counts for one variant according to
    different clinical traits.
    female, male, familial case, sporadic case, early onset, late onset
    Exemple: what is the frequency of the alternative allele within the subset of
    the cohort that has a familial case of aneurysm? How many people does this represent?
    '''
    frequencies = [0,5]
    counts = [15]
    infofields = [('AF_whole', 'AC_whole')]
    return frequencies, counts, infofields

if __name__ == "__main__":
    print("This script is intended to be used as a module.")
