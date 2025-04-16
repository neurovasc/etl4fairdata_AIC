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
import pandas as pd
import xxxvcfheaderinfo as vhi

# Logger 
loggy = logging.getLogger(__name__)

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
    # VCF HEADER: file format
    headervcf = vhi.fileformat # starts with the file format
    # VCF HEADER: contigs
    contig = open(contigfile, 'r').read()
    headervcf += contig
    # VCF HEADER: reference
    headervcf += vhi.reference
    # VCF HEADER: info field for the new phenotypic traits to add as frequencies and counts
    for i in info:
        try:
            headervcf += vhi.info_headerchunk[i]
        except KeyError as e:
            print(f"KeyError: {i} not found in dictionary for INFO field of vcf header.")
            print(f"KeyError: {e} not found in dictionary .")
            sys.exit(1)
    # VCF HEADER: annotations in the INFO field of the original OG input VCF
    annotation = open(infofile, 'r').read().split("\n")
    tokeep = [] #['##INFO=<ID=CADD_RAW', '##INFO=<ID=PHRED', '##INFO=<ID=gnomad_AF', '##INFO=<ID=CSQ']
    #
    for a in annotation:
        if any(keep in a for keep in tokeep):
            headervcf += a+"\n"
    # VCF HEADER: final line in VCF header
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
    #csq = 'CSQ=nan'
    gnomad_AF = 'gnomad_AF=nan'
    gnomad_AF_AFR = 'gnomad_AF_AFR=nan'
    gnomad_AF_EAS = 'gnomad_AF_EAS=nan'
    gnomad_AF_NFE = 'gnomad_AF_NFE=nan'
    #cadd_raw = 'CADD_RAW=nan'
    #cadd_phred = 'PHRED=nan'

    for field in fields:
        #if 'CSQ' in field:
        #    csq = field
        #if 'gnomad_AF=' in field or 'AF=' in field:
        #    gnomad_AF = field
        if 'AF_AFR=' in field:
            gnomad_AF_AFR = field
        if 'AF_EAS=' in field:
            gnomad_AF_EAS = field
        if 'AF_NFE=' in field:
            gnomad_AF_NFE = field
        #if 'CADD_RAW=' in field:
        #    cadd_raw = field
        #if 'PHRED=' in field:
        #    cadd_phred = field
    #reshaped = ";".join([gnomad_AF, gnomad_AF_AFR, gnomad_AF_EAS, gnomad_AF_NFE, cadd_raw, cadd_phred, csq])
    reshaped = ";".join([])
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
    I insist on ALLELES being counted, not individual genotypes. Info about homozygous or heterozygous is not 
    computed for now.
    '''
    frequencies = ['nan']*len(group)
    counts = ['nan']*len(group)
    infofields = get_infofieldslabels(group)
    #
    for i, g in enumerate(group):
        # OK CHECK 04152025 syntheticican2
        if 'whole' == g: # whole population in the OC vcf file
            af, af_hom, af_het, ac, ac_hom, ac_het = process_genotypes(df['genotypes'])
            frequencies[i] = {'allz': af, 'homoz': af_hom, 'heteroz': af_het} 
            counts[i] = {'allz': ac, 'homoz': ac_hom, 'heteroz': ac_het}
        # OK CHECK 04152025 syntheticican2
        if "female" == g: # female population in the OC vcf file
            # filter based on clinical/phenotypical trait 
            femalesex_list = ['F', 'female', 'femme']
            female_df = df[df['sex'].isin(femalesex_list)]
            # process 
            af, af_hom, af_het, ac, ac_hom, ac_het = process_genotypes(female_df['genotypes'])
            frequencies[i] = {'allz': af, 'homoz': af_hom, 'heteroz': af_het} 
            counts[i] = {'allz': ac, 'homoz': ac_hom, 'heteroz': ac_het}
        # OK CHECK 04152025 syntheticican2
        if "male" == g: # male population in the OC vcf file
            # filter based on clinical/phenotypical trait
            malesex_list = ['M', 'male', 'homme', 'H']
            male_df = df[df['sex'].isin(malesex_list)]
            # process
            af, af_hom, af_het, ac, ac_hom, ac_het = process_genotypes(male_df['genotypes'])
            frequencies[i] = {'allz': af, 'homoz': af_hom, 'heteroz': af_het} 
            counts[i] = {'allz': ac, 'homoz': ac_hom, 'heteroz': ac_het}
        # OK CHECK 04152025 syntheticican2
        if "earlyonset" in g: # early onset population in the OC vcf file
            # filter based on clinical/phenotypical trait
            earlyonset_df = df[(df['firstDiagnosisAge'] < 35)]
            # process
            af, af_hom, af_het, ac, ac_hom, ac_het = process_genotypes(earlyonset_df['genotypes'])
            frequencies[i] = {'allz': af, 'homoz': af_hom, 'heteroz': af_het} 
            counts[i] = {'allz': ac, 'homoz': ac_hom, 'heteroz': ac_het}
        # OK CHECK 04152025 syntheticican2
        if "obese" in g:
            # filter based on clinical/phenotypical trait
            df['bodyMassIndex'] = pd.to_numeric(df['bodyMassIndex'], errors='coerce')
            obese_df = df[(df['bodyMassIndex'] >= 30) & (df['bodyMassIndex'] < 100)]
            # process
            af, af_hom, af_het, ac, ac_hom, ac_het = process_genotypes(obese_df['genotypes'])
            frequencies[i] = {'allz': af, 'homoz': af_hom, 'heteroz': af_het} 
            counts[i] = {'allz': ac, 'homoz': ac_hom, 'heteroz': ac_het}
        # OK CHECK 04152025 syntheticican2
        if "familialcase" in g:
            # filter based on clinical/phenotypical trait
            familial_df = df[(df['familialCase'] == True)]
            # process
            af, af_hom, af_het, ac, ac_hom, ac_het = process_genotypes(familial_df['genotypes'])
            frequencies[i] = {'allz': af, 'homoz': af_hom, 'heteroz': af_het}
            counts[i] = {'allz': ac, 'homoz': ac_hom, 'heteroz': ac_het}
        # OK CHECK 04152025 syntheticican2
        if "aht" in g:
            # filter based on clinical/phenotypical trait
            hta_df = df[(df['medicalHistory-arterialHypertension'] == True)]
            # process
            af, af_hom, af_het, ac, ac_hom, ac_het = process_genotypes(hta_df['genotypes'])
            frequencies[i] = {'allz': af, 'homoz': af_hom, 'heteroz': af_het}
            counts[i] = {'allz': ac, 'homoz': ac_hom, 'heteroz': ac_het}
        # OK CHECK 04152025 syntheticican2
        if "diabetes" in g:
            # filter based on clinical/phenotypical trait
            diabetes_df = df[(df['medicalHistory-diabetes'] == True)]
            # process
            af, af_hom, af_het, ac, ac_hom, ac_het = process_genotypes(diabetes_df['genotypes'])
            frequencies[i] = {'allz': af, 'homoz': af_hom, 'heteroz': af_het}
            counts[i] = {'allz': ac, 'homoz': ac_hom, 'heteroz': ac_het}
        # OK CHECK 04152025 syntheticican2
        if "neversmoked" in g:
            # filter based on clinical/phenotypical trait
            neversmoked_df = df[(df['lifestyle-tobacco_neversmoked'] == True)]
            # process
            af, af_hom, af_het, ac, ac_hom, ac_het = process_genotypes(neversmoked_df['genotypes'])
            frequencies[i] = {'allz': af, 'homoz': af_hom, 'heteroz': af_het}
            counts[i] = {'allz': ac, 'homoz': ac_hom, 'heteroz': ac_het}
        # OK CHECK 04152025 syntheticican2
        if "currentlysmoking" in g:
            # filter based on clinical/phenotypical trait
            currentlysmoking_df = df[(df['lifestyle-tobacco_currentlysmoking'] == True)]
            # process
            af, af_hom, af_het, ac, ac_hom, ac_het = process_genotypes(currentlysmoking_df['genotypes'])
            frequencies[i] = {'allz': af, 'homoz': af_hom, 'heteroz': af_het}
            counts[i] = {'allz': ac, 'homoz': ac_hom, 'heteroz': ac_het}
        # OK CHECK 04152025 syntheticican2
        if "multipleica" in g:
            # filter based on clinical/phenotypical trait
            multipleica_df = df[(df['multipleAneurysms'] == True)]
            # process
            af, af_hom, af_het, ac, ac_hom, ac_het = process_genotypes(multipleica_df['genotypes'])
            frequencies[i] = {'allz': af, 'homoz': af_hom, 'heteroz': af_het}
            counts[i] = {'allz': ac, 'homoz': ac_hom, 'heteroz': ac_het}
        if "treatment" in g:
            # filter based on clinical/phenotypical trait
            treatment_df = df[df['aic-1-treatment'] | df['aic-2-treatment'] | \
                              df['aic-3-treatment'] | df['aic-4-treatment'] | df['aic-5-treatment']] 
            # process
            af, af_hom, af_het, ac, ac_hom, ac_het = process_genotypes(treatment_df['genotypes'])
            frequencies[i] = {'allz': af, 'homoz': af_hom, 'heteroz': af_het}
            counts[i] = {'allz': ac, 'homoz': ac_hom, 'heteroz': ac_het}
        if "treatmentcoils" in g:
            # filter based on clinical/phenotypical trait
            treatmentcoils_df = df[(df['aic-1-treatmentType'] == 'Coils') | (df['aic-2-treatmentType'] == 'Coils') | \
                                    (df['aic-3-treatmentType'] == 'Coils') | (df['aic-4-treatmentType'] == 'Coils') | \
                                    (df['aic-5-treatmentType'] == 'Coils')]
            # process
            af, af_hom, af_het, ac, ac_hom, ac_het = process_genotypes(treatmentcoils_df['genotypes'])
            frequencies[i] = {'allz': af, 'homoz': af_hom, 'heteroz': af_het}
            counts[i] = {'allz': ac, 'homoz': ac_hom, 'heteroz': ac_het}
        #
        verbose = 0
        if verbose:
            print(f"Group: {g}, \
                frequencies: {frequencies[i]}, \
                counts: {counts[i]}")
    return frequencies, counts, infofields
#
def get_infofieldslabels(group):
    infofields = []
    for g in group:
        try:
            # here we split the vcf header info fields corresponding to the group
            # of interest, and extract only the labels for the frequencies and counts
            # to appear in the line of the variant in the VCF file
            labels = vhi.info_headerchunk[g].split('\n')
            l1 = labels[0].split('ID=')[1].split(',')[0] # allele frequency
            l2 = labels[1].split('ID=')[1].split(',')[0] # allele count
            l3 = labels[2].split('ID=')[1].split(',')[0] # allele count in homozygous states
            l4 = labels[3].split('ID=')[1].split(',')[0] # allele count in heterozygous states
            infofields.append((l1, l2, l3, l4))   
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
def process_genotypes(genotypes):
    ''' Processes a genotype string to count/compute
    - the number of alternative alleles 
    - the number of alternative alleles in homozygous state
    - the number of alternative alleles in heterozygous state
    - the frequency of the alternative allele
    - the frequency of the alternative allele in homozygous state
    - the frequency of the alternative allele in heterozygous state
    '''
    #
    genotypeslist = genotypes.tolist()
    #
    for i in range(len(genotypeslist)):
        g = genotypeslist[i]
        newg = g.replace('|', '/') # if phased genotypes
        genotypeslist[i] = newg
    #
    # TOTAL ALTERNATE ALLELE COUNT AND FREQUENCY
    totalAlleleCount = len(genotypeslist)*2 
    totalReferenceAlleleCount = sum(item.count('0') for item in genotypeslist)
    totalAlternateAlleleCount = sum(item.count('1') for item in genotypeslist)
    #
    totalAlternateAlleleFrequency = calcfreq({'0': totalReferenceAlleleCount, '1': totalAlternateAlleleCount}) # needs to be FLOAT, as described in the vcf header
    totalAlternateAlleleCount = f"{totalAlternateAlleleCount}" # needs to be INTEGER, as described in the vcf header
    #
    # ALTERNATE ALLELE COUNT AND FREQUENCY IN HOMOZYGOUS GENOTYPES
    # What is the alternate allele frequency among homozygous individuals? 
    homozygousGenotypeslist = [x for x in genotypeslist if x not in {'0/1', '1/0', './.'}]
    homozygousReferenceAlleleCount = sum(item.count('0') for item in homozygousGenotypeslist)
    homozygousAlternateAlleleCount =  sum(item.count('1') for item in homozygousGenotypeslist)
    homozygousAlternateAlleleFrequency = calcfreq({'0': homozygousReferenceAlleleCount, '1': homozygousAlternateAlleleCount})
    homozygousAlternateAlleleCount = f"{homozygousAlternateAlleleCount}" # allele counts, in homozygous status, meaning we get the number of individuals by /2 
    #
    # ALTERNATE ALLELE COUNT AND FREQUENCY IN HETEROZYGOUS GENOTYPES
    # What is the alternate allele frequency among heterozygous individuals? Always 0.5 of course
    # We only want to know the count anyway.
    heterozygousGenotypeslist = [x for x in genotypeslist if x not in {'0/0', '1/1', './.'}]
    heterozygousReferenceAlleleCount = sum(item.count('0') for item in heterozygousGenotypeslist)
    heterozygousAlternateAlleleCount =  heterozygousReferenceAlleleCount # sum(item.count('1') for item in heterozygousGenotypeslist)
    heterozygousAlternateAlleleFrequency = 0.5 #calcfreq({'0': heterozygousReferenceAlleleCount, '1': heterozygousAlternateAlleleCount})
    heterozygousReferenceAlleleCount = f"{heterozygousReferenceAlleleCount}" # allele counts, in homozygous status, meaning we get the number of individuals by /2
    #
    verbose = 0
    if verbose:
        print("Total Alleles:", len(genotypeslist)*2)
        print("totalReferenceAlleleCount:", totalReferenceAlleleCount)
        print("totalAlternateAlleleCount:",  totalAlternateAlleleCount)
        print("alternateAlleleFrequency:", totalAlternateAlleleFrequency)
        print("alternateAlleleCount", totalAlternateAlleleCount)
        print("homozygousAlternateAlleleCount:", homozygousAlternateAlleleCount)
        print("homozygousAlternateAlleleFrequency:", homozygousAlternateAlleleFrequency)
        print("heterozygousAlternateAlleleCount:", heterozygousAlternateAlleleCount)
        print("heterozygousAlternateAlleleFrequency:", heterozygousAlternateAlleleFrequency)
    #
    return totalAlternateAlleleFrequency, homozygousAlternateAlleleFrequency, heterozygousAlternateAlleleFrequency, totalAlternateAlleleCount, homozygousAlternateAlleleCount, heterozygousAlternateAlleleCount
#
def add_age_of_onset(df):
    ''' Add the age of onset to the dataframe
    Warning: this is quite time consuming, the most time
    consuming part of adding annotations to the VCF file.
        '''
    if 'firstDiagnosisAge' in df.columns.tolist():
        df['age of onset'] = df['firstDiagnosisAge']
        return df
    else:
        df['date de naissance'] = pd.to_datetime(df['date de naissance'], dayfirst=True)
        df['date du 1er diagnostic'] = pd.to_datetime(df['date du 1er diagnostic'], errors='coerce', dayfirst=True)
        df['age of onset'] = df.apply(lambda row: (row['date du 1er diagnostic'].year - row['date de naissance'].year 
                                            - ((row['date du 1er diagnostic'].month, row['date du 1er diagnostic'].day) 
                                                < (row['date de naissance'].month, row['date de naissance'].day))) 
                                if pd.notna(row['date du 1er diagnostic']) else pd.NA, axis=1)
        return df

#
def check_values_equal(val1, val2):
    '''
    '''
    if val1 != val2:
        print(f"Values are not equal: {val1} != {val2}")
        sys.exit(1)
    print("Values are equal.")
#
if __name__ == "__main__":
    print("This script is intended to be used as a module.")

# color text
def color_truefalse(text):
    """Returns text wrapped in ANSI color codes."""
    color_code = 31 if text == "True" else 32
    return f"\033[{color_code}m{text}\033[0m"