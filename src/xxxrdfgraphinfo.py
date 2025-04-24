from rdflib import Namespace, Literal, BNode
import pandas as pd

# namespace
AGGRVAR = Namespace('https://ican.univ-nantes.io/variants-kg-schema/')
ICAN = Namespace('http://ican.ressource.org/')
FALDO = Namespace('http://biohackathon.org/resource/faldo#')
GENO = Namespace('http://purl.obolibrary.org/obo/GENO_')
HPO = Namespace('http://purl.obolibrary.org/obo/HP_')
LINKML = Namespace('https://w3id.org/linkml/')
MONDO = Namespace('http://purl.obolibrary.org/obo/MONDO_')
NCIT = Namespace('http://purl.obolibrary.org/obo/NCIT_')
OWL = Namespace('http://www.w3.org/2002/07/owl#')
RDF = Namespace('http://www.w3.org/1999/02/22-rdf-syntax-ns#')
RDFS = Namespace('http://www.w3.org/2000/01/rdf-schema#')
OWL = Namespace('https://www.w3.org/TR/2004/REC-owl-semantics-20040210/#owl_')
SIO = Namespace('http://semanticscience.org/resource/SIO_')
SKOS = Namespace('http://www.w3.org/2004/02/skos/core#')
SO = Namespace('http://purl.obolibrary.org/obo/SO_')
XSD = Namespace('http://www.w3.org/2001/XMLSchema#')
RO = Namespace('http://purl.obolibrary.org/obo/RO_')
HP = Namespace('http://purl.obolibrary.org/obo/HP_')
ORPHA = Namespace('http://www.orpha.net/ORDO/Orphanet_')

namespace = {'faldo' : FALDO, 'aggrvar' : AGGRVAR, 'aic' : ICAN,
             'geno' : GENO, 'hpo' : HPO, 'linkml' : LINKML, 'mondo' : MONDO,
             'ncit' : NCIT, 'owl' : OWL, 'rdf' : RDF, 'rdfs' : RDFS,
             'sio' : SIO, 'skos' : SKOS, 'so' : SO, 'xsd' : XSD, 'ro' : RO,
             'hp' : HP, 'orpha' : ORPHA}
#
# Building nodes
#
def build_variant_node(chromosome, position, reference, alternate):
    ''' Returns variant internal id node, avoid creating empty node
    '''
    return f"viid/{chromosome}_{position}_{reference}_{alternate}"
#
def build_variant_region_node(chromosome, position, reference, alternate):
    ''' Returns variant Region (faldo) node, avoid creating empty node
    Coordinates are inclusive [start, end]
    '''
    position_start, position_end = get_startend_coordinates(position, reference, alternate)
    return f"Region/{chromosome}_{position_start}_{position_end}_{reference}_{alternate}"
#
def build_variant_HGVSid_node(chromosome, position, reference, alternate):
    ''' Returns variant hgvsid node, avoid creating empty node
    '''
    return f"HGVSid/{chromosome}_{position}_{reference}_{alternate}"
#
def build_variant_refallele_n(chromosome, position, reference, alternate):
    '''
    '''
    return f"refallele/{chromosome}_{position}_{reference}_{alternate}"
#
def build_variant_altallele_n(chromosome, position, reference, alternate):
    '''
    '''
    return f"altallele/{chromosome}_{position}_{reference}_{alternate}"
#
def build_region_exactpositions_nandv(chromosome, position, reference, alternate):
    '''
    '''
    position_start, position_end = get_startend_coordinates(position, reference, alternate)
    start_n = f"RegionStart/{chromosome}_{position_start}"
    end_n = f"RegionEnd/{chromosome}_{position_end}"
    start_v = position_start
    end_v = position_end
    return start_n, end_n, start_v, end_v
#
def build_referencesequence_n(chromosome):
    '''
    '''
    return f"RegionSequence/{chromosome}"
#
def build_cohort_n():
    ''' This just builds the name of the cohort node.
    The Cohort's description with attributes and Ontology terms
    describing the type of patients it contains is build
    with a separate script, in a separate graph (future works)
    # TODO: Add the cohort description in semantic ontology terms
    '''
    return f"Cohort/iCAN/cases"
#
def build_subpopulation_n(phenotype):
    '''
    '''
    return f"Subpopulation/{phenotype}"
#
def build_phenotype_n(phenotype):
    '''
    '''
    phenocodes = {'female': NCIT['C16576'], 
                  'male' : NCIT['C20197'],
                  'aht': HP['0000822'], 
                  'obese': NCIT['C159658'],
                  'currentlysmoking': NCIT['C17934'],
                  'familial': ORPHA['231160']}
    #print(f"Phenotype: {phenotype}")
    #print(f"Phenotype code: {phenocodes[phenotype]}")
    #print(f"Phenotype codes: {phenocodes}")
    phenotype_code = phenocodes[phenotype]
    return phenotype_code

#
# Building strings / values
#
def build_variant_HGVSid(chromosome, position, reference, alternate):
    ''' Build HGVSid for hg38
    '''
    chromosome = chromosome.strip('chr')
    chridandmore = pd.read_csv('datasets/hg38_ids_report_sequences.tsv', sep='\t') # hardcoded
    # source : https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/
    refseqids = dict(zip(chridandmore['Sequence name'], chridandmore['RefSeq seq accession']))

    HGVSid = pd.NA
    if len(reference) == 1 and len(alternate) == 1:
        # Substitution: g.123A>G
        HGVSid = f"{refseqids[chromosome]}:g.{position}{reference}>{alternate}"
    elif len(reference) > 1 and alternate == reference[0]:
        # Deletion: g.123_125del
        start = position + 1
        end = position + len(reference) - 1
        HGVSid = f"{refseqids[chromosome]}:g.{start}_{end}del"
    elif reference == reference[0] and len(alternate) > 1:
        # Insertion: g.123_124insGTA (inserts between pos and pos+1)
        inserted_seq = alternate[1:]
        start = position
        end = position + 1
        HGVSid = f"{refseqids[chromosome]}:g.{start}_{end}ins{inserted_seq}"
    else:
        raise ValueError("Unsupported variant type or inconsistent REF/ALT format.")
    return f"HGVSid:{HGVSid}"
#
def build_variant_refallele(chromosome, position, reference, alternate):
    return reference
#
def build_variant_altallele(chromosome, position, reference, alternate):
    return alternate
#
def build_referencesequence_genbanid(chromosome):
    chromosome = chromosome.strip('chr')
    chridandmore = pd.read_csv('datasets/hg38_ids_report_sequences.tsv', sep='\t') # hardcoded
    # source : https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/
    genbankid = dict(zip(chridandmore['Sequence name'], chridandmore['GenBank seq accession']))
    return genbankid[chromosome]
#
def build_referencesequence_refseqid(chromosome):
    chromosome = chromosome.strip('chr')
    chridandmore = pd.read_csv('datasets/hg38_ids_report_sequences.tsv', sep='\t') # hardcoded
    # source : https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/
    refseqid = dict(zip(chridandmore['Sequence name'], chridandmore['RefSeq seq accession']))
    return refseqid[chromosome] 
#
# Other utilitary functions
def get_startend_coordinates(position, reference, alternate):
    refallelelen = len(reference)
    altallelelen = len(alternate)
    span = max(refallelelen, altallelelen) - 1
    position_start = position
    position_end = position + span
    return position_start, position_end
#
# Info fields
def get_info_fields_labels_cohortAFAC(info):
    ''' From infos in the INFO field of the VCF, extract the labels
    In this function, only the AF and AC fields from within the cohort itself
    is processed. These can be considered as *internal* annotations
    inherent to the cohort and non obtainable from other external databases
    '''
    infofields = info.split(';')
    array_of_observations = []
    for field in infofields:
        if field.startswith('AF_') or field.startswith('AC_'): # internal annotations start with these
            text = field.split('=')[0]
            value = field.split('=')[1]
            label = text.split('_')[1]
            try:
                zygosity = text.split('_')[2]
            except:
                zygosity = 'unspecified'
            observedFeature = 'feature:'+label+'-'+'zygosity:'+zygosity
            observation = { 'observedFeature': observedFeature,
                    'phenotype': label, 
                    'zygosity': zygosity}
            #
            if field.startswith('AF_'):
                matching_elements = [i for i in array_of_observations if i['observedFeature'] == observedFeature]
                if len(matching_elements) == 0:
                    observation['frequency'] =  value
                    array_of_observations.append(observation)
                elif len(matching_elements) == 1:
                    observation['frequency'] = value
                    array_of_observations[array_of_observations.index(matching_elements[0])]['frequency'] = value
                else:
                    print("Exiting because more than two matching observations?")
                    exit()
            #
            if field.startswith('AC_'):
                matching_elements = [i for i in array_of_observations if i['observedFeature'] == observedFeature]
                if len(matching_elements) == 0:
                    observation['count'] = value
                    array_of_observations.append(observation)
                elif len(matching_elements) == 1:
                    observation['count'] = value
                    array_of_observations[array_of_observations.index(matching_elements[0])]['count'] = value
                else:
                    print("Exiting because more than two matching observations?")
                    exit()

    return array_of_observations
#
def build_observation_n(variant_iid, observation):
    '''
    '''
    #
    return f"{variant_iid}/ObservedFeature/{observation}"
#
def build_observation_frequency_n(variant_iid, observation):
    '''
    '''
    return f"{variant_iid}/ObservedFeature/{observation}/Frequency"
#
def build_observation_count_n(variant_iid, observation):
    '''
    '''
    return f"{variant_iid}/ObservedFeature/{observation}/Count"

#
# Building the whole graph
def build_rdfgraph(g, df):
    ''' Follwing schema: 
    schema: https://docs.google.com/drawings/d/1xfawlZxZgUYMsIuDHQgSAZnKI3FV0lTmR7JQKXJ3v58
    '''
    for index, row in df.iterrows():
        # Variant chromosome, position, ref, alt, info
        chromosome = row['chromosome'].strip('chr')
        position = int(row['position'])
        reference = row['reference']
        alternate = row['alternate']
        info = row['info']
        ############################################################
        # Variant and its immediate attributes excuding sio:001403
        # Variant HGVSid, Region, Reference Allele and Alternate Allele
        #
        variant_iid = build_variant_node(chromosome, position, reference, alternate) # variant node
        variant_region_n = build_variant_region_node(chromosome, position, reference, alternate) # region node 
        variant_hgvsid_n = build_variant_HGVSid_node(chromosome, position, reference, alternate) # hgvsid node
        variant_refallele_n = build_variant_refallele_n(chromosome, position, reference, alternate) # refallele node
        variant_altallele_n = build_variant_altallele_n(chromosome, position, reference, alternate) # altallele node
        #
        # schema: light pink
        g.add((ICAN[variant_iid], RDF.type, SO["0001059"])) # type
        g.add((ICAN[variant_iid], RDF.type, GENO["0000660"])) # type
        g.add((ICAN[variant_iid], FALDO['location'], ICAN[variant_region_n])) # attribute 1
        g.add((ICAN[variant_iid], SIO['000671'], ICAN[variant_hgvsid_n])) # attribute 2
        g.add((ICAN[variant_iid], GENO['0000385'], ICAN[variant_refallele_n])) # attribute 3
        g.add((ICAN[variant_iid], GENO['0000382'], ICAN[variant_altallele_n])) # attribute 4
        #
        ############################################################
        # Alternate and Reference alleles
        #
        variant_refallele = build_variant_refallele(chromosome, position, reference, alternate) # refallele value
        variant_altallele = build_variant_altallele(chromosome, position, reference, alternate) # altallele value
        #
        # schema: light yellow
        g.add((ICAN[variant_refallele_n], RDF.type, GENO["0000036"])) # refallele type
        g.add((ICAN[variant_altallele_n], RDF.type, GENO["0000002"])) # altallele type
        g.add((ICAN[variant_refallele_n], SIO['000300'], Literal(variant_refallele, datatype=XSD.string)))
        g.add((ICAN[variant_altallele_n], SIO['000300'], Literal(variant_altallele, datatype=XSD.string)))
        #
        ############################################################
        # Variant HGVSid is an identifiers
        #
        variant_hgvsid_v = build_variant_HGVSid(chromosome, position, reference, alternate)
        #
        # schema: light blue
        g.add((ICAN[variant_hgvsid_n], RDF.type, SIO['000675']))
        g.add((ICAN[variant_hgvsid_n], SIO['000300'], Literal(variant_hgvsid_v, datatype=XSD.string)))
        #
        ############################################################
        # Variant Region and its attributes
        #
        region_exactposition_start_n, region_exactposition_end_n, region_exactposition_start, region_exactposition_end = build_region_exactpositions_nandv(chromosome, position, reference, alternate)
        region_referencesequence_n = build_referencesequence_n(chromosome)
        region_referencesequence_genbankid = build_referencesequence_genbanid(chromosome)
        region_referencesequence_refseqid = build_referencesequence_refseqid(chromosome)
        #
        # schema: light green
        g.add((ICAN[variant_region_n], RDF.type, FALDO['Region']))
        g.add((ICAN[variant_region_n], FALDO['begin'], ICAN[region_exactposition_start_n]))
        g.add((ICAN[variant_region_n], FALDO['end'], ICAN[region_exactposition_end_n]))
        g.add((ICAN[variant_region_n], FALDO['reference'], ICAN[region_referencesequence_n]))
        g.add((ICAN[region_exactposition_start_n], RDF.type, FALDO['ExactPosition']))
        g.add((ICAN[region_exactposition_end_n], RDF.type, FALDO['ExactPosition']))
        #
        # schema: light purple
        g.add((ICAN[region_referencesequence_n], RDF.type, SO['0000353']))
        g.add((ICAN[region_referencesequence_n], RDF.label, Literal(chromosome, datatype=XSD.string))) # label 1 or 22
        g.add((ICAN[region_referencesequence_n], SO['000300'], Literal(region_referencesequence_genbankid, datatype=XSD.string))) # ena sequenceid
        g.add((ICAN[region_referencesequence_n], OWL.sameAs, Literal(region_referencesequence_refseqid, datatype=XSD.string))) # refseqid
        #
        # schema: light cyan
        g.add((ICAN[region_exactposition_start_n], FALDO['position'], Literal(region_exactposition_start, datatype=XSD.integer)))
        g.add((ICAN[region_exactposition_end_n], FALDO['position'], Literal(region_exactposition_end, datatype=XSD.integer)))        
        #
        ##############################################################
        # Linking variant to Cohort and Subpopulations' Observations
        #
        # Variant to Cohort
        cohort = build_cohort_n() # Details can be build in another script, separately
        #
        # schema: light orange
        g.add((ICAN[variant_iid], SIO['001403'], ICAN[cohort])) # variant is associated with Cohort
        g.add((ICAN[cohort], RDF.type, NCIT['C61512'])) # ICAN is a cohort
        #
        # Observations in info field of VCF
        # An observation is a combination of a variant, a phenotype (or several phenotypes) and a zygosity
        # It can have an AF and an AC, or just an AC 
        #
        infofields = get_info_fields_labels_cohortAFAC(info)
        observed_features = list(set(i['observedFeature'] for i in infofields))

        for observation in observed_features:
                # Matching observation 
                matching_observation = [i for i in infofields if i['observedFeature'] == observation]
                # An observation is a combination of (a variant), a phenotype and a zygosity
                #print(f"Building graph for {variant_hgvsid_v} within {observation} .")
                if len(matching_observation) > 1:
                    print("More than one observation for this variant+phenotype+zygosity combo.")
                    print(observation)
                    print(variant_iid)
                    print(matching_observation)
                    print("This is not expected.")
                    print("Exiting.")
                    exit()
                else:
                    matching_observation = matching_observation[0]
                #
                phenotype = matching_observation['phenotype'] # obese --> value
                zygosity = matching_observation['zygosity'] # heterozygous --> value
                observation_n = build_observation_n(variant_iid, observation) # variant_iid / feature:obese - zygosity:het
                observationFrequency_n = build_observation_frequency_n(variant_iid, observation)
                observationCount_n = build_observation_count_n(variant_iid, observation)
                #
                # Variant to Observation 
                # Observation to its attributes: count, frequency, zygosity and Subpopulation
                # schema: quite obvious red
                #
                try:
                    frequency = float(matching_observation['frequency']) # --> value
                except:
                    frequency = 0
                try:
                    count = int(matching_observation['count']) # --> value
                except:
                    count = 0
                # Do not build an observation if count or frequency are 0
                # It will save time and space
                # OR: because an observation does not always have count and frequecy
                if count > 5 and frequency > 0:
                    g.add((ICAN[variant_iid], SIO['001403'], ICAN[observation_n])) # variant is associated with Observation
                    g.add((ICAN[observation_n], RDF.type, SIO['000649'])) # Observation is an sio:000649 'information_processing' class
                    #
                    # Frequency and Counts are not always both present for observations.
                    #
                    # Frequency
                    # schema: obviously vibrant blue
                    try:
                        #frequency = matching_observation['frequency'] # --> value
                        # Frequency
                        g.add((ICAN[observation_n], SIO['000900'], ICAN[observationFrequency_n])) # observation has frequency
                        g.add((ICAN[observationFrequency_n], RDF.type, SIO['001367'])) # frequency is frequency
                        g.add((ICAN[observationFrequency_n], SIO['000300'],Literal(frequency, datatype=XSD.float)))
                    except: 
                        pass
                    # Count
                    # schema: obviously vibrant blue
                    try:
                        #count = matching_observation['count'] # --> value
                        # Count
                        g.add((ICAN[observation_n], SIO['000216'], ICAN[observationCount_n])) # has count (measurement value)
                        g.add((ICAN[observationCount_n], RDF.type, SIO['000794'])) # count is count
                        g.add((ICAN[observationCount_n], SIO['000300'], Literal(count, datatype=XSD.integer)))
                    except:
                        pass
                    #
                    # Zygosity information is alsways present for an observation
                    #
                    # Zygosity
                    # schema: obviously vibrant blue
                    zygosityOntology = {'unspecified' : '0000137', 'het' : '0000135', 'hom' : '0000136'} # GENO
                    g.add((ICAN[observation_n], GENO['0000608'], GENO[zygosityOntology[zygosity]]))
                    #
                    # Subpopulation is always present for an observation
                    #
                    # Subpopulation
                    # schema: grass green, quite green
                    subpopulation = build_subpopulation_n(phenotype) # feature:X
                    # 
                    # observation is associated with population group
                    g.add((ICAN[observation_n], SIO['001403'], ICAN[subpopulation]))
                    g.add((ICAN[subpopulation], RDF.type, NCIT['C17005'])) # subpopulation is a NCIT's "Population group"
                    g.add((ICAN[subpopulation], RDF.label, Literal(phenotype, datatype=XSD.string)))
                    try:
                        phenotype_n_code = build_phenotype_n(phenotype) # feature:X
                        g.add((ICAN[subpopulation], RO['0016001'], phenotype_n_code)) # has phenotype or disease
                    except:
                        #print(f"{phenotype}: phenotype not in dictionary")
                        pass 
                else:
                    pass
        #       
    return g