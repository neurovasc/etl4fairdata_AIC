from rdflib import Namespace, Literal
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
SIO = Namespace('http://semanticscience.org/resource/SIO_')
SKOS = Namespace('http://www.w3.org/2004/02/skos/core#')
SO = Namespace('http://purl.obolibrary.org/obo/SO_')
XSD = Namespace('http://www.w3.org/2001/XMLSchema#')

namespace = {'faldo' : FALDO, 'aggrvar' : AGGRVAR, 'aic' : ICAN,
             'geno' : GENO, 'hpo' : HPO, 'linkml' : LINKML, 'mondo' : MONDO,
             'ncit' : NCIT, 'owl' : OWL, 'rdf' : RDF, 'rdfs' : RDFS,
             'sio' : SIO, 'skos' : SKOS, 'so' : SO, 'xsd' : XSD}
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
    refallelelen = len(reference)
    altallelelen = len(alternate)
    span = max(refallelelen, altallelelen) - 1
    position_start = position
    position_end = position + span
    return f"Region/{chromosome}_{position_start}_{position_end}_{reference}_{alternate}"
#
def build_variant_HGVSid_node(chromosome, position, reference, alternate):
    ''' Returns variant hgvsid node, avoid creating empty node
    '''
    return f"HGVSid/{chromosome}_{position}_{reference}_{alternate}"
#
def build_variant_refallele_n(chromosome, position, reference, alternate):
    return f"refallele/{chromosome}_{position}_{reference}_{alternate}"
#
def build_variant_altallele_n(chromosome, position, reference, alternate):
    return f"altallele/{chromosome}_{position}_{reference}_{alternate}"
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
    HGVSid = f"HGVSid:{refseqids[chromosome]}:g.{position}{reference}>{alternate}"
    return HGVSid
#
def build_variant_refallele(chromosome, position, reference, alternate):
    return reference
#
def build_variant_altallele(chromosome, position, reference, alternate):
    return alternate
#
# Building the whole graph
def build_rdfgraph(g, df):
    '''
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
        g.add((ICAN[variant_refallele_n], RDF.type, GENO["0000036"])) # refallele type
        g.add((ICAN[variant_altallele_n], RDF.type, GENO["0000002"])) # altallele type
        g.add((ICAN[variant_refallele_n], SIO['000300'], Literal(variant_refallele, datatype=XSD.string)))
        g.add((ICAN[variant_altallele_n], SIO['000300'], Literal(variant_altallele, datatype=XSD.string)))
        #
        ############################################################
        # Variant HGVSid is an identifiers
        #
        variant_hgvsid_v = build_variant_HGVSid(chromosome, position, reference, alternate)
        g.add((ICAN[variant_hgvsid_n], RDF.type, SIO['000675']))
        g.add((ICAN[variant_hgvsid_n], SIO['000300'], Literal(variant_hgvsid_v, datatype=XSD.string)))
        #
        ############################################################
        # Variant Region and its attributes
        #
        g.add((ICAN[variant_region_n], RDF.type, FALDO['Region']))
    return g