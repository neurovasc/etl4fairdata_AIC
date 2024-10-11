# author: bodrug-a, chatgpt, github Copilot
# 2024-08-01
# This script takes in a VCF file and converts it to RDF-ttl format
# Unithread processing is (should be) limited to 100 variants
# Multithreading is used for processing more than 100 variants

import os
import io
import glob
import pandas as pd
from rdflib import Graph, URIRef, Literal, Namespace
from rdflib.namespace import RDF, XSD
import argparse
import concurrent.futures
import threading



# Arguments
#
def limit_type(x):
    x = int(x)
    if x > 100:
        raise argparse.ArgumentTypeError("Maximum unithread processing is 100 variants.")
    return x
agp = argparse.ArgumentParser(description='Convert VCF to RDF-ttl')
agp.add_argument('-r', '--varrdf', type=str, help='Path to the output RDF file')
agp.add_argument('-b', '--bcf', type=str, help='Path to the BCF file') 
agp.add_argument('-l', '--limit', type=limit_type, help='Number of variants to process')
agp.add_argument('-g', '--gnomad', type=str, help='Path to the gnomad file', \
                 default='/data/gnomad/gnomad.exomes.v4.1.sites.chr1.vcf.bgz')
agp.add_argument('-d', '--dbsnp', type=str, help='Path to the dbsnp file', \
                 default='/data/dbsnp/GCF_000001405.40.gz')
agp.add_argument('-x', '--temp', type=str, help='Temp folder', \
                 default='temp/')
agp.add_argument('-t', '--thread', type=int, help='Number of threads to use', default=4)    
agp.add_argument('-k', '--chunksize', type=int, help='Chunk size for processing', default=10)

args = agp.parse_args()

stream = io.BytesIO()

# Define namespaces
# Custom
ican = Namespace("http://ican.ressource.org/")
# Bio
so = Namespace("http://purl.obolibrary.org/obo/SO_")
sio = Namespace("http://semanticscience.org/resource/SIO_")
dbsnp = Namespace("http://bio2rdf.org/dbsnp:")
geno = Namespace("http://purl.obolibrary.org/obo/GENO_")
# Semantics
rdf = Namespace("http://www.w3.org/1999/02/22-rdf-syntax-ns#")
rdfs = Namespace("http://www.w3.org/2000/01/rdf-schema#")
xsd = Namespace("http://www.w3.org/2001/XMLSchema#")

def create_rdfgraph_namespace():
# Create a Graph
    g = Graph()

    # Bind prefixes
    g.bind("ican", ican)
    g.bind("so", so)
    g.bind("sio", sio)
    g.bind("dbsnp", dbsnp)
    g.bind("geno", geno)
    g.bind("rdf", rdf)
    g.bind("rdfs", rdfs)
    g.bind("xsd", xsd)

    return g

def initiate_rdf_from_vcf(bcf, varrdf):
    '''
    Build the rdf graph from the vcf file
    '''
    # Convert vcf.gz or bcf to bcftools query file 
    queryvcf = args.temp + '/temp_vcf.query'
    os.system('bcftools query'
                " -f'%CHROM\t%POS\t%REF\t%ALT\t%INFO\n' "  +
                bcf + '>' + queryvcf)
    # Read the queryvcf file and build the graph
    # No multithreading
    if args.limit:
        g1 = create_rdfgraph_namespace()
        df = pd.read_csv(queryvcf, sep='\t', \
                        names=['chromosome', 'position', 'reference', 'alternate', 'info'], \
                        nrows=args.limit)
        rowg = build_rdfgraph(g1, df) # this functions modifies the graph
        g1.serialize(destination=varrdf, format="turtle")
    # Multithreading : sometimes not worth it
    else:
        # Clean any previously existing temp turtle files
        clean_temp_ttlfiles()
        #
        chunksize = args.chunksize
        with open(queryvcf, 'r') as file:
            lines = file.readlines()
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.thread) as executor:
            print("Submitting to concurrent executor with", args.thread, "threads.")
            futures = [executor.submit(process_queryline, lines[chunk:chunk+chunksize], chunk) for chunk in range(0, len(lines), chunksize)]
            print("Finished building futures. Waiting for results.")
            results = [future.result() for future in concurrent.futures.as_completed(futures)]
        nb_chunks = len(results)
        print(f"\nFinished processing all {nb_chunks}.")
        print("Merging all temp .ttl files...")
        count = 0
        # The turtle files should not be loaded into memory
        # but rather serialized to the final file
        temporary_merge = args.temp + '/temp_merge.ttl'
        # Empty file if it already exists
        if os.path.exists(temporary_merge):
            open(temporary_merge, 'w').close()
        # Get all the turtle files in temp
        ttlfiles = glob.glob(args.temp + '/' +str(os.path.basename(args.varrdf))+'_*_intermediate.ttl')
        # Merge all the turtle files into one, without loading them into memory all at once
        with open(temporary_merge, 'w') as f:
            for ttl in ttlfiles:
                count += 1
                print(f"\rMerging files: {count}/{nb_chunks}", end="")
                graph = Graph()
                graph.parse(ttl, format='turtle')
                f.write(graph.serialize(format='turtle'))
        # Create graph for final serialization
        g2 = create_rdfgraph_namespace() # merged graph
        g2.parse(temporary_merge, format='turtle')
        g2.serialize(varrdf, format="turtle")
        print("\nDone.")
        print(f"Checking if valid turtle file: {varrdf}.")
        if check_turtle(varrdf):
            print("Valid turtle file.")
        else:
            print("Invalid turtle file.")
#
def check_turtle(file):
    '''
    Check if the file is a valid turtle file
    '''
    if not os.path.isfile(file):
        print(f"File does not exist: {file}")
        return False
    try:
        g = Graph()
        g.parse(file, format='turtle') # this will load everything into memory
        return True
    except Exception as e:
        print(f"Error: {e}")
        return False
#
def process_queryline(linechunk, chunk):
    '''
    Process the query line from the bcftools query
    For each bcftools query line a small ttl file is created
    With the namespace and the info of the one variant
    '''
    minig = create_rdfgraph_namespace()
    output = args.temp + '/' +str(os.path.basename(args.varrdf))+'_'+str(chunk)+'_intermediate.ttl'

    df = pd.DataFrame()
    for line in linechunk:
        #print(line)
        dfl = pd.read_csv(io.StringIO(line), sep='\t', \
                        names=['chromosome', 'position', 'reference', 'alternate', 'info'])
        df = pd.concat([df, dfl])
    #     
    minig = build_rdfgraph(minig, df)
    #
    minig.serialize(destination=output, format="turtle")
    print(f"\rWriting to output files: {output}", end="")
    return minig
#
def clean_temp_ttlfiles():
    ''' Deleting intermediate turtle files created during multithreading
    '''
    temp_folder = args.temp
    temp_files = glob.glob(os.path.join(temp_folder, '*.ttl'))
    for file in temp_files:
        os.remove(file)
#
def build_rdfgraph(g, df):
    for index, row in df.iterrows():
        # variant object
        '''
        variant = Variant.from_df(row)
        '''
        #
        chromosome = row['chromosome']
        position = int(row['position'])
        reference = row['reference']
        alternate = row['alternate']
        info = row['info']
        #print(chromosome, position, reference, alternate, info)
        annotations = info.split(';')
        fq_gnomad = 'nan'
        fq_ican_AF_whole = 'nan'
        fq_ican_AF_female = 'nan'
        fq_ican_AF_male = 'nan'
        fq_ican_AF_aht = 'nan'
        fq_ican_AF_diabetes = 'nan'
        fq_ican_AF_dyslipidemia = 'nan'
        fq_ican_AF_obese = 'nan'
        fq_ican_AF_overweight = 'nan'
        fq_ican_AF_sporadic = 'nan'
        fq_ican_AF_familial = 'nan'
        fq_ican_AF_ruptured = 'nan'
        fq_ican_AF_multipleica = 'nan'
        csq = 'nan'
        caddraw = 'nan'
        caddphred = 'nan'
        for a in annotations:
            if 'gnomad_AF=' in a:
                fq_gnomad = a.split('=')[1]
            if 'AF_whole=' in a:
                fq_ican_AF_whole = a.split('=')[1]
            if 'AF_female=' in a:
                fq_ican_AF_female = a.split('=')[1]
            if 'AF_male=' in a:
                fq_ican_AF_male = a.split('=')[1]
            if 'AF_aht=' in a:
                fq_ican_AF_aht = a.split('=')[1]
            if 'AF_diabetes=' in a:
                fq_ican_AF_diabetes = a.split('=')[1]
            if 'AF_dyslipidemia=' in a:
                fq_ican_AF_dyslipidemia = a.split('=')[1]
            if 'AF_obese=' in a:
                fq_ican_AF_obese = a.split('=')[1]
            if 'AF_overweight=' in a:
                fq_ican_AF_overweight = a.split('=')[1]
            if 'AF_sporadic=' in a:
                fq_ican_AF_sporadic = a.split('=')[1]
            if 'AF_familial=' in a:
                fq_ican_AF_familial = a.split('=')[1]
            if 'AF_ruptured=' in a:
                fq_ican_AF_ruptured = a.split('=')[1]
            if 'AF_multipleica=' in a:
                fq_ican_AF_multipleica = a.split('=')[1]
            if 'CSQ=' in a :
                consequences = a.split('=')[1]
            if 'CADD_RAW=' in a:
                caddraw = a.split('=')[1]
            if 'PHRED=' in a:
                caddphred = a.split('=')[1]
                
        rsid = get_rsids_fromdbsnp(chromosome, position, reference, alternate)

        # Define the variant URI
        lilprefix = 'icanexomehg38' 
        variant_id = f"{lilprefix}-{chromosome}-{position}-{reference}-{alternate}"
        variant_reference_id = f"{lilprefix}-{chromosome}-{position}-{reference}-{alternate}"
        variant_alternate_id = f"{lilprefix}-{chromosome}-{position}-{reference}-{alternate}"
        chromosome_id = f"{lilprefix}-{chromosome}"
        position_id = f"{lilprefix}-{chromosome}-{position}"
        # variant, ref and alt
        variant_uri = ican["variantId/"+variant_id]
        variant_reference_uri = ican["variantReference/"+variant_reference_id]
        variant_alternate_uri = ican["variantAlternate/"+variant_alternate_id]
        # chromosome and position
        chromosome_uri = ican["chromosome/"+chromosome_id]
        position_uri = ican["chromosome/"+chromosome_id+"/position/"+position_id]
        # Define the variant frequency URIs
        variant_fq_gnomad_uri = ican["variantAlternate/"+'fq/gnomad/'+variant_alternate_id]
        variant_fq_ican_uri = ican["variantAlternate/"+'fq/ican/cohort/'+variant_alternate_id]
        variant_fq_ican_uri_female = ican["variantAlternate/"+'fq/ican/female/'+variant_alternate_id]
        variant_fq_ican_uri_male = ican["variantAlternate/"+'fq/ican/male/'+variant_alternate_id]
        variant_fq_ican_uri_ruptured = ican["variantAlternate/"+'fq/ican/ruptured/'+variant_alternate_id]
        variant_fq_ican_uri_multipleica = ican["variantAlternate/"+'fq/ican/multipleica/'+variant_alternate_id]
        variant_fq_ican_uri_aht = ican["variantAlternate/"+'fq/ican/aht/'+variant_alternate_id]
        variant_fq_ican_uri_diabetes = ican["variantAlternate/"+'fq/ican/diabetes/'+variant_alternate_id]
        variant_fq_ican_uri_dyslipidemia = ican["variantAlternate/"+'fq/ican/dyslipidemia/'+variant_alternate_id]
        variant_fq_ican_uri_obese = ican["variantAlternate/"+'fq/ican/obese/'+variant_alternate_id]
        variant_fq_ican_uri_overweight = ican["variantAlternate/"+'fq/ican/overweight/'+variant_alternate_id]
        variant_fq_ican_uri_sporadic = ican["variantAlternate/"+'fq/ican/sporadic/'+variant_alternate_id]
        variant_fq_ican_uri_familial = ican["variantAlternate/"+'fq/ican/familial/'+variant_alternate_id]
        # Define the variant annotation URIs
        variant_annotation_uri = ican["variantAlternate/"+'annotation/'+variant_alternate_id]
        caddraw_uri = variant_annotation_uri+'/cadd/'
        caddphred_uri = variant_annotation_uri+'/caddphred/'

        # Add variant triplets
        # Story telling

        # The variant is a so:variant
        g.add((variant_uri, RDF.type, so["0001060"]))
        # The variant sio:is_located_in the chromosome
        g.add((variant_uri, sio["SIO_000061"], chromosome_uri))
        # The chromosome has a value that is a string
        g.add((chromosome_uri, sio["SIO_000300"], Literal(chromosome, datatype=XSD.string))) # is located in
        # The variant has sio:start_position at position
        g.add((variant_uri, sio["SIO_000791"], position_uri))
        # The position has a value that is an integer
        g.add((position_uri, sio["SIO_000300"], Literal(position, datatype=XSD.integer))) # at position
        # The variant has a sio:synonymous id in dbsnp
        if rsid != 'nan':
            g.add((variant_uri, sio["SIO_000122"], dbsnp[rsid])) # is synonymous to
        # The variant sio:has_attributes that are the reference and alternate alleles
        g.add((variant_uri, sio["SIO_000223"], variant_reference_uri)) # has property
        g.add((variant_uri, sio["SIO_000223"], variant_alternate_uri)) # has property
        # The reference is a geno:reference_allele
        g.add((variant_reference_uri, RDF.type, geno["0000036"])) # reference allele
        # The reference has a value that is a string (nucleotide)
        g.add((variant_reference_uri, sio["SIO_000300"], Literal(reference, datatype=XSD.string))) # reference nucleotide
        # The reference has a source, but we don't add it for now, as we work only with a single sourcr
        # that is the ican cohort.
        #g.add((variant_reference_uri, sio["SIO_000253"], ican['cohort.ican'])) # data source
        g.add((variant_reference_uri, rdfs['label'], Literal("Reference allele for the variant "+variant_id + ' [ dbspn:'+rsid+' ]', datatype=XSD.string)))
        # The alternate is a geno:alternate_allele
        g.add((variant_alternate_uri, RDF.type, geno["0000002"])) # alternate allele
        g.add((variant_alternate_uri, sio["000300"], Literal(alternate, datatype=XSD.string))) # alternate nucleotide
        g.add((variant_alternate_uri, sio["000253"], ican['cohort.ican'])) # data source
        g.add((variant_alternate_uri, rdfs['label'], Literal("Alternate allele for the variant "+variant_id + ' [ dbspn:'+rsid+' ]', datatype=XSD.string)))
        g.add((variant_alternate_uri, sio["SIO_000900"], variant_fq_gnomad_uri)) # alternate allele frequencies
        g.add((variant_alternate_uri, sio["SIO_000900"], variant_fq_ican_uri))
        g.add((variant_alternate_uri, sio["SIO_000900"], variant_fq_ican_uri_female))
        g.add((variant_alternate_uri, sio["SIO_000900"], variant_fq_ican_uri_male))
        # Add alternate allele info with frequencies
        # gnomad
        g.add((variant_fq_gnomad_uri, RDF.type, sio["SIO_001367"]))
        g.add((variant_fq_gnomad_uri, sio["000253"], ican['cohort/gnomad']))
        g.add((variant_fq_gnomad_uri, sio["000300"], Literal(fq_gnomad, datatype=XSD.float)))
        # ican
        g.add((variant_fq_ican_uri, RDF.type, sio["SIO_001367"]))
        g.add((variant_fq_ican_uri, sio["000253"], ican['cohort/ican']))
        g.add((variant_fq_ican_uri, sio["000300"], Literal(fq_ican_AF_whole, datatype=XSD.float)))
        # ican females
        g.add((variant_fq_ican_uri_female, RDF.type, sio["SIO_001367"]))
        g.add((variant_fq_ican_uri_female, sio["000253"], ican['cohort/ican/female']))
        g.add((variant_fq_ican_uri_female, sio["000300"], Literal(fq_ican_AF_female, datatype=XSD.float)))
        g.add((variant_fq_ican_uri_female, sio["000628"], sio["SIO_010052"]))
        # ican males
        g.add((variant_fq_ican_uri_male, RDF.type, sio["SIO_001367"]))
        g.add((variant_fq_ican_uri_male, sio["000253"], ican['cohort/ican/male']))
        g.add((variant_fq_ican_uri_male, sio["000300"], Literal(fq_ican_AF_male, datatype=XSD.float)))
        g.add((variant_fq_ican_uri_male, sio["000628"], sio["SIO_010048"]))

        #

    return g
#
def calculate_fqgnomad(chromosome, position, reference, alternate):
    ''' Not used for now, as the gnomad frequencies are present in
    the INFO field of the OG VCF I am working with, and I just parse them
    '''
    gnomadfile = args.gnomad
    fq = 'nan'
    output = os.popen("bcftools query -r "+str(chromosome)+":"+str(position)+" -f'%CHROM\t%POS\t%REF\t%ALT\t%AF\n' "  +
                    gnomadfile).read()
    df = pd.read_csv(io.StringIO(output), sep='\t', names=['chromosome', 'position', 'reference', 'alternate', 'af'])
    for index, row in df.iterrows():
        if row['position'] == position and row['reference'] == reference and row['alternate'] == alternate:
            if type(row['af']) == str:
                return 'nan'
            else:
                return round(row['af'], 6)
    return fq
#
def get_rsids_fromdbsnp(chromosome, position, reference, alternate):
    ''' Apparently not that useful according to R. Blanchet
    '''
    dbsnpfile = args.dbsnp
    rsid = 'nan'
    dict = {"chr1" : "NC_000001.11", "chr2" : "NC_000002.12", "chr3" : "NC_000003.12", 
            "chr4" : "NC_000004.12", "chr5" : "NC_000005.10", "chr6" : "NC_000006.12", 
            "chr7" : "NC_000007.14", "chr8" : "NC_000008.11", "chr9" : "NC_000009.12", 
            "chr10" : "NC_000010.11", "chr11" : "NC_000011.10", "chr12" : "NC_000012.12", 
            "chr13" : "NC_000013.11", "chr14" : "NC_000014.9", "chr15" : "NC_000015.10", 
            "chr16" : "NC_000016.10", "chr17" : "NC_000017.11", "chr18" : "NC_000018.10", 
            "chr19" : "NC_000019.10", "chr20" : "NC_000020.11", "chr21" : "NC_000021.9", 
            "chr22" : "NC_000022.11", "chrX" : "NC_000023.10", "chrY" : "NC_000024.9"}
    output = os.popen("bcftools query -r "+str(dict[chromosome])+":"+str(position)+" -f'%CHROM\t%POS\t%REF\t%ALT\t%ID\n' "  +
                    dbsnpfile).read()
    df = pd.read_csv(io.StringIO(output), sep='\t', names=['chromosome', 'position', 'reference', 'alternate', 'rsid'])

    for index, row in df.iterrows():
        if row['position'] == position and row['reference'] == reference and row['alternate'] == alternate:
            return row['rsid']
        else:
            return rsid
    return rsid
#
def bcf2query(bcf, samples):
    ''' Takes in a bcf file and queries it with bcftools
    In a way that makes it easier to parse for rdf creation.
    The input bcf should have multiallelic positions split by allele.
    '''
    queryvcf = bcf +'.samples.query'
    if args.convert:
        os.system('bcftools query -S ' +
                samples +
                " -f'%CHROM\t%POS\t%REF\t%ALT\t%INFO\n' "  +
                bcf + '>' + queryvcf)
    return queryvcf
#
'''
class Variant:
    def __init__(self, chromosome, position, ref, alt, info):
        self.chromosome = chromosome
        self.position = position
        self.ref = ref
        self.alt = alt
        self.info = info
 
    @classmethod
    def get_frequencies_from_info(self.info):
        return Variant
'''
#
if __name__ == "__main__":
    # Info about parameters
    print("Script launched with:")
    bcf = args.bcf ; print(f"\tinput - BCF: {bcf}.")
    #
    varrdf = args.varrdf
    # Build rdf from vcf aggregagte
    initiate_rdf_from_vcf(bcf, varrdf)
    print("Wrote output RDF to:")
    print(f"\toutput - ttl: {varrdf}.")
    if args.thread:
        print("Cleaning temp files...\r")
        clean_temp_ttlfiles()
        print("Done.")
