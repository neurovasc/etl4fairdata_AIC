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
so = Namespace("http://www.sequenceontology.org/miso/current_svn/term/SO:")
sio = Namespace("http://semanticscience.org/resource/")
dbsnp = Namespace("http://bio2rdf.org/dbsnp:")
geno = Namespace("http://example.org/geno/GENO_")
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
        rowg = build_rdfgraph(g1, df) # this functions modified the graph
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
        g.parse(file, format='turtle')
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
        chromosome = row['chromosome']
        position = int(row['position'])
        reference = row['reference']
        alternate = row['alternate']
        info = row['info']
        #print(chromosome, position, reference, alternate, info)
        fq_gnomad = calculate_fqgnomad(chromosome, position, reference, alternate)
        fq_ican_AF_whole = info.split(';')[0].split('=')[1]
        rsid = get_rsids_fromdbsnp(chromosome, position, reference, alternate)

        # Define the variant URI
        variant_id = f"ican_hg38-{chromosome}-{position}-{reference}-{alternate}"
        variant_reference_id = f"ican_hg38-{chromosome}-{position}-{reference}-{alternate}"
        variant_alternate_id = f"ican_hg38-{chromosome}-{position}-{reference}-{alternate}"
        variant_uri = ican["variantId/"+variant_id]
        variant_reference_uri = ican["variantReference/"+variant_reference_id]
        variant_alternate_uri = ican["variantAlternate/"+variant_alternate_id]
        # Define the variant frequency URIs
        variant_fq_gnomad_uri = ican["variantAlternate/"+'fq/gnomad/'+variant_alternate_id]
        variant_fq_ican_uri = ican["variantAlternate/"+'fq/ican/'+variant_alternate_id]

        # Add variant triplets
        g.add((variant_uri, RDF.type, so["0001060"]))
        g.add((variant_uri, sio["SIO_000061"], Literal(chromosome, datatype=XSD.string))) # is located in
        g.add((variant_uri, sio["SIO_000791"], Literal(position, datatype=XSD.integer))) # at position
        g.add((variant_uri, sio["SIO_000122"], dbsnp[rsid])) # is synonymous to
        g.add((variant_uri, sio["SIO_000223"], variant_reference_uri)) # has property
        g.add((variant_uri, sio["SIO_000223"], variant_alternate_uri)) # has property

        # Add reference description triplets
        g.add((variant_reference_uri, RDF.type, geno["0000152"])) # reference allele
        g.add((variant_reference_uri, sio["SIO_000300"], Literal(reference, datatype=XSD.string))) # reference nucleotide
        g.add((variant_reference_uri, sio["SIO_000253"], ican['cohort.ican'])) # data source
        g.add((variant_reference_uri, rdfs['label'], Literal("Reference allele for the variant "+variant_id + ' [ dbspn:'+rsid+' ]', datatype=XSD.string)))
        # Add alternate description triplets
        g.add((variant_alternate_uri, RDF.type, geno["0000476"])) # alternate allele
        g.add((variant_alternate_uri, sio["000300"], Literal(alternate, datatype=XSD.string))) # alternate nucleotide
        g.add((variant_alternate_uri, sio["000253"], ican['cohort.ican'])) # data source
        g.add((variant_alternate_uri, rdfs['label'], Literal("Alternate allele for the variant "+variant_id + ' [ dbspn:'+rsid+' ]', datatype=XSD.string)))
        g.add((variant_alternate_uri, sio["SIO_000900"], variant_fq_gnomad_uri)) # alternate allele frequencies
        g.add((variant_alternate_uri, sio["SIO_000900"], variant_fq_ican_uri))
        # Add alternate allele info with frequencies
        # gnomad
        g.add((variant_fq_gnomad_uri, RDF.type, sio["SIO_001367"]))
        g.add((variant_fq_gnomad_uri, sio["000253"], ican['cohort/gnomad']))
        g.add((variant_fq_gnomad_uri, sio["000300"], Literal(fq_gnomad, datatype=XSD.float)))
        # ican
        g.add((variant_fq_ican_uri, RDF.type, sio["SIO_001367"]))
        g.add((variant_fq_ican_uri, sio["000253"], ican['cohort/ican']))
        g.add((variant_fq_ican_uri, sio["000300"], Literal(fq_ican_AF_whole, datatype=XSD.float)))
    return g
#
def calculate_fqgnomad(chromosome, position, reference, alternate):
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
def calculate_fqican(genotypes):
    alternateallele = genotypes.count('1')
    referenceallele = genotypes.count('0')
    if (alternateallele+referenceallele) == 0:
        return 'nan'
    else:
        alternate_fq = alternateallele/(alternateallele+referenceallele)
        return round(alternate_fq, 6)
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
