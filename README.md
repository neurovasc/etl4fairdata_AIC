# etl4fairdata_AIC
This repository contains the code for the etl4fairdata_AIC snakemake pipeline. 
The pipeline extracts and transforms data from different sources and prepares it to be loaded in a FAIR genomic ecosystem.
The pipeline needs:
- a csv containing clinical/phenotypical data, extracted from medical databases for example
- a VCF file containing individual genotypes
The pipeline main products are:
- an aggregated vcf file where individual genotypes are aggregated into frequencies and counts, relative to different clinical/phenotypical traits
- a ttl file representing the vcf according to a data schema

# Enviroment
 conda export --from-history > ~/Devlopment/etl4fairdata_AIC/conda/environement-etl4fair.yml
 https://www.youtube.com/watch?v=r9PWnEmz_tc
 snakemake --dag targets | dot -Tpng > dag.png
 python scripts starting with xxx are not meant to be executed, but are ressources for other scripts because I like fragmenting my code a lot.

# TODOs
- For variants that do not have a gnomad frequency, should the frequency be set to 0 
in the rdf file or simply not included? Now it appears as "NaN"^^xsd:float
[x 13092024] Add logger and verbose argument for vcfaggregate3rdf.py 
- Add custom check function that prevents the usage of both --limit and (--chunk or --threads)
- Add custom adaptation of chunksize according to the number of variants in the vcf file
[DONE, december 24] - Maybe add function in vcfaggregate3rdf.py that checks if the number of variants
is the same in the vcf and in the final rdf -> it is also a sanity check in the snakemake pipeline
- Issue with clean up function -> I want to keep the intermediate merged file, but it is deleted. why?
- Add tests
[x 13092024] - Multithread computeallelefrequency.py
[DONE, november 24] - Lower the complexity of computeallelefrequency.py : loop over genotyped df as little
as possible
- Build rdf in vcfaggregate2rdf in a prettier way: now it is spaghetti code with rdflib
[DONE, january 24] - Add homozygous and heterozygous counts
- create a sturdier test dataset

# What happens after the Snakemake pipeline ?
Setting up a fuseki server for sparql queries
>fuseki-server --file=/home/bodrug-a/Devlopment/etl4fairdata_AIC/data-deliverable/aicdataset-QCed.VEP.AFctrls.GND.CADD.aggregate.ttl /ican
Please visit: https://10-54-1-83.gcp.glicid.fr/

# About input files 
QCed.VEP.AFctrls.GND.CADD.bcf 
- has positions with genotypes that are all reference like
- has samples that are not in the phenotype file
- has samples that are a merge of two samples in the phenotype file
signaled like sample1_sample2. These need to be eliminated as this is
contamination.
- Quality filtered variants
- VEP annotation
- CADD scores
- exomes
extraction_GAIA_ICAN_26-09-2023.csv

 # Authors & Licence
 Alexandrina Bodrug while at Institut du Thorax
[MIT LICENCE](LICENCE)