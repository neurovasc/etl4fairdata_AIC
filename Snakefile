import os
import glob

########################################################################################
# datasets
datasetdir = "datasets/"
datasets = {"test" : "test", 
            "synthetic1" : "synthetic1", 
            "aic" : "aic", 
            "syntheticican2" : "syntheticican2",
            "synican3" : "synican3"}
###########################
dataset = datasets["synican3"]# change this to the dataset you want to work with
###########################
input_dir = datasetdir + dataset + "/" +  dataset + "-data-input"
inputcsv = glob.glob(input_dir + "/*.csv")[0]
inputvcf = glob.glob(input_dir + "/*.vcf.gz")[0]
intermediate_dir = datasetdir + dataset + "/" +  dataset + "-data-intermediate"
deliverable_dir = datasetdir + dataset + "/" +  dataset + "-data-deliverable"


rule setup:
    run:
        os.makedirs(intermediate_dir, exist_ok=True)
        os.makedirs(deliverable_dir, exist_ok=True)
        print(f"Directory holding data: {input_dir}")
        print(f"Files in {input_dir}: ", os.listdir(input_dir))
        print(f"Setting up directories: {intermediate_dir} and  {deliverable_dir}")
        print(f"Input csv files: {inputcsv}")
        print(f"Input vcf files: {inputvcf}")
#
# extract biosample ids from genotype vcf file
rule sample_selection_vcfids:
    input:
        vcf = inputvcf
    output:
        list1 = intermediate_dir+"/biosample-ids-in-vcf.lst"
    shell:
        "bcftools query -l {input.vcf} | tr '_' '\n'> {output.list1}"
#
# 1) extract biosample ids from phenotype csv file using the list of biosample 
# ids from the genotype vcf file
# 2) write the output to a new phenotype csv file
# 3) write the list of biosample ids from the genotype vcf 
# file that also have a phenotype in the phenotype csv file (for rule sample_selection_withgeno)
rule sample_selection_withpheno:
    input:
        csv = inputcsv,
        list1 = intermediate_dir+"/biosample-ids-in-vcf.lst"
    output:
        csv2 = intermediate_dir + "/genopheno-pheno-sampled.csv",
        list2 = intermediate_dir+"/biosample-ids-in-vcf-with-pheno-in-csv.lst"
    shell:
        """ 
        echo python src/sampleselection.py \
            -s {input.list1} \
            -c {input.csv} \
            -o {output.list2} \
            -f {output.csv2};
        python src/sampleselection.py \
            -s {input.list1} \
            -c {input.csv} \
            -o {output.list2} \
            -f {output.csv2}
        """
#
# extract the biosamples in the genotype vcf file that also have a phenotype in the phenotype csv file
rule sample_selection_withgeno:
    input:
        vcf = inputvcf,
        list2 = intermediate_dir+"/biosample-ids-in-vcf-with-pheno-in-csv.lst"
    output:
        vcf = intermediate_dir + "/genopheno-geno-sampled.vcf.gz"
    shell:
        "bcftools view -S {input.list2} {input.vcf} -c1:alt1 -O b -o {output.vcf}"
#
# Reorder samples in phenotype file to match the order in the VCF file
# Delete duplicate rows in phenotype file
rule reorder_samples_in_phenotypecsv:
    input:
        csv = intermediate_dir + "/genopheno-pheno-sampled.csv",
        list2 = intermediate_dir+"/biosample-ids-in-vcf-with-pheno-in-csv.lst"
    output:
        csv = intermediate_dir + "/genopheno-pheno-sampled-reordered.csv",
    run:
        import pandas as pd
        import csv
        def reorder_csv(phenotypes, samples, output):
            df = pd.read_csv(phenotypes)
            samples_df = pd.read_csv(samples, header=None, names=['identifier'])
            samples_list = samples_df['identifier'].tolist()
            try:
                identifier_column = 'N°ADN IRT 1' 
            except:
                pass
            try:
                identifier_column = 'biosampleId' 
            except:
                pass
            sorted_df = df.set_index(identifier_column).loc[samples_list].reset_index()
            sorted_df = sorted_df.drop_duplicates(subset=[identifier_column])
            return sorted_df
        df = reorder_csv(input.csv, input.list2, output.csv)
        df.to_csv(output.csv, index=False, quoting=csv.QUOTE_ALL)
#
# Sanity check: there is the same amount of samples in the VCF and the phenotype file
# The samples are in the same order in both files.
# The same samples are present in both files.
rule sanity_check_samples:
    input:
        csv = intermediate_dir + "/genopheno-pheno-sampled-reordered.csv",
        vcf = intermediate_dir + "/genopheno-geno-sampled.vcf.gz"
    output:
        validation = intermediate_dir + "/sanity_check_samples_complete.txt"
    run:
        import pandas as pd
        import subprocess
        def check_samples(phenotypes, vcf, output_file):
            df = pd.read_csv(phenotypes, header=0)
            vcf_samples = subprocess.check_output("bcftools query -l " + vcf, shell=True)
            vcf_samples = vcf_samples.decode('utf-8').split("\n")
            vcf_samples.pop() # Remove the last empty element
            
            # Check if the number of samples is the same
            if len(df) != len(vcf_samples):
                raise ValueError("The number of samples in the phenotype file and the VCF file is different.")
            # Check if the samples are in the same order
            try:
                thesamples = df['N°ADN IRT 1']
            except:
                thesamples = df['biosampleId']
            if not all(thesamples == vcf_samples):
                raise ValueError("The samples are not in the same order in the phenotype file and the VCF file.")
            # Write a success message to the output file
            with open(output_file, "w") as f:
                f.write("OK : Sanity check passed successfully: same number of samples, in the same order.\n")
                f.write(f"{input.csv}: {len(df)} biosample ids \n")
                f.write(f"{input.vcf}: {len(vcf_samples)} biosample ids \n")

        check_samples(input.csv, input.vcf, output.validation)
#
# Extract minimal genotype information from the VCF file
# Before extraction: eliminate positions with no alternative alleles (--min-ac=1)
# Before extraction: reformat the bcf so that one single allele per line appears 
# BEWARE WARNING WARNING: It will do funky things with multiallelic sites
# Not a problem for primary (1) allele frequency calculation
# https://stackoverflow.com/questions/73611363/better-splitting-of-mutliallelic-sites-then-bcftools-norm-m-any
rule vcf2querygenotype:
    input:
        vcf = intermediate_dir + "/genopheno-geno-sampled.vcf.gz"
    output:
        query = intermediate_dir + "/genopheno-geno-sampled-querygenotype.tsv",
        contigs = intermediate_dir + "/genopheno-geno-sampled-contigs.txt",
        info = intermediate_dir + "/genopheno-geno-sampled-info.txt"
    shell:
        """
        bcftools norm -a {input.vcf} | bcftools view --min-ac=1 | \
        bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT,]\t%INFO\n' > {output.query} && \
        bcftools view -h {input.vcf} | grep '^##contig'  > {output.contigs} && \
        bcftools view -h {input.vcf} | grep '^##INFO'  > {output.info}
        """
#
# Build frequencies per variant 
# Each variant has genotype information, I slap it into the phenotype file, and calculate the 
# frequencies I want to use for the new VCF aggregate
# Frequencies should have meaningful names and definitions, preferentially using ontologies
rule computeallelefrequencies:
    input:
        query = intermediate_dir + "/genopheno-geno-sampled-querygenotype.tsv",
        clinical = intermediate_dir + "/genopheno-pheno-sampled-reordered.csv",
        sequences = intermediate_dir + "/genopheno-geno-sampled-contigs.txt",
        info = intermediate_dir + "/genopheno-geno-sampled-info.txt"
    output:
        aggregatevcf = deliverable_dir+"/"+datasets[dataset]+"-genotypes.aggregate.vcf.gz"
    shell:
        "python3 src/computeallelefrequencies.py -g {input.query} -c {input.clinical} -o {output.aggregatevcf} -s {input.sequences} -i {input.info}"
#
# Sanity check: is the bgziped VCF file valid?
rule sanity_check_vcfvalidity:
    input:
        vcf = deliverable_dir+"/"+datasets[dataset]+"-genotypes.aggregate.vcf.gz"
    output:
        validation = intermediate_dir + "/sanity_check_vcfvalidity_complete.txt"
    run:
        import subprocess
        def check_vcf_validity(vcf, output_file):
            try:
                subprocess.check_output("bcftools view " + vcf + " > /dev/null", shell=True)
                subprocess.check_output("vcf-validator " + vcf + " > /dev/null", shell=True)
            except subprocess.CalledProcessError as e:
                raise ValueError("The bgziped VCF file is not valid.")
            
            with open(output_file, "w") as f:
                f.write("OK : Sanity check passed successfully: aggregated vcf is a valid vcf.\n")

        check_vcf_validity(input.vcf, output.validation)

#
# Build RDF from the VCF aggregate
# This step uses a RDF schema for variant information
# inspired by faldo, disgenet, swat4hls paper, etc
# This script uses gnomad file and dbsnp file, by default the paths are those
# on my PC pp-irs1-4071ylt but you can change them with -g and -d parameters
rule vcfaggregate2rdf:
    input:
        vcf = deliverable_dir+"/"+datasets[dataset]+"-genotypes.aggregate.vcf.gz"
    params:
        limit = 100
    output:
        rdf=deliverable_dir+"/"+datasets[dataset]+"-genotypes.aggregate.ttl"
    shell:
        "python3 src/vcfaggregate2rdf_v2.py -v {input.vcf} -r {output.rdf}"

#snakemake --rulegraph | dot -Tpdf > dag.pdf