import os
# # # #
#
#
# Extract biosample IDs from the VCF file
# Some biosamples in the VCF are a merge of two biosamples in the phenotype file
rule get_biosamples_in_vcf:
    input:
        bcf="data-input/QCed.VEP.AFctrls.GND.CADD.bcf"
    output:
        list="data-intermediate/biosamplesinvcf.lst"
    shell:
        "bcftools query -l {input} | tr '_' '\n'> {output}"
# # # #
#
#
# Extract biosample IDs from the phenotype file
# That have a genomic analysis (ID present in VCF)
# Warning: some biosamples in the phenotype file are merged in the VCF
# The merged samples are not ADN 1 and ADN 2, but samples from different rows
rule sample_selection:
    input:
        code="src/sampleselection.py",
        csv="data-input/extraction_GAIA_ICAN_26-09-2023.csv",
        vcfsamples="data-intermediate/biosamplesinvcf.lst"
    output:
        csvsamples="data-intermediate/aicdataset-samplelist.lst",
        csvfilt="data-intermediate/aicdataset-extraction_GAIA_ICAN_26-09-2023.csv"
    shell: 
        "python3 {input.code} -c {input.csv} -s {input.vcfsamples} -o {output.csvsamples} -f {output.csvfilt}"
# # # #
#
#
# Extract samples with phenotype information and AIC or related from the VCF file
# Samples are ordered in the order of the file aicdataset-samplelist.lst
rule vcf_sample_selection:
    input:
        samples="data-intermediate/aicdataset-samplelist.lst",
        bcf="data-input/QCed.VEP.AFctrls.GND.CADD.bcf"
    output:
        bcf="data-intermediate/aicdataset-QCed.VEP.AFctrls.GND.CADD.bcf"
    shell:
        "bcftools view -S {input.samples} {input.bcf} -O b -o {output.bcf}"
# # # #
#
#
# Reorder samples in phenotype file to match the order in the VCF file
# Delete duplicate rows in phenotype file
rule reorder_samples_in_phenotypefile:
    input:
        csv="data-intermediate/aicdataset-extraction_GAIA_ICAN_26-09-2023.csv",
        samples="data-intermediate/aicdataset-samplelist.lst"
    output:
        csv="data-intermediate/aicdataset-extraction_GAIA_ICAN_26-09-2023.reordered.csv"
    run:
        import pandas as pd
        def reorder_csv(phenotypes, samples, output):
            # Read the CSV file
            df = pd.read_csv(phenotypes)
            
            # Read the CSV samples list
            samples_df = pd.read_csv(samples, header=None, names=['identifier'])
            samples_list = samples_df['identifier'].tolist()

            # Identifiers
            identifier_column = 'N°ADN IRT 1' 
            
            # Sort phenotype file according to the order of identifiers in samples file
            sorted_df = df.set_index(identifier_column).loc[samples_list].reset_index()
            sorted_df = sorted_df.drop_duplicates(subset=[identifier_column])

            return sorted_df

        df = reorder_csv(input.csv, input.samples, output.csv)
        df.to_csv(output.csv, index=False)
# # # #
#
#
# Sanity check: there is the same amount of samples in the VCF and the phenotype file
# The samples are in the same order in both files.
# The same samples are present in both files.
rule sanity_check_samples:
    input:
        csv="data-intermediate/aicdataset-extraction_GAIA_ICAN_26-09-2023.reordered.csv",
        bcf="data-intermediate/aicdataset-QCed.VEP.AFctrls.GND.CADD.bcf"
    run:
        import pandas as pd
        import subprocess
        def check_samples(phenotypes, bcf):
            # Read the CSV file
            df = pd.read_csv(phenotypes, header=0)
            #print("phenotype df len", len(df))
            # Read the VCF file
            bcf_samples = subprocess.check_output("bcftools query -l " + bcf, shell=True)
            bcf_samples = bcf_samples.decode('utf-8').split("\n")
            bcf_samples.pop()
            #print("bcf samples len:", len(bcf_samples))
            #print("bcf samples:", bcf_samples)
            
            # Check if the number of samples is the same
            if len(df) != len(bcf_samples):
                raise ValueError("The number of samples in the phenotype file and the VCF file is different.")
            # Check if the samples are in the same order
            if not all(df['N°ADN IRT 1'] == bcf_samples):
                raise ValueError("The samples are not in the same order in the phenotype file and the VCF file.")
            return True

        check_samples(input.csv, input.bcf)
# # # #
#
#
# Extract minimal genotype information from the VCF file
# Before extraction: eliminate positions with no alternative alleles (--min-ac=1)
# Before extraction: reformat the bcf so that one single allele per line appears 
# BEWARE WARNING WARNING: It will do funky things with multiallelic sites
# Not a problem for primary (1) allele frequency calculation
# https://stackoverflow.com/questions/73611363/better-splitting-of-mutliallelic-sites-then-bcftools-norm-m-any
rule vcf2querygenotype:
    input:
        bcf="data-intermediate/aicdataset-QCed.VEP.AFctrls.GND.CADD.bcf"
    output:
        query="data-intermediate/aicdataset-querygenotype.tsv"
    shell:
        "bcftools norm -a {input.bcf} | bcftools view --min-ac=1 | \
        bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT,]\n' > {output.query}"
# # # #
#
#
# Build frequencies per variant 
# Each variant has genotype information, I slap it into the phenotype file, and calculate the 
# frequencies I want to use for the new VCF aggregate
# Frequencies should have meaningful names and definitions, preferentially using ontologies
rule computeallelefrequencies:
    input:
        code="src/computeallelefrequencies.py",
        query="data-intermediate/aicdataset-querygenotype.tsv",
        csv="data-intermediate/aicdataset-extraction_GAIA_ICAN_26-09-2023.csv"
    output:
        gtfreq="data-intermediate/aicdataset-genotypefrequencies.tsv"
    shell:
        "python3 {input.code} -q {input.query} -c {input.csv} -o {output.gtfreq}"


