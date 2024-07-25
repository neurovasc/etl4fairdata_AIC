# Extract biosample IDs from the VCF file
# Some biosamples in the VCF are a merge of two biosamples in the phenotype file
rule get_biosamples_in_vcf:
    input:
        bcf="data-input/QCed.VEP.AFctrls.GND.CADD.bcf"
    output:
        list="data-intermediate/biosamplesinvcf.lst"
    shell:
        "bcftools query -l {input} | tr '_' '\n'> {output}"
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
# Extract samples with phenotype information and AIC or related from the VCF file
rule vcf_sample_selection:
    input:
        samples="data-intermediate/aicdataset-samplelist.lst",
        bcf="data-input/QCed.VEP.AFctrls.GND.CADD.bcf"
    output:
        bcf="data-intermediate/aicdataset-QCed.VEP.AFctrls.GND.CADD.bcf"
    shell:
        "bcftools view -S {input.samples} {input.bcf} -O b -o {output.bcf}"
# Extract minimal genotype information from the VCF file
rule vcf2querygenotype:
    input:
        bcf="data-intermediate/aicdataset-QCed.VEP.AFctrls.GND.CADD.bcf"
    output:
        query="data-intermediate/aicdataset-querygenotype.tsv"
    shell:
        "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' {input.bcf} > {output.query}"


