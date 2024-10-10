import random

# Function to generate a random genotype
def generate_genotype():
    return random.choice(['0/0', '0/1', '1/1'])

# File paths
input_vcf = "fake.vcf"  # Path to the input VCF file
output_vcf = "fakegenotyped.vcf"  # Path to save the output VCF file

# Open the input VCF file
with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
    for line in infile:
        # If the line starts with '##', it's part of the header, so copy it as is
        if line.startswith('##'):
            outfile.write(line)
        # If the line starts with '#CHROM', it's the column header, add FORMAT and sample names
        elif line.startswith('#CHROM'):
            # Add "FORMAT" column (for GT) and maintain the existing samples
            outfile.write(line.strip() + "\tFORMAT" + "\t" + "\t".join([f"SAMPLE{i}" for i in range(1, 95)]) + "\n")
        else:
            # For each data line, we assume the VCF has 8 fields before genotypes
            fields = line.strip().split("\t")
            
            # We need to add the genotype column
            format_field = "GT"
            
            # Generate random genotypes for 94 samples
            genotypes = "\t".join([generate_genotype() for _ in range(94)])
            
            # Write the line with added genotype information
            outfile.write("\t".join(fields[:8]) + f"\t{format_field}\t{genotypes}\n")
