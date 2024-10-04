# Header for the VCF file
#
fileformat = ("##fileformat=VCFv4.2\n")
reference = ('##reference=file:///LAB-DATA/BiRD/resources/species/human/cng.fr/hs38me/dragen_cnv//reference.bin\n')
columns = (
   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
) # no FORMAT column needed, because no individual level data
whole = (
    '##INFO=<ID=AF_whole,Number=A,Type=Float,Description="Allele Frequency in AIC population">\n'
    '##INFO=<ID=AC_whole,Number=A,Type=Integer,Description="Allele Count in AIC population">\n'
    )
male = (
    '##INFO=<ID=AF_male,Number=A,Type=Float,Description="Allele Frequency of AIC males">\n'
    '##INFO=<ID=AC_male,Number=A,Type=Integer,Description="Allele Count of AIC males">\n'
)
female = (
    '##INFO=<ID=AF_female,Number=A,Type=Float,Description="Allele Frequency of AIC females">\n'
    '##INFO=<ID=AC_female,Number=A,Type=Integer,Description="Allele Count of AIC females">\n'
)
familialcase = (
    '##INFO=<ID=AF_familial,Number=A,Type=Float,Description="Allele Frequency of AIC familial cases">\n'
    '##INFO=<ID=AC_familial,Number=A,Type=Integer,Description="Allele Count of AIC familail cases">\n'
)
uncertaincase = (
    '##INFO=<ID=AF_uncertain,Number=A,Type=Float,Description="Allele Frequency of AIC uncertain cases, neither familial nor sporadic">\n'
    '##INFO=<ID=AC_uncertain,Number=A,Type=Integer,Description="Allele Count of AIC uncertain cases, neither familial nor sporadic">\n'
)
sporadiccase = (
    '##INFO=<ID=AF_sporadic,Number=A,Type=Float,Description="Allele Frequency of AIC sporadic cases">\n'
    '##INFO=<ID=AC_sporadic,Number=A,Type=Integer,Description="Allele Count of AIC sporadic cases">\n'
)
earlyonset = (
    '##INFO=<ID=AF_earlyonset,Number=A,Type=Float,Description="Allele Frequency of AIC early onset cases">\n'
    '##INFO=<ID=AC_earlyonset,Number=A,Type=Integer,Description="Allele Count of AIC early onset cases">\n'
)
lateonset = (
    '##INFO=<ID=AF_lateonset,Number=A,Type=Float,Description="Allele Frequency of AIC late onset cases">\n'
    '##INFO=<ID=AC_lateonset,Number=A,Type=Integer,Description="Allele Count of AIC late onset cases">\n'
)
info_headerchunk = { 
    'whole' : whole,
    'male' : male, 
    'female' : female,
    'familialcase' : familialcase,
    'uncertaincase' : uncertaincase,
    'sporadiccase' : sporadiccase, 
    'earlyonset' : earlyonset, 
    'lateonset' : lateonset 
}
#