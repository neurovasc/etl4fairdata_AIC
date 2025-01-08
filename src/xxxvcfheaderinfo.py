# Header for the VCF file
#
fileformat = ("##fileformat=VCFv4.2\n")
reference = ('##reference=file:///LAB-DATA/BiRD/resources/species/human/cng.fr/hs38me/dragen_cnv//reference.bin\n')
columns = (
   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
) # no FORMAT column needed, because no individual level data
whole = (
    '##INFO=<ID=AF_whole,Number=A,Type=Float,Description="Alt Allele Frequency in the cohort (ICAN)">\n'
    '##INFO=<ID=AC_whole,Number=A,Type=Integer,Description="Alt Allele Count in the cohort (ICAN)">\n'
    '##INFO=<ID=AC_whole_hom,Number=1,Type=Integer,Description="Alt Allele Count in the cohort (ICAN) in homozygous genotypes">\n'
    '##INFO=<ID=AC_whole_het,Number=1,Type=Integer,Description="Alt Allele Count in the cohort (ICAN) in heterozygous genotypes">\n'
    )
male = (
    '##INFO=<ID=AF_male,Number=A,Type=Float,Description="Alt Allele Frequency within individuals with male sex [NCIT:00000]">\n'
    '##INFO=<ID=AC_male,Number=A,Type=Integer,Description="Alt Allele Count in individuals with male sex [NCIT:00000]">\n'
    '##INFO=<ID=AC_male_hom,Number=1,Type=Integer,Description="Alt Allele Count in individuals with male sex [NCIT:00000] with homozygous genotypes">\n'
    '##INFO=<ID=AC_male_het,Number=1,Type=Integer,Description="Alt Allele Count in individuals with male sex [NCIT:00000] with heterozygous genotypes">\n'
    )
female = (
    '##INFO=<ID=AF_female,Number=A,Type=Float,Description="Alt Allele Frequency within individuals with female sex [NCIT:00000]">\n'
    '##INFO=<ID=AC_female,Number=A,Type=Integer,Description="Alt Allele Count in individuals with female sex [NCIT:00000]">\n'
    '##INFO=<ID=AC_female_hom,Number=1,Type=Integer,Description="Alt Allele Count in individuals with female sex [NCIT:00000] with homozygous genotypes">\n'
    '##INFO=<ID=AC_female_het,Number=1,Type=Integer,Description="Alt Allele Count in individuals with female sex [NCIT:00000] with heterozygous genotypes">\n'
    )

earlyonset = (
    '##INFO=<ID=AF_earlyonset,Number=A,Type=Float,Description="Alt Allele Frequency within individuals with early onset [NCIT:00000]">\n'
    '##INFO=<ID=AC_earlyonset,Number=A,Type=Integer,Description="Alt Allele Count in individuals with early onset [NCIT:00000]">\n'
    '##INFO=<ID=AC_earlyonset_hom,Number=1,Type=Integer,Description="Alt Allele Count in individuals with early onset [NCIT:00000] with homozygous genotypes">\n'
    '##INFO=<ID=AC_earlyonset_het,Number=1,Type=Integer,Description="Alt Allele Count in individuals with early onset [NCIT:00000] with heterozygous genotypes">\n'
)
obese = (
    '##INFO=<ID=AF_obese,Number=A,Type=Float,Description="Alt Allele Frequency within individuals with obesity class I & II & III [NCIT:00000]">\n'
    '##INFO=<ID=AC_obese,Number=A,Type=Integer,Description="Alt Allele Count in individuals with obesity class I & II & III [NCIT:00000]">\n'
    '##INFO=<ID=AC_obese_hom,Number=1,Type=Integer,Description="Alt Allele Count in individuals with obesity class I & II & III [NCIT:00000] with homozygous genotypes">\n'
    '##INFO=<ID=AC_obese_het,Number=1,Type=Integer,Description="Alt Allele Count in individuals with obesity class I & II & III [NCIT:00000] with heterozygous genotypes">\n'
)
overweight = (
    '##INFO=<ID=AF_overweight,Number=A,Type=Float,Description="Alt Allele Frequency within individuals that are overweight [NCIT:00000]">\n'
    '##INFO=<ID=AC_overweight,Number=A,Type=Integer,Description="Alt Allele Count in individuals that are overweight [NCIT:00000]">\n'
    '##INFO=<ID=AC_overweight_hom,Number=1,Type=Integer,Description="Alt Allele Count in individuals that are overweight [NCIT:00000] with homozygous genotypes">\n'
    '##INFO=<ID=AC_overweight_het,Number=1,Type=Integer,Description="Alt Allele Count in individuals that are overweight [NCIT:00000] with heterozygous genotypes">\n'
)
familialcase = (
    '##INFO=<ID=AF_familial,Number=A,Type=Float,Description="Alt Allele Frequency within individuals with at least one 1st degree relative carrying an IC aneurysm [NCIT:00000]">\n'
    '##INFO=<ID=AC_familial,Number=A,Type=Integer,Description="Alt Allele Count in individuals with at least one 1st degree relative carrying an IC aneurysm [NCIT:00000]">\n'
)
sporadiccase = (
    '##INFO=<ID=AF_sporadic,Number=A,Type=Float,Description="Alt Allele Frequency within individuals with no relatives carrying an IC aneurysm [NCIT:00000]">\n'
    '##INFO=<ID=AC_sporadic,Number=A,Type=Integer,Description="Alt Allele Count in individuals with no relatives carrying an IC aneurysm [NCIT:00000]">\n'
)
discoveryincidental = (
    '##INFO=<ID=AF_discoveryincidental,Number=A,Type=Float,Description="Alt Allele Frequency within individuals with an incidental finding of an IC aneurysm [NCIT:00000]">\n'
    '##INFO=<ID=AC_discoveryincidental,Number=A,Type=Integer,Description="Alt Allele Count in individuals with an incidental finding of an IC aneurysm [NCIT:00000]">\n'
)
discoveryfamilialscreening = (
    '##INFO=<ID=AF_discoveryfamilialscreening,Number=A,Type=Float,Description="Alt Allele Frequency within individuals whose IC aneurysms were found through familial screening[NCIT:00000]">\n'
    '##INFO=<ID=AC_discoveryfamilialscreening,Number=A,Type=Integer,Description="Alt Allele Count in individuals whose IC aneurysms were found through familial screening [NCIT:00000]">\n'
)
discoveryruptured = (
    '##INFO=<ID=AF_discoveryruptured,Number=A,Type=Float,Description="Alt Allele Frequency within individuals whose IC aneurysms were found as a consequence of its rupture [NCIT:00000]">\n'
    '##INFO=<ID=AC_discoveryruptured,Number=A,Type=Integer,Description="Alt Allele Count in individuals whose IC aneurysms were found as a consequence of its rupture [NCIT:00000]">\n'
) 
discoveryischemic = (
    '##INFO=<ID=AF_discoveryischemic,Number=A,Type=Float,Description="Alt Allele Frequency within individuals whose IC aneurysms were found as a consequence of compressif or ischemic aneurysm symptoms (headaches, diplopia) [NCIT:00000]">\n'
    '##INFO=<ID=AC_discoveryischemic,Number=A,Type=Integer,Description="Alt Allele Count in individuals whose IC aneurysms was found as a consequence of compressif or ischemic aneurysm symptoms (headaches, diplopia) [NCIT:00000]">\n'
)
ruptured =(
    '##INFO=<ID=AF_ruptured,Number=A,Type=Float,Description="Alt Allele Frequency within individuals with a ruptured IC aneurysm [NCIT:00000]">\n'
    '##INFO=<ID=AC_ruptured,Number=A,Type=Integer,Description="Alt Allele Count in individuals with a ruptured IC aneurysm [NCIT:00000]">\n'
)
multipleica = (
    '##INFO=<ID=AF_multipleica,Number=A,Type=Float,Description="Alt Allele Frequency within individuals with multiple IC aneurysms [NCIT:00000]">\n'
    '##INFO=<ID=AC_multipleica,Number=A,Type=Integer,Description="Alt Allele Count in individuals with multiple IC aneurysms [NCIT:00000]">\n'
)
aht = (
    '##INFO=<ID=AF_aht,Number=A,Type=Float,Description="Alt Allele Frequency within individuals with arterial hypertension [NCIT:00000]">\n'
    '##INFO=<ID=AC_aht,Number=A,Type=Integer,Description="Alt Allele Count in individuals with arterial hypertension [NCIT:00000]">\n'
)
diabetes =(
    '##INFO=<ID=AF_diabetes,Number=A,Type=Float,Description="Alt Allele Frequency within individuals with diabetes [NCIT:00000]">\n'
    '##INFO=<ID=AC_diabetes,Number=A,Type=Integer,Description="Alt Allele Count in individuals with diabetes [NCIT:00000]">\n'
)
dyslipidemia = (
    '##INFO=<ID=AF_dyslipidemia,Number=A,Type=Float,Description="Alt Allele Frequency within individuals with dyslipidemia [NCIT:00000]">\n'
    '##INFO=<ID=AC_dyslipidemia,Number=A,Type=Integer,Description="Alt Allele Count in individuals with dyslipidemia [NCIT:00000]">\n'
)
#
info_headerchunk = { 
    'whole' : whole,
    'male' : male, 
    'female' : female,
    'earlyonset' : earlyonset,
    'obese' : obese,
    'overweight' : overweight,
    'familialcase' : familialcase,
    'sporadiccase' : sporadiccase,
    'discoveryincidental' : discoveryincidental,
    'discoveryfamilialscreening' : discoveryfamilialscreening,
    'discoveryischemic' : discoveryischemic,
    'discoveryruptured' : discoveryruptured,
    'ruptured' : ruptured,
    'multipleica' : multipleica,
    'aht' : aht,
    'diabetes' : diabetes,
    'dyslipidemia' : dyslipidemia
}
#