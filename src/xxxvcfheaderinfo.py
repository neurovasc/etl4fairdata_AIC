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
    '##INFO=<ID=AC_familial_hom,Number=1,Type=Integer,Description="Alt Allele Count in individuals have familial forms of ICA [NCIT:00000] with homozygous genotypes">\n'
    '##INFO=<ID=AC_familial_het,Number=1,Type=Integer,Description="Alt Allele Count in individuals have familial forms of ICA [NCIT:00000] with heterozygous genotypes">\n'
)
sporadiccase = (
    '##INFO=<ID=AF_sporadic,Number=A,Type=Float,Description="Alt Allele Frequency within individuals with no relatives carrying an IC aneurysm [NCIT:00000]">\n'
    '##INFO=<ID=AC_sporadic,Number=A,Type=Integer,Description="Alt Allele Count in individuals with no relatives carrying an IC aneurysm [NCIT:00000]">\n'
    '##INFO=<ID=AC_sporadic_hom,Number=1,Type=Integer,Description="Alt Allele Count in individuals have sporadic forms of ICA [NCIT:00000] with homozygous genotypes">\n'
    '##INFO=<ID=AC_sporadic_het,Number=1,Type=Integer,Description="Alt Allele Count in individuals have sporadic forms of ICA [NCIT:00000] with heterozygous genotypes">\n'
)
multipleica = (
    '##INFO=<ID=AF_multipleica,Number=A,Type=Float,Description="Alt Allele Frequency within individuals with multiple ica [NCIT:00000]">\n'
    '##INFO=<ID=AC_multipleica,Number=A,Type=Integer,Description="Alt Allele Count in individuals with multiple ica [NCIT:00000]">\n'
    '##INFO=<ID=AC_multipleica_hom,Number=1,Type=Integer,Description="Alt Allele Count in individuals with multiple ica [NCIT:00000] with homozygous genotypes">\n'
    '##INFO=<ID=AC_multipleica_het,Number=1,Type=Integer,Description="Alt Allele Count in individuals with multiple ica [NCIT:00000] with heterozygous genotypes">\n'
)
aht = (
    '##INFO=<ID=AF_aht,Number=A,Type=Float,Description="Alt Allele Frequency within individuals with arterial hypertension [NCIT:00000]">\n'
    '##INFO=<ID=AC_aht,Number=A,Type=Integer,Description="Alt Allele Count in individuals with arterial hypertension [NCIT:00000]">\n'
    '##INFO=<ID=AC_aht_hom,Number=1,Type=Integer,Description="Alt Allele Count in individuals have arterial hypertension [NCIT:00000] with homozygous genotypes">\n'
    '##INFO=<ID=AC_aht_het,Number=1,Type=Integer,Description="Alt Allele Count in individuals have arterial hypertension [NCIT:00000] with heterozygous genotypes">\n'    
)
diabetes =(
    '##INFO=<ID=AF_diabetes,Number=A,Type=Float,Description="Alt Allele Frequency within individuals with diabetes [NCIT:00000]">\n'
    '##INFO=<ID=AC_diabetes,Number=A,Type=Integer,Description="Alt Allele Count in individuals with diabetes [NCIT:00000]">\n'
    '##INFO=<ID=AC_diabetes_hom,Number=1,Type=Integer,Description="Alt Allele Count in individuals that have diabetes [NCIT:00000] with homozygous genotypes">\n'
    '##INFO=<ID=AC_diabetes_het,Number=1,Type=Integer,Description="Alt Allele Count in individuals that have diabetes [NCIT:00000] with heterozygous genotypes">\n'  
)
neversmoked = (
    '##INFO=<ID=AF_neversmoked,Number=A,Type=Float,Description="Alt Allele Frequency within individuals that never smoked [NCIT:00000]">\n'
    '##INFO=<ID=AC_neversmoked,Number=A,Type=Integer,Description="Alt Allele Count in individuals that never smoked [NCIT:00000]">\n'
    '##INFO=<ID=AC_neversmoked_hom,Number=1,Type=Integer,Description="Alt Allele Count in individuals that never smoked [NCIT:00000] with homozygous genotypes">\n'
    '##INFO=<ID=AC_neversmoked_het,Number=1,Type=Integer,Description="Alt Allele Count in individuals that never smoked [NCIT:00000] with heterozygous genotypes">\n'  

)
currentlysmoking = (
    '##INFO=<ID=AF_currentlysmoking,Number=A,Type=Float,Description="Alt Allele Frequency within individuals that are currently smoking [NCIT:00000]">\n'
    '##INFO=<ID=AC_currentlysmoking,Number=A,Type=Integer,Description="Alt Allele Count in individuals that are currently smoking [NCIT:00000]">\n'
    '##INFO=<ID=AC_currentlysmoking_hom,Number=1,Type=Integer,Description="Alt Allele Count in individuals that are currently smoking [NCIT:00000] with homozygous genotypes">\n'
    '##INFO=<ID=AC_currentlysmoking_het,Number=1,Type=Integer,Description="Alt Allele Count in individuals that are currently smoking [NCIT:00000] with heterozygous genotypes">\n'  
)
treatment = (
    '##INFO=<ID=AF_treatment,Number=A,Type=Float,Description="Alt Allele Frequency within individuals that had an ica treatment [NCIT:00000]">\n'
    '##INFO=<ID=AC_treatment,Number=A,Type=Integer,Description="Alt Allele Count in individuals that had an ica treatment [NCIT:00000]">\n'
    '##INFO=<ID=AC_treatment_hom,Number=1,Type=Integer,Description="Alt Allele Count in individuals that had an ica treatment [NCIT:00000] with homozygous genotypes">\n'
    '##INFO=<ID=AC_treatment_het,Number=1,Type=Integer,Description="Alt Allele Count in individuals that had an ica treatment [NCIT:00000] with heterozygous genotypes">\n'  
)
treatmentcoils = (
    '##INFO=<ID=AF_treatmentcoils,Number=A,Type=Float,Description="Alt Allele Frequency within individuals that had an ica treatment with coils [NCIT:00000]">\n'
    '##INFO=<ID=AC_treatmentcoils,Number=A,Type=Integer,Description="Alt Allele Count in individuals that had an ica treatment with coils [NCIT:00000]">\n'
    '##INFO=<ID=AC_treatmentcoils_hom,Number=1,Type=Integer,Description="Alt Allele Count in individuals that had an ica treatment with coils [NCIT:00000] with homozygous genotypes">\n'
    '##INFO=<ID=AC_treatmentcoils_het,Number=1,Type=Integer,Description="Alt Allele Count in individuals that had an ica treatment with coils [NCIT:00000] with heterozygous genotypes">\n'  
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
    'multipleica' : multipleica,
    'aht' : aht,
    'diabetes' : diabetes,
    'neversmoked' : neversmoked,
    'currentlysmoking' : currentlysmoking,
    'treatment' : treatment, 
    'treatmentcoils' : treatmentcoils
}
#