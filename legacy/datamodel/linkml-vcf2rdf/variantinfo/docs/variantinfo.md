
# variants-kg


**metamodel version:** 1.7.0

**version:** None


Schema for describing sequence alterations (SO:0001059) and their relationships to other elements of a variant catalog format.
The ontology aims at describing genomic variations found in patients carrying intracranial aneurysms.
The ontology and schema is organized arount variants, their annotations relevant to the study of 
intracranial aneurysms, and phenotype frequencies and counts relevant to the study of intracranial aneurysms.


### Classes

 * [AlternateAllele](AlternateAllele.md) - Represents the alternate allele (geno:0000002).
 * [AssociatedCharacteristic](AssociatedCharacteristic.md) - Represents a phenotype, clinical trait or grouping associated with a sequence alteration.
 * [Container](Container.md)
 * [Count](Count.md) - Represents the count of alternative alleles in individuals with a certain zygosity.
 * [Frequency](Frequency.md) - Represents the frequency of an allele in individuals with a certain zygosity.
 * [ReferenceAllele](ReferenceAllele.md) - Represents the reference allele (geno:0000036).
 * [SequenceAlteration](SequenceAlteration.md) - A representation of a sequence alteration (so:0001059).
 * [VariationSite](VariationSite.md) - Represents the location of a sequence alteration.
 * [VariationSiteBegin](VariationSiteBegin.md) - Represents the beginning of the location of a sequence alteration.
 * [VariationSiteEnd](VariationSiteEnd.md) - Represents the end of the location of a sequence alteration.
 * [VariationSiteReference](VariationSiteReference.md) - Represents the reference sequence (contig, sequence, chromosome).
 * [Zygosity](Zygosity.md) - Represents the zygosity of an associated characteristic.

### Mixins


### Slots

 * [➞value](alternateAllele__value.md) - The value of the alternate allele.
 * [➞has_identifier](associatedCharacteristic__has_identifier.md) - A unique identifier for the associated characteristic.
 * [➞has_zygosity](associatedCharacteristic__has_zygosity.md) - The zygocity of the associated characteristic.
 * [➞label](associatedCharacteristic__label.md) - A human-readable label for the associated characteristic.
 * [➞variantcatalog](container__variantcatalog.md)
 * [➞value](count__value.md) - The count value.
 * [➞value](frequency__value.md) - The frequency value.
 * [➞value](referenceAllele__value.md) - The value of the reference allele.
 * [➞has_alternate_allele](sequenceAlteration__has_alternate_allele.md) - Links the sequence alteration to its alternate allele (geno:0000002) using the property geno:0000382.
 * [➞has_identifier](sequenceAlteration__has_identifier.md) - A unique identifier for the sequence alteration.
 * [➞has_location](sequenceAlteration__has_location.md) - Links the sequence alteration to its location using the property faldo:Region.
 * [➞has_reference_allele](sequenceAlteration__has_reference_allele.md) - Links the sequence alteration to its reference allele (geno:0000036) using the property geno:0000385.
 * [➞is_associated_with](sequenceAlteration__is_associated_with.md) - Links the sequence alteration to a phenotype (female), clinical trait (diabetes) or grouping (whole cohort).
 * [➞position](variationSiteBegin__position.md) - The position of the beginning of the location of the sequence alteration.
 * [➞position](variationSiteEnd__position.md) - The position of the end of the location of the sequence alteration.
 * [➞label](variationSiteReference__label.md) - A human-readable label for the reference sequence.
 * [➞sameAs](variationSiteReference__sameAs.md) - The value of the reference sequence in another database (ncbi sequence i.e. https://www.ncbi.nlm.nih.gov/nuccore/NC_000022.11).
 * [➞value](variationSiteReference__value.md) - The value of the reference sequence (ena sequence i.e. https://www.ebi.ac.uk/ena/browser/view/CM000684.2).
 * [➞begins_at](variationSite__begins_at.md) - The beginning of the location of the sequence alteration.
 * [➞ends_at](variationSite__ends_at.md) - The end of the location of the sequence alteration.
 * [➞has_reference](variationSite__has_reference.md) - The reference sequence (contig, sequence, chromosome).
 * [➞has_frequency](zygosity__has_frequency.md) - Frequency of the allele in indibiduals with a certain zygosity.
 * [➞has_measurement_value](zygosity__has_measurement_value.md) - Count of alternative alleles in individuals with a certain zygosity.
 * [➞type](zygosity__type.md) - The type of zygosity.

### Enums

 * [ZygosityType](ZygosityType.md) - Enumeration of zygosity types.

### Subsets


### Types


#### Built in

 * **Bool**
 * **Curie**
 * **Decimal**
 * **ElementIdentifier**
 * **NCName**
 * **NodeIdentifier**
 * **URI**
 * **URIorCURIE**
 * **XSDDate**
 * **XSDDateTime**
 * **XSDTime**
 * **float**
 * **int**
 * **str**

#### Defined

 * [Boolean](types/Boolean.md)  (**Bool**)  - A binary (true or false) value
 * [Curie](types/Curie.md)  (**Curie**)  - a compact URI
 * [Date](types/Date.md)  (**XSDDate**)  - a date (year, month and day) in an idealized calendar
 * [DateOrDatetime](types/DateOrDatetime.md)  (**str**)  - Either a date or a datetime
 * [Datetime](types/Datetime.md)  (**XSDDateTime**)  - The combination of a date and time
 * [Decimal](types/Decimal.md)  (**Decimal**)  - A real number with arbitrary precision that conforms to the xsd:decimal specification
 * [Double](types/Double.md)  (**float**)  - A real number that conforms to the xsd:double specification
 * [Float](types/Float.md)  (**float**)  - A real number that conforms to the xsd:float specification
 * [Integer](types/Integer.md)  (**int**)  - An integer
 * [Jsonpath](types/Jsonpath.md)  (**str**)  - A string encoding a JSON Path. The value of the string MUST conform to JSON Point syntax and SHOULD dereference to zero or more valid objects within the current instance document when encoded in tree form.
 * [Jsonpointer](types/Jsonpointer.md)  (**str**)  - A string encoding a JSON Pointer. The value of the string MUST conform to JSON Point syntax and SHOULD dereference to a valid object within the current instance document when encoded in tree form.
 * [Ncname](types/Ncname.md)  (**NCName**)  - Prefix part of CURIE
 * [Nodeidentifier](types/Nodeidentifier.md)  (**NodeIdentifier**)  - A URI, CURIE or BNODE that represents a node in a model.
 * [Objectidentifier](types/Objectidentifier.md)  (**ElementIdentifier**)  - A URI or CURIE that represents an object in the model.
 * [Sparqlpath](types/Sparqlpath.md)  (**str**)  - A string encoding a SPARQL Property Path. The value of the string MUST conform to SPARQL syntax and SHOULD dereference to zero or more valid objects within the current instance document when encoded as RDF.
 * [String](types/String.md)  (**str**)  - A character string
 * [Time](types/Time.md)  (**XSDTime**)  - A time object represents a (local) time of day, independent of any particular day
 * [Uri](types/Uri.md)  (**URI**)  - a complete URI
 * [Uriorcurie](types/Uriorcurie.md)  (**URIorCURIE**)  - a URI or a CURIE
