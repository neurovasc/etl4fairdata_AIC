# variants-kg

Schema for describing sequence alterations (SO:0001059) and their relationships to other elements of a variant catalog format.
The ontology aims at describing genomic variations found in patients carrying intracranial aneurysms.
The ontology and schema is organized arount variants, their annotations relevant to the study of 
intracranial aneurysms, and phenotype frequencies and counts relevant to the study of intracranial aneurysms. 

URI: https://ican.univ-nantes.io/variants-kg

Name: variants-kg



## Classes

| Class | Description |
| --- | --- |
| [AlternateAllele](AlternateAllele.md) | Represents the alternate allele (geno:0000002). |
| [AssociatedCharacteristic](AssociatedCharacteristic.md) | Represents a phenotype, clinical trait or grouping associated with a sequence alteration. |
| [Container](Container.md) | None |
| [Count](Count.md) | Represents the count of alternative alleles in individuals with a certain zygosity. |
| [Frequency](Frequency.md) | Represents the frequency of an allele in individuals with a certain zygosity. |
| [ReferenceAllele](ReferenceAllele.md) | Represents the reference allele (geno:0000036). |
| [SequenceAlteration](SequenceAlteration.md) | A representation of a sequence alteration (so:0001059). |
| [VariationSite](VariationSite.md) | Represents the location of a sequence alteration. |
| [VariationSiteBegin](VariationSiteBegin.md) | Represents the beginning of the location of a sequence alteration. |
| [VariationSiteEnd](VariationSiteEnd.md) | Represents the end of the location of a sequence alteration. |
| [VariationSiteReference](VariationSiteReference.md) | Represents the reference sequence (contig, sequence, chromosome). |
| [Zygosity](Zygosity.md) | Represents the zygosity of an associated characteristic. |



## Slots

| Slot | Description |
| --- | --- |
| [begins_at](begins_at.md) | The beginning of the location of the sequence alteration |
| [ends_at](ends_at.md) | The end of the location of the sequence alteration |
| [has_alternate_allele](has_alternate_allele.md) | Links the sequence alteration to its alternate allele (geno:0000002) using th... |
| [has_frequency](has_frequency.md) | Frequency of the allele in indibiduals with a certain zygosity |
| [has_identifier](has_identifier.md) | A unique identifier for the sequence alteration |
| [has_location](has_location.md) | Links the sequence alteration to its location using the property faldo:Region |
| [has_measurement_value](has_measurement_value.md) | Count of alternative alleles in individuals with a certain zygosity |
| [has_reference](has_reference.md) | The reference sequence (contig, sequence, chromosome) |
| [has_reference_allele](has_reference_allele.md) | Links the sequence alteration to its reference allele (geno:0000036) using th... |
| [has_zygosity](has_zygosity.md) | The zygocity of the associated characteristic |
| [is_associated_with](is_associated_with.md) | Links the sequence alteration to a phenotype (female), clinical trait (diabet... |
| [label](label.md) | A human-readable label for the reference sequence |
| [position](position.md) | The position of the beginning of the location of the sequence alteration |
| [sameAs](sameAs.md) | The value of the reference sequence in another database (ncbi sequence i |
| [type](type.md) | The type of zygosity |
| [value](value.md) | The value of the reference allele |
| [variantcatalog](variantcatalog.md) |  |


## Enumerations

| Enumeration | Description |
| --- | --- |
| [ZygosityType](ZygosityType.md) | Enumeration of zygosity types |


## Types

| Type | Description |
| --- | --- |
| [Boolean](Boolean.md) | A binary (true or false) value |
| [Curie](Curie.md) | a compact URI |
| [Date](Date.md) | a date (year, month and day) in an idealized calendar |
| [DateOrDatetime](DateOrDatetime.md) | Either a date or a datetime |
| [Datetime](Datetime.md) | The combination of a date and time |
| [Decimal](Decimal.md) | A real number with arbitrary precision that conforms to the xsd:decimal speci... |
| [Double](Double.md) | A real number that conforms to the xsd:double specification |
| [Float](Float.md) | A real number that conforms to the xsd:float specification |
| [Integer](Integer.md) | An integer |
| [Jsonpath](Jsonpath.md) | A string encoding a JSON Path |
| [Jsonpointer](Jsonpointer.md) | A string encoding a JSON Pointer |
| [Ncname](Ncname.md) | Prefix part of CURIE |
| [Nodeidentifier](Nodeidentifier.md) | A URI, CURIE or BNODE that represents a node in a model |
| [Objectidentifier](Objectidentifier.md) | A URI or CURIE that represents an object in the model |
| [Sparqlpath](Sparqlpath.md) | A string encoding a SPARQL Property Path |
| [String](String.md) | A character string |
| [Time](Time.md) | A time object represents a (local) time of day, independent of any particular... |
| [Uri](Uri.md) | a complete URI |
| [Uriorcurie](Uriorcurie.md) | a URI or a CURIE |


## Subsets

| Subset | Description |
| --- | --- |
