
# Class: SequenceAlteration

A representation of a sequence alteration (so:0001059).

URI: [https://ican.univ-nantes.io/variants-kg/SequenceAlteration](https://ican.univ-nantes.io/variants-kg/SequenceAlteration)


[![img](https://yuml.me/diagram/nofunky;dir:TB/class/[VariationSite],[AssociatedCharacteristic]<is_associated_with%201..*-++[SequenceAlteration&#124;has_identifier:string],[VariationSite]<has_location%201..1-++[SequenceAlteration],[AlternateAllele]<has_alternate_allele%201..1-++[SequenceAlteration],[ReferenceAllele]<has_reference_allele%201..1-++[SequenceAlteration],[Container]++-%20variantcatalog%200..*>[SequenceAlteration],[ReferenceAllele],[Container],[AssociatedCharacteristic],[AlternateAllele])](https://yuml.me/diagram/nofunky;dir:TB/class/[VariationSite],[AssociatedCharacteristic]<is_associated_with%201..*-++[SequenceAlteration&#124;has_identifier:string],[VariationSite]<has_location%201..1-++[SequenceAlteration],[AlternateAllele]<has_alternate_allele%201..1-++[SequenceAlteration],[ReferenceAllele]<has_reference_allele%201..1-++[SequenceAlteration],[Container]++-%20variantcatalog%200..*>[SequenceAlteration],[ReferenceAllele],[Container],[AssociatedCharacteristic],[AlternateAllele])

## Referenced by Class

 *  **None** *[➞variantcatalog](container__variantcatalog.md)*  <sub>0..\*</sub>  **[SequenceAlteration](SequenceAlteration.md)**

## Attributes


### Own

 * [➞has_identifier](sequenceAlteration__has_identifier.md)  <sub>1..1</sub>
     * Description: A unique identifier for the sequence alteration.
     * Range: [String](types/String.md)
 * [➞has_reference_allele](sequenceAlteration__has_reference_allele.md)  <sub>1..1</sub>
     * Description: Links the sequence alteration to its reference allele (geno:0000036) using the property geno:0000385.
     * Range: [ReferenceAllele](ReferenceAllele.md)
 * [➞has_alternate_allele](sequenceAlteration__has_alternate_allele.md)  <sub>1..1</sub>
     * Description: Links the sequence alteration to its alternate allele (geno:0000002) using the property geno:0000382.
     * Range: [AlternateAllele](AlternateAllele.md)
 * [➞has_location](sequenceAlteration__has_location.md)  <sub>1..1</sub>
     * Description: Links the sequence alteration to its location using the property faldo:Region.
     * Range: [VariationSite](VariationSite.md)
 * [➞is_associated_with](sequenceAlteration__is_associated_with.md)  <sub>1..\*</sub>
     * Description: Links the sequence alteration to a phenotype (female), clinical trait (diabetes) or grouping (whole cohort).
     * Range: [AssociatedCharacteristic](AssociatedCharacteristic.md)

## Other properties

|  |  |  |
| --- | --- | --- |
| **Mappings:** | | so:0001059 |