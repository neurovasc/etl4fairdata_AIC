
# Class: AssociatedCharacteristic

Represents a phenotype, clinical trait or grouping associated with a sequence alteration.

URI: [https://ican.univ-nantes.io/variants-kg/AssociatedCharacteristic](https://ican.univ-nantes.io/variants-kg/AssociatedCharacteristic)


[![img](https://yuml.me/diagram/nofunky;dir:TB/class/[Zygosity],[Zygosity]<has_zygosity%200..1-++[AssociatedCharacteristic&#124;has_identifier:string%20%3F;label:string%20%3F],[SequenceAlteration]++-%20is_associated_with%201..*>[AssociatedCharacteristic],[SequenceAlteration])](https://yuml.me/diagram/nofunky;dir:TB/class/[Zygosity],[Zygosity]<has_zygosity%200..1-++[AssociatedCharacteristic&#124;has_identifier:string%20%3F;label:string%20%3F],[SequenceAlteration]++-%20is_associated_with%201..*>[AssociatedCharacteristic],[SequenceAlteration])

## Referenced by Class

 *  **None** *[➞is_associated_with](sequenceAlteration__is_associated_with.md)*  <sub>1..\*</sub>  **[AssociatedCharacteristic](AssociatedCharacteristic.md)**

## Attributes


### Own

 * [➞has_identifier](associatedCharacteristic__has_identifier.md)  <sub>0..1</sub>
     * Description: A unique identifier for the associated characteristic.
     * Range: [String](types/String.md)
 * [➞label](associatedCharacteristic__label.md)  <sub>0..1</sub>
     * Description: A human-readable label for the associated characteristic.
     * Range: [String](types/String.md)
 * [➞has_zygosity](associatedCharacteristic__has_zygosity.md)  <sub>0..1</sub>
     * Description: The zygocity of the associated characteristic.
     * Range: [Zygosity](Zygosity.md)
