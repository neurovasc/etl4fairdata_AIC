
# Class: Zygosity

Represents the zygosity of an associated characteristic.

URI: [https://ican.univ-nantes.io/variants-kg/Zygosity](https://ican.univ-nantes.io/variants-kg/Zygosity)


[![img](https://yuml.me/diagram/nofunky;dir:TB/class/[Frequency]<has_frequency%201..1-++[Zygosity&#124;type:ZygosityType],[Count]<has_measurement_value%201..1-++[Zygosity],[AssociatedCharacteristic]++-%20has_zygosity%200..1>[Zygosity],[Frequency],[Count],[AssociatedCharacteristic])](https://yuml.me/diagram/nofunky;dir:TB/class/[Frequency]<has_frequency%201..1-++[Zygosity&#124;type:ZygosityType],[Count]<has_measurement_value%201..1-++[Zygosity],[AssociatedCharacteristic]++-%20has_zygosity%200..1>[Zygosity],[Frequency],[Count],[AssociatedCharacteristic])

## Referenced by Class

 *  **None** *[➞has_zygosity](associatedCharacteristic__has_zygosity.md)*  <sub>0..1</sub>  **[Zygosity](Zygosity.md)**

## Attributes


### Own

 * [➞type](zygosity__type.md)  <sub>1..1</sub>
     * Description: The type of zygosity.
     * Range: [ZygosityType](ZygosityType.md)
 * [➞has_measurement_value](zygosity__has_measurement_value.md)  <sub>1..1</sub>
     * Description: Count of alternative alleles in individuals with a certain zygosity.
     * Range: [Count](Count.md)
 * [➞has_frequency](zygosity__has_frequency.md)  <sub>1..1</sub>
     * Description: Frequency of the allele in indibiduals with a certain zygosity.
     * Range: [Frequency](Frequency.md)

## Other properties

|  |  |  |
| --- | --- | --- |
| **Mappings:** | | geno:0000133 |