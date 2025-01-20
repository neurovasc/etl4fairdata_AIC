
# Class: VariationSite

Represents the location of a sequence alteration.

URI: [https://ican.univ-nantes.io/variants-kg/VariationSite](https://ican.univ-nantes.io/variants-kg/VariationSite)


[![img](https://yuml.me/diagram/nofunky;dir:TB/class/[VariationSiteReference],[VariationSiteEnd],[VariationSiteBegin],[VariationSiteReference]<has_reference%201..1-++[VariationSite],[VariationSiteEnd]<ends_at%201..1-++[VariationSite],[VariationSiteBegin]<begins_at%201..1-++[VariationSite],[SequenceAlteration]++-%20has_location%201..1>[VariationSite],[SequenceAlteration])](https://yuml.me/diagram/nofunky;dir:TB/class/[VariationSiteReference],[VariationSiteEnd],[VariationSiteBegin],[VariationSiteReference]<has_reference%201..1-++[VariationSite],[VariationSiteEnd]<ends_at%201..1-++[VariationSite],[VariationSiteBegin]<begins_at%201..1-++[VariationSite],[SequenceAlteration]++-%20has_location%201..1>[VariationSite],[SequenceAlteration])

## Referenced by Class

 *  **None** *[➞has_location](sequenceAlteration__has_location.md)*  <sub>1..1</sub>  **[VariationSite](VariationSite.md)**

## Attributes


### Own

 * [➞begins_at](variationSite__begins_at.md)  <sub>1..1</sub>
     * Description: The beginning of the location of the sequence alteration.
     * Range: [VariationSiteBegin](VariationSiteBegin.md)
 * [➞ends_at](variationSite__ends_at.md)  <sub>1..1</sub>
     * Description: The end of the location of the sequence alteration.
     * Range: [VariationSiteEnd](VariationSiteEnd.md)
 * [➞has_reference](variationSite__has_reference.md)  <sub>1..1</sub>
     * Description: The reference sequence (contig, sequence, chromosome).
     * Range: [VariationSiteReference](VariationSiteReference.md)

## Other properties

|  |  |  |
| --- | --- | --- |
| **Mappings:** | | faldo:Region |