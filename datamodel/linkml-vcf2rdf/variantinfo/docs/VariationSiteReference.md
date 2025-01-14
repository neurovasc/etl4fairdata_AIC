
# Class: VariationSiteReference

Represents the reference sequence (contig, sequence, chromosome).

URI: [https://ican.univ-nantes.io/variants-kg/VariationSiteReference](https://ican.univ-nantes.io/variants-kg/VariationSiteReference)


[![img](https://yuml.me/diagram/nofunky;dir:TB/class/[VariationSite]++-%20has_reference%201..1>[VariationSiteReference&#124;value:string;label:string;sameAs:string],[VariationSite])](https://yuml.me/diagram/nofunky;dir:TB/class/[VariationSite]++-%20has_reference%201..1>[VariationSiteReference&#124;value:string;label:string;sameAs:string],[VariationSite])

## Referenced by Class

 *  **None** *[➞has_reference](variationSite__has_reference.md)*  <sub>1..1</sub>  **[VariationSiteReference](VariationSiteReference.md)**

## Attributes


### Own

 * [➞value](variationSiteReference__value.md)  <sub>1..1</sub>
     * Description: The value of the reference sequence (ena sequence i.e. https://www.ebi.ac.uk/ena/browser/view/CM000684.2).
     * Range: [String](types/String.md)
 * [➞label](variationSiteReference__label.md)  <sub>1..1</sub>
     * Description: A human-readable label for the reference sequence.
     * Range: [String](types/String.md)
 * [➞sameAs](variationSiteReference__sameAs.md)  <sub>1..1</sub>
     * Description: The value of the reference sequence in another database (ncbi sequence i.e. https://www.ncbi.nlm.nih.gov/nuccore/NC_000022.11).
     * Range: [String](types/String.md)

## Other properties

|  |  |  |
| --- | --- | --- |
| **Mappings:** | | so:0000353 |