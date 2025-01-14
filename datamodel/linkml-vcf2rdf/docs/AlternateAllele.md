

# Class: AlternateAllele


_Represents the alternate allele (geno:0000002)._





URI: [geno:0000002](http://purl.obolibrary.org/obo/GENO_0000002)






```mermaid
 classDiagram
    class AlternateAllele
    click AlternateAllele href "../AlternateAllele"
      AlternateAllele : value
        
      
```




<!-- no inheritance hierarchy -->


## Slots

| Name | Cardinality and Range | Description | Inheritance |
| ---  | --- | --- | --- |
| [value](value.md) | 1 <br/> [String](String.md) | The value of the alternate allele | direct |





## Usages

| used by | used in | type | used |
| ---  | --- | --- | --- |
| [SequenceAlteration](SequenceAlteration.md) | [has_alternate_allele](has_alternate_allele.md) | range | [AlternateAllele](AlternateAllele.md) |






## Identifier and Mapping Information







### Schema Source


* from schema: https://ican.univ-nantes.io/variants-kg




## Mappings

| Mapping Type | Mapped Value |
| ---  | ---  |
| self | geno:0000002 |
| native | https://ican.univ-nantes.io/variants-kg/:AlternateAllele |







## LinkML Source

<!-- TODO: investigate https://stackoverflow.com/questions/37606292/how-to-create-tabbed-code-blocks-in-mkdocs-or-sphinx -->

### Direct

<details>
```yaml
name: AlternateAllele
description: Represents the alternate allele (geno:0000002).
from_schema: https://ican.univ-nantes.io/variants-kg
attributes:
  value:
    name: value
    description: The value of the alternate allele.
    from_schema: https://ican.univ-nantes.io/variants-kg
    slot_uri: sio:000300
    domain_of:
    - ReferenceAllele
    - AlternateAllele
    - VariationSiteReference
    - Frequency
    - Count
    range: string
    required: true
class_uri: geno:0000002

```
</details>

### Induced

<details>
```yaml
name: AlternateAllele
description: Represents the alternate allele (geno:0000002).
from_schema: https://ican.univ-nantes.io/variants-kg
attributes:
  value:
    name: value
    description: The value of the alternate allele.
    from_schema: https://ican.univ-nantes.io/variants-kg
    slot_uri: sio:000300
    alias: value
    owner: AlternateAllele
    domain_of:
    - ReferenceAllele
    - AlternateAllele
    - VariationSiteReference
    - Frequency
    - Count
    range: string
    required: true
class_uri: geno:0000002

```
</details>