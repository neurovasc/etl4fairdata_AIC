

# Class: Count


_Represents the count of alternative alleles in individuals with a certain zygosity._





URI: [https://ican.univ-nantes.io/variants-kg/:Count](https://ican.univ-nantes.io/variants-kg/:Count)






```mermaid
 classDiagram
    class Count
    click Count href "../Count"
      Count : value
        
      
```




<!-- no inheritance hierarchy -->


## Slots

| Name | Cardinality and Range | Description | Inheritance |
| ---  | --- | --- | --- |
| [value](value.md) | 1 <br/> [Integer](Integer.md) | The count value | direct |





## Usages

| used by | used in | type | used |
| ---  | --- | --- | --- |
| [Zygosity](Zygosity.md) | [has_measurement_value](has_measurement_value.md) | range | [Count](Count.md) |






## Identifier and Mapping Information







### Schema Source


* from schema: https://ican.univ-nantes.io/variants-kg




## Mappings

| Mapping Type | Mapped Value |
| ---  | ---  |
| self | https://ican.univ-nantes.io/variants-kg/:Count |
| native | https://ican.univ-nantes.io/variants-kg/:Count |







## LinkML Source

<!-- TODO: investigate https://stackoverflow.com/questions/37606292/how-to-create-tabbed-code-blocks-in-mkdocs-or-sphinx -->

### Direct

<details>
```yaml
name: Count
description: Represents the count of alternative alleles in individuals with a certain
  zygosity.
from_schema: https://ican.univ-nantes.io/variants-kg
attributes:
  value:
    name: value
    description: The count value.
    from_schema: https://ican.univ-nantes.io/variants-kg
    slot_uri: sio:000300
    domain_of:
    - ReferenceAllele
    - AlternateAllele
    - VariationSiteReference
    - Frequency
    - Count
    range: integer
    required: true

```
</details>

### Induced

<details>
```yaml
name: Count
description: Represents the count of alternative alleles in individuals with a certain
  zygosity.
from_schema: https://ican.univ-nantes.io/variants-kg
attributes:
  value:
    name: value
    description: The count value.
    from_schema: https://ican.univ-nantes.io/variants-kg
    slot_uri: sio:000300
    alias: value
    owner: Count
    domain_of:
    - ReferenceAllele
    - AlternateAllele
    - VariationSiteReference
    - Frequency
    - Count
    range: integer
    required: true

```
</details>