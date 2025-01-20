

# Class: VariationSiteEnd


_Represents the end of the location of a sequence alteration._





URI: [faldo:ExactPosition](http://biohackathon.org/resource/faldo#ExactPosition)






```mermaid
 classDiagram
    class VariationSiteEnd
    click VariationSiteEnd href "../VariationSiteEnd"
      VariationSiteEnd : position
        
      
```




<!-- no inheritance hierarchy -->


## Slots

| Name | Cardinality and Range | Description | Inheritance |
| ---  | --- | --- | --- |
| [position](position.md) | 1 <br/> [Integer](Integer.md) | The position of the end of the location of the sequence alteration | direct |





## Usages

| used by | used in | type | used |
| ---  | --- | --- | --- |
| [VariationSite](VariationSite.md) | [ends_at](ends_at.md) | range | [VariationSiteEnd](VariationSiteEnd.md) |






## Identifier and Mapping Information







### Schema Source


* from schema: https://ican.univ-nantes.io/variants-kg




## Mappings

| Mapping Type | Mapped Value |
| ---  | ---  |
| self | faldo:ExactPosition |
| native | https://ican.univ-nantes.io/variants-kg/:VariationSiteEnd |







## LinkML Source

<!-- TODO: investigate https://stackoverflow.com/questions/37606292/how-to-create-tabbed-code-blocks-in-mkdocs-or-sphinx -->

### Direct

<details>
```yaml
name: VariationSiteEnd
description: Represents the end of the location of a sequence alteration.
from_schema: https://ican.univ-nantes.io/variants-kg
attributes:
  position:
    name: position
    description: The position of the end of the location of the sequence alteration.
    from_schema: https://ican.univ-nantes.io/variants-kg
    slot_uri: faldo:position
    domain_of:
    - VariationSiteBegin
    - VariationSiteEnd
    range: integer
    required: true
class_uri: faldo:ExactPosition

```
</details>

### Induced

<details>
```yaml
name: VariationSiteEnd
description: Represents the end of the location of a sequence alteration.
from_schema: https://ican.univ-nantes.io/variants-kg
attributes:
  position:
    name: position
    description: The position of the end of the location of the sequence alteration.
    from_schema: https://ican.univ-nantes.io/variants-kg
    slot_uri: faldo:position
    alias: position
    owner: VariationSiteEnd
    domain_of:
    - VariationSiteBegin
    - VariationSiteEnd
    range: integer
    required: true
class_uri: faldo:ExactPosition

```
</details>