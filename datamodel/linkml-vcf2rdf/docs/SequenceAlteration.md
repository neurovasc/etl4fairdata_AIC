

# Class: SequenceAlteration


_A representation of a sequence alteration (so:0001059)._





URI: [so:0001059](http://purl.obolibrary.org/obo/SO_0001059)






```mermaid
 classDiagram
    class SequenceAlteration
    click SequenceAlteration href "../SequenceAlteration"
      SequenceAlteration : has_alternate_allele
        
          
    
    
    SequenceAlteration --> "1" AlternateAllele : has_alternate_allele
    click AlternateAllele href "../AlternateAllele"

        
      SequenceAlteration : has_identifier
        
      SequenceAlteration : has_location
        
          
    
    
    SequenceAlteration --> "1" VariationSite : has_location
    click VariationSite href "../VariationSite"

        
      SequenceAlteration : has_reference_allele
        
          
    
    
    SequenceAlteration --> "1" ReferenceAllele : has_reference_allele
    click ReferenceAllele href "../ReferenceAllele"

        
      SequenceAlteration : is_associated_with
        
          
    
    
    SequenceAlteration --> "1..*" AssociatedCharacteristic : is_associated_with
    click AssociatedCharacteristic href "../AssociatedCharacteristic"

        
      
```




<!-- no inheritance hierarchy -->


## Slots

| Name | Cardinality and Range | Description | Inheritance |
| ---  | --- | --- | --- |
| [has_identifier](has_identifier.md) | 1 <br/> [String](String.md) | A unique identifier for the sequence alteration | direct |
| [has_reference_allele](has_reference_allele.md) | 1 <br/> [ReferenceAllele](ReferenceAllele.md) | Links the sequence alteration to its reference allele (geno:0000036) using th... | direct |
| [has_alternate_allele](has_alternate_allele.md) | 1 <br/> [AlternateAllele](AlternateAllele.md) | Links the sequence alteration to its alternate allele (geno:0000002) using th... | direct |
| [has_location](has_location.md) | 1 <br/> [VariationSite](VariationSite.md) | Links the sequence alteration to its location using the property faldo:Region | direct |
| [is_associated_with](is_associated_with.md) | 1..* <br/> [AssociatedCharacteristic](AssociatedCharacteristic.md) | Links the sequence alteration to a phenotype (female), clinical trait (diabet... | direct |





## Usages

| used by | used in | type | used |
| ---  | --- | --- | --- |
| [Container](Container.md) | [variantcatalog](variantcatalog.md) | range | [SequenceAlteration](SequenceAlteration.md) |






## Identifier and Mapping Information







### Schema Source


* from schema: https://ican.univ-nantes.io/variants-kg




## Mappings

| Mapping Type | Mapped Value |
| ---  | ---  |
| self | so:0001059 |
| native | https://ican.univ-nantes.io/variants-kg/:SequenceAlteration |







## LinkML Source

<!-- TODO: investigate https://stackoverflow.com/questions/37606292/how-to-create-tabbed-code-blocks-in-mkdocs-or-sphinx -->

### Direct

<details>
```yaml
name: SequenceAlteration
description: A representation of a sequence alteration (so:0001059).
from_schema: https://ican.univ-nantes.io/variants-kg
attributes:
  has_identifier:
    name: has_identifier
    description: A unique identifier for the sequence alteration.
    from_schema: https://ican.univ-nantes.io/variants-kg
    rank: 1000
    slot_uri: sio:000300
    domain_of:
    - SequenceAlteration
    - AssociatedCharacteristic
    range: string
    required: true
  has_reference_allele:
    name: has_reference_allele
    description: Links the sequence alteration to its reference allele (geno:0000036)
      using the property geno:0000385.
    from_schema: https://ican.univ-nantes.io/variants-kg
    rank: 1000
    slot_uri: geno:0000385
    domain_of:
    - SequenceAlteration
    range: ReferenceAllele
    required: true
  has_alternate_allele:
    name: has_alternate_allele
    description: Links the sequence alteration to its alternate allele (geno:0000002)
      using the property geno:0000382.
    from_schema: https://ican.univ-nantes.io/variants-kg
    rank: 1000
    slot_uri: geno:0000382
    domain_of:
    - SequenceAlteration
    range: AlternateAllele
    required: true
  has_location:
    name: has_location
    description: Links the sequence alteration to its location using the property
      faldo:Region.
    from_schema: https://ican.univ-nantes.io/variants-kg
    rank: 1000
    slot_uri: faldo:location
    domain_of:
    - SequenceAlteration
    range: VariationSite
    required: true
  is_associated_with:
    name: is_associated_with
    description: Links the sequence alteration to a phenotype (female), clinical trait
      (diabetes) or grouping (whole cohort).
    from_schema: https://ican.univ-nantes.io/variants-kg
    rank: 1000
    slot_uri: sio:001403
    domain_of:
    - SequenceAlteration
    range: AssociatedCharacteristic
    required: true
    multivalued: true
class_uri: so:0001059

```
</details>

### Induced

<details>
```yaml
name: SequenceAlteration
description: A representation of a sequence alteration (so:0001059).
from_schema: https://ican.univ-nantes.io/variants-kg
attributes:
  has_identifier:
    name: has_identifier
    description: A unique identifier for the sequence alteration.
    from_schema: https://ican.univ-nantes.io/variants-kg
    rank: 1000
    slot_uri: sio:000300
    alias: has_identifier
    owner: SequenceAlteration
    domain_of:
    - SequenceAlteration
    - AssociatedCharacteristic
    range: string
    required: true
  has_reference_allele:
    name: has_reference_allele
    description: Links the sequence alteration to its reference allele (geno:0000036)
      using the property geno:0000385.
    from_schema: https://ican.univ-nantes.io/variants-kg
    rank: 1000
    slot_uri: geno:0000385
    alias: has_reference_allele
    owner: SequenceAlteration
    domain_of:
    - SequenceAlteration
    range: ReferenceAllele
    required: true
  has_alternate_allele:
    name: has_alternate_allele
    description: Links the sequence alteration to its alternate allele (geno:0000002)
      using the property geno:0000382.
    from_schema: https://ican.univ-nantes.io/variants-kg
    rank: 1000
    slot_uri: geno:0000382
    alias: has_alternate_allele
    owner: SequenceAlteration
    domain_of:
    - SequenceAlteration
    range: AlternateAllele
    required: true
  has_location:
    name: has_location
    description: Links the sequence alteration to its location using the property
      faldo:Region.
    from_schema: https://ican.univ-nantes.io/variants-kg
    rank: 1000
    slot_uri: faldo:location
    alias: has_location
    owner: SequenceAlteration
    domain_of:
    - SequenceAlteration
    range: VariationSite
    required: true
  is_associated_with:
    name: is_associated_with
    description: Links the sequence alteration to a phenotype (female), clinical trait
      (diabetes) or grouping (whole cohort).
    from_schema: https://ican.univ-nantes.io/variants-kg
    rank: 1000
    slot_uri: sio:001403
    alias: is_associated_with
    owner: SequenceAlteration
    domain_of:
    - SequenceAlteration
    range: AssociatedCharacteristic
    required: true
    multivalued: true
class_uri: so:0001059

```
</details>