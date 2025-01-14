

# Slot: is_associated_with


_Links the sequence alteration to a phenotype (female), clinical trait (diabetes) or grouping (whole cohort)._





URI: [sio:001403](http://semanticscience.org/resource/001403)



<!-- no inheritance hierarchy -->





## Applicable Classes

| Name | Description | Modifies Slot |
| --- | --- | --- |
| [SequenceAlteration](SequenceAlteration.md) | A representation of a sequence alteration (so:0001059) |  no  |







## Properties

* Range: [AssociatedCharacteristic](AssociatedCharacteristic.md)

* Multivalued: True

* Required: True





## Identifier and Mapping Information







### Schema Source


* from schema: https://ican.univ-nantes.io/variants-kg




## Mappings

| Mapping Type | Mapped Value |
| ---  | ---  |
| self | sio:001403 |
| native | https://ican.univ-nantes.io/variants-kg/:is_associated_with |




## LinkML Source

<details>
```yaml
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

```
</details>