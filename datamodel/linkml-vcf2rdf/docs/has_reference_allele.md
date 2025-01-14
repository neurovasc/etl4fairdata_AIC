

# Slot: has_reference_allele


_Links the sequence alteration to its reference allele (geno:0000036) using the property geno:0000385._





URI: [geno:0000385](http://purl.obolibrary.org/obo/GENO_0000385)



<!-- no inheritance hierarchy -->





## Applicable Classes

| Name | Description | Modifies Slot |
| --- | --- | --- |
| [SequenceAlteration](SequenceAlteration.md) | A representation of a sequence alteration (so:0001059) |  no  |







## Properties

* Range: [ReferenceAllele](ReferenceAllele.md)

* Required: True





## Identifier and Mapping Information







### Schema Source


* from schema: https://ican.univ-nantes.io/variants-kg




## Mappings

| Mapping Type | Mapped Value |
| ---  | ---  |
| self | geno:0000385 |
| native | https://ican.univ-nantes.io/variants-kg/:has_reference_allele |




## LinkML Source

<details>
```yaml
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

```
</details>