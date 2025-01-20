

# Slot: has_alternate_allele


_Links the sequence alteration to its alternate allele (geno:0000002) using the property geno:0000382._





URI: [geno:0000382](http://purl.obolibrary.org/obo/GENO_0000382)



<!-- no inheritance hierarchy -->





## Applicable Classes

| Name | Description | Modifies Slot |
| --- | --- | --- |
| [SequenceAlteration](SequenceAlteration.md) | A representation of a sequence alteration (so:0001059) |  no  |







## Properties

* Range: [AlternateAllele](AlternateAllele.md)

* Required: True





## Identifier and Mapping Information







### Schema Source


* from schema: https://ican.univ-nantes.io/variants-kg




## Mappings

| Mapping Type | Mapped Value |
| ---  | ---  |
| self | geno:0000382 |
| native | https://ican.univ-nantes.io/variants-kg/:has_alternate_allele |




## LinkML Source

<details>
```yaml
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

```
</details>