

# Slot: sameAs


_The value of the reference sequence in another database (ncbi sequence i.e. https://www.ncbi.nlm.nih.gov/nuccore/NC_000022.11)._





URI: [owl:sameAs](http://www.w3.org/2002/07/owl#sameAs)



<!-- no inheritance hierarchy -->





## Applicable Classes

| Name | Description | Modifies Slot |
| --- | --- | --- |
| [VariationSiteReference](VariationSiteReference.md) | Represents the reference sequence (contig, sequence, chromosome) |  no  |







## Properties

* Range: [String](String.md)

* Required: True





## Identifier and Mapping Information







### Schema Source


* from schema: https://ican.univ-nantes.io/variants-kg




## Mappings

| Mapping Type | Mapped Value |
| ---  | ---  |
| self | owl:sameAs |
| native | https://ican.univ-nantes.io/variants-kg/:sameAs |




## LinkML Source

<details>
```yaml
name: sameAs
description: The value of the reference sequence in another database (ncbi sequence
  i.e. https://www.ncbi.nlm.nih.gov/nuccore/NC_000022.11).
from_schema: https://ican.univ-nantes.io/variants-kg
rank: 1000
slot_uri: owl:sameAs
alias: sameAs
owner: VariationSiteReference
domain_of:
- VariationSiteReference
range: string
required: true

```
</details>