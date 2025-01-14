

# Slot: has_location


_Links the sequence alteration to its location using the property faldo:Region._





URI: [faldo:location](http://biohackathon.org/resource/faldo#location)



<!-- no inheritance hierarchy -->





## Applicable Classes

| Name | Description | Modifies Slot |
| --- | --- | --- |
| [SequenceAlteration](SequenceAlteration.md) | A representation of a sequence alteration (so:0001059) |  no  |







## Properties

* Range: [VariationSite](VariationSite.md)

* Required: True





## Identifier and Mapping Information







### Schema Source


* from schema: https://ican.univ-nantes.io/variants-kg




## Mappings

| Mapping Type | Mapped Value |
| ---  | ---  |
| self | faldo:location |
| native | https://ican.univ-nantes.io/variants-kg/:has_location |




## LinkML Source

<details>
```yaml
name: has_location
description: Links the sequence alteration to its location using the property faldo:Region.
from_schema: https://ican.univ-nantes.io/variants-kg
rank: 1000
slot_uri: faldo:location
alias: has_location
owner: SequenceAlteration
domain_of:
- SequenceAlteration
range: VariationSite
required: true

```
</details>