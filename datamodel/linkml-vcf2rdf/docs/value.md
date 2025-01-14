

# Slot: value



URI: [https://ican.univ-nantes.io/variants-kg/:value](https://ican.univ-nantes.io/variants-kg/:value)



<!-- no inheritance hierarchy -->





## Applicable Classes

| Name | Description | Modifies Slot |
| --- | --- | --- |
| [AlternateAllele](AlternateAllele.md) | Represents the alternate allele (geno:0000002) |  no  |
| [Count](Count.md) | Represents the count of alternative alleles in individuals with a certain zyg... |  no  |
| [Frequency](Frequency.md) | Represents the frequency of an allele in individuals with a certain zygosity |  no  |
| [ReferenceAllele](ReferenceAllele.md) | Represents the reference allele (geno:0000036) |  no  |
| [VariationSiteReference](VariationSiteReference.md) | Represents the reference sequence (contig, sequence, chromosome) |  no  |







## Properties

* Range: [String](String.md)





## Identifier and Mapping Information








## Mappings

| Mapping Type | Mapped Value |
| ---  | ---  |
| self | https://ican.univ-nantes.io/variants-kg/:value |
| native | https://ican.univ-nantes.io/variants-kg/:value |




## LinkML Source

<details>
```yaml
name: value
alias: value
domain_of:
- ReferenceAllele
- AlternateAllele
- VariationSiteReference
- Frequency
- Count
range: string

```
</details>