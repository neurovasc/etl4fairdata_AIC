

# Class: Container



URI: [https://ican.univ-nantes.io/variants-kg/:Container](https://ican.univ-nantes.io/variants-kg/:Container)






```mermaid
 classDiagram
    class Container
    click Container href "../Container"
      Container : variantcatalog
        
          
    
    
    Container --> "*" SequenceAlteration : variantcatalog
    click SequenceAlteration href "../SequenceAlteration"

        
      
```




<!-- no inheritance hierarchy -->


## Slots

| Name | Cardinality and Range | Description | Inheritance |
| ---  | --- | --- | --- |
| [variantcatalog](variantcatalog.md) | * <br/> [SequenceAlteration](SequenceAlteration.md) |  | direct |









## Identifier and Mapping Information







### Schema Source


* from schema: https://ican.univ-nantes.io/variants-kg




## Mappings

| Mapping Type | Mapped Value |
| ---  | ---  |
| self | https://ican.univ-nantes.io/variants-kg/:Container |
| native | https://ican.univ-nantes.io/variants-kg/:Container |







## LinkML Source

<!-- TODO: investigate https://stackoverflow.com/questions/37606292/how-to-create-tabbed-code-blocks-in-mkdocs-or-sphinx -->

### Direct

<details>
```yaml
name: Container
from_schema: https://ican.univ-nantes.io/variants-kg
attributes:
  variantcatalog:
    name: variantcatalog
    from_schema: https://ican.univ-nantes.io/variants-kg
    rank: 1000
    domain_of:
    - Container
    range: SequenceAlteration
    multivalued: true
    inlined: true
    inlined_as_list: true

```
</details>

### Induced

<details>
```yaml
name: Container
from_schema: https://ican.univ-nantes.io/variants-kg
attributes:
  variantcatalog:
    name: variantcatalog
    from_schema: https://ican.univ-nantes.io/variants-kg
    rank: 1000
    alias: variantcatalog
    owner: Container
    domain_of:
    - Container
    range: SequenceAlteration
    multivalued: true
    inlined: true
    inlined_as_list: true

```
</details>