# Auto generated from variantinfo.yaml by pythongen.py version: 0.0.1
# Generation date: 2025-01-14T14:41:41
# Schema: variants-kg
#
# id: https://ican.univ-nantes.io/variants-kg
# description: Schema for describing sequence alterations (SO:0001059) and their relationships to other elements of a variant catalog format.
#   The ontology aims at describing genomic variations found in patients carrying intracranial aneurysms.
#   The ontology and schema is organized arount variants, their annotations relevant to the study of
#   intracranial aneurysms, and phenotype frequencies and counts relevant to the study of intracranial aneurysms.
# license: https://creativecommons.org/publicdomain/zero/1.0/

import dataclasses
import re
from dataclasses import dataclass
from datetime import (
    date,
    datetime,
    time
)
from typing import (
    Any,
    ClassVar,
    Dict,
    List,
    Optional,
    Union
)

from jsonasobj2 import (
    JsonObj,
    as_dict
)
from linkml_runtime.linkml_model.meta import (
    EnumDefinition,
    PermissibleValue,
    PvFormulaOptions
)
from linkml_runtime.utils.curienamespace import CurieNamespace
from linkml_runtime.utils.dataclass_extensions_376 import dataclasses_init_fn_with_kwargs
from linkml_runtime.utils.enumerations import EnumDefinitionImpl
from linkml_runtime.utils.formatutils import (
    camelcase,
    sfx,
    underscore
)
from linkml_runtime.utils.metamodelcore import (
    bnode,
    empty_dict,
    empty_list
)
from linkml_runtime.utils.slot import Slot
from linkml_runtime.utils.yamlutils import (
    YAMLRoot,
    extended_float,
    extended_int,
    extended_str
)
from rdflib import (
    Namespace,
    URIRef
)

from linkml_runtime.linkml_model.types import Float, Integer, String

metamodel_version = "1.7.0"
version = None

# Overwrite dataclasses _init_fn to add **kwargs in __init__
dataclasses._init_fn = dataclasses_init_fn_with_kwargs

# Namespaces
AIC = CurieNamespace('aic', 'https://example.org/aiccohortvariants#')
FALDO = CurieNamespace('faldo', 'http://biohackathon.org/resource/faldo#')
GENO = CurieNamespace('geno', 'http://purl.obolibrary.org/obo/GENO_')
HPO = CurieNamespace('hpo', 'http://purl.obolibrary.org/obo/HP_')
LINKML = CurieNamespace('linkml', 'https://w3id.org/linkml/')
MONDO = CurieNamespace('mondo', 'http://purl.obolibrary.org/obo/MONDO_')
NCIT = CurieNamespace('ncit', 'http://purl.obolibrary.org/obo/NCIT_')
OWL = CurieNamespace('owl', 'http://www.w3.org/2002/07/owl#')
RDFS = CurieNamespace('rdfs', 'https://www.w3.org/2000/01/rdf-schema#')
SIO = CurieNamespace('sio', 'http://semanticscience.org/resource/')
SO = CurieNamespace('so', 'http://purl.obolibrary.org/obo/SO_')
DEFAULT_ = CurieNamespace('', 'https://ican.univ-nantes.io/variants-kg/')


# Types

# Class references



@dataclass(repr=False)
class SequenceAlteration(YAMLRoot):
    """
    A representation of a sequence alteration (so:0001059).
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = SO["0001059"]
    class_class_curie: ClassVar[str] = "so:0001059"
    class_name: ClassVar[str] = "SequenceAlteration"
    class_model_uri: ClassVar[URIRef] = URIRef("https://ican.univ-nantes.io/variants-kg/SequenceAlteration")

    has_identifier: str = None
    has_reference_allele: Union[dict, "ReferenceAllele"] = None
    has_alternate_allele: Union[dict, "AlternateAllele"] = None
    has_location: Union[dict, "VariationSite"] = None
    is_associated_with: Union[Union[dict, "AssociatedCharacteristic"], List[Union[dict, "AssociatedCharacteristic"]]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.has_identifier):
            self.MissingRequiredField("has_identifier")
        if not isinstance(self.has_identifier, str):
            self.has_identifier = str(self.has_identifier)

        if self._is_empty(self.has_reference_allele):
            self.MissingRequiredField("has_reference_allele")
        if not isinstance(self.has_reference_allele, ReferenceAllele):
            self.has_reference_allele = ReferenceAllele(**as_dict(self.has_reference_allele))

        if self._is_empty(self.has_alternate_allele):
            self.MissingRequiredField("has_alternate_allele")
        if not isinstance(self.has_alternate_allele, AlternateAllele):
            self.has_alternate_allele = AlternateAllele(**as_dict(self.has_alternate_allele))

        if self._is_empty(self.has_location):
            self.MissingRequiredField("has_location")
        if not isinstance(self.has_location, VariationSite):
            self.has_location = VariationSite(**as_dict(self.has_location))

        if self._is_empty(self.is_associated_with):
            self.MissingRequiredField("is_associated_with")
        if not isinstance(self.is_associated_with, list):
            self.is_associated_with = [self.is_associated_with] if self.is_associated_with is not None else []
        self.is_associated_with = [v if isinstance(v, AssociatedCharacteristic) else AssociatedCharacteristic(**as_dict(v)) for v in self.is_associated_with]

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class ReferenceAllele(YAMLRoot):
    """
    Represents the reference allele (geno:0000036).
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GENO["0000036"]
    class_class_curie: ClassVar[str] = "geno:0000036"
    class_name: ClassVar[str] = "ReferenceAllele"
    class_model_uri: ClassVar[URIRef] = URIRef("https://ican.univ-nantes.io/variants-kg/ReferenceAllele")

    value: str = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.value):
            self.MissingRequiredField("value")
        if not isinstance(self.value, str):
            self.value = str(self.value)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class AlternateAllele(YAMLRoot):
    """
    Represents the alternate allele (geno:0000002).
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GENO["0000002"]
    class_class_curie: ClassVar[str] = "geno:0000002"
    class_name: ClassVar[str] = "AlternateAllele"
    class_model_uri: ClassVar[URIRef] = URIRef("https://ican.univ-nantes.io/variants-kg/AlternateAllele")

    value: str = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.value):
            self.MissingRequiredField("value")
        if not isinstance(self.value, str):
            self.value = str(self.value)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class VariationSite(YAMLRoot):
    """
    Represents the location of a sequence alteration.
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = FALDO["Region"]
    class_class_curie: ClassVar[str] = "faldo:Region"
    class_name: ClassVar[str] = "VariationSite"
    class_model_uri: ClassVar[URIRef] = URIRef("https://ican.univ-nantes.io/variants-kg/VariationSite")

    begins_at: Union[dict, "VariationSiteBegin"] = None
    ends_at: Union[dict, "VariationSiteEnd"] = None
    has_reference: Union[dict, "VariationSiteReference"] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.begins_at):
            self.MissingRequiredField("begins_at")
        if not isinstance(self.begins_at, VariationSiteBegin):
            self.begins_at = VariationSiteBegin(**as_dict(self.begins_at))

        if self._is_empty(self.ends_at):
            self.MissingRequiredField("ends_at")
        if not isinstance(self.ends_at, VariationSiteEnd):
            self.ends_at = VariationSiteEnd(**as_dict(self.ends_at))

        if self._is_empty(self.has_reference):
            self.MissingRequiredField("has_reference")
        if not isinstance(self.has_reference, VariationSiteReference):
            self.has_reference = VariationSiteReference(**as_dict(self.has_reference))

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class VariationSiteBegin(YAMLRoot):
    """
    Represents the beginning of the location of a sequence alteration.
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = FALDO["ExactPosition"]
    class_class_curie: ClassVar[str] = "faldo:ExactPosition"
    class_name: ClassVar[str] = "VariationSiteBegin"
    class_model_uri: ClassVar[URIRef] = URIRef("https://ican.univ-nantes.io/variants-kg/VariationSiteBegin")

    position: int = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.position):
            self.MissingRequiredField("position")
        if not isinstance(self.position, int):
            self.position = int(self.position)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class VariationSiteEnd(YAMLRoot):
    """
    Represents the end of the location of a sequence alteration.
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = FALDO["ExactPosition"]
    class_class_curie: ClassVar[str] = "faldo:ExactPosition"
    class_name: ClassVar[str] = "VariationSiteEnd"
    class_model_uri: ClassVar[URIRef] = URIRef("https://ican.univ-nantes.io/variants-kg/VariationSiteEnd")

    position: int = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.position):
            self.MissingRequiredField("position")
        if not isinstance(self.position, int):
            self.position = int(self.position)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class VariationSiteReference(YAMLRoot):
    """
    Represents the reference sequence (contig, sequence, chromosome).
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = SO["0000353"]
    class_class_curie: ClassVar[str] = "so:0000353"
    class_name: ClassVar[str] = "VariationSiteReference"
    class_model_uri: ClassVar[URIRef] = URIRef("https://ican.univ-nantes.io/variants-kg/VariationSiteReference")

    value: str = None
    label: str = None
    sameAs: str = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.value):
            self.MissingRequiredField("value")
        if not isinstance(self.value, str):
            self.value = str(self.value)

        if self._is_empty(self.label):
            self.MissingRequiredField("label")
        if not isinstance(self.label, str):
            self.label = str(self.label)

        if self._is_empty(self.sameAs):
            self.MissingRequiredField("sameAs")
        if not isinstance(self.sameAs, str):
            self.sameAs = str(self.sameAs)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class AssociatedCharacteristic(YAMLRoot):
    """
    Represents a phenotype, clinical trait or grouping associated with a sequence alteration.
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = URIRef("https://ican.univ-nantes.io/variants-kg/AssociatedCharacteristic")
    class_class_curie: ClassVar[str] = None
    class_name: ClassVar[str] = "AssociatedCharacteristic"
    class_model_uri: ClassVar[URIRef] = URIRef("https://ican.univ-nantes.io/variants-kg/AssociatedCharacteristic")

    has_identifier: Optional[str] = None
    label: Optional[str] = None
    has_zygosity: Optional[Union[dict, "Zygosity"]] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self.has_identifier is not None and not isinstance(self.has_identifier, str):
            self.has_identifier = str(self.has_identifier)

        if self.label is not None and not isinstance(self.label, str):
            self.label = str(self.label)

        if self.has_zygosity is not None and not isinstance(self.has_zygosity, Zygosity):
            self.has_zygosity = Zygosity(**as_dict(self.has_zygosity))

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class Zygosity(YAMLRoot):
    """
    Represents the zygosity of an associated characteristic.
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = GENO["0000133"]
    class_class_curie: ClassVar[str] = "geno:0000133"
    class_name: ClassVar[str] = "Zygosity"
    class_model_uri: ClassVar[URIRef] = URIRef("https://ican.univ-nantes.io/variants-kg/Zygosity")

    type: Union[str, "ZygosityType"] = None
    has_measurement_value: Union[dict, "Count"] = None
    has_frequency: Union[dict, "Frequency"] = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.type):
            self.MissingRequiredField("type")
        if not isinstance(self.type, ZygosityType):
            self.type = ZygosityType(self.type)

        if self._is_empty(self.has_measurement_value):
            self.MissingRequiredField("has_measurement_value")
        if not isinstance(self.has_measurement_value, Count):
            self.has_measurement_value = Count(**as_dict(self.has_measurement_value))

        if self._is_empty(self.has_frequency):
            self.MissingRequiredField("has_frequency")
        if not isinstance(self.has_frequency, Frequency):
            self.has_frequency = Frequency(**as_dict(self.has_frequency))

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class Frequency(YAMLRoot):
    """
    Represents the frequency of an allele in individuals with a certain zygosity.
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = URIRef("https://ican.univ-nantes.io/variants-kg/Frequency")
    class_class_curie: ClassVar[str] = None
    class_name: ClassVar[str] = "Frequency"
    class_model_uri: ClassVar[URIRef] = URIRef("https://ican.univ-nantes.io/variants-kg/Frequency")

    value: float = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.value):
            self.MissingRequiredField("value")
        if not isinstance(self.value, float):
            self.value = float(self.value)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class Count(YAMLRoot):
    """
    Represents the count of alternative alleles in individuals with a certain zygosity.
    """
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = URIRef("https://ican.univ-nantes.io/variants-kg/Count")
    class_class_curie: ClassVar[str] = None
    class_name: ClassVar[str] = "Count"
    class_model_uri: ClassVar[URIRef] = URIRef("https://ican.univ-nantes.io/variants-kg/Count")

    value: int = None

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if self._is_empty(self.value):
            self.MissingRequiredField("value")
        if not isinstance(self.value, int):
            self.value = int(self.value)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class Container(YAMLRoot):
    _inherited_slots: ClassVar[List[str]] = []

    class_class_uri: ClassVar[URIRef] = URIRef("https://ican.univ-nantes.io/variants-kg/Container")
    class_class_curie: ClassVar[str] = None
    class_name: ClassVar[str] = "Container"
    class_model_uri: ClassVar[URIRef] = URIRef("https://ican.univ-nantes.io/variants-kg/Container")

    variantcatalog: Optional[Union[Union[dict, SequenceAlteration], List[Union[dict, SequenceAlteration]]]] = empty_list()

    def __post_init__(self, *_: List[str], **kwargs: Dict[str, Any]):
        if not isinstance(self.variantcatalog, list):
            self.variantcatalog = [self.variantcatalog] if self.variantcatalog is not None else []
        self.variantcatalog = [v if isinstance(v, SequenceAlteration) else SequenceAlteration(**as_dict(v)) for v in self.variantcatalog]

        super().__post_init__(**kwargs)


# Enumerations
class ZygosityType(EnumDefinitionImpl):
    """
    Enumeration of zygosity types.
    """
    Homozygous = PermissibleValue(
        text="Homozygous",
        description="Homozygous site.",
        meaning=GENO["0000136"])
    Heterozygous = PermissibleValue(
        text="Heterozygous",
        description="Heterozygous site.",
        meaning=GENO["0000135"])
    Disomic = PermissibleValue(
        text="Disomic",
        description="Disomic site.",
        meaning=GENO["0000391"])

    _defn = EnumDefinition(
        name="ZygosityType",
        description="Enumeration of zygosity types.",
    )

# Slots
class slots:
    pass

slots.sequenceAlteration__has_identifier = Slot(uri=SIO['000300'], name="sequenceAlteration__has_identifier", curie=SIO.curie('000300'),
                   model_uri=DEFAULT_.sequenceAlteration__has_identifier, domain=None, range=str)

slots.sequenceAlteration__has_reference_allele = Slot(uri=GENO['0000385'], name="sequenceAlteration__has_reference_allele", curie=GENO.curie('0000385'),
                   model_uri=DEFAULT_.sequenceAlteration__has_reference_allele, domain=None, range=Union[dict, ReferenceAllele])

slots.sequenceAlteration__has_alternate_allele = Slot(uri=GENO['0000382'], name="sequenceAlteration__has_alternate_allele", curie=GENO.curie('0000382'),
                   model_uri=DEFAULT_.sequenceAlteration__has_alternate_allele, domain=None, range=Union[dict, AlternateAllele])

slots.sequenceAlteration__has_location = Slot(uri=FALDO.location, name="sequenceAlteration__has_location", curie=FALDO.curie('location'),
                   model_uri=DEFAULT_.sequenceAlteration__has_location, domain=None, range=Union[dict, VariationSite])

slots.sequenceAlteration__is_associated_with = Slot(uri=SIO['001403'], name="sequenceAlteration__is_associated_with", curie=SIO.curie('001403'),
                   model_uri=DEFAULT_.sequenceAlteration__is_associated_with, domain=None, range=Union[Union[dict, AssociatedCharacteristic], List[Union[dict, AssociatedCharacteristic]]])

slots.referenceAllele__value = Slot(uri=SIO['000300'], name="referenceAllele__value", curie=SIO.curie('000300'),
                   model_uri=DEFAULT_.referenceAllele__value, domain=None, range=str)

slots.alternateAllele__value = Slot(uri=SIO['000300'], name="alternateAllele__value", curie=SIO.curie('000300'),
                   model_uri=DEFAULT_.alternateAllele__value, domain=None, range=str)

slots.variationSite__begins_at = Slot(uri=FALDO.begin, name="variationSite__begins_at", curie=FALDO.curie('begin'),
                   model_uri=DEFAULT_.variationSite__begins_at, domain=None, range=Union[dict, VariationSiteBegin])

slots.variationSite__ends_at = Slot(uri=FALDO.end, name="variationSite__ends_at", curie=FALDO.curie('end'),
                   model_uri=DEFAULT_.variationSite__ends_at, domain=None, range=Union[dict, VariationSiteEnd])

slots.variationSite__has_reference = Slot(uri=FALDO.reference, name="variationSite__has_reference", curie=FALDO.curie('reference'),
                   model_uri=DEFAULT_.variationSite__has_reference, domain=None, range=Union[dict, VariationSiteReference])

slots.variationSiteBegin__position = Slot(uri=FALDO.position, name="variationSiteBegin__position", curie=FALDO.curie('position'),
                   model_uri=DEFAULT_.variationSiteBegin__position, domain=None, range=int)

slots.variationSiteEnd__position = Slot(uri=FALDO.position, name="variationSiteEnd__position", curie=FALDO.curie('position'),
                   model_uri=DEFAULT_.variationSiteEnd__position, domain=None, range=int)

slots.variationSiteReference__value = Slot(uri=DEFAULT_.value, name="variationSiteReference__value", curie=DEFAULT_.curie('value'),
                   model_uri=DEFAULT_.variationSiteReference__value, domain=None, range=str)

slots.variationSiteReference__label = Slot(uri=RDFS.label, name="variationSiteReference__label", curie=RDFS.curie('label'),
                   model_uri=DEFAULT_.variationSiteReference__label, domain=None, range=str)

slots.variationSiteReference__sameAs = Slot(uri=OWL.sameAs, name="variationSiteReference__sameAs", curie=OWL.curie('sameAs'),
                   model_uri=DEFAULT_.variationSiteReference__sameAs, domain=None, range=str)

slots.associatedCharacteristic__has_identifier = Slot(uri=SIO['000300'], name="associatedCharacteristic__has_identifier", curie=SIO.curie('000300'),
                   model_uri=DEFAULT_.associatedCharacteristic__has_identifier, domain=None, range=Optional[str])

slots.associatedCharacteristic__label = Slot(uri=RDFS.label, name="associatedCharacteristic__label", curie=RDFS.curie('label'),
                   model_uri=DEFAULT_.associatedCharacteristic__label, domain=None, range=Optional[str])

slots.associatedCharacteristic__has_zygosity = Slot(uri=GENO['0000608'], name="associatedCharacteristic__has_zygosity", curie=GENO.curie('0000608'),
                   model_uri=DEFAULT_.associatedCharacteristic__has_zygosity, domain=None, range=Optional[Union[dict, Zygosity]])

slots.zygosity__type = Slot(uri=RDFS.type, name="zygosity__type", curie=RDFS.curie('type'),
                   model_uri=DEFAULT_.zygosity__type, domain=None, range=Union[str, "ZygosityType"])

slots.zygosity__has_measurement_value = Slot(uri=SIO['000216'], name="zygosity__has_measurement_value", curie=SIO.curie('000216'),
                   model_uri=DEFAULT_.zygosity__has_measurement_value, domain=None, range=Union[dict, Count])

slots.zygosity__has_frequency = Slot(uri=SIO['000900'], name="zygosity__has_frequency", curie=SIO.curie('000900'),
                   model_uri=DEFAULT_.zygosity__has_frequency, domain=None, range=Union[dict, Frequency])

slots.frequency__value = Slot(uri=SIO['000300'], name="frequency__value", curie=SIO.curie('000300'),
                   model_uri=DEFAULT_.frequency__value, domain=None, range=float)

slots.count__value = Slot(uri=SIO['000300'], name="count__value", curie=SIO.curie('000300'),
                   model_uri=DEFAULT_.count__value, domain=None, range=int)

slots.container__variantcatalog = Slot(uri=DEFAULT_.variantcatalog, name="container__variantcatalog", curie=DEFAULT_.curie('variantcatalog'),
                   model_uri=DEFAULT_.container__variantcatalog, domain=None, range=Optional[Union[Union[dict, SequenceAlteration], List[Union[dict, SequenceAlteration]]]])