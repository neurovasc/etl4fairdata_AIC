-- # Class: "SequenceAlteration" Description: "A representation of a sequence alteration (so:0001059)."
--     * Slot: id Description: 
--     * Slot: has_identifier Description: A unique identifier for the sequence alteration.
--     * Slot: Container_id Description: Autocreated FK slot
--     * Slot: has_reference_allele_id Description: Links the sequence alteration to its reference allele (geno:0000036) using the property geno:0000385.
--     * Slot: has_alternate_allele_id Description: Links the sequence alteration to its alternate allele (geno:0000002) using the property geno:0000382.
--     * Slot: has_location_id Description: Links the sequence alteration to its location using the property faldo:Region.
-- # Class: "ReferenceAllele" Description: "Represents the reference allele (geno:0000036)."
--     * Slot: id Description: 
--     * Slot: value Description: The value of the reference allele.
-- # Class: "AlternateAllele" Description: "Represents the alternate allele (geno:0000002)."
--     * Slot: id Description: 
--     * Slot: value Description: The value of the alternate allele.
-- # Class: "VariationSite" Description: "Represents the location of a sequence alteration."
--     * Slot: id Description: 
--     * Slot: begins_at_id Description: The beginning of the location of the sequence alteration.
--     * Slot: ends_at_id Description: The end of the location of the sequence alteration.
--     * Slot: has_reference_id Description: The reference sequence (contig, sequence, chromosome).
-- # Class: "VariationSiteBegin" Description: "Represents the beginning of the location of a sequence alteration."
--     * Slot: id Description: 
--     * Slot: position Description: The position of the beginning of the location of the sequence alteration.
-- # Class: "VariationSiteEnd" Description: "Represents the end of the location of a sequence alteration."
--     * Slot: id Description: 
--     * Slot: position Description: The position of the end of the location of the sequence alteration.
-- # Class: "VariationSiteReference" Description: "Represents the reference sequence (contig, sequence, chromosome)."
--     * Slot: id Description: 
--     * Slot: value Description: The value of the reference sequence (ena sequence i.e. https://www.ebi.ac.uk/ena/browser/view/CM000684.2).
--     * Slot: label Description: A human-readable label for the reference sequence.
--     * Slot: sameAs Description: The value of the reference sequence in another database (ncbi sequence i.e. https://www.ncbi.nlm.nih.gov/nuccore/NC_000022.11).
-- # Class: "AssociatedCharacteristic" Description: "Represents a phenotype, clinical trait or grouping associated with a sequence alteration."
--     * Slot: id Description: 
--     * Slot: has_identifier Description: A unique identifier for the associated characteristic.
--     * Slot: label Description: A human-readable label for the associated characteristic.
--     * Slot: has_zygosity_id Description: The zygocity of the associated characteristic.
-- # Class: "Zygosity" Description: "Represents the zygosity of an associated characteristic."
--     * Slot: id Description: 
--     * Slot: type Description: The type of zygosity.
--     * Slot: has_measurement_value_id Description: Count of alternative alleles in individuals with a certain zygosity.
--     * Slot: has_frequency_id Description: Frequency of the allele in indibiduals with a certain zygosity.
-- # Class: "Frequency" Description: "Represents the frequency of an allele in individuals with a certain zygosity."
--     * Slot: id Description: 
--     * Slot: value Description: The frequency value.
-- # Class: "Count" Description: "Represents the count of alternative alleles in individuals with a certain zygosity."
--     * Slot: id Description: 
--     * Slot: value Description: The count value.
-- # Class: "Container" Description: ""
--     * Slot: id Description: 
-- # Class: "SequenceAlteration_is_associated_with" Description: ""
--     * Slot: SequenceAlteration_id Description: Autocreated FK slot
--     * Slot: is_associated_with_id Description: Links the sequence alteration to a phenotype (female), clinical trait (diabetes) or grouping (whole cohort).

CREATE TABLE "ReferenceAllele" (
	id INTEGER NOT NULL, 
	value TEXT NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "AlternateAllele" (
	id INTEGER NOT NULL, 
	value TEXT NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "VariationSiteBegin" (
	id INTEGER NOT NULL, 
	position INTEGER NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "VariationSiteEnd" (
	id INTEGER NOT NULL, 
	position INTEGER NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "VariationSiteReference" (
	id INTEGER NOT NULL, 
	value TEXT NOT NULL, 
	label TEXT NOT NULL, 
	"sameAs" TEXT NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "Frequency" (
	id INTEGER NOT NULL, 
	value FLOAT NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "Count" (
	id INTEGER NOT NULL, 
	value INTEGER NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "Container" (
	id INTEGER NOT NULL, 
	PRIMARY KEY (id)
);
CREATE TABLE "VariationSite" (
	id INTEGER NOT NULL, 
	begins_at_id INTEGER NOT NULL, 
	ends_at_id INTEGER NOT NULL, 
	has_reference_id INTEGER NOT NULL, 
	PRIMARY KEY (id), 
	FOREIGN KEY(begins_at_id) REFERENCES "VariationSiteBegin" (id), 
	FOREIGN KEY(ends_at_id) REFERENCES "VariationSiteEnd" (id), 
	FOREIGN KEY(has_reference_id) REFERENCES "VariationSiteReference" (id)
);
CREATE TABLE "Zygosity" (
	id INTEGER NOT NULL, 
	type VARCHAR(12) NOT NULL, 
	has_measurement_value_id INTEGER NOT NULL, 
	has_frequency_id INTEGER NOT NULL, 
	PRIMARY KEY (id), 
	FOREIGN KEY(has_measurement_value_id) REFERENCES "Count" (id), 
	FOREIGN KEY(has_frequency_id) REFERENCES "Frequency" (id)
);
CREATE TABLE "SequenceAlteration" (
	id INTEGER NOT NULL, 
	has_identifier TEXT NOT NULL, 
	"Container_id" INTEGER, 
	has_reference_allele_id INTEGER NOT NULL, 
	has_alternate_allele_id INTEGER NOT NULL, 
	has_location_id INTEGER NOT NULL, 
	PRIMARY KEY (id), 
	FOREIGN KEY("Container_id") REFERENCES "Container" (id), 
	FOREIGN KEY(has_reference_allele_id) REFERENCES "ReferenceAllele" (id), 
	FOREIGN KEY(has_alternate_allele_id) REFERENCES "AlternateAllele" (id), 
	FOREIGN KEY(has_location_id) REFERENCES "VariationSite" (id)
);
CREATE TABLE "AssociatedCharacteristic" (
	id INTEGER NOT NULL, 
	has_identifier TEXT, 
	label TEXT, 
	has_zygosity_id INTEGER, 
	PRIMARY KEY (id), 
	FOREIGN KEY(has_zygosity_id) REFERENCES "Zygosity" (id)
);
CREATE TABLE "SequenceAlteration_is_associated_with" (
	"SequenceAlteration_id" INTEGER, 
	is_associated_with_id INTEGER NOT NULL, 
	PRIMARY KEY ("SequenceAlteration_id", is_associated_with_id), 
	FOREIGN KEY("SequenceAlteration_id") REFERENCES "SequenceAlteration" (id), 
	FOREIGN KEY(is_associated_with_id) REFERENCES "AssociatedCharacteristic" (id)
);