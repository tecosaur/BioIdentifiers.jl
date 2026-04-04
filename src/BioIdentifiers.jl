# SPDX-FileCopyrightText: © 2025 TEC <contact@tecosaur.net>
# SPDX-License-Identifier: MPL-2.0

"""
    BioIdentifiers

Structured and validated types for biological and biomedical identifiers.

This module provides type-safe representations of identifiers commonly used in
biological and biomedical databases, extending the FastIdentifiers framework.
Each identifier type validates input format, provides canonical string
representations, and includes persistent URLs where applicable.

# Supported Identifiers

**Protein & Structure**
- [`AFDB`](@ref): AlphaFold Database protein structure predictions
- [`PDB`](@ref): Protein Data Bank macromolecular structures
- [`UniProt`](@ref): Universal Protein Resource sequences
- [`UniRef`](@ref): UniProt Reference Clusters
- [`IntAct`](@ref): Molecular interaction database
- [`InterPro`](@ref): Protein family/domain classification
- [`Pfam`](@ref): Protein family HMM profiles
- [`PXD`](@ref): ProteomeXchange dataset identifiers

**Genomics & Genetics**
- [`ENSG`](@ref EnsemblIdentifier), [`ENST`](@ref EnsemblIdentifier), [`ENSP`](@ref EnsemblIdentifier), [`ENSE`](@ref EnsemblIdentifier), [`ENSR`](@ref EnsemblIdentifier), [`ENSF`](@ref EnsemblIdentifier), [`ENSFM`](@ref EnsemblIdentifier): Ensembl genome annotations
- [`NCBIGene`](@ref): NCBI Entrez Gene database
- [`RefSeq`](@ref): NCBI Reference Sequence database
- [`HGNC`](@ref): HUGO Gene Nomenclature Committee
- [`INSDC`](@ref): International Nucleotide Sequence Database Collaboration
- [`OMIM`](@ref): Online Mendelian Inheritance in Man

**Variation & Clinical**
- [`CA`](@ref): ClinGen Allele Registry
- [`ClinVar`](@ref): Clinical genomic variant database
- [`dbSNP`](@ref): Single Nucleotide Polymorphism database
- [`dbVar`](@ref): Genomic structural variation database
- [`GWAS`](@ref): NHGRI-EBI GWAS Catalog studies

**Expression & Functional Genomics**
- [`ArrayExpress`](@ref): Functional genomics data archive
- [`GEO`](@ref): Gene Expression Omnibus

**Studies & Samples**
- [`BioProject`](@ref): Biological research projects
- [`BioSample`](@ref): Biological sample metadata
- [`ClinicalTrials`](@ref): ClinicalTrials.gov registry
- [`dbGaP`](@ref): Database of Genotypes and Phenotypes
- [`EGA`](@ref): European Genome-phenome Archive
- [`RRID`](@ref): Research Resource Identifiers
- [`SRA`](@ref): Sequence Read Archive

**Ontologies & Controlled Vocabularies**
- [`CL`](@ref): Cell Ontology
- [`DOID`](@ref): Disease Ontology
- [`ECO`](@ref): Evidence & Conclusion Ontology
- [`EFO`](@ref): Experimental Factor Ontology
- [`GO`](@ref): Gene Ontology
- [`HPO`](@ref): Human Phenotype Ontology
- [`MeSH`](@ref): Medical Subject Headings
- [`MONDO`](@ref): Monarch Disease Ontology
- [`MP`](@ref): Mammalian Phenotype Ontology
- [`PATO`](@ref): Phenotype And Trait Ontology
- [`SO`](@ref): Sequence Ontology
- [`UBERON`](@ref): Uber-anatomy Ontology

**Chemical & Metabolic**
- [`ChEBI`](@ref): Chemical Entities of Biological Interest
- [`ChEMBL`](@ref): Bioactive compound database
- [`DrugBank`](@ref): Drug and pharmaceutical database
- [`HMDB`](@ref): Human Metabolome Database
- [`KEGG`](@ref): Kyoto Encyclopedia of Genes and Genomes
- [`MetaboLights`](@ref): Metabolomics study archive
- [`PubChem`](@ref): Chemical compound/substance/assay database

**Networks & Interactions**
- [`BioGRID`](@ref): Biological interaction datasets
- [`Reactome`](@ref): Curated biological pathways
- [`WikiPathways`](@ref): Community pathway database

**Cell Lines & Model Organisms**
- [`Cellosaurus`](@ref): Cell line registry
- [`FlyBase`](@ref): Drosophila gene database
- [`MGI`](@ref): Mouse Genome Informatics
- [`NCBITaxon`](@ref): NCBI Taxonomy database
- [`SGD`](@ref): Saccharomyces Genome Database
- [`WormBase`](@ref): Caenorhabditis gene database
"""
module BioIdentifiers

using FastIdentifiers
import FastIdentifiers: purlprefix, shortcode
using FastIdentifiers: @defid

FastIdentifiers.@reexport

export BioIdentifier, AFDB, ArrayExpress, BioGRID, BioProject, BioSample, CA,
    CL, Cellosaurus, ChEBI, ChEMBL, ClinicalTrials, ClinVar, dbGaP, dbSNP,
    dbVar, DOID, DrugBank, ECO, EFO, EGA, ENSG, ENST, ENSP, ENSE, ENSR, ENSF,
    ENSFM, FlyBase, GEO, GO, GWAS, HGNC, HMDB, HPO, INSDC, IntAct, InterPro,
    KEGG, MeSH, MetaboLights, MGI, MONDO, MP, NCBIGene, NCBITaxon, OMIM, PATO,
    PDB, Pfam, PubChem, PXD, Reactome, RefSeq, RRID, SGD, SO, SRA, UBERON,
    UniProt, UniRef, WikiPathways, WormBase

"""
    BioIdentifier <: AbstractIdentifier

An abstract type representing a biological or biomedical identifier.

See also: `AbstractIdentifier`.
"""
abstract type BioIdentifier <: AbstractIdentifier end

# Local helpers for hand-written types

"""
    chopprefixes(str, prefixes...) -> (matched::Bool, rest::SubString)

Case-insensitive prefix stripping.  All prefixes must be lowercase ASCII.
"""
function chopprefixes(str::AbstractString, prefixes::AbstractString...)
    s = SubString(string(str))
    matched = false
    for prefix in prefixes
        ncodeunits(s) >= ncodeunits(prefix) || continue
        ok = true
        for i in 1:ncodeunits(prefix)
            (codeunit(s, i) | 0x20) == codeunit(prefix, i) || (ok = false; break)
        end
        ok || continue
        s = SubString(s, ncodeunits(prefix) + 1)
        matched = true
    end
    matched, s
end

# ──────────────────────────────────────────────────────────────────
# @defid types
# ──────────────────────────────────────────────────────────────────

# ArrayExpress

"""
    ArrayExpress <: BioIdentifier

ArrayExpress functional genomics experiment identifier.

ArrayExpress is EBI's functional genomics archive storing gene expression data.
Identifiers can be either native EBI experiments (`E-<platform>-<number>`) or
imported GEO experiments (`E-GEOD-<number>`). Parsing may throw a `MalformedIdentifier`
if the format is invalid.

# Examples

```julia
julia> parse(ArrayExpress, "E-MTAB-1234")
ArrayExpress:E-MTAB-1234

julia> parse(ArrayExpress, "E-GEOD-5678")
ArrayExpress:E-GEOD-5678
```
"""
@defid(ArrayExpress <: BioIdentifier,
       ("E-", :platform(letters(4:5, upper=true)), "-", :id(digits(UInt32))),
       prefix="", purlprefix="https://www.ebi.ac.uk/biostudies/arrayexpress/studies/")

# BioGRID

"""
    BioGRID <: BioIdentifier

A [BioGRID](https://thebiogrid.org/) identifier for molecular interaction records.

The Biological General Repository for Interaction Datasets (BioGRID) curates
genetic and protein interactions from the primary literature across all major
model organisms. Identifiers are positive integers.

# Examples

```julia
julia> parse(BioGRID, "123456")
BioGRID:123456

julia> parse(BioGRID, "https://thebiogrid.org/789012")
BioGRID:789012
```
"""
@defid(BioGRID <: BioIdentifier,
       :id(digits(UInt32)),
       prefix="", purlprefix="https://thebiogrid.org/")

# BioProject

"""
    BioProject <: BioIdentifier

BioProject umbrella project identifier for multi-omics submissions.

BioProject provides a single place to find links to diverse data types that may be
scattered across multiple databases. Identifiers follow format PRJNA/PRJEB/PRJDB<number>
for NCBI/ENA/DDBJ respectively. Parsing may throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(BioProject, "PRJNA123456")
BioProject:PRJNA123456

julia> parse(BioProject, "PRJEB789012")
BioProject:PRJEB789012
```
"""
@defid(BioProject <: BioIdentifier,
       (:prefix(choice("PRJNA", "PRJEB", "PRJDB")), :id(digits(UInt32))),
       prefix="", purlprefix="https://www.ncbi.nlm.nih.gov/bioproject/")

# BioSample

"""
    BioSample <: BioIdentifier

BioSample unique accession for biological samples.

BioSample database contains descriptions of biological source materials used in
experimental assays. Identifiers follow format SAMN/SAMEA/SAMD<number> for
NCBI/EBI/DDBJ respectively. Parsing may throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(BioSample, "SAMN12345678")
BioSample:SAMN12345678

julia> parse(BioSample, "SAMEA9876543")
BioSample:SAMEA9876543
```
"""
@defid(BioSample <: BioIdentifier,
       (:prefix(choice("SAMN", "SAMEA", "SAMD")), :id(digits(UInt32))),
       prefix="", purlprefix="https://www.ncbi.nlm.nih.gov/biosample/")

# CA

"""
    CA <: BioIdentifier

ClinGen Allele Registry globally stable allele/variant key.

The ClinGen Allele Registry provides stable identifiers for variants allowing
for the aggregation of information about variants across multiple resources.
Identifiers follow format `CA<number>`. Parsing may throw a `MalformedIdentifier`
if the format is invalid.

# Examples

```julia
julia> parse(CA, "CA123456")
CA:123456

julia> parse(CA, "https://reg.clinicalgenome.org/allele/CA789012")
CA:789012
```
"""
@defid(CA <: BioIdentifier,
       ("CA", :id(digits(UInt32))),
       prefix="", purlprefix="https://reg.clinicalgenome.org/allele/")

# CL

"""
    CL <: BioIdentifier

A [Cell Ontology](https://obofoundry.org/ontology/cl.html) term identifier for cross-species cell types.

Part of the [OBO Foundry](https://obofoundry.org/), the Cell Ontology represents
in vivo and in vitro cell types across animals, plants, fungi and other
organisms. It is used for cell type annotation in single-cell RNA-seq
experiments and phenotype ontologies. Identifiers follow the format `CL:`
followed by 7 zero-padded digits.

# Examples

```julia
julia> parse(CL, "CL:0000001")
CL:0000001

julia> parse(CL, "0000001")
CL:0000001
```
"""
@defid(CL <: BioIdentifier,
       :id(digits(7, pad=7)),
       prefix="CL:", purlprefix="https://purl.obolibrary.org/obo/CL_")

# Cellosaurus

"""
    Cellosaurus <: BioIdentifier

A [Cellosaurus](https://www.cellosaurus.org/) cell line registry identifier.

Maintained by the [SIB Swiss Institute of Bioinformatics](https://www.sib.swiss/),
Cellosaurus is the most comprehensive knowledge resource on cell lines,
providing information on origin, cross-contamination, synonyms, and
references. Identifiers follow the format `CVCL_` followed by 1–7 digits.

# Examples

```julia
julia> parse(Cellosaurus, "CVCL_0001")
Cellosaurus:CVCL_1

julia> parse(Cellosaurus, "CVCL_1234567")
Cellosaurus:CVCL_1234567
```
"""
@defid(Cellosaurus <: BioIdentifier,
       ("CVCL_", :id(digits(1:7))),
       prefix="", purlprefix="https://www.cellosaurus.org/")

# ChEBI

"""
    ChEBI <: BioIdentifier

A [ChEBI](https://www.ebi.ac.uk/chebi/) identifier for chemical entities of biological interest.

Maintained by the [European Bioinformatics Institute](https://www.ebi.ac.uk/)
(EMBL-EBI), ChEBI is a freely available ontology of molecular entities focused
on small chemical compounds. It provides systematic nomenclature, ontological
classification, and links to reactions and pathways. Identifiers follow the
format `CHEBI:` followed by a positive integer.

# Examples

```julia
julia> parse(ChEBI, "CHEBI:15377")
ChEBI:15377

julia> parse(ChEBI, "15377")
ChEBI:15377
```
"""
@defid(ChEBI <: BioIdentifier,
       :id(digits(UInt32)),
       prefix="CHEBI:", purlprefix="https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:")

# ChEMBL

"""
    ChEMBL <: BioIdentifier

A [ChEMBL](https://www.ebi.ac.uk/chembl/) identifier for bioactive compounds, targets, and assays.

Maintained by [EMBL-EBI](https://www.ebi.ac.uk/), ChEMBL is a manually curated
database of bioactive molecules with drug-like properties, widely used in drug
discovery and chemical biology. Identifiers follow the format `CHEMBL` followed
by a positive integer.

# Examples

```julia
julia> parse(ChEMBL, "CHEMBL25")
ChEMBL:25

julia> parse(ChEMBL, "CHEMBL1075272")
ChEMBL:1075272
```
"""
@defid(ChEMBL <: BioIdentifier,
       ("CHEMBL", :id(digits(UInt32))),
       prefix="", purlprefix="https://www.ebi.ac.uk/chembl/compound_report_card/")

# ClinVar

"""
    ClinVar <: BioIdentifier

ClinVar aggregated clinical variant archive identifier.

ClinVar aggregates information about genomic variation and its relationship to
human health. Identifiers follow formats VCV/RCV/SCV<number>.<version> for
variant/variant-condition/submitter records respectively. Parsing may throw a
`MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(ClinVar, "VCV000123456.1")
ClinVar:VCV000123456.1

julia> parse(ClinVar, "RCV000789012.2")
ClinVar:RCV000789012.2
```
"""
@defid(ClinVar <: BioIdentifier,
       (choice(seq("VCV", :id(digits(UInt32, pad=9))),
               seq("RCV", :id(digits(1:8, pad=8))),
               seq("SCV", :id(digits(1:8, pad=8)))),
        optional(".", :ver(digits(1:5)))),
       prefix="", purlprefix="https://www.ncbi.nlm.nih.gov/clinvar/")

# ClinicalTrials

"""
    ClinicalTrials <: BioIdentifier

A [ClinicalTrials.gov](https://clinicaltrials.gov/) registry identifier for clinical studies.

Maintained by the [U.S. National Library of Medicine](https://www.nlm.nih.gov/),
ClinicalTrials.gov is the largest registry of clinical studies worldwide,
covering interventional and observational studies across all diseases and
conditions. Identifiers follow the format `NCT` followed by 8 zero-padded
digits.

# Examples

```julia
julia> parse(ClinicalTrials, "NCT01234567")
ClinicalTrials:NCT01234567

julia> parse(ClinicalTrials, "https://clinicaltrials.gov/ct2/show/NCT01234567")
ClinicalTrials:NCT01234567
```
"""
@defid(ClinicalTrials <: BioIdentifier,
       ("NCT", :id(digits(8, pad=8))),
       prefix="", purlprefix="https://clinicaltrials.gov/ct2/show/")

# dbGaP

"""
    dbGaP <: BioIdentifier

dbGaP genotype-phenotype study archive identifier.

The database of Genotypes and Phenotypes (dbGaP) archives and distributes the data
and results from studies investigating the interaction of genotype and phenotype.
Identifiers follow format `phs<6 digits>.v<version>.p<participant_set>`. Parsing may
throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(dbGaP, "phs000001.v1.p1")
dbGaP:phs000001.v1.p1

julia> parse(dbGaP, "phs123456.v2.p3")
dbGaP:phs123456.v2.p3
```
"""
@defid(dbGaP <: BioIdentifier,
       ("phs", :study(digits(6, pad=6)), ".v", :ver(digits(1:2)),
        ".p", :participant(digits(1:2))),
       prefix="", purlprefix="https://www.ncbi.nlm.nih.gov/gap/")

# dbSNP

"""
    dbSNP <: BioIdentifier

dbSNP reference SNP cluster identifier.

The Single Nucleotide Polymorphism database (dbSNP) contains human single nucleotide
variations, microsatellites, and small-scale insertions and deletions. Identifiers
follow format `rs<number>`. Parsing may throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(dbSNP, "rs123456")
dbSNP:rs123456

julia> parse(dbSNP, "https://www.ncbi.nlm.nih.gov/snp/rs789012")
dbSNP:rs789012
```
"""
@defid(dbSNP <: BioIdentifier,
       ("rs", :id(digits(UInt32))),
       prefix="", purlprefix="https://www.ncbi.nlm.nih.gov/snp/")

# dbVar

"""
    dbVar <: BioIdentifier

dbVar structural variant archive identifier.

The database of genomic structural variation (dbVar) contains insertions, deletions,
duplications, inversions, mobile elements, translocations, and complex chromosomal
rearrangements. Identifiers follow formats nsv/nssv/nstd<number> for study variant/
submitted variant/study records respectively. Parsing may throw a `MalformedIdentifier`
if the format is invalid.

# Examples

```julia
julia> parse(dbVar, "nsv123456")
dbVar:nsv123456

julia> parse(dbVar, "nssv789012")
dbVar:nssv789012
```
"""
@defid(dbVar <: BioIdentifier,
       (:prefix(choice("nsv", "nssv", "nstd")), :id(digits(UInt32))),
       prefix="", purlprefix="https://www.ncbi.nlm.nih.gov/dbvar/")

# DOID

"""
    DOID <: BioIdentifier

A [Disease Ontology](https://disease-ontology.org/) (DO) identifier for human diseases.

Part of the [OBO Foundry](https://obofoundry.org/), the Human Disease Ontology
provides standardised, reusable descriptions of human disease terms that
integrate medical vocabularies (MeSH, ICD, SNOMED CT). Identifiers follow
the format `DOID:` followed by 7 zero-padded digits.

# Examples

```julia
julia> parse(DOID, "DOID:0001816")
DOID:0001816

julia> parse(DOID, "0001816")
DOID:0001816
```
"""
@defid(DOID <: BioIdentifier,
       :id(digits(7, pad=7)),
       prefix="DOID:", purlprefix="https://purl.obolibrary.org/obo/DOID_")

# ECO

"""
    ECO <: BioIdentifier

Evidence & Conclusion Ontology identifier.

The Evidence & Conclusion Ontology (ECO) describes types of scientific evidence
within the realm of biological research. ECO terms are used in UniProt, GO,
and ClinVar annotations to describe the evidence supporting biological
assertions. Parsing may throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(ECO, "ECO:0000313")
ECO:0000313

julia> parse(ECO, "0000269")
ECO:0000269
```
"""
@defid(ECO <: BioIdentifier,
       :id(digits(7, pad=7)),
       prefix="ECO:", purlprefix="https://purl.obolibrary.org/obo/ECO_")

# EFO

"""
    EFO <: BioIdentifier

Experimental Factor Ontology identifier for traits and experimental variables.

The Experimental Factor Ontology provides a systematic description of many experimental
variables available in EBI databases. These variables are used in genomic studies
including phenotypes, biological/experimental factors and measurement end points.
Identifiers follow format `EFO:<number>`. Parsing may throw a `MalformedIdentifier`
if the format is invalid.

# Examples

```julia
julia> parse(EFO, "EFO:0000001")
EFO:0000001

julia> parse(EFO, "0000001")
EFO:0000001
```
"""
@defid(EFO <: BioIdentifier,
       :id(digits(7, pad=7)),
       prefix="EFO:", purlprefix="https://www.ebi.ac.uk/ols/ontologies/efo/terms?short_form=EFO_")

# EGA

"""
    EGA <: BioIdentifier

An [EGA](https://ega-archive.org/) (European Genome-phenome Archive) controlled-access data identifier.

Operated by [EMBL-EBI](https://www.ebi.ac.uk/) and [CRG](https://www.crg.eu/),
the EGA provides permanent archiving and controlled sharing of personally
identifiable genetic and phenotypic data from biomedical research. Identifiers
use prefix `EGAD` (dataset) or `EGAS` (study) followed by 11 zero-padded
digits.

# Examples

```julia
julia> parse(EGA, "EGAD00001000001")
EGA:EGAD00001000001

julia> parse(EGA, "EGAS00001000002")
EGA:EGAS00001000002
```
"""
@defid(EGA <: BioIdentifier,
       (:prefix(choice("EGAD", "EGAS")), :id(digits(11, pad=11))),
       prefix="", purlprefix="https://ega-archive.org/")

# GEO

"""
    GEO <: BioIdentifier

NCBI Gene Expression Omnibus identifier.

The Gene Expression Omnibus is a public functional genomics data repository supporting
MIAME-compliant data submissions. Identifiers follow formats GPL/GSM/GSE/GDS<number>
for platform/sample/series/dataset records respectively. Parsing may throw a
`MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(GEO, "GSE123456")
GEO:GSE123456

julia> parse(GEO, "GPL570")
GEO:GPL570
```
"""
@defid(GEO <: BioIdentifier,
       (:prefix(choice("GPL", "GSM", "GSE", "GDS")), :id(digits(UInt32))),
       prefix="", purlprefix="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=")

# GO

"""
    GO <: BioIdentifier

A [Gene Ontology](http://geneontology.org/) (GO) term identifier.

The Gene Ontology, maintained by the [GO Consortium](http://geneontology.org/),
is the most widely used structured vocabulary in genomics. It covers three
domains — cellular component, molecular function, and biological process — and
is used for functional annotation across all sequenced organisms. Identifiers
follow the format `GO:` followed by 7 zero-padded digits.

# Examples

```julia
julia> parse(GO, "GO:0006915")
GO:0006915

julia> parse(GO, "0006915")
GO:0006915
```
"""
@defid(GO <: BioIdentifier,
       :id(digits(7, pad=7)),
       prefix="GO:", purlprefix="https://purl.obolibrary.org/obo/GO_")

# HGNC

"""
    HGNC <: BioIdentifier

HUGO Gene Nomenclature Committee authoritative human gene identifier.

The HGNC provides unique symbols and names for human loci, including protein coding
genes, ncRNA genes and pseudogenes, to allow unambiguous scientific communication.
Identifiers follow format `HGNC:<number>`. Parsing may throw a `MalformedIdentifier`
if the format is invalid.

# Examples

```julia
julia> parse(HGNC, "HGNC:5")
HGNC:5

julia> parse(HGNC, "5")
HGNC:5
```
"""
@defid(HGNC <: BioIdentifier,
       :id(digits(1:5)),
       prefix="HGNC:", purlprefix="https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:")

# HMDB

"""
    HMDB <: BioIdentifier

A [Human Metabolome Database](https://hmdb.ca/) (HMDB) entry identifier.

Maintained by the [Wishart Lab](http://www.wishartlab.com/) at the University
of Alberta, HMDB contains detailed information about small molecule metabolites
found in the human body, supporting metabolomics, clinical chemistry, and
biomarker discovery. Identifiers follow the format `HMDB` followed by 7
zero-padded digits.

# Examples

```julia
julia> parse(HMDB, "HMDB0000001")
HMDB:0000001

julia> parse(HMDB, "https://hmdb.ca/metabolites/HMDB0123456")
HMDB:0123456
```
"""
@defid(HMDB <: BioIdentifier,
       ("HMDB", :id(digits(7, pad=7))),
       prefix="", purlprefix="https://hmdb.ca/metabolites/")

# HPO

"""
    HPO <: BioIdentifier

A [Human Phenotype Ontology](https://hpo.jax.org/) (HPO) term identifier.

Developed by the [Monarch Initiative](https://monarchinitiative.org/) and
part of the [OBO Foundry](https://obofoundry.org/), HPO provides a
standardised vocabulary of phenotypic abnormalities observed in human disease.
It is the standard for phenotype-driven diagnosis, rare disease databases
(e.g. Orphanet, OMIM), and clinical genomics. Identifiers follow the format
`HP:` followed by 7 zero-padded digits.

# Examples

```julia
julia> parse(HPO, "HP:0000001")
HPO:0000001

julia> parse(HPO, "0000001")
HPO:0000001
```
"""
@defid(HPO <: BioIdentifier,
       :id(digits(7, pad=7)),
       prefix="HP:", purlprefix="https://purl.obolibrary.org/obo/HP_")

# IntAct

"""
    IntAct <: BioIdentifier

IntAct protein interaction record identifier (EBI).

IntAct provides a freely available, open source database system and analysis tools
for molecular interaction data. All interactions are derived from literature curation
or direct user submissions. Identifiers follow format `EBI-<number>`. Parsing may
throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(IntAct, "EBI-123456")
IntAct:EBI-123456

julia> parse(IntAct, "https://www.ebi.ac.uk/intact/interaction/EBI-789012")
IntAct:EBI-789012
```
"""
@defid(IntAct <: BioIdentifier,
       ("EBI-", :id(digits(UInt32))),
       prefix="", purlprefix="https://www.ebi.ac.uk/intact/interaction/")

# KEGG

"""
    KEGG <: BioIdentifier

A [KEGG](https://www.kegg.jp/) (Kyoto Encyclopedia of Genes and Genomes) pathway map identifier.

KEGG is an integrated database resource for systems-level understanding of
biological processes, maintained by [Kanehisa Laboratories](https://www.kanehisa.jp/).
Pathway map identifiers consist of a 3–4 letter organism code (e.g. `hsa` for
human, `mmu` for mouse) followed by 5 zero-padded digits.

# Examples

```julia
julia> parse(KEGG, "hsa00010")
KEGG:hsa00010

julia> parse(KEGG, "https://www.kegg.jp/entry/mmu04210")
KEGG:mmu04210
```
"""
@defid(KEGG <: BioIdentifier,
       (:species(letters(3:4, lower=true)), :id(digits(5, pad=5))),
       prefix="", purlprefix="https://www.kegg.jp/entry/")

# MetaboLights

"""
    MetaboLights <: BioIdentifier

A [MetaboLights](https://www.ebi.ac.uk/metabolights/) metabolomics study archive identifier.

Hosted by [EMBL-EBI](https://www.ebi.ac.uk/), MetaboLights is a database for
metabolomics experiments and derived information, serving as a recommended
repository for metabolomics data submission by many journals. Identifiers
follow the format `MTBLS` followed by a positive integer.

# Examples

```julia
julia> parse(MetaboLights, "MTBLS1")
MetaboLights:MTBLS1

julia> parse(MetaboLights, "https://www.ebi.ac.uk/metabolights/MTBLS123")
MetaboLights:MTBLS123
```
"""
@defid(MetaboLights <: BioIdentifier,
       ("MTBLS", :id(digits(UInt32))),
       prefix="", purlprefix="https://www.ebi.ac.uk/metabolights/")

# MONDO

"""
    MONDO <: BioIdentifier

A [MONDO](https://mondo.monarchinitiative.org/) (Monarch Disease Ontology) unified disease identifier.

Developed by the [Monarch Initiative](https://monarchinitiative.org/) and part
of the [OBO Foundry](https://obofoundry.org/), MONDO harmonises disease
definitions by merging multiple disease resources (OMIM, Orphanet, EFO, DOID,
NCIt) into a single coherent ontology. Identifiers follow the format `MONDO:`
followed by 7 zero-padded digits.

# Examples

```julia
julia> parse(MONDO, "MONDO:0000001")
MONDO:0000001

julia> parse(MONDO, "0000001")
MONDO:0000001
```
"""
@defid(MONDO <: BioIdentifier,
       :id(digits(7, pad=7)),
       prefix="MONDO:", purlprefix="https://purl.obolibrary.org/obo/MONDO_")

# NCBIGene

"""
    NCBIGene <: BioIdentifier

NCBI Entrez Gene numeric primary key.

Entrez Gene provides a unified query environment for genes defined by sequence
and/or in Map Viewer. Records include nomenclature, Reference Sequences,
maps, pathways, variations, phenotypes, and links to genome, phenotype, and
locus-specific resources. Identifiers are positive integers (≤ 10 digits).
Parsing may throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(NCBIGene, "7157")
NCBIGene:7157

julia> parse(NCBIGene, "https://www.ncbi.nlm.nih.gov/gene/5594")
NCBIGene:5594
```
"""
@defid(NCBIGene <: BioIdentifier,
       :id(digits(UInt32)),
       prefix="", purlprefix="https://www.ncbi.nlm.nih.gov/gene/")

# NCBITaxon

"""
    NCBITaxon <: BioIdentifier

NCBI Taxonomy species/strain identifier.

The NCBI Taxonomy database contains the names of all organisms that are represented
in the genetic databases with at least one nucleotide or protein sequence. The
taxonomy database includes phylogenetic and taxonomic information for the organisms.
Identifiers are positive integers. Parsing may throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(NCBITaxon, "9606")
NCBITaxon:9606

julia> parse(NCBITaxon, "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=10090")
NCBITaxon:10090
```
"""
@defid(NCBITaxon <: BioIdentifier,
       :id(digits(UInt32)),
       prefix="", purlprefix="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=")

# OMIM

"""
    OMIM <: BioIdentifier

An [OMIM](https://omim.org/) (Online Mendelian Inheritance in Man) entry identifier for genetic disorders and genes.

Maintained by [Johns Hopkins University](https://www.jhu.edu/), OMIM is a
continuously updated catalog of human genes and genetic disorders, with
particular focus on the molecular relationship between genetic variation and
phenotypic expression. Identifiers consist of 6 digits, optionally followed
by a `.` and a 4-digit allelic variant suffix.

# Examples

```julia
julia> parse(OMIM, "600918")
OMIM:600918

julia> parse(OMIM, "613659.0001")
OMIM:613659.0001
```
"""
@defid(OMIM <: BioIdentifier,
       (:id(digits(6, pad=6)), optional(".", :suffix(digits(4, pad=4)))),
       prefix="", purlprefix="https://omim.org/entry/")

# PDB

"""
    PDB <: BioIdentifier

A [Protein Data Bank](https://www.rcsb.org/) (PDB) structure identifier.

Managed by the [Worldwide Protein Data Bank](https://www.wwpdb.org/) (wwPDB),
the PDB is the single global archive of experimentally determined 3D structures
of biological macromolecules. Identifiers are 4-character alphanumeric codes
(digits and uppercase letters), with the first character always a digit.

# Examples

```julia
julia> parse(PDB, "1A3N")
PDB:1A3N

julia> parse(PDB, "7QDH")
PDB:7QDH
```
"""
@defid(PDB <: BioIdentifier,
       :id(alphnum(4, upper=true)),
       prefix="", purlprefix="https://www.ebi.ac.uk/pdbe/entry/pdb/")

# PubChem

"""
    PubChem <: BioIdentifier

A [PubChem](https://pubchem.ncbi.nlm.nih.gov/) identifier for chemical compounds, substances, and bioassays.

Operated by [NCBI](https://www.ncbi.nlm.nih.gov/), PubChem is the world's
largest freely accessible chemical database, providing information on chemical
structures, properties, biological activities, and safety data. Identifiers
use a prefix indicating the record type — `CID` (compound), `SID`
(substance), or `AID` (bioassay) followed by a positive integer.

# Examples

```julia
julia> parse(PubChem, "CID2244")
PubChem:CID2244

julia> parse(PubChem, "SID12345678")
PubChem:SID12345678
```
"""
@defid(PubChem <: BioIdentifier,
       choice(seq("CID", :id(digits(UInt32))),
              seq("SID", :id(digits(UInt32))),
              seq("AID", :id(digits(UInt32)))))

function FastIdentifiers.purl(pc::PubChem)
    sc = shortcode(pc)
    type = if startswith(sc, "SID")
        "substance"
    elseif startswith(sc, "AID")
        "assay"
    else
        "compound"
    end
    "https://pubchem.ncbi.nlm.nih.gov/" * type * "/" * sc[4:end]
end

# PXD

"""
    PXD <: BioIdentifier

A [ProteomeXchange](http://www.proteomexchange.org/) dataset identifier for mass spectrometry proteomics data.

ProteomeXchange is a consortium providing a single point of submission for
proteomics datasets to partner repositories (PRIDE, PeptideAtlas, MassIVE,
jPOST). Both `PXD` (original) and `RPXD` (reprocessed) prefixes are
accepted during parsing; the canonical output uses `PXD`. Identifiers
consist of the prefix followed by 6 zero-padded digits.

# Examples

```julia
julia> parse(PXD, "PXD000001")
PXD:PXD000001

julia> parse(PXD, "RPXD123456")
PXD:PXD123456
```
"""
@defid(PXD <: BioIdentifier,
       (skip(choice("PXD", "RPXD"), print="PXD"), :id(digits(6, pad=6))),
       prefix="", purlprefix="https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=")

# Reactome

"""
    Reactome <: BioIdentifier

A [Reactome](https://reactome.org/) curated biological pathway identifier.

Reactome is a free, open-source, peer-reviewed pathway knowledgebase maintained
by [OICR](https://oicr.on.ca/), [NYU Langone](https://nyulangone.org/), and
[EMBL-EBI](https://www.ebi.ac.uk/). Identifiers follow the format
`R-<species>-<number>`, where the 3-letter species code indicates the organism
(e.g. `HSA` for human, `MMU` for mouse).

# Examples

```julia
julia> parse(Reactome, "R-HSA-199420")
Reactome:R-HSA-199420

julia> parse(Reactome, "https://reactome.org/PathwayBrowser/#/R-MMU-123456")
Reactome:R-MMU-123456
```
"""
@defid(Reactome <: BioIdentifier,
       ("R-", :species(letters(3, upper=true)), "-", :id(digits(UInt32))),
       prefix="", purlprefix="https://reactome.org/PathwayBrowser/#/")

# RefSeq

"""
    RefSeq <: BioIdentifier

RefSeq curated reference sequence identifier.

The Reference Sequence database is a comprehensive, integrated, non-redundant,
well-annotated set of reference sequences including genomic, transcript, and
protein sequences. Identifiers follow formats like `NM_<number>.<version>`,
`NP_<number>.<version>`, etc. Parsing may throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(RefSeq, "NM_000546.6")
RefSeq:NM_000546.6

julia> parse(RefSeq, "NP_000537.3")
RefSeq:NP_000537.3
```
"""
@defid(RefSeq <: BioIdentifier,
       (:prefix(choice("NM", "NP", "XM", "XP", "NR", "XR", "NG", "NT", "NW", "NC")),
        "_", :id(digits(UInt32, pad=6)), ".", :ver(digits(1:2))),
       prefix="", purlprefix="https://www.ncbi.nlm.nih.gov/nuccore/")

# RRID

"""
    RRID <: BioIdentifier

A [Research Resource Identifier](https://scicrunch.org/resources) (RRID) for research materials and tools.

Managed by [SciCrunch](https://scicrunch.org/) and promoted by the
[RRID Initiative](https://www.rrids.org/), RRIDs provide persistent unique
identifiers for antibodies, cell lines, organisms, plasmids, software tools,
and other research resources. They are designed to be cited in the methods
sections of research articles for reproducibility. Identifiers follow the
format `RRID:<provider>_<id>`, where the provider code indicates the source
registry (e.g. `AB` for antibodies, `SCR` for software).

# Examples

```julia
julia> parse(RRID, "RRID:AB_123456")
RRID:AB_123456

julia> parse(RRID, "https://scicrunch.org/resolver/RRID:AB_123456")
RRID:AB_123456
```
"""
@defid(RRID <: BioIdentifier,
       (:prefix(charset(2:12, 'A':'Z', 'a':'z', upper=true)),
        "_", :id(charset(3:12, 'A':'Z', 'a':'z', '0':'9', upper=true))),
       prefix="RRID:", purlprefix="https://scicrunch.org/resolver/RRID:")

# SO

"""
    SO <: BioIdentifier

Sequence Ontology identifier.

The Sequence Ontology is a structured controlled vocabulary for the parts of
a biological sequence. SO terms are used extensively in variant annotation
tools like VEP and SnpEff to describe sequence features and variant effects.
Parsing may throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(SO, "SO:0001234")
SO:0001234

julia> parse(SO, "0000704")
SO:0000704
```
"""
@defid(SO <: BioIdentifier,
       :id(digits(7, pad=7)),
       prefix="SO:", purlprefix="https://purl.obolibrary.org/obo/SO_")

# SRA

"""
    SRA <: BioIdentifier

An [SRA](https://www.ncbi.nlm.nih.gov/sra) (Sequence Read Archive) identifier for sequencing data.

Operated by [NCBI](https://www.ncbi.nlm.nih.gov/), the Sequence Read Archive
is the largest public repository of high-throughput sequencing data. Identifiers
use a prefix indicating the record type — `SRR` (run), `SRX` (experiment),
`SRS` (sample), or `SRP` (study) followed by digits.

# Examples

```julia
julia> parse(SRA, "SRR123456")
SRA:SRR123456

julia> parse(SRA, "SRX789012")
SRA:SRX789012
```
"""
@defid(SRA <: BioIdentifier,
       (:prefix(choice("SRR", "SRX", "SRS", "SRP")), :id(digits(1:18, pad=6))),
       prefix="", purlprefix="https://www.ncbi.nlm.nih.gov/sra/")

# UniProt (must precede UniRef and AFDB which embed it)

"""
    UniProt <: BioIdentifier

A [UniProt](https://www.uniprot.org/) (Universal Protein Resource) accession for protein sequences.

Maintained by the [UniProt Consortium](https://www.uniprot.org/help/about)
(EMBL-EBI, SIB, PIR), UniProt is the most comprehensive resource for protein
sequence and functional annotation. Accessions come in two formats: 6-character
Swiss-Prot entries (e.g. `P04637`) starting with O/P/Q, and 10-character
TrEMBL entries (e.g. `A0A0G2JMC8`) starting with A–N or R–Z.

# Examples

```julia
julia> parse(UniProt, "P04637")
UniProt:P04637

julia> parse(UniProt, "A0A0G2JMC8")
UniProt:A0A0G2JMC8

julia> parse(UniProt, "https://www.uniprot.org/uniprot/Q9Y6K5")
UniProt:Q9Y6K5
```
"""
@defid(UniProt <: BioIdentifier,
       choice(seq(choice("O", "P", "Q"), :d1(digits(1)), :mid(alphnum(3, upper=true)), :d2(digits(1))),
              seq(charset(1, 'A':'N', 'R':'Z', upper=true), :d1(digits(1)), :mid(alphnum(8, upper=true)))),
       prefix="", purlprefix="https://www.uniprot.org/uniprot/")

# UniRef

"""
    UniRef <: BioIdentifier

UniProt Reference Clusters identifier.

UniRef provides clustered sets of sequences from the UniProt Knowledgebase
at 50%, 90%, and 100% sequence identity. These clusters are widely used
by AlphaFold, MMseqs2, and other bioinformatics tools to reduce redundancy.
Identifiers follow format UniRef<level>_<UniProt> where level is 50, 90, or 100.
Parsing may throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(UniRef, "UniRef50_Q9Y6K5")
UniRef:UniRef50_Q9Y6K5

julia> parse(UniRef, "UniRef90_P69905")
UniRef:UniRef90_P69905
```
"""
@defid(UniRef <: BioIdentifier,
       (:cluster(choice("UniRef50", "UniRef90", "UniRef100")),
        "_", :member(embed(UniProt))),
       prefix="", purlprefix="https://www.uniprot.org/uniref/")

# AFDB

"""
    AFDB <: BioIdentifier

AlphaFold Database identifier for predicted 3-D protein models.

AlphaFold DB identifiers reference AI-predicted protein structures from DeepMind/EMBL-EBI.
They follow the format `AF-<UniProt>-F<number>` where UniProt is a valid UniProt accession
and the F-number represents the fragment/model number. Parsing may throw a `MalformedIdentifier`
if the format is invalid.

# Examples

```julia
julia> parse(AFDB, "AF-P69905-F1")
AFDB:AF-P69905-F1

julia> parse(AFDB, "https://alphafold.com/entry/AF-Q9Y6K5-F1")
AFDB:AF-Q9Y6K5-F1
```
"""
@defid(AFDB <: BioIdentifier,
       ("AF-", :uniprot(embed(UniProt)), "-F", :fragment(digits(1:3))),
       prefix="", purlprefix="https://alphafold.com/entry/")

# WikiPathways

"""
    WikiPathways <: BioIdentifier

A [WikiPathways](https://www.wikipathways.org/) community-curated biological pathway identifier.

WikiPathways is an open, collaborative platform for biological pathway curation,
widely used in pathway enrichment analysis tools such as
[Cytoscape](https://cytoscape.org/) and [GSEA](https://www.gsea-msigdb.org/).
Identifiers follow the format `WP` followed by 1–5 digits.

# Examples

```julia
julia> parse(WikiPathways, "WP554")
WikiPathways:WP554

julia> purl(parse(WikiPathways, "WP554"))
"https://www.wikipathways.org/instance/WP554"
```
"""
@defid(WikiPathways <: BioIdentifier,
       ("WP", :id(digits(1:5))),
       prefix="", purlprefix="https://www.wikipathways.org/instance/")

# ──────────────────────────────────────────────────────────────────
# Hand-written types
# ──────────────────────────────────────────────────────────────────

# Ensembl

@defid(GenericEnsembl <: BioIdentifier,
       (choice(:kind,
           :Gene         => "ENSG",
           :Transcript   => "ENST",
           :Protein      => "ENSP",
           :Exon         => "ENSE",
           :Regulatory   => "ENSR",
           :Family       => "ENSF",
           :FamilyMember => "ENSFM"),
        :number(digits(11, pad=11)),
        optional(".", :version(digits(max=65535)))),
       purlprefix="https://www.ensembl.org/id/")

"""
    EnsemblIdentifier{T} <: BioIdentifier

An [Ensembl](https://www.ensembl.org/) stable identifier for genome features.

Maintained by [EMBL-EBI](https://www.ebi.ac.uk/) and the
[Wellcome Sanger Institute](https://www.sanger.ac.uk/), Ensembl provides
reference genome annotations for vertebrates and model organisms. Each
feature type has a distinct prefix: `ENSG` (gene), `ENST` (transcript),
`ENSP` (protein), `ENSE` (exon), `ENSR` (regulatory), `ENSF` (family),
`ENSFM` (family member). Identifiers consist of the prefix followed by
11 zero-padded digits, with an optional `.version` suffix.

# Examples

```julia
julia> parse(ENSG, "ENSG00000141510")
ENSG:ENSG00000141510

julia> parse(ENST, "ENST00000269305.9")
ENST:ENST00000269305.9

julia> purl(parse(ENSG, "ENSG00000141510"))
"https://www.ensembl.org/id/ENSG00000141510"
```
"""
struct EnsemblIdentifier{T} <: BioIdentifier
    inner::GenericEnsembl
end

const ENSG  = EnsemblIdentifier{:Gene}
const ENST  = EnsemblIdentifier{:Transcript}
const ENSP  = EnsemblIdentifier{:Protein}
const ENSE  = EnsemblIdentifier{:Exon}
const ENSR  = EnsemblIdentifier{:Regulatory}
const ENSF  = EnsemblIdentifier{:Family}
const ENSFM = EnsemblIdentifier{:FamilyMember}

EnsemblIdentifier{T}(number::Integer) where {T} =
    EnsemblIdentifier{T}(GenericEnsembl(T, number))
EnsemblIdentifier{T}(number::Integer, version::Integer) where {T} =
    EnsemblIdentifier{T}(GenericEnsembl(T, number, version))

function Base.parse(::Type{EnsemblIdentifier{kind}}, id::AbstractString) where {kind}
    inner = parse(GenericEnsembl, id)
    inner.kind == kind ||
        throw(MalformedIdentifier{EnsemblIdentifier{kind}}(id, 1, "expected $(kind) prefix"))
    EnsemblIdentifier{kind}(inner)
end

function Base.parse(::Type{EnsemblIdentifier}, id::AbstractString)
    inner = parse(GenericEnsembl, id)
    EnsemblIdentifier{inner.kind}(inner)
end

function Base.tryparse(::Type{EnsemblIdentifier{kind}}, id::AbstractString) where {kind}
    inner = tryparse(GenericEnsembl, id)
    isnothing(inner) && return nothing
    inner.kind == kind || return nothing
    EnsemblIdentifier{kind}(inner)
end

function Base.tryparse(::Type{EnsemblIdentifier}, id::AbstractString)
    inner = tryparse(GenericEnsembl, id)
    isnothing(inner) && return nothing
    EnsemblIdentifier{inner.kind}(inner)
end

function Base.show(io::IO, @nospecialize(ens::EnsemblIdentifier))
    if get(io, :limit, false) === true
        if get(io, :typeinfo, Nothing) != typeof(ens)
            print(io, nameof(typeof(ens)), ':')
        end
        shortcode(io, ens)
    else
        show(io, typeof(ens))
        print(io, '(', ens.inner.number)
        isnothing(ens.inner.version) || print(io, ", ", ens.inner.version)
        print(io, ')')
    end
end

Base.isless(@nospecialize(a::EnsemblIdentifier), @nospecialize(b::EnsemblIdentifier)) = isless(a.inner, b.inner)
Base.write(io::IO, @nospecialize(id::EnsemblIdentifier)) = write(io, id.inner)
Base.string(@nospecialize(id::EnsemblIdentifier)) = shortcode(id)
shortcode(io::IO, @nospecialize(id::EnsemblIdentifier)) = shortcode(io, id.inner)
shortcode(@nospecialize(id::EnsemblIdentifier)) = shortcode(id.inner)
purlprefix(@nospecialize(_::Type{<:EnsemblIdentifier})) = purlprefix(GenericEnsembl)
idcode(@nospecialize(ens::EnsemblIdentifier)) = ens.inner.number

# INSDC

"""
    INSDC <: BioIdentifier

International Nucleotide Sequence Database Collaboration accession.

INSDC accessions provide stable identifiers for nucleotide sequences across
the three major databases: GenBank (NCBI), ENA (EBI), and DDBJ. Identifiers
follow format <prefix><digits> with optional version (.version). Parsing may
throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(INSDC, "AY123456.1")
INSDC:AY123456.1

julia> parse(INSDC, "AB123456")
INSDC:AB123456
```
"""
struct INSDC <: BioIdentifier
    prefix::SubString{String}
    number::UInt64
    version::Union{UInt16,Nothing}
end

const INSDC_PREFIXES = [
    "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z",
    "AA", "AB", "AC", "AD", "AE", "AF", "AG", "AH", "AI", "AJ", "AK", "AL", "AM", "AN", "AO", "AP", "AQ", "AR", "AS", "AT", "AU", "AV", "AW", "AX", "AY", "AZ",
    "BA", "BB", "BC", "BD", "BE", "BF", "BG", "BH", "BI", "BJ", "BK", "BL", "BM", "BN", "BO", "BP", "BQ", "BR", "BS", "BT", "BU", "BV", "BW", "BX", "BY", "BZ"
]

function Base.parse(::Type{INSDC}, id::AbstractString)
    id = SubString(string(id))
    _, id = chopprefixes(id, "https://", "http://", "www.ncbi.nlm.nih.gov/nuccore/", "ncbi.nlm.nih.gov/nuccore/", "www.ebi.ac.uk/ena/browser/view/", "ebi.ac.uk/ena/browser/view/")
    # Extract version if present
    base_id, version = if contains(id, '.')
        parts = split(id, '.', limit=2)
        dotpos = findfirst('.', id)
        version_str = parts[2]
        isempty(version_str) && throw(MalformedIdentifier{INSDC}(id, dotpos + 1, "version cannot be empty after period"))
        parsed_version = tryparse(UInt16, version_str)
        isnothing(parsed_version) && throw(MalformedIdentifier{INSDC}(id, dotpos + 1, "must contain only digits"))
        parts[1], parsed_version
    else
        id, nothing
    end
    # Find matching prefix
    prefix_match = nothing
    for prefix in INSDC_PREFIXES
        if startswith(base_id, prefix)
            remainder = base_id[length(prefix)+1:end]
            if all(isdigit, remainder) && !isempty(remainder)
                prefix_match = prefix
                break
            end
        end
    end
    prefix_match === nothing && throw(MalformedIdentifier{INSDC}(id, 1, "must start with valid INSDC prefix followed by digits"))
    remainder = base_id[length(prefix_match)+1:end]
    num = tryparse(UInt64, remainder)
    isnothing(num) && throw(MalformedIdentifier{INSDC}(id, length(prefix_match) + 1, "must contain only digits"))
    INSDC(SubString(prefix_match), num, version)
end

function Base.tryparse(::Type{INSDC}, id::AbstractString)
    try parse(INSDC, id) catch e; e isa MalformedIdentifier || rethrow(); nothing end
end

purlprefix(::Type{INSDC}) = "https://www.ncbi.nlm.nih.gov/nuccore/"

function shortcode(insdc::INSDC)
    base = "$(insdc.prefix)$(insdc.number)"
    isnothing(insdc.version) ? base : "$(base).$(insdc.version)"
end

idcode(insdc::INSDC) = insdc.number

function Base.show(io::IO, insdc::INSDC)
    show(io, INSDC)
    show(io, (insdc.prefix, insdc.number, insdc.version))
end

function Base.print(io::IO, insdc::INSDC)
    get(io, :limit, false) === true && get(io, :compact, false) === true ||
        print(io, "INSDC:")
    print(io, shortcode(insdc))
end


## DrugBank

"""
    DrugBank <: BioIdentifier

A [DrugBank](https://go.drugbank.com/) identifier for drug and pharmaceutical compound entries.

DrugBank combines detailed pharmacological data with drug target information,
providing a resource for drug discovery, pharmacology research, and clinical
decision support. Identifiers follow the format `DB` followed by 5 zero-padded
digits.

# Examples

```julia
julia> parse(DrugBank, "DB00945")
DrugBank:DB00945

julia> purl(parse(DrugBank, "DB00945"))
"https://go.drugbank.com/drugs/DB00945"
```
"""
@defid(DrugBank <: BioIdentifier,
       ("DB", :id(digits(5, pad=5))),
       prefix="", purlprefix="https://go.drugbank.com/drugs/")


## GWAS

"""
    GWAS <: BioIdentifier

NHGRI-EBI GWAS Catalog study accession. The GWAS Catalog provides a curated
collection of published genome-wide association studies, their SNP-trait
associations, and supporting metadata. Identifiers follow the format `GCST`
followed by 1–8 digits.

# Examples

```julia
julia> parse(GWAS, "GCST90000001")
GWAS:GCST90000001

julia> purl(parse(GWAS, "GCST90000001"))
"https://www.ebi.ac.uk/gwas/studies/GCST90000001"
```
"""
@defid(GWAS <: BioIdentifier,
       ("GCST", :id(digits(1:8))),
       prefix="", purlprefix="https://www.ebi.ac.uk/gwas/studies/")


## MP

"""
    MP <: BioIdentifier

Mammalian Phenotype Ontology (MP) term identifier. The MP ontology provides a
structured vocabulary for annotating mammalian phenotypes observed in gene
knockout experiments, QTL studies, and spontaneous mutations. It is a core
ontology used by the Mouse Genome Database (MGD) and the International Mouse
Phenotyping Consortium (IMPC). Identifiers follow the format `MP:` followed
by 7 zero-padded digits.

# Examples

```julia
julia> parse(MP, "MP:0000001")
MP:0000001

julia> purl(parse(MP, "MP:0000001"))
"https://purl.obolibrary.org/obo/MP_0000001"
```
"""
@defid(MP <: BioIdentifier,
       :id(digits(7, pad=7)),
       prefix="MP:", purlprefix="https://purl.obolibrary.org/obo/MP_")


## PATO

"""
    PATO <: BioIdentifier

Phenotype And Trait Ontology (PATO) term identifier. PATO provides a
structured vocabulary for describing phenotypic qualities such as shape,
size, colour, and mass. It is used as a building block by other phenotype
ontologies (including HPO and MP) to compose precise phenotype descriptions.
Identifiers follow the format `PATO:` followed by 7 zero-padded digits.

# Examples

```julia
julia> parse(PATO, "PATO:0000001")
PATO:0000001

julia> purl(parse(PATO, "PATO:0000001"))
"https://purl.obolibrary.org/obo/PATO_0000001"
```
"""
@defid(PATO <: BioIdentifier,
       :id(digits(7, pad=7)),
       prefix="PATO:", purlprefix="https://purl.obolibrary.org/obo/PATO_")


## UBERON

"""
    UBERON <: BioIdentifier

Uber-anatomy ontology (UBERON) cross-species anatomical structure identifier.
UBERON is an integrated ontology covering anatomical structures in animals,
bridging species-specific anatomy ontologies to enable comparative analyses.
It is widely used in multi-species phenotype annotation and gene expression
studies. Identifiers follow the format `UBERON:` followed by 7 zero-padded
digits.

# Examples

```julia
julia> parse(UBERON, "UBERON:0000001")
UBERON:0000001

julia> purl(parse(UBERON, "UBERON:0000001"))
"https://purl.obolibrary.org/obo/UBERON_0000001"
```
"""
@defid(UBERON <: BioIdentifier,
       :id(digits(7, pad=7)),
       prefix="UBERON:", purlprefix="https://purl.obolibrary.org/obo/UBERON_")


## MeSH

"""
    MeSH <: BioIdentifier

Medical Subject Headings descriptor or supplementary concept identifier.

MeSH is the NLM's controlled vocabulary for indexing PubMed and MEDLINE
literature. Descriptor IDs use a `D` prefix followed by 6-9 digits;
supplementary concept IDs use a `C` prefix.

# Examples

```julia
julia> parse(MeSH, "D001249")
MeSH:D001249

julia> parse(MeSH, "C000657245")
MeSH:C657245
```
"""
@defid(MeSH <: BioIdentifier,
       (choice("D", "C"), :id(digits(6:9, pad=6))),
       prefix="", purlprefix="https://meshb.nlm.nih.gov/record/ui?ui=")


## InterPro

"""
    InterPro <: BioIdentifier

An [InterPro](https://www.ebi.ac.uk/interpro/) protein family and domain identifier.

Maintained by [EMBL-EBI](https://www.ebi.ac.uk/), InterPro integrates protein
signatures from Pfam, PRINTS, PROSITE, SMART, and other member databases into
a single classification resource. Identifiers have the format `IPR` followed
by 6 zero-padded digits.

# Examples

```julia
julia> parse(InterPro, "IPR011009")
InterPro:IPR011009

julia> purl(parse(InterPro, "IPR011009"))
"https://www.ebi.ac.uk/interpro/entry/InterPro/IPR011009"
```
"""
@defid(InterPro <: BioIdentifier,
       ("IPR", :id(digits(6, pad=6))),
       prefix="", purlprefix="https://www.ebi.ac.uk/interpro/entry/InterPro/")


## Pfam

"""
    Pfam <: BioIdentifier

A [Pfam](https://www.ebi.ac.uk/interpro/entry/pfam/) protein family identifier.

Originally developed at the Sanger Institute and now hosted by
[EMBL-EBI](https://www.ebi.ac.uk/) as part of [`InterPro`](@ref), Pfam
classifies protein sequences into families and domains using hidden Markov
model profiles. Identifiers have the format `PF` followed by 5 zero-padded
digits.

# Examples

```julia
julia> parse(Pfam, "PF00069")
Pfam:PF00069

julia> purl(parse(Pfam, "PF00069"))
"https://www.ebi.ac.uk/interpro/entry/pfam/PF00069"
```
"""
@defid(Pfam <: BioIdentifier,
       ("PF", :id(digits(5, pad=5))),
       prefix="", purlprefix="https://www.ebi.ac.uk/interpro/entry/pfam/")


## MGI (Mouse Genome Informatics)

"""
    MGI <: BioIdentifier

Mouse Genome Informatics gene identifier. MGI is the authoritative
resource for the laboratory mouse genome, maintained by the Jackson
Laboratory. A founding member of the Alliance of Genome Resources.

Format: `MGI:` followed by 1–7 digits.

# Examples

```julia
julia> parse(MGI, "MGI:96417")
MGI:96417

julia> purl(parse(MGI, "MGI:96417"))
"https://www.informatics.jax.org/marker/MGI:96417"
```
"""
@defid(MGI <: BioIdentifier,
       :id(digits(1:7)),
       prefix="MGI:", purlprefix="https://www.informatics.jax.org/marker/MGI:")


## SGD (Saccharomyces Genome Database)

"""
    SGD <: BioIdentifier

Saccharomyces Genome Database gene identifier. SGD is the authoritative
resource for the budding yeast *Saccharomyces cerevisiae* genome. A
founding member of the Alliance of Genome Resources.

Format: `S` followed by 9 zero-padded digits (database identifier form).

Note: SGD also uses systematic ORF names (e.g. `YAL001C`) in publications,
but this type represents the stable database identifier.

# Examples

```julia
julia> parse(SGD, "S000002493")
SGD:S000002493
```
"""
@defid(SGD <: BioIdentifier,
       ("S", :id(digits(9, pad=9))),
       prefix="", purlprefix="https://www.yeastgenome.org/locus/")


## FlyBase

"""
    FlyBase <: BioIdentifier

A [FlyBase](https://flybase.org/) gene identifier for *Drosophila* species.

FlyBase is the primary database for *Drosophila* genetics and genomics,
and a founding member of the [Alliance of Genome Resources](https://www.alliancegenome.org/).
Gene identifiers have the format `FBgn` followed by 7 zero-padded digits.

# Examples

```julia
julia> parse(FlyBase, "FBgn0031701")
FlyBase:FBgn0031701

julia> purl(parse(FlyBase, "FBgn0031701"))
"https://flybase.org/reports/FBgn0031701"
```
"""
@defid(FlyBase <: BioIdentifier,
       ("FBgn", :id(digits(7, pad=7))),
       prefix="", purlprefix="https://flybase.org/reports/")


## WormBase

"""
    WormBase <: BioIdentifier

A [WormBase](https://wormbase.org/) gene identifier for *Caenorhabditis elegans* and related nematodes.

WormBase is the central resource for nematode biology, and a founding member
of the [Alliance of Genome Resources](https://www.alliancegenome.org/). Gene
identifiers have the format `WBGene` followed by 8 zero-padded digits.

# Examples

```julia
julia> parse(WormBase, "WBGene00006763")
WormBase:WBGene00006763

julia> purl(parse(WormBase, "WBGene00006763"))
"https://wormbase.org/species/c_elegans/gene/WBGene00006763"
```
"""
@defid(WormBase <: BioIdentifier,
       ("WBGene", :id(digits(8, pad=8))),
       prefix="", purlprefix="https://wormbase.org/species/c_elegans/gene/")

end
