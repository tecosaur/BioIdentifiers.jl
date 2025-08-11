# SPDX-FileCopyrightText: © 2025 TEC <contact@tecosaur.net>
# SPDX-License-Identifier: MPL-2.0

"""
    BioIdentifiers

Structured and validated types for biological and biomedical identifiers.

This module provides type-safe representations of identifiers commonly used in
biological and biomedical databases, extending the DigitalIdentifiersBase framework.
Each identifier type validates input format, provides canonical string representations,
and includes persistent URLs where applicable.

# Supported Identifiers

**Protein & Structure**
- `AFDB`: AlphaFold Database protein structure predictions
- `UniProt`: Universal Protein Resource database
- `IntAct`: Molecular interaction database

**Genomics & Genetics**
- Ensembl (`ENSG`, `ENST`, `ENSP`, `ENSE`, `ENSR`, `ENSF`, `ENSFM`), Ensembl genome annotations
- `NCBIGene`: NCBI Gene database
- `RefSeq`: NCBI Reference Sequence database
- `HGNC`: HUGO Gene Nomenclature Committee
- `OMIM`: Online Mendelian Inheritance in Man

**Variation & Clinical**
- `ClinVar`: Clinical genomic variant database
- `dbSNP`: Single Nucleotide Polymorphism database
- `dbVar`: Genomic structural variation database

**Expression & Functional Genomics**
- `ArrayExpress`: Functional genomics data archive
- `GEO`: Gene Expression Omnibus

**Studies & Samples**
- `BioProject`: Biological research projects
- `BioSample`: Biological sample metadata
- `SRA`: Sequence Read Archive
- `EGA`: European Genome-phenome Archive
- `dbGaP`: Database of Genotypes and Phenotypes
- `ClinicalTrials`: Clinical trial registry

**Ontologies & Controlled Vocabularies**
- `GO`: Gene Ontology
- `HPO`: Human Phenotype Ontology
- `MONDO`: Monarch Disease Ontology
- `DOID`: Disease Ontology
- `CL`: Cell Ontology
- `EFO`: Experimental Factor Ontology

**Chemical & Metabolic**
- `ChEBI`: Chemical Entities of Biological Interest
- `HMDB`: Human Metabolome Database
- `KEGG`: Kyoto Encyclopedia of Genes and Genomes
- `MetaboLights`: Metabolomics data repository

**Networks & Interactions**
- `BioGRID`: Biological General Repository for Interaction Datasets
- `Reactome`: Pathway database

**Cell Lines & Model Organisms**
- `Cellosaurus`: Cell line database
- `NCBITaxon`: NCBI Taxonomy database

**Proteomics**
- `PXD`: ProteomeXchange dataset identifiers

**Resources & Reagents**
- `RRID`: Research Resource Identifiers

**Alleles & Variants**
- `CA`: ClinGen Allele Registry
- `VIAF`: Virtual International Authority File (when used for biological entities)

Each identifier can be constructed from its canonical string form and provides
methods for validation, canonical output, and URL generation where applicable.
"""
module BioIdentifiers

import DigitalIdentifiersBase: AcademicIdentifier, MalformedIdentifier, parseid, parsefor, chopprefixes, shortcode, purl, idcode, purlprefix

export BioIdentifier, AFDB, ArrayExpress, BioGRID, BioProject, BioSample, CA, CAS, CL, Cellosaurus, ChEBI, ChEMBL, ClinVar, ClinicalTrials, dbGaP, dbSNP, dbVar, DOID, DrugBank, ECO, EFO, EGA, ENSG, ENST, ENSP, ENSE, ENSR, ENSF, ENSFM, GO, GEO, GWAS, HGNC, HMDB, HPO, ICD10, INSDC, IntAct, KEGG, LOINC, MP, MetaboLights, MONDO, NCBIGene, NCBITaxon, OMIM, PATO, PDB, PubChem, PXD, Reactome, RefSeq, RRID, SO, SRA, UBERON, UniProt, UniRef, WikiPathways

abstract type BioIdentifier <: AcademicIdentifier end

abstract type OboIdentifier <: BioIdentifier end

macro numericid(name::Symbol, kws...)
    namelower = Symbol(lowercase(string(name)))
    urls, prefix, inttype, digits, maxdigits = String[], "", UInt64, 0, 0
    for kw in kws
        Meta.isexpr(kw, :(=), 2) || throw(ArgumentError("Expected keyword argument for numericid macro, got: $kw"))
        key, val = kw.args
        if key == :url
            if val isa String
                push!(urls, val)
            elseif Meta.isexpr(val, :vect)
                urls = Vector{String}(val.args)
            else
                throw(ArgumentError("Invalid url value: $val"))
            end
        elseif key == :prefix
            prefix = val::String
        elseif key == :digits
            digits = val::Int
        elseif key == :inttype
            inttype = val::Symbol
        else
            throw(ArgumentError("Unknown keyword argument for numericid macro: $(kw.args[1])"))
        end
    end
    rmprefixes = String[]
    if !isempty(urls)
        push!(rmprefixes, "https://", "http://")
        any(u -> startswith(u, "www."), urls) && push!(rmprefixes, "www.")
        for url in urls
            push!(rmprefixes, lowercase(chopprefix(url, "www.")))
        end
    end
    isempty(prefix) || push!(rmprefixes, lowercase(prefix))
    output = Expr[]
    typearg = Expr(:(::), Expr(:curly, :Type, name))
    push!(output, :(function parseid(::Type{$name}, id::SubString)
        _, id = chopprefixes(id, $(rmprefixes...))
        num = parsefor($(name), $inttype, id)
        num isa $inttype || return num
        $(name)(num)
    end))
    isempty(urls) || push!(output, :(purlprefix(::Type{$name}) = $("https://$(first(urls))")))
    isempty(prefix) || if iszero(digits)
        push!(output, :(shortcode($namelower::$name) = string($prefix, $namelower.id)))
    else
        push!(output, :(shortcode($namelower::$name) = string($prefix, lpad($namelower.id, $digits, '0'))))
    end
    push!(output, :(idcode($namelower::$name) = $namelower.id))
    push!(output, :(Base.show(io::IO, $namelower::$name) = (show(io, $name); print(io, '(', $namelower.id, ')'))))
    Expr(:block, output...) |> esc
end


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
struct ArrayExpress <: BioIdentifier
    platform::SubString{String}
    number::UInt32
end

function parseid(::Type{ArrayExpress}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "www.ebi.ac.uk/biostudies/arrayexpress/studies/", "ebi.ac.uk/biostudies/arrayexpress/studies/")
    startswith(id, "E-") || return MalformedIdentifier{ArrayExpress}(id, "must start with 'E-'")
    parts = split(id, '-')
    length(parts) == 3 || return MalformedIdentifier{ArrayExpress}(id, "must follow format E-<platform>-<number>")
    platform, number = parts[2], parts[3]
    (platform == "GEOD" || (4 <= length(platform) <= 5 && all(c -> isuppercase(c) && isletter(c), platform))) || return MalformedIdentifier{ArrayExpress}(id, "platform code must be 4-5 uppercase letters or 'GEOD'")
    num = parsefor(ArrayExpress, UInt32, number)
    num isa UInt32 || return num
    ArrayExpress(SubString(platform), num)
end

purlprefix(::Type{ArrayExpress}) = "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/"
shortcode(ae::ArrayExpress) = "E-$(ae.platform)-$(ae.number)"
idcode(ae::ArrayExpress) = ae.number

Base.show(io::IO, ae::ArrayExpress) = (show(io, ArrayExpress); show(io, (ae.platform, ae.number)))


# BioGRID

"""
    BioGRID <: BioIdentifier

BioGRID molecular interaction database identifier.

BioGRID is a biomedical interaction repository with data compiled through comprehensive
curation efforts. Identifiers are simple positive integers. Parsing may throw a
`MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(BioGRID, "123456")
BioGRID:123456

julia> parse(BioGRID, "https://thebiogrid.org/789012")
BioGRID:789012
```
"""
struct BioGRID <: BioIdentifier
    id::UInt64
end

@numericid BioGRID url = ["www.thebiogrid.org/", "wiki.thebiogrid.org/"] inttype = UInt64


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
struct BioProject <: BioIdentifier
    database::UInt8
    number::UInt64
end

const BIOPROJECT_PREFIXES = ("prjna", "prjeb", "prjdb")

function parseid(::Type{BioProject}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "www.ncbi.nlm.nih.gov/bioproject/", "ncbi.nlm.nih.gov/bioproject/")
    for (i, prefix) in enumerate(BIOPROJECT_PREFIXES)
        found, remainder = chopprefixes(id, prefix)
        if found
            number = parsefor(BioProject, UInt64, remainder)
            number isa UInt64 || return number
            return BioProject(UInt8(i), number)
        end
    end
    MalformedIdentifier{BioProject}(id, "must start with PRJNA, PRJEB, or PRJDB")
end

purlprefix(::Type{BioProject}) = "https://www.ncbi.nlm.nih.gov/bioproject/"
shortcode(bp::BioProject) = string(uppercase(BIOPROJECT_PREFIXES[bp.database]), bp.number)
idcode(bp::BioProject) = bp.number

Base.show(io::IO, bp::BioProject) = (show(io, BioProject); show(io, (uppercase(BIOPROJECT_PREFIXES[bp.database]), bp.number)))


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
struct BioSample <: BioIdentifier
    database::UInt8
    number::UInt64
end

const BIOSAMPLE_PREFIXES = ("samn", "samea", "samd")

function parseid(::Type{BioSample}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "www.ncbi.nlm.nih.gov/biosample/", "ncbi.nlm.nih.gov/biosample/")
    for (i, prefix) in enumerate(BIOSAMPLE_PREFIXES)
        found, remainder = chopprefixes(id, prefix)
        if found
            number = parsefor(BioSample, UInt64, remainder)
            number isa UInt64 || return number
            return BioSample(UInt8(i), number)
        end
    end
    return MalformedIdentifier{BioSample}(id, "must start with SAMN, SAMEA, or SAMD")
end

purlprefix(::Type{BioSample}) = "https://www.ncbi.nlm.nih.gov/biosample/"
shortcode(bs::BioSample) = string(uppercase(BIOSAMPLE_PREFIXES[bs.database]), bs.number)
idcode(bs::BioSample) = bs.number

Base.show(io::IO, bs::BioSample) = (show(io, BioSample); show(io, (uppercase(BIOSAMPLE_PREFIXES[bs.database]), bs.number)))


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
CA:CA123456

julia> parse(CA, "https://reg.clinicalgenome.org/allele/CA789012")
CA:CA789012
```
"""
struct CA <: BioIdentifier
    id::UInt64
end

@numericid CA url = "reg.clinicalgenome.org/allele/" prefix = "CA" inttype = UInt64




# CL

"""
    CL <: BioIdentifier

Cell Ontology term identifier for cross-species cell types.

The Cell Ontology is designed to represent the in vivo and in vitro cells of animals,
plants, fungi and other organisms. Identifiers follow format `CL:<7 digits>`.
Parsing may throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(CL, "CL:0000001")
CL:CL:0000001

julia> parse(CL, "0000001")
CL:CL:0000001
```
"""
struct CL <: BioIdentifier
    id::UInt32
end

function parseid(::Type{CL}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "purl.obolibrary.org/obo/CL_", "CL:")
    all(isdigit, id) && length(id) == 7 || return MalformedIdentifier{CL}(id, "must be 7 digits (optionally prefixed with CL:)")
    num = parsefor(CL, UInt32, id)
    num isa UInt32 || return num
    CL(num)
end

function purl(cl::CL)
    digits = lpad(cl.id, 7, '0')
    "https://purl.obolibrary.org/obo/CL_$(digits)"
end

shortcode(cl::CL) = string("CL:", lpad(cl.id, 7, '0'))
idcode(cl::CL) = cl.id


# Cellosaurus

"""
    Cellosaurus <: BioIdentifier

Cellosaurus authoritative cell line registry identifier.

Cellosaurus is a knowledge resource on cell lines providing extensive information
including origin, synonyms, and references. Identifiers follow format `CVCL_<4-7 digits>`.
Parsing may throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(Cellosaurus, "CVCL_0001")
Cellosaurus:CVCL_0001

julia> parse(Cellosaurus, "CVCL_1234567")
Cellosaurus:CVCL_1234567
```
"""
struct Cellosaurus <: BioIdentifier
    id::UInt32
end

@numericid Cellosaurus url = "www.cellosaurus.org/" prefix = "CVCL_" digits = 4 inttype = UInt32





# ChEBI

"""
    ChEBI <: BioIdentifier

Chemical Entities of Biological Interest (ChEBI) identifier.

ChEBI is a dictionary of molecular entities focused on small chemical compounds.
It provides an ontological classification and unambiguous nomenclature for chemical
entities. Identifiers follow format `CHEBI:<number>`. Parsing may throw a
`MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(ChEBI, "CHEBI:15377")
ChEBI:CHEBI:15377

julia> parse(ChEBI, "15377")
ChEBI:CHEBI:15377
```
"""
struct ChEBI <: BioIdentifier
    id::UInt32
end

@numericid ChEBI url = "www.ebi.ac.uk/chebi/searchId.do?chebiId=" prefix = "CHEBI:" inttype = UInt32

function purl(chebi::ChEBI)
    "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:$(chebi.id)"
end




# ChEMBL

"""
    ChEMBL <: BioIdentifier

ChEMBL database identifier for bioactive compounds, targets, and assays.

ChEMBL is a manually curated database of bioactive molecules with drug-like
properties. Identifiers follow format CHEMBL<number> for molecules, targets,
and assays. Parsing may throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(ChEMBL, "CHEMBL25")
ChEMBL:CHEMBL25

julia> parse(ChEMBL, "CHEMBL1075272")
ChEMBL:CHEMBL1075272
```
"""
struct ChEMBL <: BioIdentifier
    id::UInt64
end

@numericid ChEMBL url = ["www.ebi.ac.uk/chembl/compound_report_card/", "www.ebi.ac.uk/chembl/target_report_card/", "www.ebi.ac.uk/chembl/assay_report_card/"] prefix = "CHEMBL" inttype = UInt64

function parseid(::Type{ChEMBL}, id::SubString)
    _, id = chopprefixes(id, "www.ebi.ac.uk/chembl/compound_report_card/", "www.ebi.ac.uk/chembl/target_report_card/", "www.ebi.ac.uk/chembl/assay_report_card/", "ebi.ac.uk/chembl/compound_report_card/", "ebi.ac.uk/chembl/target_report_card/", "ebi.ac.uk/chembl/assay_report_card/", "CHEMBL")
    remainder = startswith(id, "CHEMBL") ? id[7:end] : id
    all(isdigit, remainder) || return MalformedIdentifier{ChEMBL}(id, "must contain only digits after 'CHEMBL'")
    isempty(remainder) && return MalformedIdentifier{ChEMBL}(id, "must have a number after 'CHEMBL'")
    num = parsefor(ChEMBL, UInt64, remainder)
    num isa UInt64 || return num
    num == 0 && return MalformedIdentifier{ChEMBL}(id, "number must be positive")
    ChEMBL(num)
end



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
struct ClinVar <: BioIdentifier
    type::UInt8
    number::UInt64
    version::UInt16
end

const CLINVAR_PREFIXES = ("vcv", "rcv", "scv")

function parseid(::Type{ClinVar}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "www.ncbi.nlm.nih.gov/clinvar/", "ncbi.nlm.nih.gov/clinvar/")
    for (i, prefix) in enumerate(CLINVAR_PREFIXES)
        found, remainder = chopprefixes(id, prefix)
        if found
            number_str, version = if '.' in remainder
                parts = split(remainder, '.', limit=2)
                version_str = parts[2]
                isempty(version_str) && return MalformedIdentifier{ClinVar}(id, "version cannot be empty after period")
                parsed_version = parsefor(ClinVar, UInt16, version_str)
                parsed_version isa UInt16 || return parsed_version
                parts[1], parsed_version
            else
                remainder, UInt16(0)
            end
            number = parsefor(ClinVar, UInt64, number_str)
            number isa UInt64 || return number
            return ClinVar(UInt8(i), number, version)
        end
    end
    return MalformedIdentifier{ClinVar}(id, "must start with VCV, RCV, or SCV")
end

purlprefix(::Type{ClinVar}) = "https://www.ncbi.nlm.nih.gov/clinvar/"
function shortcode(cv::ClinVar)
    prefix = uppercase(CLINVAR_PREFIXES[cv.type])
    # VCV uses 9 digits, RCV/SCV can use 8 or 9 digits
    pad_length = cv.type == 1 ? 9 : 8  # VCV=1, RCV=2, SCV=3
    if cv.version == 0
        string(prefix, lpad(cv.number, pad_length, '0'))
    else
        string(prefix, lpad(cv.number, pad_length, '0'), '.', cv.version)
    end
end
idcode(cv::ClinVar) = cv.number

Base.show(io::IO, cv::ClinVar) = (show(io, ClinVar); show(io, (uppercase(CLINVAR_PREFIXES[cv.type]), cv.number, cv.version)))


# ClinicalTrials

"""
    ClinicalTrials <: BioIdentifier

ClinicalTrials.gov registry identifier for interventional/observational studies.

ClinicalTrials.gov is a database of privately and publicly funded clinical studies
conducted around the world. Identifiers follow format `NCT<8 digits>`. Parsing may
throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(ClinicalTrials, "NCT01234567")
ClinicalTrials:NCT01234567

julia> parse(ClinicalTrials, "https://clinicaltrials.gov/ct2/show/NCT01234567")
ClinicalTrials:NCT01234567
```
"""
struct ClinicalTrials <: BioIdentifier
    id::UInt32
end

@numericid ClinicalTrials url = "clinicaltrials.gov/ct2/show/" prefix = "NCT" digits = 8 inttype = UInt32

function parseid(::Type{ClinicalTrials}, id::SubString)
    _, id = chopprefixes(id, "clinicaltrials.gov/ct2/show/", "clinicaltrials.gov/show/")
    found, remainder = chopprefixes(id, "nct")
    found || return MalformedIdentifier{ClinicalTrials}(id, "must start with 'NCT'")
    length(remainder) == 8 && all(isdigit, remainder) || return MalformedIdentifier{ClinicalTrials}(id, "must have exactly 8 digits after 'NCT'")
    num = parsefor(ClinicalTrials, UInt32, remainder)
    num isa UInt32 || return num
    ClinicalTrials(num)
end



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
struct dbGaP <: BioIdentifier
    study::UInt32
    version::UInt16
    participant::UInt16
end

function parseid(::Type{dbGaP}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=", "ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=")
    found, remainder = chopprefixes(id, "phs")
    found || return MalformedIdentifier{dbGaP}(id, "must start with 'phs'")
    parts = split(remainder, '.')
    length(parts) == 3 || return MalformedIdentifier{dbGaP}(id, "must follow format phs<number>.v<version>.p<participant>")
    study_str = parts[1]
    length(study_str) == 6 && all(isdigit, study_str) || return MalformedIdentifier{dbGaP}(id, "study number must be exactly 6 digits")
    study = parsefor(dbGaP, UInt32, study_str)
    study isa UInt32 || return study
    version_part = parts[2]
    version_found, version_remainder = chopprefixes(version_part, "v")
    version_found || return MalformedIdentifier{dbGaP}(id, "version part must start with 'v'")
    version = parsefor(dbGaP, UInt16, version_remainder)
    version isa UInt16 || return version
    participant_part = parts[3]
    participant_found, participant_remainder = chopprefixes(participant_part, "p")
    participant_found || return MalformedIdentifier{dbGaP}(id, "participant part must start with 'p'")
    participant = parsefor(dbGaP, UInt16, participant_remainder)
    participant isa UInt16 || return participant
    dbGaP(study, version, participant)
end

function purl(gap::dbGaP)
    id = shortcode(gap)
    "https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=$(id)"
end

shortcode(dbgap::dbGaP) = "phs$(lpad(dbgap.study, 6, '0')).v$(dbgap.version).p$(dbgap.participant)"
idcode(dbgap::dbGaP) = dbgap.study

Base.show(io::IO, dbgap::dbGaP) = (show(io, dbGaP); show(io, (dbgap.study, dbgap.version, dbgap.participant)))


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
struct dbSNP <: BioIdentifier
    id::UInt32
end

@numericid dbSNP url = "www.ncbi.nlm.nih.gov/snp/" prefix = "rs" inttype = UInt32





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
struct dbVar <: BioIdentifier
    type::UInt8
    number::UInt32
end

const DBVAR_PREFIXES = ("nsv", "nssv", "nstd")

function parseid(::Type{dbVar}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "www.ncbi.nlm.nih.gov/dbvar/", "ncbi.nlm.nih.gov/dbvar/")
    for (i, prefix) in enumerate(DBVAR_PREFIXES)
        found, remainder = chopprefixes(id, prefix)
        if found
            number = parsefor(dbVar, UInt32, remainder)
            number isa UInt32 || return number
            return dbVar(UInt8(i), number)
        end
    end
    return MalformedIdentifier{dbVar}(id, "must start with nsv, nssv, or nstd")
end

purlprefix(::Type{dbVar}) = "https://www.ncbi.nlm.nih.gov/dbvar/"
shortcode(dbvar::dbVar) = string(uppercase(DBVAR_PREFIXES[dbvar.type]), dbvar.number)
idcode(dbvar::dbVar) = dbvar.number

Base.show(io::IO, dbvar::dbVar) = (show(io, dbVar); show(io, (uppercase(DBVAR_PREFIXES[dbvar.type]), dbvar.number)))


# DOID

"""
    DOID <: BioIdentifier

Disease Ontology identifier for human diseases.

The Human Disease Ontology provides a standardized ontology for human disease with
the purpose of providing the biomedical community with consistent, reusable and
sustainable descriptions of human disease terms. Identifiers follow format `DOID:<number>`.
Parsing may throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(DOID, "DOID:0001816")
DOID:DOID:0001816

julia> parse(DOID, "0001816")
DOID:DOID:0001816
```
"""
struct DOID <: BioIdentifier
    id::UInt32
end

@numericid DOID url = "purl.obolibrary.org/obo/DOID_" prefix = "DOID:" digits = 7 inttype = UInt32

function purl(doid::DOID)
    digits = lpad(doid.id, 7, '0')
    "https://purl.obolibrary.org/obo/DOID_$(digits)"
end



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
EFO:EFO:0000001

julia> parse(EFO, "0000001")
EFO:EFO:0000001
```
"""
struct EFO <: BioIdentifier
    id::UInt32
end

@numericid EFO url = "www.ebi.ac.uk/ols/ontologies/efo/terms?short_form=" prefix = "EFO:" digits = 7 inttype = UInt32

function purl(efo::EFO)
    "https://www.ebi.ac.uk/ols/ontologies/efo/terms?short_form=EFO:$(lpad(efo.id, 7, '0'))"
end



# EGA

"""
    EGA <: BioIdentifier

European Genome-phenome Archive controlled-access identifier.

The European Genome-phenome Archive provides a service for permanent archiving and
sharing of personally identifiable genetic and phenotypic data resulting from
biomedical research projects. Identifiers follow formats EGAD/EGAS<11 digits> for
dataset/study records respectively. Parsing may throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(EGA, "EGAD00001000001")
EGA:EGAD00001000001

julia> parse(EGA, "EGAS00001000002")
EGA:EGAS00001000002
```
"""
struct EGA <: BioIdentifier
    type::UInt8
    number::UInt64
end

const EGA_PREFIXES = ("egad", "egas")

function parseid(::Type{EGA}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "ega-archive.org/datasets/", "ega-archive.org/studies/")
    for (i, prefix) in enumerate(EGA_PREFIXES)
        found, remainder = chopprefixes(id, prefix)
        if found
            length(remainder) == 11 && all(isdigit, remainder) || return MalformedIdentifier{EGA}(id, "must have exactly 11 digits after prefix")
            num = parsefor(EGA, UInt64, remainder)
            num isa UInt64 || return num
            return EGA(UInt8(i), num)
        end
    end
    return MalformedIdentifier{EGA}(id, "must start with EGAD or EGAS")
end

function purl(ega::EGA)
    id = shortcode(ega)
    base_url = ega.type == 1 ? "https://ega-archive.org/datasets/" : "https://ega-archive.org/studies/"
    "$(base_url)$(id)"
end

shortcode(ega::EGA) = string(uppercase(EGA_PREFIXES[ega.type]), lpad(ega.number, 11, '0'))
idcode(ega::EGA) = ega.number

Base.show(io::IO, ega::EGA) = (show(io, EGA); show(io, (uppercase(EGA_PREFIXES[ega.type]), ega.number)))


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
ECO:ECO:0000313

julia> parse(ECO, "0000269")
ECO:ECO:0000269
```
"""
struct ECO <: BioIdentifier
    id::UInt32
end

@numericid ECO url = "purl.obolibrary.org/obo/ECO_" prefix = "ECO:" digits = 7 inttype = UInt32

function parseid(::Type{ECO}, id::SubString)
    _, id = chopprefixes(id, "purl.obolibrary.org/obo/ECO_", "ECO:")
    length(id) == 7 && all(isdigit, id) || return MalformedIdentifier{ECO}(id, "must be 7 digits (optionally prefixed with ECO:)")
    num = parsefor(ECO, UInt32, id)
    num isa UInt32 || return num
    ECO(num)
end

function purl(eco::ECO)
    digits = lpad(eco.id, 7, '0')
    "https://purl.obolibrary.org/obo/ECO_$(digits)"
end




# Ensembl

"""
    EnsemblIdentifier{T} <: BioIdentifier

Ensembl stable identifier for genome features.

Ensembl creates, integrates and distributes reference datasets and analysis tools
that enable genomics research. Identifiers follow various formats for different
feature types, optionally with version suffix. Parsing may throw a `MalformedIdentifier`
if the format is invalid.

Type parameter `T` indicates the feature type:
- `:Gene` (ENSG) - Gene identifiers
- `:Transcript` (ENST) - Transcript identifiers
- `:Protein` (ENSP) - Protein/peptide identifiers
- `:Exon` (ENSE) - Exon identifiers
- `:Regulatory` (ENSR) - Regulatory feature identifiers
- `:Family` (ENSF) - Protein family identifiers
- `:FamilyMember` (ENSFM) - Protein family member identifiers

# Examples

```julia
julia> parse(ENSG, "ENSG00000141510")
ENSG:ENSG00000141510

julia> parse(ENST, "ENST00000269305.v4")
ENST:ENST00000269305.4

julia> parse(ENSE, "ENSE00001234567")
ENSE:ENSE00001234567
```
"""
struct EnsemblIdentifier{T} <: BioIdentifier
    number::UInt64
    version::Union{UInt16,Nothing}
end

# Type aliases for specific Ensembl identifier types
const ENSG = EnsemblIdentifier{:Gene}
const ENST = EnsemblIdentifier{:Transcript}
const ENSP = EnsemblIdentifier{:Protein}
const ENSE = EnsemblIdentifier{:Exon}
const ENSR = EnsemblIdentifier{:Regulatory}
const ENSF = EnsemblIdentifier{:Family}
const ENSFM = EnsemblIdentifier{:FamilyMember}
# Parallel tuples for type symbols and prefixes
const ENSEMBL_TYPES = (:Gene, :Transcript, :Protein, :Exon, :Regulatory, :Family, :FamilyMember)
const ENSEMBL_PREFIXES = ("ensg", "enst", "ensp", "ense", "ensr", "ensf", "ensfm")

function parseid(::Type{EnsemblIdentifier{kind}}, id::SubString) where {kind}
    isempty(id) && return MalformedIdentifier{EnsemblIdentifier{kind}}(id, "cannot be empty")
    _, id = chopprefixes(id, "https://", "http://", "www.ensembl.org/id/", "ensembl.org/id/")
    version = nothing
    original_id = id
    if contains(id, ".v")
        return MalformedIdentifier{EnsemblIdentifier{kind}}(original_id, "version should use dot notation (.1) not v notation (v1)")
    elseif contains(id, ".")
        parts = split(id, ".")
        length(parts) == 2 || return MalformedIdentifier{EnsemblIdentifier{kind}}(original_id, "invalid version format")
        id, version_str = parts
        version_num = parsefor(EnsemblIdentifier{kind}, UInt16, version_str)
        version_num isa UInt16 || return version_num
        version = version_num
    end
    type_index = findfirst(==(kind), ENSEMBL_TYPES)
    type_index === nothing && error("Unknown Ensembl type: $kind")
    expected_prefix = ENSEMBL_PREFIXES[type_index]
    if startswith(lowercase(id), expected_prefix)
        remainder = @view id[length(expected_prefix)+1:end]
    else
        return MalformedIdentifier{EnsemblIdentifier{kind}}(original_id, "must start with $(uppercase(expected_prefix))")
    end
    length(remainder) == 11 && all(isdigit, remainder) ||
        return MalformedIdentifier{EnsemblIdentifier{kind}}(original_id, "must have exactly 11 digits after prefix")
    num = parsefor(EnsemblIdentifier{kind}, UInt64, remainder)
    num isa UInt64 || return num
    EnsemblIdentifier{kind}(num, version)
end

function parseid(::Type{EnsemblIdentifier}, id::SubString)
    isempty(id) && return MalformedIdentifier{EnsemblIdentifier}(id, "cannot be empty")
    _, id = chopprefixes(id, "https://", "http://", "www.ensembl.org/id/", "ensembl.org/id/")
    version = nothing
    original_id = id
    if contains(id, ".v")
        return MalformedIdentifier{EnsemblIdentifier}(original_id, "version should use dot notation (.1) not v notation (v1)")
    elseif contains(id, ".")
        parts = split(id, ".")
        length(parts) == 2 || return MalformedIdentifier{EnsemblIdentifier}(original_id, "invalid version format")
        id, version_str = parts
        version_num = parsefor(EnsemblIdentifier, UInt16, version_str)
        version_num isa UInt16 || return version_num
        version = version_num
    end
    for (i, prefix) in enumerate(ENSEMBL_PREFIXES)
        if startswith(lowercase(id), prefix)
            remainder = @view id[length(prefix)+1:end]
            length(remainder) == 11 && all(isdigit, remainder) ||
                return MalformedIdentifier{EnsemblIdentifier}(original_id, "must have exactly 11 digits after prefix")
            num = parsefor(EnsemblIdentifier, UInt64, remainder)
            num isa UInt64 || return num
            return EnsemblIdentifier{ENSEMBL_TYPES[i]}(num, version)
        end
    end
    MalformedIdentifier{EnsemblIdentifier}(original_id, "must start with a valid Ensembl prefix (ENSG, ENST, ENSP, ENSE, ENSR, ENSF, ENSFM)")
end

purlprefix(::Type{<:EnsemblIdentifier}) = "https://www.ensembl.org/id/"
function shortcode(ens::EnsemblIdentifier{T}) where {T}
    type_index = findfirst(==(T), ENSEMBL_TYPES)
    prefix = uppercase(ENSEMBL_PREFIXES[type_index])
    base = string(prefix, lpad(ens.number, 11, '0'))
    if isnothing(ens.version)
        base
    else
        string(base, '.', ens.version)
    end
end

idcode(ens::EnsemblIdentifier) = ens.number

function Base.show(io::IO, ens::EnsemblIdentifier{T}) where {T}
    type_index = findfirst(==(T), ENSEMBL_TYPES)
    prefix = uppercase(ENSEMBL_PREFIXES[type_index])
    show(io, EnsemblIdentifier{T})
    print(io, '(', ens.number)
    isnothing(ens.version) || print(io, ", ", ens.version)
    print(io, ')')
end


# GO

"""
    GO <: BioIdentifier

Gene Ontology term identifier.

The Gene Ontology provides a framework for the model of biology through a structured
controlled vocabulary. It covers three domains: cellular component, molecular function,
and biological process. Identifiers follow format `GO:<7 digits>`. Parsing may throw
a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(GO, "GO:0006915")
GO:GO:0006915

julia> parse(GO, "0006915")
GO:GO:0006915
```
"""
struct GO <: BioIdentifier
    id::UInt32
end

@numericid GO url = "purl.obolibrary.org/obo/GO_" prefix = "GO:" digits = 7 inttype = UInt32

function parseid(::Type{GO}, id::SubString)
    _, id = chopprefixes(id, "purl.obolibrary.org/obo/GO_", "GO:")
    all(isdigit, id) && length(id) == 7 || return MalformedIdentifier{GO}(id, "must be 7 digits (optionally prefixed with GO:)")
    num = parsefor(GO, UInt32, id)
    num isa UInt32 || return num
    GO(num)
end

function purl(go::GO)
    digits = lpad(go.id, 7, '0')
    "https://purl.obolibrary.org/obo/GO_$(digits)"
end



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
struct GEO <: BioIdentifier
    type::UInt8
    number::UInt64
end

const GEO_PREFIXES = ("gpl", "gsm", "gse", "gds")

function parseid(::Type{GEO}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", "ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=")
    for (i, prefix) in enumerate(GEO_PREFIXES)
        found, remainder = chopprefixes(id, prefix)
        if found
            number = parsefor(GEO, UInt64, remainder)
            number isa UInt64 || return number
            return GEO(UInt8(i), number)
        end
    end
    return MalformedIdentifier{GEO}(id, "must start with GPL, GSM, GSE, or GDS")
end

function purl(geo::GEO)
    id = shortcode(geo)
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$(id)"
end

shortcode(geo::GEO) = string(uppercase(GEO_PREFIXES[geo.type]), geo.number)
idcode(geo::GEO) = geo.number

Base.show(io::IO, geo::GEO) = (show(io, GEO); show(io, (uppercase(GEO_PREFIXES[geo.type]), geo.number)))


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
HGNC:HGNC:5

julia> parse(HGNC, "5")
HGNC:HGNC:5
```
"""
struct HGNC <: BioIdentifier
    id::UInt32
end

@numericid HGNC url = "www.genenames.org/data/gene-symbol-report/#!/hgnc_id/" prefix = "HGNC:" inttype = UInt32




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

function parseid(::Type{INSDC}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "www.ncbi.nlm.nih.gov/nuccore/", "ncbi.nlm.nih.gov/nuccore/", "www.ebi.ac.uk/ena/browser/view/", "ebi.ac.uk/ena/browser/view/")

    # Extract version if present
    base_id, version = if contains(id, '.')
        parts = split(id, '.', limit=2)
        version_str = parts[2]
        isempty(version_str) && return MalformedIdentifier{INSDC}(id, "version cannot be empty after period")
        parsed_version = parsefor(INSDC, UInt16, version_str)
        parsed_version isa UInt16 || return parsed_version
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

    prefix_match === nothing && return MalformedIdentifier{INSDC}(id, "must start with valid INSDC prefix followed by digits")

    remainder = base_id[length(prefix_match)+1:end]
    num = parsefor(INSDC, UInt64, remainder)
    num isa UInt64 || return num

    INSDC(SubString(prefix_match), num, version)
end

purlprefix(::Type{INSDC}) = "https://www.ncbi.nlm.nih.gov/nuccore/"

function shortcode(insdc::INSDC)
    base = "$(insdc.prefix)$(insdc.number)"
    isnothing(insdc.version) ? base : "$(base).$(insdc.version)"
end

idcode(insdc::INSDC) = insdc.number

Base.show(io::IO, insdc::INSDC) = (show(io, INSDC); show(io, (insdc.prefix, insdc.number, insdc.version)))


# HMDB

"""
    HMDB <: BioIdentifier

Human Metabolome Database entry identifier.

The Human Metabolome Database contains detailed information about small molecule
metabolites found in the human body. It is intended to be used for applications
in metabolomics, clinical chemistry, and general biomarker discovery. Identifiers
follow format `HMDB<7 digits>`. Parsing may throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(HMDB, "HMDB0000001")
HMDB:HMDB0000001

julia> parse(HMDB, "https://hmdb.ca/metabolites/HMDB0123456")
HMDB:HMDB0123456
```
"""
struct HMDB <: BioIdentifier
    id::UInt32
end

@numericid HMDB url = "hmdb.ca/metabolites/" prefix = "HMDB" digits = 7 inttype = UInt32

function parseid(::Type{HMDB}, id::SubString)
    _, id = chopprefixes(id, "hmdb.ca/metabolites/")
    found, remainder = chopprefixes(id, "hmdb")
    found || return MalformedIdentifier{HMDB}(id, "must start with 'HMDB'")
    length(remainder) == 7 && all(isdigit, remainder) || return MalformedIdentifier{HMDB}(id, "must have exactly 7 digits after 'HMDB'")
    num = parsefor(HMDB, UInt32, remainder)
    num isa UInt32 || return num
    HMDB(num)
end



# HPO

"""
    HPO <: BioIdentifier

Human Phenotype Ontology term identifier.

The Human Phenotype Ontology provides a standardized vocabulary of phenotypic
abnormalities encountered in human disease. Each term describes a phenotypic
feature of human disease. Identifiers follow format `HP:<7 digits>`. Parsing may
throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(HPO, "HP:0000001")
HPO:HP:0000001

julia> parse(HPO, "0000001")
HPO:HP:0000001
```
"""
struct HPO <: BioIdentifier
    id::UInt32
end

function parseid(::Type{HPO}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "purl.obolibrary.org/obo/HP_", "HP:")
    all(isdigit, id) && length(id) == 7 || return MalformedIdentifier{HPO}(id, "must be 7 digits (optionally prefixed with HP:)")
    num = parsefor(HPO, UInt32, id)
    num isa UInt32 || return num
    HPO(num)
end

function purl(hpo::HPO)
    digits = lpad(hpo.id, 7, '0')
    "https://purl.obolibrary.org/obo/HP_$(digits)"
end

shortcode(hpo::HPO) = string("HP:", lpad(hpo.id, 7, '0'))
idcode(hpo::HPO) = hpo.id


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
struct IntAct <: BioIdentifier
    id::UInt32
end

@numericid IntAct url = "www.ebi.ac.uk/intact/interaction/" prefix = "EBI-" inttype = UInt32





# KEGG

"""
    KEGG <: BioIdentifier

KEGG pathway map identifier.

KEGG is a database resource for understanding high-level functions and utilities
of the biological system from molecular-level information. Pathway maps provide
graphical representations of cellular processes. Identifiers follow format
`<3-4 letters><5 digits>`. Parsing may throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(KEGG, "hsa00010")
KEGG:hsa00010

julia> parse(KEGG, "https://www.kegg.jp/entry/mmu04210")
KEGG:mmu04210
```
"""
struct KEGG <: BioIdentifier
    species::UInt32  # Encoded 3-4 letter species code
    pathway::UInt32  # 5-digit pathway number
end

function parseid(::Type{KEGG}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "www.kegg.jp/entry/", "kegg.jp/entry/")
    (8 <= length(id) <= 9) || return MalformedIdentifier{KEGG}(id, "must be 3-4 lowercase letters followed by 5 digits")
    letter_part = id[1:end-5]
    digit_part = id[end-4:end]
    (3 <= length(letter_part) <= 4) && all(c -> islowercase(c) && isletter(c), letter_part) || return MalformedIdentifier{KEGG}(id, "first part must be 3-4 lowercase letters")
    length(digit_part) == 5 && all(isdigit, digit_part) || return MalformedIdentifier{KEGG}(id, "last part must be exactly 5 digits")
    # Encode 3-4 letter species code as UInt32 using ASCII
    species_code = UInt32(0)
    for (i, c) in enumerate(codeunits(letter_part))
        species_code |= UInt32(c) << (8 * (4 - i))  # Pack ASCII values into 8-bit slots
    end
    pathway = parsefor(KEGG, UInt32, digit_part)
    pathway isa UInt32 || return pathway
    KEGG(species_code, pathway)
end

purlprefix(::Type{KEGG}) = "https://www.kegg.jp/entry/"
function shortcode(kegg::KEGG)
    # Decode species code from UInt32
    species_bytes = Vector{UInt8}(undef, 4)
    temp = kegg.species
    len = 0
    for i in 1:4
        byte = (temp >> (8 * (4 - i))) & 0xff
        byte == 0 && break
        species_bytes[i] = byte
        len = i
    end
    string(String(resize!(species_bytes, len)), lpad(kegg.pathway, 5, '0'))
end
idcode(kegg::KEGG) = kegg.pathway

Base.show(io::IO, kegg::KEGG) = (show(io, KEGG); show(io, (kegg.species, kegg.pathway)))


# MetaboLights

"""
    MetaboLights <: BioIdentifier

MetaboLights metabolomics study archive identifier.

MetaboLights is a database for metabolomics experiments and derived information.
The database supports the public dissemination of metabolomics data in support
of the scientific endeavor. Identifiers follow format `MTBLS<number>`. Parsing may
throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(MetaboLights, "MTBLS1")
MetaboLights:MTBLS1

julia> parse(MetaboLights, "https://www.ebi.ac.uk/metabolights/MTBLS123")
MetaboLights:MTBLS123
```
"""
struct MetaboLights <: BioIdentifier
    id::UInt32
end

@numericid MetaboLights url = "www.ebi.ac.uk/metabolights/" prefix = "MTBLS" inttype = UInt32





# MONDO

"""
    MONDO <: BioIdentifier

Monarch Disease Ontology unified disease identifier.

The Monarch Disease Ontology aims to harmonize disease definitions across the world.
It is a semi-automatically constructed ontology that merges in multiple disease
resources to yield a coherent merged ontology. Identifiers follow format `MONDO:<number>`.
Parsing may throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(MONDO, "MONDO:0000001")
MONDO:MONDO:0000001

julia> parse(MONDO, "0000001")
MONDO:MONDO:0000001
```
"""
struct MONDO <: BioIdentifier
    id::UInt32
end

@numericid MONDO url = "purl.obolibrary.org/obo/MONDO_" prefix = "MONDO:" digits = 7 inttype = UInt32

function purl(mondo::MONDO)
    digits = lpad(mondo.id, 7, '0')
    "https://purl.obolibrary.org/obo/MONDO_$(digits)"
end



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
struct NCBIGene <: BioIdentifier
    id::UInt32
end

@numericid NCBIGene url = "www.ncbi.nlm.nih.gov/gene/" inttype = UInt32

function parseid(::Type{NCBIGene}, id::SubString)
    _, id = chopprefixes(id, "www.ncbi.nlm.nih.gov/gene/", "ncbi.nlm.nih.gov/gene/")
    length(id) <= 10 && all(isdigit, id) && !isempty(id) || return MalformedIdentifier{NCBIGene}(id, "must be a positive integer with ≤ 10 digits")
    number = parsefor(NCBIGene, UInt32, id)
    number isa UInt32 || return number
    number == 0 && return MalformedIdentifier{NCBIGene}(id, "must be a positive integer")
    NCBIGene(number)
end



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
struct NCBITaxon <: BioIdentifier
    id::UInt32
end

@numericid NCBITaxon url = ["www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=", "taxon:"] inttype = UInt32



# OMIM

"""
    OMIM <: BioIdentifier

Online Mendelian Inheritance in Man disorder/gene entry identifier.

OMIM is a continuously updated catalog of human genes and genetic disorders and
traits, with a particular focus on the molecular relationship between genetic
variation and phenotypic expression. Identifiers follow format `<6 digits>`
optionally with `.<4 digits>` suffix. Parsing may throw a `MalformedIdentifier`
if the format is invalid.

# Examples

```julia
julia> parse(OMIM, "600918")
OMIM:600918

julia> parse(OMIM, "613659.0001")
OMIM:613659.0001
```
"""
struct OMIM <: BioIdentifier
    main::UInt32
    suffix::Union{UInt16,Nothing}
end

function parseid(::Type{OMIM}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "omim.org/entry/")
    suffix = nothing
    if contains(id, ".")
        parts = split(id, ".")
        length(parts) == 2 || return MalformedIdentifier{OMIM}(id, "invalid format")
        id, suffix_str = parts
        length(suffix_str) == 4 && all(isdigit, suffix_str) || return MalformedIdentifier{OMIM}(id, "suffix must be exactly 4 digits")
        suffix_num = parsefor(OMIM, UInt16, suffix_str)
        suffix_num isa UInt16 || return suffix_num
        suffix = suffix_num
    end
    length(id) == 6 && all(isdigit, id) || return MalformedIdentifier{OMIM}(id, "main identifier must be exactly 6 digits")
    main_num = parsefor(OMIM, UInt32, id)
    main_num isa UInt32 || return main_num
    OMIM(main_num, suffix)
end

purlprefix(::Type{OMIM}) = "https://omim.org/entry/"
shortcode(omim::OMIM) = isnothing(omim.suffix) ? string(lpad(omim.main, 6, '0')) : "$(lpad(omim.main, 6, '0')).$(lpad(omim.suffix, 4, '0'))"
idcode(omim::OMIM) = omim.main

Base.show(io::IO, omim::OMIM) = (show(io, OMIM); show(io, (omim.main, omim.suffix)))


# PDB

"""
    PDB <: BioIdentifier

Protein Data Bank structure identifier.

The Protein Data Bank is the single worldwide archive of structural data
of biological macromolecules. Identifiers are 4-character alphanumeric codes
(digits and uppercase letters). Parsing may throw a `MalformedIdentifier`
if the format is invalid.

# Examples

```julia
julia> parse(PDB, "1A3N")
PDB:1A3N

julia> parse(PDB, "7QDH")
PDB:7QDH
```
"""
struct PDB <: BioIdentifier
    id::SubString{String}
end

function parseid(::Type{PDB}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "www.ebi.ac.uk/pdbe/entry/pdb/", "ebi.ac.uk/pdbe/entry/pdb/", "www.rcsb.org/structure/", "rcsb.org/structure/")
    length(id) == 4 || return MalformedIdentifier{PDB}(id, "must be exactly 4 characters")
    all(c -> isdigit(c) || (isuppercase(c) && isletter(c)), id) || return MalformedIdentifier{PDB}(id, "must contain only digits and uppercase letters")
    PDB(id)
end

purlprefix(::Type{PDB}) = "https://www.ebi.ac.uk/pdbe/entry/pdb/"
shortcode(pdb::PDB) = pdb.id
idcode(pdb::PDB) = pdb.id

Base.show(io::IO, pdb::PDB) = (show(io, PDB); print(io, ':', pdb.id))


# PubChem

"""
    PubChem <: BioIdentifier

PubChem database identifier for chemical compounds, substances, and assays.

PubChem is the largest freely accessible chemical database, providing information
on chemical structures, properties, and biological activities. Identifiers
follow formats CID<number>, SID<number>, or AID<number> for compounds,
substances, and assays respectively. Parsing may throw a `MalformedIdentifier`
if the format is invalid.

# Examples

```julia
julia> parse(PubChem, "2244")
PubChem:CID:2244

julia> parse(PubChem, "SID12345678")
PubChem:SID:12345678
```
"""
struct PubChem <: BioIdentifier
    type::UInt8  # 1=CID, 2=SID, 3=AID
    number::UInt64
end

const PUBCHEM_TYPES = ("CID", "SID", "AID")

function parseid(::Type{PubChem}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "pubchem.ncbi.nlm.nih.gov/compound/", "pubchem.ncbi.nlm.nih.gov/substance/", "pubchem.ncbi.nlm.nih.gov/assay/")

    for (i, prefix) in enumerate(PUBCHEM_TYPES)
        if startswith(id, prefix)
            remainder = id[length(prefix)+1:end]
            all(isdigit, remainder) || return MalformedIdentifier{PubChem}(id, "$(prefix) must contain only digits after prefix")
            isempty(remainder) && return MalformedIdentifier{PubChem}(id, "must have a number after $(prefix)")
            num = parsefor(PubChem, UInt64, remainder)
            num isa UInt64 || return num
            num == 0 && return MalformedIdentifier{PubChem}(id, "number must be positive")
            return PubChem(UInt8(i), num)
        end
    end

    # Default to CID (just digits)
    all(isdigit, id) || return MalformedIdentifier{PubChem}(id, "CID must contain only digits")
    isempty(id) && return MalformedIdentifier{PubChem}(id, "cannot be empty")
    num = parsefor(PubChem, UInt64, id)
    num isa UInt64 || return num
    num == 0 && return MalformedIdentifier{PubChem}(id, "number must be positive")
    PubChem(UInt8(1), num)
end

function purl(pubchem::PubChem)
    type_name = lowercase(PUBCHEM_TYPES[pubchem.type])
    "https://pubchem.ncbi.nlm.nih.gov/$(type_name)/$(pubchem.number)"
end

function shortcode(pubchem::PubChem)
    pubchem.type == 1 ? string(pubchem.number) : "$(PUBCHEM_TYPES[pubchem.type])$(pubchem.number)"
end

idcode(pubchem::PubChem) = pubchem.number

Base.show(io::IO, pubchem::PubChem) = (show(io, PubChem); print(io, ':', PUBCHEM_TYPES[pubchem.type], ':', pubchem.number))


# PXD

"""
    PXD <: BioIdentifier

ProteomeXchange public mass spectrometry dataset identifier.

ProteomeXchange provides a single point of submission of mass spectrometry
proteomics data to the main existing proteomics repositories, and a unified
identifier space. Identifiers follow formats PXD/RPXD<6 digits> for original/
re-processed datasets respectively. Parsing may throw a `MalformedIdentifier`
if the format is invalid.

# Examples

```julia
julia> parse(PXD, "PXD000001")
PXD:PXD000001

julia> parse(PXD, "RPXD123456")
PXD:RPXD123456
```
"""
struct PXD <: BioIdentifier
    reprocessed::Bool  # false=PXD, true=RPXD
    number::UInt32
end

function parseid(::Type{PXD}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "proteomecentral.proteomexchange.org/cgi/GetDataset?ID=")
    found_rpxd, remainder_rpxd = chopprefixes(id, "rpxd")
    if found_rpxd
        length(remainder_rpxd) == 6 && all(isdigit, remainder_rpxd) || return MalformedIdentifier{PXD}(id, "must have exactly 6 digits after prefix")
        num = parsefor(PXD, UInt32, remainder_rpxd)
        num isa UInt32 || return num
        return PXD(true, num)
    end
    found_pxd, remainder_pxd = chopprefixes(id, "pxd")
    if found_pxd
        length(remainder_pxd) == 6 && all(isdigit, remainder_pxd) || return MalformedIdentifier{PXD}(id, "must have exactly 6 digits after prefix")
        num = parsefor(PXD, UInt32, remainder_pxd)
        num isa UInt32 || return num
        return PXD(false, num)
    end
    return MalformedIdentifier{PXD}(id, "must start with PXD or RPXD")
end

function purl(pxd::PXD)
    id = shortcode(pxd)
    "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=$(id)"
end

shortcode(pxd::PXD) = string(pxd.reprocessed ? "RPXD" : "PXD", lpad(pxd.number, 6, '0'))
idcode(pxd::PXD) = pxd.number

Base.show(io::IO, pxd::PXD) = (show(io, PXD); show(io, (pxd.reprocessed ? "RPXD" : "PXD", pxd.number)))


# Reactome

"""
    Reactome <: BioIdentifier

Reactome manually curated pathway identifier.

Reactome is a free, open-source, curated and peer-reviewed pathway database.
The goal is to provide intuitive bioinformatics tools for the visualization,
interpretation and analysis of pathway knowledge. Identifiers follow format
`R-<3 letters>-<number>`. Parsing may throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(Reactome, "R-HSA-199420")
Reactome:R-HSA-199420

julia> parse(Reactome, "https://reactome.org/PathwayBrowser/#/R-MMU-123456")
Reactome:R-MMU-123456
```
"""
struct Reactome <: BioIdentifier
    species::UInt32  # Encoded 3-letter species code (24 bits)
    number::UInt32
end

function parseid(::Type{Reactome}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "reactome.org/PathwayBrowser/#/")
    found, remainder = chopprefixes(id, "r-")
    found || return MalformedIdentifier{Reactome}(id, "must start with 'R-'")
    parts = split(remainder, '-')
    length(parts) == 2 || return MalformedIdentifier{Reactome}(id, "must follow format R-<species>-<number>")
    species_str = SubString(string(parts[1]))
    number_str = parts[2]
    length(species_str) == 3 || return MalformedIdentifier{Reactome}(id, "species code must be exactly 3 characters")

    # Only uppercase if necessary
    species_final = if all(c -> c ∈ 'A':'Z', species_str)
        species_str
    elseif all(c -> c ∈ 'a':'z', species_str)
        SubString(uppercase(species_str))
    elseif all(c -> c ∈ 'A':'Z' || c ∈ 'a':'z', species_str)
        SubString(uppercase(species_str))
    else
        return MalformedIdentifier{Reactome}(id, "species code must contain only letters")
    end

    number = parsefor(Reactome, UInt32, number_str)
    number isa UInt32 || return number
    # Encode 3-letter species code into UInt32 (8 bits per char)
    species_code = UInt32(codeunit(species_final, 1)) << 16 | UInt32(codeunit(species_final, 2)) << 8 | UInt32(codeunit(species_final, 3))
    Reactome(species_code, number)
end

function purl(reactome::Reactome)
    id = shortcode(reactome)
    "https://reactome.org/PathwayBrowser/#/$(id)"
end

function shortcode(reactome::Reactome)
    # Decode species from UInt32
    c1 = Char((reactome.species >> 16) & 0xff)
    c2 = Char((reactome.species >> 8) & 0xff)
    c3 = Char(reactome.species & 0xff)
    "R-$(c1)$(c2)$(c3)-$(reactome.number)"
end
idcode(reactome::Reactome) = reactome.number

Base.show(io::IO, reactome::Reactome) = (show(io, Reactome); show(io, (reactome.species, reactome.number)))


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
struct RefSeq <: BioIdentifier
    width::UInt8
    prefix::UInt8
    version::UInt16
    number::UInt32
end

const REFSEQ_PREFIXES = ("nm", "np", "xm", "xp", "nr", "xr", "ng", "nt", "nw", "nc")

function parseid(::Type{RefSeq}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "www.ncbi.nlm.nih.gov/nuccore/", "ncbi.nlm.nih.gov/nuccore/")
    parts = split(id, '.')
    length(parts) == 2 || return MalformedIdentifier{RefSeq}(id, "must include version number after dot")
    main_part, version_str = parts
    version = parsefor(RefSeq, UInt16, version_str)
    version isa UInt16 || return version
    underscore_pos = findfirst('_', main_part)
    isnothing(underscore_pos) && return MalformedIdentifier{RefSeq}(id, "must contain underscore separator")
    prefix_str = main_part[1:underscore_pos-1]
    digits = main_part[underscore_pos+1:end]
    length(prefix_str) >= 2 && all(isletter, prefix_str) || return MalformedIdentifier{RefSeq}(id, "prefix must be letters followed by underscore")
    w = ncodeunits(digits)
    num = parsefor(RefSeq, UInt32, digits)
    num isa UInt32 || return num
    prefix_match = findfirst(==(prefix_str), REFSEQ_PREFIXES)
    if isnothing(prefix_match)
        return MalformedIdentifier{RefSeq}(id, "unsupported RefSeq prefix: $(prefix_str)")
    end
    # Accept only 6/8/9; infer correct width for irregular input
    if w ∉ (6, 8, 9)
        w = if num < 1_000_000; 6
            elseif num < 100_000_000; 8
            else 9 end
    end
    RefSeq(UInt8(w), UInt8(prefix_match), version, num)
end

purlprefix(::Type{RefSeq}) = "https://www.ncbi.nlm.nih.gov/nuccore/"
shortcode(refseq::RefSeq) = string(uppercase(REFSEQ_PREFIXES[refseq.prefix]), '_', lpad(refseq.number, refseq.width, '0'), '.', refseq.version)
idcode(refseq::RefSeq) = refseq.number

Base.show(io::IO, refseq::RefSeq) = (show(io, RefSeq); show(io, (uppercase(REFSEQ_PREFIXES[refseq.prefix]), refseq.number, refseq.width, refseq.version)))


# RRID

"""
    RRID <: BioIdentifier

Research Resource Identifier.

Research Resource Identifiers are persistent unique identifiers for research
resources. RRIDs are designed to be used in the methods sections of research
articles to unambiguously identify the research resources used. Identifiers
follow format `RRID:<provider>:<id>`. Parsing may throw a `MalformedIdentifier`
if the format is invalid.

# Examples

```julia
julia> parse(RRID, "RRID:IMSR_JAX:000664")
RRID:RRID:IMSR_JAX:000664

julia> parse(RRID, "https://scicrunch.org/resolver/RRID:AB_123456")
RRID:RRID:AB_123456
```
"""
struct RRID <: BioIdentifier
    provider::String
    resource_id::String
end

function parseid(::Type{RRID}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "scicrunch.org/resolver/", "identifiers.org/rrid/")

    has_rrid, id = chopprefixes(id, "rrid:")

    colon_pos = findfirst(':', id)
    colon_pos === nothing && return MalformedIdentifier{RRID}(id, "must contain provider and resource separated by colon")

    provider = @view id[1:colon_pos-1]
    resource_id = @view id[colon_pos+1:end]

    isempty(provider) && return MalformedIdentifier{RRID}(id, "provider cannot be empty")
    isempty(resource_id) && return MalformedIdentifier{RRID}(id, "resource ID cannot be empty")

    all(c -> isalnum(c) || c == '_', provider) || return MalformedIdentifier{RRID}(id, "provider must contain only letters, numbers, and underscore")
    all(c -> isalnum(c) || c == '_', resource_id) || return MalformedIdentifier{RRID}(id, "resource ID must contain only letters, numbers, and underscore")

    RRID(string(provider), string(resource_id))
end

function purl(rrid::RRID)
    id = shortcode(rrid)
    "https://scicrunch.org/resolver/$(id)"
end

shortcode(rrid::RRID) = "RRID:$(rrid.provider):$(rrid.resource_id)"

Base.show(io::IO, rrid::RRID) = (show(io, RRID); show(io, (rrid.provider, rrid.resource_id)))


# SRA

"""
    SRA <: BioIdentifier

Sequence Read Archive identifier.

The Sequence Read Archive is a bioinformatics database that provides a public
repository for DNA sequencing data, especially the read data generated by
next-generation sequencing technologies. Identifiers follow formats SRR/SRX/SRS/SRP<number>
for run/experiment/sample/study records respectively. Parsing may throw a
`MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(SRA, "SRR123456")
SRA:SRR123456

julia> parse(SRA, "SRX789012")
SRA:SRX789012
```
"""
struct SRA <: BioIdentifier
    type::UInt8
    number::UInt64
end

const SRA_PREFIXES = ("srr", "srx", "srs", "srp")

function parseid(::Type{SRA}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "www.ncbi.nlm.nih.gov/sra/", "ncbi.nlm.nih.gov/sra/")
    for (i, prefix) in enumerate(SRA_PREFIXES)
        found, remainder = chopprefixes(id, prefix)
        if found
            number = parsefor(SRA, UInt64, remainder)
            number isa UInt64 || return number
            return SRA(UInt8(i), number)
        end
    end
    return MalformedIdentifier{SRA}(id, "must start with SRR, SRX, SRS, or SRP")
end

purlprefix(::Type{SRA}) = "https://www.ncbi.nlm.nih.gov/sra/"
shortcode(sra::SRA) = string(uppercase(SRA_PREFIXES[sra.type]), lpad(sra.number, max(6, ndigits(sra.number)), '0'))
idcode(sra::SRA) = sra.number

Base.show(io::IO, sra::SRA) = (show(io, SRA); show(io, (uppercase(SRA_PREFIXES[sra.type]), sra.number)))


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
SO:SO:0001234

julia> parse(SO, "0000704")
SO:SO:0000704
```
"""
struct SO <: BioIdentifier
    id::UInt32
end

@numericid SO url = "purl.obolibrary.org/obo/SO_" prefix = "SO:" digits = 7 inttype = UInt32

function parseid(::Type{SO}, id::SubString)
    _, id = chopprefixes(id, "purl.obolibrary.org/obo/SO_", "SO:")
    length(id) == 7 && all(isdigit, id) || return MalformedIdentifier{SO}(id, "must be 7 digits (optionally prefixed with SO:)")
    num = parsefor(SO, UInt32, id)
    num isa UInt32 || return num
    SO(num)
end

function purl(so::SO)
    digits = lpad(so.id, 7, '0')
    "https://purl.obolibrary.org/obo/SO_$(digits)"
end



# UniProt

"""
    UniProt <: BioIdentifier

UniProt canonical protein sequence identifier.

UniProt is a comprehensive, high-quality and freely accessible resource of
protein sequence and functional information. Identifiers are either 6-character
(Swiss-Prot) or 10-character (TrEMBL) alphanumeric codes following specific
patterns. Parsing may throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(UniProt, "P04637")
UniProt:P04637

julia> parse(UniProt, "https://www.uniprot.org/uniprot/Q9Y6K5")
UniProt:Q9Y6K5
```
"""
struct UniProt <: BioIdentifier
    id::SubString{String}
end

function parseid(::Type{UniProt}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "www.uniprot.org/uniprot/", "uniprot.org/uniprot/", "uniprot:")
    if length(id) == 6
        (id[1] in "OPQ") && isdigit(id[2]) && all(c -> isdigit(c) || (isuppercase(c) && isletter(c)), id[3:5]) && isdigit(id[6]) || return MalformedIdentifier{UniProt}(id, "6-character format must match [OPQ][0-9][A-Z0-9]{3}[0-9]")
    elseif length(id) == 10
        c1 = id[1]
        (c1 ∈ 'A':'N' || c1 ∈ 'R':'Z' || c1 ∈ 'a':'n' || c1 ∈ 'r':'z') || return MalformedIdentifier{UniProt}(id, "first character must be A-N or R-Z")
        isdigit(id[2]) || return MalformedIdentifier{UniProt}(id, "second character must be a digit")
        for i in 3:10
            c = id[i]
            (isdigit(c) || c ∈ 'A':'Z' || c ∈ 'a':'z') || return MalformedIdentifier{UniProt}(id, "characters 3-10 must be letters or digits")
        end
    else
        return MalformedIdentifier{UniProt}(id, "must be 6 or 10 alphanumeric characters")
    end
    UniProt(id)
end

purlprefix(::Type{UniProt}) = "https://www.uniprot.org/uniprot/"
shortcode(up::UniProt) = String(up.id)

Base.show(io::IO, up::UniProt) = (show(io, UniProt); show(io, (up.id,)))



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
struct UniRef <: BioIdentifier
    level::UInt8  # 50, 90, or 100
    uniprot::SubString{String}
end

function parseid(::Type{UniRef}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "www.uniprot.org/uniref/", "uniprot.org/uniref/")
    startswith(id, "UniRef") || return MalformedIdentifier{UniRef}(id, "must start with 'UniRef'")
    remainder = id[7:end]  # Skip "UniRef"
    parts = split(remainder, '_', limit=2)
    length(parts) == 2 || return MalformedIdentifier{UniRef}(id, "must follow format UniRef<level>_<UniProt>")
    level_str, uniprot_str = parts[1], parts[2]

    # Validate clustering level
    level_str ∈ ("50", "90", "100") || return MalformedIdentifier{UniRef}(id, "clustering level must be 50, 90, or 100")
    level = parse(UInt8, level_str)

    # Validate UniProt ID format (basic check)
    isempty(uniprot_str) && return MalformedIdentifier{UniRef}(id, "UniProt ID cannot be empty")
    (length(uniprot_str) == 6 || length(uniprot_str) == 10) || return MalformedIdentifier{UniRef}(id, "UniProt ID must be 6 or 10 characters")
    all(c -> isdigit(c) || (isuppercase(c) && isletter(c)), uniprot_str) || return MalformedIdentifier{UniRef}(id, "UniProt ID must contain only uppercase alphanumeric characters")

    UniRef(level, SubString(uniprot_str))
end

purlprefix(::Type{UniRef}) = "https://www.uniprot.org/uniref/"
shortcode(uniref::UniRef) = "UniRef$(uniref.level)_$(uniref.uniprot)"
idcode(uniref::UniRef) = uniref.uniprot

Base.show(io::IO, uniref::UniRef) = (show(io, UniRef); show(io, (uniref.level, uniref.uniprot)))


# WikiPathways

"""
    WikiPathways <: BioIdentifier

WikiPathways community pathway database identifier.

WikiPathways is a community-curated pathway database providing biological
pathway information. Identifiers follow format WP<digits> and are commonly
used in pathway enrichment analysis tools like Cytoscape and GSEA.
Parsing may throw a `MalformedIdentifier` if the format is invalid.

# Examples

```julia
julia> parse(WikiPathways, "WP554")
WikiPathways:WP554

julia> parse(WikiPathways, "WP1234")
WikiPathways:WP1234
```
"""
struct WikiPathways <: BioIdentifier
    id::UInt32
end

@numericid WikiPathways url = "www.wikipathways.org/instance/" prefix = "WP" inttype = UInt32





# AFDP

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
struct AFDB <: BioIdentifier
    uniprot::UniProt
    fragment::UInt8
end

function parseid(::Type{AFDB}, id::SubString)
    _, id = chopprefixes(id, "https://", "http://", "alphafold.com/entry/", "alphafold.ebi.ac.uk/entry/")
    startswith(id, "AF-") || return MalformedIdentifier{AFDB}(id, "must start with 'AF-'")
    parts = split(id[4:end], '-')
    length(parts) == 2 && startswith(parts[2], 'F') || return MalformedIdentifier{AFDB}(id, "must follow format AF-<UniProt>-F<number>")
    uniprot_part = parts[1]
    fragment_str = parts[2][2:end]
    # Validate using actual UniProt parser instead of loose checks
    uniprot_result = parseid(UniProt, SubString(uniprot_part))
    uniprot_result isa UniProt || return MalformedIdentifier{AFDB}(id, "invalid UniProt ID: $(uniprot_part)")
    fragment = parsefor(AFDB, UInt8, fragment_str)  # Use UInt8 for fragment numbers ≤ 255
    fragment isa UInt8 || return fragment
    AFDB(uniprot_result, fragment)
end

purlprefix(::Type{AFDB}) = "https://alphafold.com/entry/"
shortcode(afdb::AFDB) = "AF-$(shortcode(afdb.uniprot))-F$(afdb.fragment)"
idcode(afdb::AFDB) = afdb.fragment

Base.show(io::IO, afdb::AFDB) = (show(io, AFDB); show(io, (shortcode(afdb.uniprot), afdb.fragment)))

end
