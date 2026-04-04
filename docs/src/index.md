# BioIdentifiers.jl

Structured and validated types for biological and biomedical identifiers.

Built on [FastIdentifiers.jl](https://code.tecosaur.net/tec/FastIdentifiers), each
type provides parsing (from shortcodes, prefixed forms, and URLs), canonical string
output, persistent URL generation, and round-trip guarantees.

## Index

```@index
```

## Identifier Types

### Protein & Structure

```@docs
AFDB
PDB
UniProt
UniRef
IntAct
InterPro
Pfam
PXD
```

### Genomics & Genetics

Ensembl identifiers (`ENSG`, `ENST`, `ENSP`, `ENSE`, `ENSR`, `ENSF`, `ENSFM`)
are typed aliases for `EnsemblIdentifier{T}` — see `?ENSG` for details.

```@docs
NCBIGene
RefSeq
HGNC
INSDC
OMIM
```

### Variation & Clinical

```@docs
CA
ClinVar
dbSNP
dbVar
GWAS
```

### Expression & Functional Genomics

```@docs
ArrayExpress
GEO
```

### Studies & Samples

```@docs
BioProject
BioSample
ClinicalTrials
dbGaP
EGA
RRID
SRA
```

### Ontologies & Controlled Vocabularies

```@docs
CL
DOID
ECO
EFO
GO
HPO
MeSH
MONDO
MP
PATO
SO
UBERON
```

### Chemical & Metabolic

```@docs
ChEBI
ChEMBL
DrugBank
HMDB
KEGG
MetaboLights
PubChem
```

### Networks & Interactions

```@docs
BioGRID
Reactome
WikiPathways
```

### Cell Lines & Model Organisms

```@docs
Cellosaurus
FlyBase
MGI
NCBITaxon
SGD
WormBase
```

## Abstract Types

```@docs
BioIdentifier
```
