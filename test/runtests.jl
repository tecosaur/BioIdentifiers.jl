# SPDX-FileCopyrightText: © 2025 TEC <contact@tecosaur.net>
# SPDX-License-Identifier: MPL-2.0

using BioIdentifiers
using BioIdentifiers: MalformedIdentifier, shortcode, purl, purlprefix, idcode
using Test

# -- Test helpers --

"""Test shortcode, string, and purl for each (input, shortcode, string, purl) tuple.
Also verifies parse round-trips through shortcode, string, and (when URL
parsing is supported via `purlprefix`) the PURL."""
function check_formats(T, examples::Vector{<:Tuple})
    has_purlprefix = !isnothing(purlprefix(T))
    for (input, sc, str, p) in examples
        id = parse(T, input)
        @test shortcode(id) == sc
        @test string(id) == str
        if isnothing(p)
            @test purl(id) === nothing
        else
            @test purl(id) == p
        end
        # Round-trip through each form
        @test parse(T, sc) == id
        @test parse(T, str) == id
        has_purlprefix && !isnothing(p) && @test parse(T, p) == id
    end
end

"""Test that each input string fails to parse."""
function check_malformed(T, inputs)
    for input in inputs
        @test tryparse(T, input) === nothing
    end
end

# -- Identifier tests (alphabetical) --

@testset "AFDB" begin
    check_formats(AFDB, [
        ("AF-P12345-F1",  "AF-P12345-F1",  "AF-P12345-F1",  "https://alphafold.com/entry/AF-P12345-F1"),
        ("AF-Q9UHC1-F255","AF-Q9UHC1-F255","AF-Q9UHC1-F255","https://alphafold.com/entry/AF-Q9UHC1-F255"),
    ])
    check_malformed(AFDB, ["", "AF-", "P12345"])
end

@testset "ArrayExpress" begin
    check_formats(ArrayExpress, [
        ("E-MTAB-1234",  "E-MTAB-1234",  "E-MTAB-1234",  "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-1234"),
        ("E-GEOD-5678",  "E-GEOD-5678",  "E-GEOD-5678",  "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-GEOD-5678"),
        ("E-TABM-12345", "E-TABM-12345", "E-TABM-12345", "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-TABM-12345"),
    ])
    check_malformed(ArrayExpress, ["", "invalid", "E-AB-1"])
end

@testset "BioGRID" begin
    check_formats(BioGRID, [
        ("123456",    "123456",    "123456",    "https://thebiogrid.org/123456"),
        ("1",         "1",         "1",         "https://thebiogrid.org/1"),
        ("999999999", "999999999", "999999999", "https://thebiogrid.org/999999999"),
    ])
    check_malformed(BioGRID, ["", "abc"])
end

@testset "BioProject" begin
    check_formats(BioProject, [
        ("PRJNA123456", "PRJNA123456", "PRJNA123456", "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA123456"),
        ("PRJEB123456", "PRJEB123456", "PRJEB123456", "https://www.ncbi.nlm.nih.gov/bioproject/PRJEB123456"),
        ("PRJDB123456", "PRJDB123456", "PRJDB123456", "https://www.ncbi.nlm.nih.gov/bioproject/PRJDB123456"),
    ])
    check_malformed(BioProject, ["", "PRJ123"])
end

@testset "BioSample" begin
    check_formats(BioSample, [
        ("SAMN12345678",  "SAMN12345678",  "SAMN12345678",  "https://www.ncbi.nlm.nih.gov/biosample/SAMN12345678"),
        ("SAMEA12345678", "SAMEA12345678", "SAMEA12345678", "https://www.ncbi.nlm.nih.gov/biosample/SAMEA12345678"),
        ("SAMD12345678",  "SAMD12345678",  "SAMD12345678",  "https://www.ncbi.nlm.nih.gov/biosample/SAMD12345678"),
    ])
    check_malformed(BioSample, ["SAM123"])
end

@testset "CA" begin
    check_formats(CA, [
        ("CA123456",    "CA123456",    "CA123456",    "https://reg.clinicalgenome.org/allele/CA123456"),
        ("CA1",         "CA1",         "CA1",         "https://reg.clinicalgenome.org/allele/CA1"),
        ("CA123456789", "CA123456789", "CA123456789", "https://reg.clinicalgenome.org/allele/CA123456789"),
    ])
    check_malformed(CA, [""])
end

@testset "Cellosaurus" begin
    check_formats(Cellosaurus, [
        ("CVCL_0001",   "CVCL_1",       "CVCL_1",       "https://www.cellosaurus.org/CVCL_1"),
        ("CVCL_1234567","CVCL_1234567", "CVCL_1234567", "https://www.cellosaurus.org/CVCL_1234567"),
    ])
    check_malformed(Cellosaurus, ["CVCL"])
end

@testset "ChEBI" begin
    check_formats(ChEBI, [
        ("CHEBI:15377", "15377", "CHEBI:15377", "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:15377"),
        ("CHEBI:1",     "1",     "CHEBI:1",     "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:1"),
        ("15377",       "15377", "CHEBI:15377", "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:15377"),
    ])
    check_malformed(ChEBI, [""])
end

@testset "ChEMBL" begin
    check_formats(ChEMBL, [
        ("CHEMBL25",      "CHEMBL25",      "CHEMBL25",      "https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL25"),
        ("CHEMBL1075272", "CHEMBL1075272", "CHEMBL1075272", "https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL1075272"),
    ])
    check_malformed(ChEMBL, [""])
end

@testset "CL" begin
    check_formats(CL, [
        ("CL:0000001", "0000001", "CL:0000001", "https://purl.obolibrary.org/obo/CL_0000001"),
        ("CL:9999999", "9999999", "CL:9999999", "https://purl.obolibrary.org/obo/CL_9999999"),
        ("0000001",    "0000001", "CL:0000001", "https://purl.obolibrary.org/obo/CL_0000001"),
    ])
    check_malformed(CL, ["CL:123"])
end

@testset "ClinicalTrials" begin
    check_formats(ClinicalTrials, [
        ("NCT01234567", "NCT01234567", "NCT01234567", "https://clinicaltrials.gov/ct2/show/NCT01234567"),
        ("NCT00000001", "NCT00000001", "NCT00000001", "https://clinicaltrials.gov/ct2/show/NCT00000001"),
    ])
    check_malformed(ClinicalTrials, ["NCT123"])
end

@testset "ClinVar" begin
    check_formats(ClinVar, [
        ("VCV000123456",   "VCV000123456",   "VCV000123456",   "https://www.ncbi.nlm.nih.gov/clinvar/VCV000123456"),
        ("RCV00012345",    "RCV00012345",    "RCV00012345",    "https://www.ncbi.nlm.nih.gov/clinvar/RCV00012345"),
        ("SCV00012345",    "SCV00012345",    "SCV00012345",    "https://www.ncbi.nlm.nih.gov/clinvar/SCV00012345"),
    ])
    # Version suffix
    id = parse(ClinVar, "VCV000123456.1")
    @test shortcode(id) == "VCV000123456.1"
    @test parse(ClinVar, shortcode(id)) == id
    check_malformed(ClinVar, ["", "VCV"])
end

@testset "dbGaP" begin
    check_formats(dbGaP, [
        ("phs000001.v1.p1",  "phs000001.v1.p1",  "phs000001.v1.p1",  "https://www.ncbi.nlm.nih.gov/gap/phs000001.v1.p1"),
        ("phs123456.v99.p99","phs123456.v99.p99","phs123456.v99.p99","https://www.ncbi.nlm.nih.gov/gap/phs123456.v99.p99"),
    ])
    # Property access
    id = parse(dbGaP, "phs000001.v1.p1")
    @test id.study == 1
    @test id.ver == 1
    @test id.participant == 1
    check_malformed(dbGaP, [""])
end

@testset "dbSNP" begin
    check_formats(dbSNP, [
        ("rs123456", "rs123456", "rs123456", "https://www.ncbi.nlm.nih.gov/snp/rs123456"),
        ("rs1",      "rs1",      "rs1",      "https://www.ncbi.nlm.nih.gov/snp/rs1"),
    ])
    check_malformed(dbSNP, [""])
end

@testset "dbVar" begin
    check_formats(dbVar, [
        ("nsv123456",  "nsv123456",  "nsv123456",  "https://www.ncbi.nlm.nih.gov/dbvar/nsv123456"),
        ("nssv123456", "nssv123456", "nssv123456", "https://www.ncbi.nlm.nih.gov/dbvar/nssv123456"),
        ("nstd123456", "nstd123456", "nstd123456", "https://www.ncbi.nlm.nih.gov/dbvar/nstd123456"),
    ])
    check_malformed(dbVar, ["ns123"])
end

@testset "DOID" begin
    check_formats(DOID, [
        ("DOID:0001816", "0001816", "DOID:0001816", "https://purl.obolibrary.org/obo/DOID_0001816"),
        ("DOID:9999999", "9999999", "DOID:9999999", "https://purl.obolibrary.org/obo/DOID_9999999"),
        ("0001816",      "0001816", "DOID:0001816", "https://purl.obolibrary.org/obo/DOID_0001816"),
    ])
    check_malformed(DOID, ["DOID:123"])
end

@testset "DrugBank" begin
    check_formats(DrugBank, [
        ("DB00001", "DB00001", "DB00001", "https://go.drugbank.com/drugs/DB00001"),
        ("DB99999", "DB99999", "DB99999", "https://go.drugbank.com/drugs/DB99999"),
    ])
    check_malformed(DrugBank, [""])
end

@testset "ECO" begin
    check_formats(ECO, [
        ("ECO:0000313", "0000313", "ECO:0000313", "https://purl.obolibrary.org/obo/ECO_0000313"),
        ("ECO:0000269", "0000269", "ECO:0000269", "https://purl.obolibrary.org/obo/ECO_0000269"),
        ("0000313",     "0000313", "ECO:0000313", "https://purl.obolibrary.org/obo/ECO_0000313"),
    ])
end

@testset "EFO" begin
    check_formats(EFO, [
        ("EFO:0000001", "0000001", "EFO:0000001", "https://www.ebi.ac.uk/ols/ontologies/efo/terms?short_form=EFO_0000001"),
        ("EFO:9999999", "9999999", "EFO:9999999", "https://www.ebi.ac.uk/ols/ontologies/efo/terms?short_form=EFO_9999999"),
        ("0000001",     "0000001", "EFO:0000001", "https://www.ebi.ac.uk/ols/ontologies/efo/terms?short_form=EFO_0000001"),
    ])
end

@testset "EGA" begin
    check_formats(EGA, [
        ("EGAD00001000001", "EGAD00001000001", "EGAD00001000001", "https://ega-archive.org/EGAD00001000001"),
        ("EGAS00001000002", "EGAS00001000002", "EGAS00001000002", "https://ega-archive.org/EGAS00001000002"),
    ])
    check_malformed(EGA, ["EGA123"])
end

@testset "Ensembl" begin
    @testset "ENSG" begin
        id = parse(ENSG, "ENSG00000141510")
        @test shortcode(id) == "ENSG00000141510"
        @test string(id) == "ENSG00000141510"
        @test purl(id) == "https://www.ensembl.org/id/ENSG00000141510"
        @test parse(ENSG, shortcode(id)) == id
        @test parse(ENSG, purl(id)) == id
    end
    @testset "ENST with version" begin
        id = parse(ENST, "ENST00000269305.9")
        @test shortcode(id) == "ENST00000269305.9"
        @test parse(ENST, shortcode(id)) == id
    end
    @testset "Other types" begin
        for (T, ex) in [(ENSP, "ENSP00000269305"), (ENSE, "ENSE00001234567"),
                         (ENSR, "ENSR00000000001"), (ENSF, "ENSF00000000001")]
            id = parse(T, ex)
            @test parse(T, shortcode(id)) == id
        end
    end
    check_malformed(ENSG, ["", "invalid", "ENST00000269305"])
end

@testset "FlyBase" begin
    check_formats(FlyBase, [
        ("FBgn0031701", "FBgn0031701", "FBgn0031701", "https://flybase.org/reports/FBgn0031701"),
        ("FBgn0000001", "FBgn0000001", "FBgn0000001", "https://flybase.org/reports/FBgn0000001"),
    ])
    check_malformed(FlyBase, [""])
end

@testset "GEO" begin
    check_formats(GEO, [
        ("GSE123456", "GSE123456", "GSE123456", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123456"),
        ("GPL570",    "GPL570",    "GPL570",    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570"),
        ("GSM12345",  "GSM12345",  "GSM12345",  "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM12345"),
        ("GDS12345",  "GDS12345",  "GDS12345",  "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GDS12345"),
    ])
    check_malformed(GEO, ["GEX1"])
end

@testset "GO" begin
    check_formats(GO, [
        ("GO:0006915", "0006915", "GO:0006915", "https://purl.obolibrary.org/obo/GO_0006915"),
        ("GO:0000001", "0000001", "GO:0000001", "https://purl.obolibrary.org/obo/GO_0000001"),
        ("0006915",    "0006915", "GO:0006915", "https://purl.obolibrary.org/obo/GO_0006915"),
    ])
    check_malformed(GO, ["GO:123"])
end

@testset "GWAS" begin
    check_formats(GWAS, [
        ("GCST90000001", "GCST90000001", "GCST90000001", "https://www.ebi.ac.uk/gwas/studies/GCST90000001"),
        ("GCST1",        "GCST1",        "GCST1",        "https://www.ebi.ac.uk/gwas/studies/GCST1"),
    ])
end

@testset "HGNC" begin
    check_formats(HGNC, [
        ("HGNC:5",     "5",     "HGNC:5",     "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:5"),
        ("HGNC:12345", "12345", "HGNC:12345", "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:12345"),
        ("5",          "5",     "HGNC:5",     "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:5"),
    ])
    check_malformed(HGNC, [""])
end

@testset "HMDB" begin
    check_formats(HMDB, [
        ("HMDB0000001", "HMDB0000001", "HMDB0000001", "https://hmdb.ca/metabolites/HMDB0000001"),
        ("HMDB1234567", "HMDB1234567", "HMDB1234567", "https://hmdb.ca/metabolites/HMDB1234567"),
    ])
    check_malformed(HMDB, ["HMDB"])
end

@testset "HPO" begin
    check_formats(HPO, [
        ("HP:0000001", "0000001", "HP:0000001", "https://purl.obolibrary.org/obo/HP_0000001"),
        ("HP:9999999", "9999999", "HP:9999999", "https://purl.obolibrary.org/obo/HP_9999999"),
        ("0000001",    "0000001", "HP:0000001", "https://purl.obolibrary.org/obo/HP_0000001"),
    ])
end

@testset "INSDC" begin
    id = parse(INSDC, "AB123456.1")
    @test shortcode(id) == "AB123456.1"
    @test string(id) == "INSDC:AB123456.1"
    @test purl(id) == "https://www.ncbi.nlm.nih.gov/nuccore/AB123456.1"
    @test parse(INSDC, shortcode(id)) == id
    # Without version
    id2 = parse(INSDC, "U12345")
    @test shortcode(id2) == "U12345"
    @test parse(INSDC, shortcode(id2)) == id2
    check_malformed(INSDC, ["", "invalid"])
end

@testset "IntAct" begin
    check_formats(IntAct, [
        ("EBI-123456", "EBI-123456", "EBI-123456", "https://www.ebi.ac.uk/intact/interaction/EBI-123456"),
        ("EBI-1",      "EBI-1",      "EBI-1",      "https://www.ebi.ac.uk/intact/interaction/EBI-1"),
    ])
end

@testset "InterPro" begin
    check_formats(InterPro, [
        ("IPR011009", "IPR011009", "IPR011009", "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR011009"),
        ("IPR000001", "IPR000001", "IPR000001", "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR000001"),
    ])
    check_malformed(InterPro, [""])
end

@testset "KEGG" begin
    check_formats(KEGG, [
        ("hsa00010", "hsa00010", "hsa00010", "https://www.kegg.jp/entry/hsa00010"),
        ("mmu04210", "mmu04210", "mmu04210", "https://www.kegg.jp/entry/mmu04210"),
    ])
    # Property access
    id = parse(KEGG, "hsa00010")
    @test id.species == "hsa"
    @test id.id == 10
    check_malformed(KEGG, [""])
end

@testset "MeSH" begin
    check_formats(MeSH, [
        ("D001249",    "D001249",    "D001249",    "https://meshb.nlm.nih.gov/record/ui?ui=D001249"),
        ("C000657245", "C657245", "C657245", "https://meshb.nlm.nih.gov/record/ui?ui=C657245"),
    ])
    check_malformed(MeSH, [""])
end

@testset "MetaboLights" begin
    check_formats(MetaboLights, [
        ("MTBLS1",     "MTBLS1",     "MTBLS1",     "https://www.ebi.ac.uk/metabolights/MTBLS1"),
        ("MTBLS12345", "MTBLS12345", "MTBLS12345", "https://www.ebi.ac.uk/metabolights/MTBLS12345"),
    ])
end

@testset "MGI" begin
    check_formats(MGI, [
        ("MGI:96417",   "96417",   "MGI:96417",   "https://www.informatics.jax.org/marker/MGI:96417"),
        ("MGI:1",       "1",       "MGI:1",       "https://www.informatics.jax.org/marker/MGI:1"),
        ("96417",       "96417",   "MGI:96417",   "https://www.informatics.jax.org/marker/MGI:96417"),
    ])
    check_malformed(MGI, [""])
end

@testset "MONDO" begin
    check_formats(MONDO, [
        ("MONDO:0000001", "0000001", "MONDO:0000001", "https://purl.obolibrary.org/obo/MONDO_0000001"),
        ("MONDO:9999999", "9999999", "MONDO:9999999", "https://purl.obolibrary.org/obo/MONDO_9999999"),
        ("0000001",       "0000001", "MONDO:0000001", "https://purl.obolibrary.org/obo/MONDO_0000001"),
    ])
end

@testset "MP" begin
    check_formats(MP, [
        ("MP:0000001", "0000001", "MP:0000001", "https://purl.obolibrary.org/obo/MP_0000001"),
        ("MP:9999999", "9999999", "MP:9999999", "https://purl.obolibrary.org/obo/MP_9999999"),
    ])
end

@testset "NCBIGene" begin
    check_formats(NCBIGene, [
        ("7157",  "7157",  "7157",  "https://www.ncbi.nlm.nih.gov/gene/7157"),
        ("1",     "1",     "1",     "https://www.ncbi.nlm.nih.gov/gene/1"),
    ])
end

@testset "NCBITaxon" begin
    check_formats(NCBITaxon, [
        ("9606",  "9606",  "9606",  "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9606"),
        ("10090", "10090", "10090", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=10090"),
    ])
end

@testset "OMIM" begin
    check_formats(OMIM, [
        ("600918", "600918", "600918", "https://omim.org/entry/600918"),
        ("100100", "100100", "100100", "https://omim.org/entry/100100"),
    ])
    # With allelic variant suffix
    id = parse(OMIM, "613659.0001")
    @test shortcode(id) == "613659.0001"
    @test parse(OMIM, shortcode(id)) == id
    check_malformed(OMIM, ["12345"])
end

@testset "PATO" begin
    check_formats(PATO, [
        ("PATO:0000001", "0000001", "PATO:0000001", "https://purl.obolibrary.org/obo/PATO_0000001"),
        ("PATO:9999999", "9999999", "PATO:9999999", "https://purl.obolibrary.org/obo/PATO_9999999"),
    ])
end

@testset "PDB" begin
    check_formats(PDB, [
        ("1A3N", "1A3N", "1A3N", "https://www.ebi.ac.uk/pdbe/entry/pdb/1A3N"),
        ("7QDH", "7QDH", "7QDH", "https://www.ebi.ac.uk/pdbe/entry/pdb/7QDH"),
    ])
    check_malformed(PDB, ["ABC"])
end

@testset "Pfam" begin
    check_formats(Pfam, [
        ("PF00069", "PF00069", "PF00069", "https://www.ebi.ac.uk/interpro/entry/pfam/PF00069"),
        ("PF00001", "PF00001", "PF00001", "https://www.ebi.ac.uk/interpro/entry/pfam/PF00001"),
    ])
    check_malformed(Pfam, [""])
end

@testset "PubChem" begin
    check_formats(PubChem, [
        ("CID2244",  "CID2244",  "CID2244",  "https://pubchem.ncbi.nlm.nih.gov/compound/2244"),
        ("SID12345", "SID12345", "SID12345", "https://pubchem.ncbi.nlm.nih.gov/substance/12345"),
        ("AID12345", "AID12345", "AID12345", "https://pubchem.ncbi.nlm.nih.gov/assay/12345"),
    ])
end

@testset "PXD" begin
    check_formats(PXD, [
        ("PXD000001",  "PXD000001", "PXD000001", "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD000001"),
        ("PXD123456",  "PXD123456", "PXD123456", "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD123456"),
    ])
    # RPXD prefix normalises to PXD
    @test shortcode(parse(PXD, "RPXD123456")) == "PXD123456"
end

@testset "Reactome" begin
    check_formats(Reactome, [
        ("R-HSA-199420", "R-HSA-199420", "R-HSA-199420", "https://reactome.org/PathwayBrowser/#/R-HSA-199420"),
        ("R-MMU-1",      "R-MMU-1",      "R-MMU-1",      "https://reactome.org/PathwayBrowser/#/R-MMU-1"),
    ])
    @test parse(Reactome, "R-HSA-199420").species == "HSA"
end

@testset "RefSeq" begin
    check_formats(RefSeq, [
        ("NM_000546.6",  "NM_000546.6",  "NM_000546.6",  "https://www.ncbi.nlm.nih.gov/nuccore/NM_000546.6"),
        ("NP_123456.1",  "NP_123456.1",  "NP_123456.1",  "https://www.ncbi.nlm.nih.gov/nuccore/NP_123456.1"),
        ("NC_000001.11", "NC_000001.11", "NC_000001.11", "https://www.ncbi.nlm.nih.gov/nuccore/NC_000001.11"),
    ])
    @test parse(RefSeq, "NM_000546.6").ver == 6
    check_malformed(RefSeq, [""])
end

@testset "RRID" begin
    check_formats(RRID, [
        ("RRID:AB_123456",  "AB_123456",  "RRID:AB_123456",  "https://scicrunch.org/resolver/RRID:AB_123456"),
        ("RRID:SCR_012345", "SCR_012345", "RRID:SCR_012345", "https://scicrunch.org/resolver/RRID:SCR_012345"),
    ])
    check_malformed(RRID, [""])
end

@testset "SGD" begin
    check_formats(SGD, [
        ("S000002493", "S000002493", "S000002493", "https://www.yeastgenome.org/locus/S000002493"),
        ("S000000001", "S000000001", "S000000001", "https://www.yeastgenome.org/locus/S000000001"),
    ])
    check_malformed(SGD, [""])
end

@testset "SO" begin
    check_formats(SO, [
        ("SO:0001234", "0001234", "SO:0001234", "https://purl.obolibrary.org/obo/SO_0001234"),
        ("SO:0000704", "0000704", "SO:0000704", "https://purl.obolibrary.org/obo/SO_0000704"),
        ("0001234",    "0001234", "SO:0001234", "https://purl.obolibrary.org/obo/SO_0001234"),
    ])
end

@testset "SRA" begin
    check_formats(SRA, [
        ("SRR123456", "SRR123456", "SRR123456", "https://www.ncbi.nlm.nih.gov/sra/SRR123456"),
        ("SRX123456", "SRX123456", "SRX123456", "https://www.ncbi.nlm.nih.gov/sra/SRX123456"),
        ("SRS123456", "SRS123456", "SRS123456", "https://www.ncbi.nlm.nih.gov/sra/SRS123456"),
        ("SRP123456", "SRP123456", "SRP123456", "https://www.ncbi.nlm.nih.gov/sra/SRP123456"),
    ])
end

@testset "UBERON" begin
    check_formats(UBERON, [
        ("UBERON:0000001", "0000001", "UBERON:0000001", "https://purl.obolibrary.org/obo/UBERON_0000001"),
        ("UBERON:9999999", "9999999", "UBERON:9999999", "https://purl.obolibrary.org/obo/UBERON_9999999"),
    ])
end

@testset "UniProt" begin
    check_formats(UniProt, [
        ("P04637",     "P04637",     "P04637",     "https://www.uniprot.org/uniprot/P04637"),
        ("Q9Y6K5",     "Q9Y6K5",     "Q9Y6K5",     "https://www.uniprot.org/uniprot/Q9Y6K5"),
        ("A0A0G2JMC8", "A0A0G2JMC8", "A0A0G2JMC8", "https://www.uniprot.org/uniprot/A0A0G2JMC8"),
    ])
    check_malformed(UniProt, ["123456"])
end

@testset "UniRef" begin
    check_formats(UniRef, [
        ("UniRef50_P12345",      "UniRef50_P12345",      "UniRef50_P12345",      "https://www.uniprot.org/uniref/UniRef50_P12345"),
        ("UniRef90_Q9UHC1",      "UniRef90_Q9UHC1",      "UniRef90_Q9UHC1",      "https://www.uniprot.org/uniref/UniRef90_Q9UHC1"),
        ("UniRef100_A0A0G2JMC8", "UniRef100_A0A0G2JMC8", "UniRef100_A0A0G2JMC8", "https://www.uniprot.org/uniref/UniRef100_A0A0G2JMC8"),
    ])
    check_malformed(UniRef, ["UniRef50_invalid"])
end

@testset "WikiPathways" begin
    check_formats(WikiPathways, [
        ("WP554",   "WP554",   "WP554",   "https://www.wikipathways.org/instance/WP554"),
        ("WP12345", "WP12345", "WP12345", "https://www.wikipathways.org/instance/WP12345"),
    ])
    check_malformed(WikiPathways, [""])
end

@testset "WormBase" begin
    check_formats(WormBase, [
        ("WBGene00006763", "WBGene00006763", "WBGene00006763", "https://wormbase.org/species/c_elegans/gene/WBGene00006763"),
        ("WBGene00000001", "WBGene00000001", "WBGene00000001", "https://wormbase.org/species/c_elegans/gene/WBGene00000001"),
    ])
    check_malformed(WormBase, [""])
end
