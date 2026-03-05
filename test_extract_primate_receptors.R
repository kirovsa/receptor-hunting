# Unit tests for extract_primate_receptors.R
#
# These tests exercise all pure (non-network) helper functions and mock the
# network-dependent parts so no real HTTP requests are made.
#
# Run with:
#   Rscript -e "testthat::test_file('test_extract_primate_receptors.R')"
# or from within an R session:
#   testthat::test_file("test_extract_primate_receptors.R")

library(testthat)
library(mockery)

# Source the module under test.
# sys.nframe() == 0 is FALSE when sourced, so the CLI block is skipped.
source("extract_primate_receptors.R", local = FALSE)

# ---------------------------------------------------------------------------
# next_link
# ---------------------------------------------------------------------------
test_that("next_link returns NULL for empty string", {
  expect_null(next_link(""))
})

test_that("next_link returns NULL when no 'next' rel is present", {
  header <- '<https://example.com/page1>; rel="prev"'
  expect_null(next_link(header))
})

test_that("next_link returns URL when 'next' rel is present", {
  header <- '<https://rest.uniprot.org/uniprotkb/search?cursor=abc>; rel="next"'
  expect_equal(next_link(header),
               "https://rest.uniprot.org/uniprotkb/search?cursor=abc")
})

test_that("next_link picks 'next' from multiple rels", {
  header <- paste0(
    '<https://example.com/prev>; rel="prev", ',
    '<https://example.com/next>; rel="next"'
  )
  expect_equal(next_link(header), "https://example.com/next")
})

test_that("next_link returns NULL for NULL input", {
  expect_null(next_link(NULL))
})

# ---------------------------------------------------------------------------
# parse_ensembl_ids
# ---------------------------------------------------------------------------
test_that("parse_ensembl_ids returns empty vector for empty string", {
  expect_equal(parse_ensembl_ids(""), character(0))
})

test_that("parse_ensembl_ids returns empty vector for NULL", {
  expect_equal(parse_ensembl_ids(NULL), character(0))
})

test_that("parse_ensembl_ids parses a single ENST ID", {
  expect_equal(parse_ensembl_ids("ENST00000351508;"),
               "ENST00000351508")
})

test_that("parse_ensembl_ids strips isoform annotation", {
  expect_equal(parse_ensembl_ids("ENST00000351508 [isoform_1];"),
               "ENST00000351508")
})

test_that("parse_ensembl_ids parses multiple IDs", {
  ids <- parse_ensembl_ids("ENST00000111111 [isoform_1];ENST00000222222;")
  expect_equal(ids, c("ENST00000111111", "ENST00000222222"))
})

test_that("parse_ensembl_ids accepts ENSG IDs", {
  ids <- parse_ensembl_ids("ENSG00000012345;ENSG00000067890 [v2];")
  expect_equal(ids, c("ENSG00000012345", "ENSG00000067890"))
})

test_that("parse_ensembl_ids accepts mixed ENSG and ENST IDs", {
  ids <- parse_ensembl_ids("ENST00000111111;ENSG00000012345;")
  expect_true("ENST00000111111" %in% ids)
  expect_true("ENSG00000012345" %in% ids)
})

test_that("parse_ensembl_ids deduplicates IDs", {
  ids <- parse_ensembl_ids("ENST00000111111;ENST00000111111;")
  expect_equal(ids, "ENST00000111111")
})

test_that("parse_ensembl_ids ignores non-Ensembl tokens", {
  ids <- parse_ensembl_ids("SOMEOTHER00000;ENST00000111111;")
  expect_equal(ids, "ENST00000111111")
})

# ---------------------------------------------------------------------------
# extract_lineage_ids
# ---------------------------------------------------------------------------
sample_xml <- '<?xml version="1.0" ?>
<TaxaSet>
  <Taxon>
    <TaxId>9606</TaxId>
    <LineageEx>
      <Taxon><TaxId>131567</TaxId></Taxon>
      <Taxon><TaxId>2759</TaxId></Taxon>
      <Taxon><TaxId>9443</TaxId></Taxon>
      <Taxon><TaxId>9526</TaxId></Taxon>
    </LineageEx>
  </Taxon>
</TaxaSet>'

test_that("extract_lineage_ids extracts all taxon IDs", {
  ids <- extract_lineage_ids(sample_xml)
  expect_true(9606L %in% ids)
  expect_true(9443L %in% ids)  # Primates
  expect_true(2759L %in% ids)  # Eukaryota
})

test_that("extract_lineage_ids returns integer(0) for empty string", {
  expect_equal(extract_lineage_ids(""), integer(0))
})

# ---------------------------------------------------------------------------
# is_primate_only
# ---------------------------------------------------------------------------
make_homology <- function(taxon_id, species) {
  list(target = list(
    taxon_id = taxon_id,
    species  = species,
    id       = paste0("ENSG_FAKE_", taxon_id)
  ))
}

test_that("is_primate_only returns FALSE for empty orthologues", {
  expect_false(is_primate_only(list()))
})

test_that("is_primate_only returns TRUE when all orthologues are primates", {
  orthos <- list(
    make_homology(9598L, "pan_troglodytes"),
    make_homology(9544L, "macaca_mulatta")
  )
  stub(is_primate_only, "is_primate_taxon", TRUE)
  expect_true(is_primate_only(orthos))
})

test_that("is_primate_only returns FALSE when a non-primate orthologue is present", {
  orthos <- list(
    make_homology(9598L,  "pan_troglodytes"),
    make_homology(10090L, "mus_musculus")
  )
  stub(is_primate_only, "is_primate_taxon", function(taxon_id) taxon_id != 10090L)
  expect_false(is_primate_only(orthos))
})

test_that("is_primate_only returns FALSE when taxon_id is missing", {
  orthos <- list(list(target = list(species = "unknown_species")))
  expect_false(is_primate_only(orthos))
})

# ---------------------------------------------------------------------------
# build_output_row
# ---------------------------------------------------------------------------
make_entry <- function() {
  list(
    "Entry"          = "P12345",
    "Gene Names"     = "RECEPTOR1 REC1",
    "Protein names"  = "Receptor protein 1"
  )
}

make_orthologues <- function() {
  list(
    list(target = list(species = "pan_troglodytes", id = "ENSPTRG00000001")),
    list(target = list(species = "macaca_mulatta",  id = "ENSMMUG00000001"))
  )
}

test_that("build_output_row captures accession", {
  row <- build_output_row(make_entry(), make_orthologues())
  expect_equal(row[["accession"]], "P12345")
})

test_that("build_output_row captures gene_names", {
  row <- build_output_row(make_entry(), make_orthologues())
  expect_true(grepl("RECEPTOR1", row[["gene_names"]]))
})

test_that("build_output_row counts orthologues correctly", {
  row <- build_output_row(make_entry(), make_orthologues())
  expect_equal(row[["num_primate_orthologues"]], 2L)
})

test_that("build_output_row sorts primate_species alphabetically", {
  row     <- build_output_row(make_entry(), make_orthologues())
  species <- strsplit(row[["primate_species"]], "; ", fixed = TRUE)[[1]]
  expect_equal(species, sort(species))
})

test_that("build_output_row handles empty orthologues", {
  row <- build_output_row(make_entry(), list())
  expect_equal(row[["num_primate_orthologues"]], 0L)
  expect_equal(row[["primate_species"]], "")
})

# ---------------------------------------------------------------------------
# is_primate_taxon (with mocked fetch_tax_xml)
# ---------------------------------------------------------------------------
xml_with_primates <- function(taxon_id) {
  sprintf('<TaxaSet><Taxon>
    <TaxId>%d</TaxId>
    <LineageEx>
      <Taxon><TaxId>131567</TaxId></Taxon>
      <Taxon><TaxId>9443</TaxId></Taxon>
    </LineageEx>
  </Taxon></TaxaSet>', taxon_id)
}

xml_without_primates <- function(taxon_id) {
  sprintf('<TaxaSet><Taxon>
    <TaxId>%d</TaxId>
    <LineageEx>
      <Taxon><TaxId>131567</TaxId></Taxon>
      <Taxon><TaxId>7711</TaxId></Taxon>
    </LineageEx>
  </Taxon></TaxaSet>', taxon_id)
}

# Helper: clear the in-memory primate-status cache between tests
clear_primate_cache <- function() {
  rm(list = ls(.is_primate_cache), envir = .is_primate_cache)
}

test_that("is_primate_taxon returns TRUE for a primate taxon", {
  clear_primate_cache()
  stub(is_primate_taxon, "fetch_tax_xml", xml_with_primates(9598L))
  expect_true(is_primate_taxon(9598L))
})

test_that("is_primate_taxon returns FALSE for a non-primate taxon", {
  clear_primate_cache()
  stub(is_primate_taxon, "fetch_tax_xml", xml_without_primates(10090L))
  expect_false(is_primate_taxon(10090L))
})

test_that("is_primate_taxon caches results (fetch_tax_xml called only once)", {
  clear_primate_cache()
  mock_fetch <- mock(xml_with_primates(9598L))
  stub(is_primate_taxon, "fetch_tax_xml", mock_fetch)
  is_primate_taxon(9598L)
  is_primate_taxon(9598L)
  expect_called(mock_fetch, 1L)
})
