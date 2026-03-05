"""
Unit tests for extract_primate_receptors.py

These tests exercise all pure (non-network) helper functions and mock the
network-dependent parts so no real HTTP requests are made.
"""

import importlib
import sys
import os
import types
import unittest
from unittest.mock import MagicMock, patch

# ---------------------------------------------------------------------------
# Import module under test
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import extract_primate_receptors as epr


class TestNextLink(unittest.TestCase):
    """_next_link helper."""

    def test_returns_none_for_empty_string(self):
        self.assertIsNone(epr._next_link(""))

    def test_returns_none_when_no_next(self):
        header = '<https://example.com/page1>; rel="prev"'
        self.assertIsNone(epr._next_link(header))

    def test_returns_url_for_next(self):
        header = '<https://rest.uniprot.org/uniprotkb/search?cursor=abc>; rel="next"'
        self.assertEqual(
            epr._next_link(header),
            "https://rest.uniprot.org/uniprotkb/search?cursor=abc",
        )

    def test_multiple_rels_picks_next(self):
        header = (
            '<https://example.com/prev>; rel="prev", '
            '<https://example.com/next>; rel="next"'
        )
        self.assertEqual(epr._next_link(header), "https://example.com/next")


class TestParseEnsemblIds(unittest.TestCase):
    """parse_ensembl_ids helper."""

    def test_empty_string_returns_empty_list(self):
        self.assertEqual(epr.parse_ensembl_ids(""), [])

    def test_none_equivalent_empty(self):
        self.assertEqual(epr.parse_ensembl_ids(None), [])  # type: ignore[arg-type]

    def test_single_enst_id(self):
        self.assertEqual(
            epr.parse_ensembl_ids("ENST00000351508;"),
            ["ENST00000351508"],
        )

    def test_enst_with_isoform_annotation(self):
        self.assertEqual(
            epr.parse_ensembl_ids("ENST00000351508 [isoform_1];"),
            ["ENST00000351508"],
        )

    def test_multiple_ids(self):
        ids = epr.parse_ensembl_ids("ENST00000111111 [isoform_1];ENST00000222222;")
        self.assertEqual(ids, ["ENST00000111111", "ENST00000222222"])

    def test_ensg_ids_are_also_accepted(self):
        ids = epr.parse_ensembl_ids("ENSG00000012345;ENSG00000067890 [v2];")
        self.assertEqual(ids, ["ENSG00000012345", "ENSG00000067890"])

    def test_mixed_ensg_and_enst(self):
        ids = epr.parse_ensembl_ids("ENST00000111111;ENSG00000012345;")
        self.assertIn("ENST00000111111", ids)
        self.assertIn("ENSG00000012345", ids)

    def test_deduplication(self):
        ids = epr.parse_ensembl_ids("ENST00000111111;ENST00000111111;")
        self.assertEqual(ids, ["ENST00000111111"])

    def test_non_ensembl_tokens_ignored(self):
        ids = epr.parse_ensembl_ids("SOMEOTHER00000;ENST00000111111;")
        self.assertEqual(ids, ["ENST00000111111"])


class TestExtractLineageIds(unittest.TestCase):
    """_extract_lineage_ids helper."""

    SAMPLE_XML = """<?xml version="1.0" ?>
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
    </TaxaSet>"""

    def test_extracts_all_taxon_ids(self):
        ids = epr._extract_lineage_ids(self.SAMPLE_XML)
        self.assertIn(9606, ids)
        self.assertIn(9443, ids)  # Primates
        self.assertIn(2759, ids)  # Eukaryota

    def test_empty_xml_returns_empty_list(self):
        self.assertEqual(epr._extract_lineage_ids(""), [])


class TestIsPrimateOnlyFilter(unittest.TestCase):
    """is_primate_only filter logic."""

    def _make_homology(self, taxon_id: int, species: str) -> dict:
        return {"target": {"taxon_id": taxon_id, "species": species, "id": f"ENSG_FAKE_{taxon_id}"}}

    def test_empty_orthologues_returns_false(self):
        self.assertFalse(epr.is_primate_only([]))

    def test_all_primates_returns_true(self):
        orthos = [
            self._make_homology(9598, "pan_troglodytes"),
            self._make_homology(9544, "macaca_mulatta"),
        ]
        with patch.object(epr, "is_primate_taxon", return_value=True):
            self.assertTrue(epr.is_primate_only(orthos))

    def test_non_primate_present_returns_false(self):
        orthos = [
            self._make_homology(9598, "pan_troglodytes"),
            self._make_homology(10090, "mus_musculus"),  # mouse
        ]

        def mock_is_primate(taxon_id: int) -> bool:
            return taxon_id != 10090

        with patch.object(epr, "is_primate_taxon", side_effect=mock_is_primate):
            self.assertFalse(epr.is_primate_only(orthos))

    def test_missing_taxon_id_returns_false(self):
        orthos = [{"target": {"species": "unknown_species"}}]
        self.assertFalse(epr.is_primate_only(orthos))


class TestBuildOutputRow(unittest.TestCase):
    """build_output_row assembler."""

    def _entry(self):
        return {
            "Entry": "P12345",
            "Gene Names": "RECEPTOR1 REC1",
            "Protein names": "Receptor protein 1",
        }

    def _orthologues(self):
        return [
            {"target": {"species": "pan_troglodytes", "id": "ENSPTRG00000001"}},
            {"target": {"species": "macaca_mulatta", "id": "ENSMMUG00000001"}},
        ]

    def test_accession_captured(self):
        row = epr.build_output_row(self._entry(), self._orthologues())
        self.assertEqual(row["accession"], "P12345")

    def test_gene_names_captured(self):
        row = epr.build_output_row(self._entry(), self._orthologues())
        self.assertIn("RECEPTOR1", row["gene_names"])

    def test_orthologue_count(self):
        row = epr.build_output_row(self._entry(), self._orthologues())
        self.assertEqual(row["num_primate_orthologues"], 2)

    def test_species_sorted_and_joined(self):
        row = epr.build_output_row(self._entry(), self._orthologues())
        species = row["primate_species"].split("; ")
        self.assertEqual(species, sorted(species))

    def test_empty_orthologues(self):
        row = epr.build_output_row(self._entry(), [])
        self.assertEqual(row["num_primate_orthologues"], 0)
        self.assertEqual(row["primate_species"], "")


class TestIsPrimateTaxon(unittest.TestCase):
    """is_primate_taxon with mocked HTTP."""

    def setUp(self):
        epr._is_primate_cache.clear()

    def _xml_with_primates_lineage(self, taxon_id: int) -> str:
        return f"""<TaxaSet><Taxon>
            <TaxId>{taxon_id}</TaxId>
            <LineageEx>
                <Taxon><TaxId>131567</TaxId></Taxon>
                <Taxon><TaxId>9443</TaxId></Taxon>
            </LineageEx>
        </Taxon></TaxaSet>"""

    def _xml_without_primates_lineage(self, taxon_id: int) -> str:
        return f"""<TaxaSet><Taxon>
            <TaxId>{taxon_id}</TaxId>
            <LineageEx>
                <Taxon><TaxId>131567</TaxId></Taxon>
                <Taxon><TaxId>7711</TaxId></Taxon>
            </LineageEx>
        </Taxon></TaxaSet>"""

    def test_primate_taxon_returns_true(self):
        mock_resp = MagicMock()
        mock_resp.text = self._xml_with_primates_lineage(9598)
        with patch.object(epr, "_get", return_value=mock_resp):
            self.assertTrue(epr.is_primate_taxon(9598))

    def test_non_primate_taxon_returns_false(self):
        mock_resp = MagicMock()
        mock_resp.text = self._xml_without_primates_lineage(10090)
        with patch.object(epr, "_get", return_value=mock_resp):
            self.assertFalse(epr.is_primate_taxon(10090))

    def test_result_is_cached(self):
        mock_resp = MagicMock()
        mock_resp.text = self._xml_with_primates_lineage(9598)
        with patch.object(epr, "_get", return_value=mock_resp) as mock_get:
            epr.is_primate_taxon(9598)
            epr.is_primate_taxon(9598)
            # Should only call _get once (second call hits the cache)
            self.assertEqual(mock_get.call_count, 1)


if __name__ == "__main__":
    unittest.main()
