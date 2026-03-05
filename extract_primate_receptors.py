#!/usr/bin/env python3
"""
extract_primate_receptors.py
============================
Extract all **human membrane-bound receptors** that are present **only in
primates**, i.e. whose orthologues are detected exclusively in primate
species.

Workflow
--------
1. Query the UniProt REST API for Swiss-Prot (reviewed) human proteins that
   are annotated as both a *receptor* (keyword KW-0675) and located at the
   *cell membrane*.
2. For every receptor that carries an Ensembl gene cross-reference, query the
   Ensembl REST API to retrieve all orthologues.
3. Determine whether each orthologue's species belongs to the primate order
   (NCBI taxon 9443) using the NCBI Taxonomy API; results are cached in
   memory to avoid redundant network calls.
4. Keep only receptors for which *every* orthologue is from a primate species
   (no orthologue found outside Primates).
5. Write results to a TSV file and print a summary table to stdout.

Usage
-----
    python extract_primate_receptors.py [OUTPUT_FILE]

    OUTPUT_FILE defaults to ``primate_membrane_receptors.tsv``.

Requirements
------------
    pip install requests
"""

import csv
import sys
import time
import argparse
import logging
from typing import Optional

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# API base URLs
# ---------------------------------------------------------------------------
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
ENSEMBL_REST_URL = "https://rest.ensembl.org"
NCBI_TAXONOMY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# NCBI taxon ID for order Primates
PRIMATES_TAXON_ID = 9443

# Delay (seconds) between successive API calls to respect rate limits
API_CALL_DELAY = 0.2

# ---------------------------------------------------------------------------
# Shared HTTP session with automatic retries
# ---------------------------------------------------------------------------
_session = requests.Session()
_retry = Retry(
    total=5,
    backoff_factor=1.0,
    status_forcelist=[429, 500, 502, 503, 504],
    allowed_methods=["GET"],
)
_session.mount("https://", HTTPAdapter(max_retries=_retry))


def _get(url: str, params: Optional[dict] = None, headers: Optional[dict] = None) -> requests.Response:
    """Perform a GET request with rate-limit delay and raise on HTTP error."""
    time.sleep(API_CALL_DELAY)
    resp = _session.get(url, params=params, headers=headers, timeout=30)
    resp.raise_for_status()
    return resp


# ---------------------------------------------------------------------------
# Step 1 – fetch human membrane-bound receptors from UniProt
# ---------------------------------------------------------------------------
UNIPROT_QUERY = (
    "organism_id:9606 "
    "AND reviewed:true "
    "AND keyword:KW-0675 "
    'AND locations:(location:"Cell membrane")'
)

UNIPROT_FIELDS = "accession,gene_names,protein_name,xref_ensembl"
UNIPROT_FORMAT = "tsv"
UNIPROT_PAGE_SIZE = 500


def fetch_uniprot_receptors() -> list[dict]:
    """Return a list of dicts for every matching UniProt entry."""
    log.info("Querying UniProt for human membrane-bound receptors …")
    params = {
        "query": UNIPROT_QUERY,
        "format": UNIPROT_FORMAT,
        "fields": UNIPROT_FIELDS,
        "size": UNIPROT_PAGE_SIZE,
    }

    results: list[dict] = []
    url: Optional[str] = UNIPROT_SEARCH_URL

    while url:
        resp = _get(url, params=params)
        params = None  # only pass params on first request; subsequent pages use the full URL from Link header

        lines = resp.text.strip().splitlines()
        if not lines:
            break

        headers_row = lines[0].split("\t")
        for line in lines[1:]:
            if line.strip():
                cols = line.split("\t")
                # Pad with empty strings if a row has fewer columns than the header
                entry = dict(zip(headers_row, cols + [""] * (len(headers_row) - len(cols))))
                results.append(entry)

        # Follow pagination via the HTTP Link header
        url = _next_link(resp.headers.get("Link", ""))

    log.info("Found %d human membrane-bound receptor entries in UniProt.", len(results))
    return results


def _next_link(link_header: str) -> Optional[str]:
    """Extract the 'next' URL from an HTTP Link header, or return None."""
    if not link_header:
        return None
    for part in link_header.split(","):
        part = part.strip()
        if 'rel="next"' in part:
            return part.split(";")[0].strip().strip("<>")
    return None


def parse_ensembl_ids(cross_ref_field: str) -> list[str]:
    """
    Parse Ensembl IDs from a UniProt cross-reference field.

    UniProt TSV returns the ``xref_ensembl`` field (column header
    *"Ensembl transcript"*) as a semicolon-separated list of entries, e.g.::

        ENST00000351508 [isoform_1];ENST00000456789;

    The Ensembl REST API homology endpoint accepts both transcript IDs
    (ENST…) and gene IDs (ENSG…), so we extract whichever prefixes appear.
    """
    ids: list[str] = []
    if not cross_ref_field:
        return ids
    for token in cross_ref_field.split(";"):
        token = token.strip()
        if not token:
            continue
        # The bare ID is the first whitespace-delimited word in each token
        ens_id = token.split(" ")[0].strip()
        if ens_id.startswith(("ENSG", "ENST")):
            ids.append(ens_id)
    return list(dict.fromkeys(ids))  # deduplicate while preserving order


# ---------------------------------------------------------------------------
# Step 2 – fetch orthologues from Ensembl
# ---------------------------------------------------------------------------

def fetch_ensembl_orthologues(ensembl_id: str) -> list[dict]:
    """
    Return a list of orthologue records for an Ensembl gene ID.

    Each record is the dict returned by the Ensembl REST API homology
    endpoint, containing at minimum ``target.species`` and
    ``target.taxon_id``.
    """
    url = f"{ENSEMBL_REST_URL}/homology/id/{ensembl_id}"
    params = {"type": "orthologues", "content-type": "application/json"}
    headers = {"Content-Type": "application/json"}

    try:
        resp = _get(url, params=params, headers=headers)
    except requests.HTTPError as exc:
        if exc.response is not None and exc.response.status_code == 404:
            log.debug("Ensembl ID %s not found (404).", ensembl_id)
            return []
        raise

    data = resp.json()
    entries = data.get("data", [])
    if not entries:
        return []

    homologies = entries[0].get("homologies", [])
    return homologies


# ---------------------------------------------------------------------------
# Step 3 – determine primate status via NCBI taxonomy (cached)
# ---------------------------------------------------------------------------
_is_primate_cache: dict[int, bool] = {}


def is_primate_taxon(taxon_id: int) -> bool:
    """
    Return True if *taxon_id* belongs to the order Primates (taxid 9443).

    Results are cached to minimise redundant network calls.
    """
    if taxon_id in _is_primate_cache:
        return _is_primate_cache[taxon_id]

    params = {
        "db": "taxonomy",
        "id": taxon_id,
        "retmode": "xml",
    }
    try:
        resp = _get(NCBI_TAXONOMY_URL, params=params)
        result = PRIMATES_TAXON_ID in _extract_lineage_ids(resp.text)
    except Exception as exc:  # noqa: BLE001
        log.warning("Taxonomy lookup failed for taxon %d: %s", taxon_id, exc)
        result = False

    _is_primate_cache[taxon_id] = result
    return result


def _extract_lineage_ids(xml_text: str) -> list[int]:
    """
    Extract all taxon IDs from the lineage in an NCBI Taxonomy XML response.

    We avoid an XML library dependency by using simple string parsing.
    """
    import re

    return [int(m) for m in re.findall(r"<TaxId>(\d+)</TaxId>", xml_text)]


# ---------------------------------------------------------------------------
# Step 4 – apply primate-only filter
# ---------------------------------------------------------------------------

def is_primate_only(orthologues: list[dict]) -> bool:
    """
    Return True when *every* orthologue belongs to a primate species.

    A protein with **no orthologues at all** (i.e. no homology data in
    Ensembl) is excluded because we cannot confirm it is conserved in any
    primate species other than human.
    """
    if not orthologues:
        return False

    for homology in orthologues:
        target = homology.get("target", {})
        taxon_id = target.get("taxon_id")
        species = target.get("species", "unknown")

        if taxon_id is None:
            log.debug("Missing taxon_id for species '%s'; treating as non-primate.", species)
            return False

        if not is_primate_taxon(int(taxon_id)):
            return False

    return True


# ---------------------------------------------------------------------------
# Step 5 – main pipeline
# ---------------------------------------------------------------------------

def build_output_row(entry: dict, orthologues: list[dict]) -> dict:
    """Assemble a single output row from a UniProt entry and its orthologues."""
    primate_species = sorted({
        h.get("target", {}).get("species", "")
        for h in orthologues
        if h.get("target", {}).get("species")
    })
    return {
        "accession": entry.get("Entry", ""),
        "gene_names": entry.get("Gene Names", "").replace("\t", " "),
        "protein_name": entry.get("Protein names", "").replace("\t", " "),
        "num_primate_orthologues": len(orthologues),
        "primate_species": "; ".join(primate_species),
    }


OUTPUT_FIELDNAMES = [
    "accession",
    "gene_names",
    "protein_name",
    "num_primate_orthologues",
    "primate_species",
]


def run(output_path: str) -> None:
    """Full extraction pipeline."""
    # 1. Fetch candidates from UniProt
    uniprot_entries = fetch_uniprot_receptors()

    primate_only_receptors: list[dict] = []
    total = len(uniprot_entries)

    for idx, entry in enumerate(uniprot_entries, start=1):
        accession = entry.get("Entry", "?")
        cross_ref = entry.get("Ensembl transcript", "") or entry.get("Cross-reference (Ensembl)", "")
        ensembl_ids = parse_ensembl_ids(cross_ref)

        if not ensembl_ids:
            log.debug("[%d/%d] %s – no Ensembl ID, skipping.", idx, total, accession)
            continue

        log.info("[%d/%d] %s – checking %d Ensembl ID(s) …", idx, total, accession, len(ensembl_ids))

        # Collect orthologues across all Ensembl IDs for this protein.
        # Deduplicate by target gene ID to avoid double-counting orthologues
        # returned by multiple transcript queries for the same source gene.
        seen_target_ids: set[str] = set()
        all_orthologues: list[dict] = []
        for eid in ensembl_ids:
            orthologs = fetch_ensembl_orthologues(eid)
            for ortho in orthologs:
                target_id = ortho.get("target", {}).get("id", "")
                if target_id not in seen_target_ids:
                    seen_target_ids.add(target_id)
                    all_orthologues.append(ortho)

        if is_primate_only(all_orthologues):
            row = build_output_row(entry, all_orthologues)
            primate_only_receptors.append(row)
            log.info("  ✓  %s (%s) is primate-only (%d orthologues).",
                     accession, row["gene_names"], row["num_primate_orthologues"])

    # Write TSV
    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUTPUT_FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerows(primate_only_receptors)

    log.info("Done. %d primate-only membrane receptors written to '%s'.",
             len(primate_only_receptors), output_path)

    # Print summary table to stdout
    print(f"\n{'='*80}")
    print(f"Human membrane-bound receptors present only in primates: {len(primate_only_receptors)}")
    print(f"{'='*80}")
    if primate_only_receptors:
        print(f"{'Accession':<12} {'Gene':<20} {'Protein name'}")
        print("-" * 80)
        for row in primate_only_receptors:
            gene = row["gene_names"].split()[0] if row["gene_names"].split() else "—"
            protein = row["protein_name"][:45] + "…" if len(row["protein_name"]) > 45 else row["protein_name"]
            print(f"{row['accession']:<12} {gene:<20} {protein}")
    print(f"{'='*80}\n")
    print(f"Full results saved to: {output_path}")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract human membrane-bound receptors that are present only in "
            "primates, based on UniProt annotations and Ensembl orthology data."
        )
    )
    parser.add_argument(
        "output",
        nargs="?",
        default="primate_membrane_receptors.tsv",
        help="Path for the output TSV file (default: primate_membrane_receptors.tsv)",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug-level logging.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    run(args.output)
