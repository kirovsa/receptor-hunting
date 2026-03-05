# receptor-hunting

A collection of scripts for identifying biologically interesting receptor sets.

---

## `extract_primate_receptors.py`

Extracts all **human membrane-bound receptors** that are **present only in
primates**, i.e. receptors whose orthologues are detected exclusively in
primate species.

### How it works

1. **UniProt query** – fetches Swiss-Prot (reviewed) human proteins annotated
   as a *receptor* (keyword `KW-0675`) with subcellular location *Cell
   membrane* (`organism_id:9606`).
2. **Ensembl orthology** – for each receptor that has an Ensembl cross-
   reference, the [Ensembl REST API][ensembl-rest] homology endpoint is
   queried to retrieve all orthologues.
3. **Primate filter** – for every orthologue species the [NCBI Taxonomy
   API][ncbi-tax] is consulted (results are cached in-process) to check
   whether the taxon falls under order Primates (taxon ID 9443).  A receptor
   is kept only when *every* orthologue belongs to a primate species.
4. **Output** – results are written to a TSV file and a summary table is
   printed to stdout.

### Requirements

```
pip install -r requirements.txt
```

Requires Python 3.10 or later.

### Usage

```bash
# Write results to the default file (primate_membrane_receptors.tsv)
python extract_primate_receptors.py

# Specify a custom output path
python extract_primate_receptors.py results/my_output.tsv

# Enable verbose debug logging
python extract_primate_receptors.py --debug
```

### Output columns

| Column | Description |
|---|---|
| `accession` | UniProt accession number |
| `gene_names` | Gene symbol(s) from UniProt |
| `protein_name` | Protein name from UniProt |
| `num_primate_orthologues` | Number of primate orthologues found in Ensembl |
| `primate_species` | Semicolon-separated list of primate species with an orthologue |

### Running the tests

```bash
python -m unittest test_extract_primate_receptors -v
```

[ensembl-rest]: https://rest.ensembl.org
[ncbi-tax]: https://www.ncbi.nlm.nih.gov/books/NBK25500/
