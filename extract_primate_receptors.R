#!/usr/bin/env Rscript
#
# extract_primate_receptors.R
# ===========================
# Extract all human membrane-bound receptors that are present only in
# primates, i.e. whose orthologues are detected exclusively in primate
# species.
#
# Workflow
# --------
# 1. Query the UniProt REST API for Swiss-Prot (reviewed) human proteins that
#    are annotated as both a receptor (keyword KW-0675) and located at the
#    cell membrane.
# 2. For every receptor that carries an Ensembl gene cross-reference, query
#    the Ensembl REST API to retrieve all orthologues.
# 3. Determine whether each orthologue's species belongs to the primate order
#    (NCBI taxon 9443) using the NCBI Taxonomy API; results are cached in
#    memory to avoid redundant network calls.
# 4. Keep only receptors for which every orthologue is from a primate species
#    (no orthologue found outside Primates).
# 5. Write results to a TSV file and print a summary table to stdout.
#
# Usage
# -----
#     Rscript extract_primate_receptors.R [OUTPUT_FILE]
#
#     OUTPUT_FILE defaults to primate_membrane_receptors.tsv.
#
# Requirements
# ------------
#     install.packages(c("httr", "jsonlite"))

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
})

# ---------------------------------------------------------------------------
# API base URLs
# ---------------------------------------------------------------------------
UNIPROT_SEARCH_URL  <- "https://rest.uniprot.org/uniprotkb/search"
ENSEMBL_REST_URL    <- "https://rest.ensembl.org"
NCBI_TAXONOMY_URL   <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# NCBI taxon ID for order Primates
PRIMATES_TAXON_ID <- 9443L

# Delay (seconds) between successive API calls to respect rate limits
API_CALL_DELAY <- 0.2

# In-memory primate-status cache (environment used for reference semantics)
.is_primate_cache <- new.env(hash = TRUE, parent = emptyenv())

# ---------------------------------------------------------------------------
# Logging helper
# ---------------------------------------------------------------------------
log_info  <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  INFO     ", ...)
log_debug <- function(...) if (isTRUE(getOption("epr.debug"))) message(format(Sys.time(), "%H:%M:%S"), "  DEBUG    ", ...)
log_warn  <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  WARNING  ", ...)

# ---------------------------------------------------------------------------
# Shared HTTP GET with automatic retries
# ---------------------------------------------------------------------------
api_get <- function(url, query = NULL, add_headers = NULL,
                    max_retries = 5, backoff = 1.0) {
  Sys.sleep(API_CALL_DELAY)
  retry_statuses <- c(429L, 500L, 502L, 503L, 504L)
  attempt <- 0L
  repeat {
    resp <- GET(url, query = query, add_headers(.headers = add_headers),
                timeout(30))
    status <- status_code(resp)
    if (!(status %in% retry_statuses) || attempt >= max_retries) {
      if (http_error(resp)) stop_for_status(resp)
      return(resp)
    }
    attempt <- attempt + 1L
    wait <- backoff * (2 ^ (attempt - 1L))
    log_warn(sprintf("HTTP %d – retrying in %.1fs (attempt %d/%d) …",
                     status, wait, attempt, max_retries))
    Sys.sleep(wait)
  }
}

# ---------------------------------------------------------------------------
# Step 1 – fetch human membrane-bound receptors from UniProt
# ---------------------------------------------------------------------------
UNIPROT_QUERY <- paste0(
  "organism_id:9606 ",
  "AND reviewed:true ",
  "AND keyword:KW-0675 ",
  'AND locations:(location:"Cell membrane")'
)

UNIPROT_FIELDS    <- "accession,gene_names,protein_name,xref_ensembl"
UNIPROT_FORMAT    <- "tsv"
UNIPROT_PAGE_SIZE <- 500L

#' Fetch all matching UniProt entries as a data.frame.
fetch_uniprot_receptors <- function() {
  log_info("Querying UniProt for human membrane-bound receptors \u2026")
  params <- list(
    query  = UNIPROT_QUERY,
    format = UNIPROT_FORMAT,
    fields = UNIPROT_FIELDS,
    size   = UNIPROT_PAGE_SIZE
  )

  results <- list()
  url     <- UNIPROT_SEARCH_URL

  repeat {
    resp   <- api_get(url, query = params)
    params <- NULL  # subsequent pages use the full URL from Link header

    body   <- content(resp, as = "text", encoding = "UTF-8")
    lines  <- strsplit(body, "\n", fixed = TRUE)[[1]]
    lines  <- lines[nzchar(trimws(lines))]
    if (length(lines) == 0L) break

    header_cols <- strsplit(lines[[1]], "\t", fixed = TRUE)[[1]]
    for (line in lines[-1]) {
      cols  <- strsplit(line, "\t", fixed = TRUE)[[1]]
      # Pad with empty strings when a row has fewer columns than the header
      if (length(cols) < length(header_cols)) {
        cols <- c(cols, rep("", length(header_cols) - length(cols)))
      }
      entry        <- as.list(cols[seq_along(header_cols)])
      names(entry) <- header_cols
      results      <- c(results, list(entry))
    }

    next_url <- next_link(headers(resp)[["link"]])
    if (is.null(next_url)) break
    url <- next_url
  }

  log_info(sprintf("Found %d human membrane-bound receptor entries in UniProt.", length(results)))
  results
}

#' Extract the 'next' URL from an HTTP Link header, or return NULL.
next_link <- function(link_header) {
  if (is.null(link_header) || !nzchar(link_header)) return(NULL)
  parts <- strsplit(link_header, ",", fixed = TRUE)[[1]]
  for (part in parts) {
    part <- trimws(part)
    if (grepl('rel="next"', part, fixed = TRUE)) {
      url <- sub(";.*$", "", part)
      url <- trimws(url)
      url <- sub("^<", "", url)
      url <- sub(">$", "", url)
      return(url)
    }
  }
  NULL
}

#' Parse Ensembl IDs from a UniProt cross-reference field.
#'
#' UniProt TSV returns the xref_ensembl field as a semicolon-separated list,
#' e.g. "ENST00000351508 [isoform_1];ENST00000456789;".
#' Both transcript IDs (ENST...) and gene IDs (ENSG...) are returned.
parse_ensembl_ids <- function(cross_ref_field) {
  if (is.null(cross_ref_field) || !nzchar(cross_ref_field)) return(character(0))
  tokens <- strsplit(cross_ref_field, ";", fixed = TRUE)[[1]]
  ids <- character(0)
  for (token in tokens) {
    token <- trimws(token)
    if (!nzchar(token)) next
    # The bare ID is the first whitespace-delimited word
    ens_id <- strsplit(token, "\\s+")[[1]][1]
    if (grepl("^ENSG|^ENST", ens_id)) {
      ids <- c(ids, ens_id)
    }
  }
  unique(ids)
}

# ---------------------------------------------------------------------------
# Step 2 – fetch orthologues from Ensembl
# ---------------------------------------------------------------------------

#' Return a list of orthologue records for an Ensembl gene/transcript ID.
fetch_ensembl_orthologues <- function(ensembl_id) {
  url     <- paste0(ENSEMBL_REST_URL, "/homology/id/", ensembl_id)
  params  <- list(type = "orthologues", "content-type" = "application/json")
  headers <- c("Content-Type" = "application/json")

  resp <- tryCatch(
    api_get(url, query = params, add_headers = headers),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("404", msg, fixed = TRUE)) {
        log_debug(sprintf("Ensembl ID %s not found (404).", ensembl_id))
        return(NULL)
      }
      stop(e)
    }
  )

  if (is.null(resp)) return(list())

  body    <- content(resp, as = "text", encoding = "UTF-8")
  data    <- fromJSON(body, simplifyVector = FALSE)
  entries <- data[["data"]]
  if (length(entries) == 0L) return(list())

  entries[[1L]][["homologies"]] %||% list()
}

# Null-coalescing helper
`%||%` <- function(x, y) if (!is.null(x)) x else y

# ---------------------------------------------------------------------------
# Step 3 – determine primate status via NCBI taxonomy (cached)
# ---------------------------------------------------------------------------

#' Fetch the NCBI Taxonomy XML for a taxon ID and return it as a string.
#' Separated out to simplify unit testing.
fetch_tax_xml <- function(taxon_id) {
  resp <- api_get(NCBI_TAXONOMY_URL,
                  query = list(db = "taxonomy", id = taxon_id, retmode = "xml"))
  content(resp, as = "text", encoding = "UTF-8")
}

#' Return TRUE if taxon_id belongs to the order Primates (taxid 9443).
#' Results are cached to minimise redundant network calls.
is_primate_taxon <- function(taxon_id) {
  key <- as.character(taxon_id)
  if (exists(key, envir = .is_primate_cache, inherits = FALSE)) {
    return(get(key, envir = .is_primate_cache, inherits = FALSE))
  }

  result <- tryCatch({
    xml_body <- fetch_tax_xml(taxon_id)
    PRIMATES_TAXON_ID %in% extract_lineage_ids(xml_body)
  }, error = function(e) {
    log_warn(sprintf("Taxonomy lookup failed for taxon %s: %s", taxon_id, conditionMessage(e)))
    FALSE
  })

  assign(key, result, envir = .is_primate_cache)
  result
}

#' Extract all taxon IDs from the lineage in an NCBI Taxonomy XML response.
#' Uses simple pattern matching to avoid an XML library dependency.
extract_lineage_ids <- function(xml_text) {
  if (!nzchar(xml_text)) return(integer(0))
  m <- gregexpr("<TaxId>(\\d+)</TaxId>", xml_text, perl = TRUE)
  matches <- regmatches(xml_text, m)[[1]]
  as.integer(sub("<TaxId>(\\d+)</TaxId>", "\\1", matches, perl = TRUE))
}

# ---------------------------------------------------------------------------
# Step 4 – apply primate-only filter
# ---------------------------------------------------------------------------

#' Return TRUE when every orthologue belongs to a primate species.
#'
#' A protein with no orthologues at all is excluded because we cannot confirm
#' it is conserved in any primate species other than human.
is_primate_only <- function(orthologues) {
  if (length(orthologues) == 0L) return(FALSE)

  for (homology in orthologues) {
    target    <- homology[["target"]] %||% list()
    taxon_id  <- target[["taxon_id"]]
    species   <- target[["species"]] %||% "unknown"

    if (is.null(taxon_id)) {
      log_debug(sprintf("Missing taxon_id for species '%s'; treating as non-primate.", species))
      return(FALSE)
    }

    if (!is_primate_taxon(as.integer(taxon_id))) return(FALSE)
  }

  TRUE
}

# ---------------------------------------------------------------------------
# Step 5 – main pipeline
# ---------------------------------------------------------------------------

#' Assemble a single output row from a UniProt entry and its orthologues.
build_output_row <- function(entry, orthologues) {
  species_vec <- character(0)
  for (h in orthologues) {
    sp <- h[["target"]][["species"]]
    if (!is.null(sp) && nzchar(sp)) species_vec <- c(species_vec, sp)
  }
  primate_species <- paste(sort(unique(species_vec)), collapse = "; ")

  list(
    accession              = entry[["Entry"]] %||% "",
    gene_names             = gsub("\t", " ", entry[["Gene Names"]] %||% ""),
    protein_name           = gsub("\t", " ", entry[["Protein names"]] %||% ""),
    num_primate_orthologues = length(orthologues),
    primate_species        = primate_species
  )
}

OUTPUT_FIELDNAMES <- c(
  "accession", "gene_names", "protein_name",
  "num_primate_orthologues", "primate_species"
)

#' Full extraction pipeline.
run <- function(output_path) {
  # 1. Fetch candidates from UniProt
  uniprot_entries <- fetch_uniprot_receptors()

  primate_only_receptors <- list()
  total <- length(uniprot_entries)

  for (idx in seq_along(uniprot_entries)) {
    entry      <- uniprot_entries[[idx]]
    accession  <- entry[["Entry"]] %||% "?"
    cross_ref  <- entry[["Ensembl transcript"]] %||%
                  entry[["Cross-reference (Ensembl)"]] %||% ""
    ensembl_ids <- parse_ensembl_ids(cross_ref)

    if (length(ensembl_ids) == 0L) {
      log_debug(sprintf("[%d/%d] %s \u2013 no Ensembl ID, skipping.", idx, total, accession))
      next
    }

    log_info(sprintf("[%d/%d] %s \u2013 checking %d Ensembl ID(s) \u2026",
                     idx, total, accession, length(ensembl_ids)))

    # Collect orthologues across all Ensembl IDs for this protein.
    # Deduplicate by target gene ID to avoid double-counting orthologues
    # returned by multiple transcript queries for the same source gene.
    seen_target_ids <- character(0)
    all_orthologues <- list()
    for (eid in ensembl_ids) {
      orthologs <- fetch_ensembl_orthologues(eid)
      for (ortho in orthologs) {
        target_id <- ortho[["target"]][["id"]] %||% ""
        if (!(target_id %in% seen_target_ids)) {
          seen_target_ids <- c(seen_target_ids, target_id)
          all_orthologues <- c(all_orthologues, list(ortho))
        }
      }
    }

    if (is_primate_only(all_orthologues)) {
      row <- build_output_row(entry, all_orthologues)
      primate_only_receptors <- c(primate_only_receptors, list(row))
      log_info(sprintf("  \u2713  %s (%s) is primate-only (%d orthologues).",
                       accession, row[["gene_names"]], row[["num_primate_orthologues"]]))
    }
  }

  # Write TSV
  out_df <- do.call(rbind, lapply(primate_only_receptors, as.data.frame,
                                  stringsAsFactors = FALSE))
  if (is.null(out_df) || nrow(out_df) == 0L) {
    # Write header-only file when no results
    write(paste(OUTPUT_FIELDNAMES, collapse = "\t"), file = output_path)
  } else {
    out_df <- out_df[OUTPUT_FIELDNAMES]
    write.table(out_df, file = output_path, sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = TRUE)
  }

  log_info(sprintf("Done. %d primate-only membrane receptors written to '%s'.",
                   length(primate_only_receptors), output_path))

  # Print summary table to stdout
  sep <- strrep("=", 80)
  cat(sprintf("\n%s\n", sep))
  cat(sprintf("Human membrane-bound receptors present only in primates: %d\n",
              length(primate_only_receptors)))
  cat(sprintf("%s\n", sep))
  if (length(primate_only_receptors) > 0L) {
    cat(sprintf("%-12s %-20s %s\n", "Accession", "Gene", "Protein name"))
    cat(strrep("-", 80), "\n")
    for (row in primate_only_receptors) {
      gene_parts <- strsplit(trimws(row[["gene_names"]]), "\\s+")[[1]]
      gene    <- if (length(gene_parts) > 0L) gene_parts[1] else "\u2014"
      protein <- row[["protein_name"]]
      if (nchar(protein) > 45L) protein <- paste0(substr(protein, 1L, 45L), "\u2026")
      cat(sprintf("%-12s %-20s %s\n", row[["accession"]], gene, protein))
    }
  }
  cat(sprintf("%s\n\n", sep))
  cat(sprintf("Full results saved to: %s\n", output_path))

  invisible(primate_only_receptors)
}

# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------
# sys.nframe() == 0L is TRUE only when this script is invoked directly via
# Rscript (no enclosing call frames), and FALSE when source()'d from another
# script or test file – the R equivalent of Python's if __name__ == "__main__".
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)

  # Parse --debug flag
  debug_flag <- "--debug" %in% args
  args       <- args[args != "--debug"]

  if (length(args) >= 1L && args[1] %in% c("-h", "--help")) {
    cat("Usage: Rscript extract_primate_receptors.R [OUTPUT_FILE] [--debug]\n\n")
    cat("Extract human membrane-bound receptors present only in primates.\n\n")
    cat("Arguments:\n")
    cat("  OUTPUT_FILE  Path for the output TSV file\n")
    cat("               (default: primate_membrane_receptors.tsv)\n")
    cat("  --debug      Enable debug-level logging\n")
    quit(status = 0L)
  }

  output_file <- if (length(args) >= 1L) args[1] else "primate_membrane_receptors.tsv"

  if (debug_flag) options(epr.debug = TRUE)

  run(output_file)
}
