#!/usr/bin/env perl
#
# extract_primate_receptors.pl
# ============================
# Extract all **human membrane-bound receptors** that are present **only in
# primates**, i.e. whose orthologues are detected exclusively in primate
# species.
#
# Workflow
# --------
# 1. Query the UniProt REST API for Swiss-Prot (reviewed) human proteins that
#    are annotated as both a *receptor* (keyword KW-0675) and located at the
#    *cell membrane*.
# 2. For every receptor that carries an Ensembl gene/transcript cross-reference,
#    query the Ensembl REST API to retrieve all orthologues.
# 3. Determine whether each orthologue's species belongs to Primates (NCBI taxon
#    9443) using the NCBI Taxonomy API; results are cached in memory.
# 4. Keep only receptors for which *every* orthologue is from a primate species.
# 5. Write results to a TSV file and print a summary table to stdout.
#
# Usage
# -----
#   perl extract_primate_receptors.pl [OUTPUT_FILE] [--debug]
#
#   OUTPUT_FILE defaults to primate_membrane_receptors.tsv
#
# Requirements
# ------------
#   cpanm HTTP::Tiny JSON URI::Escape Getopt::Long Time::HiRes
#

use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use Time::HiRes qw(usleep);
use HTTP::Tiny;
use JSON qw(decode_json);
use URI::Escape qw(uri_escape_utf8);

# -----------------------------
# Logging
# -----------------------------
my $DEBUG = 0;

sub log_info  { my ($msg) = @_; print STDERR _ts() . "  INFO     $msg\n"; }
sub log_warn  { my ($msg) = @_; print STDERR _ts() . "  WARNING  $msg\n"; }
sub log_debug { my ($msg) = @_; return unless $DEBUG; print STDERR _ts() . "  DEBUG    $msg\n"; }

sub _ts {
    my @t = localtime();
    return sprintf("%02d:%02d:%02d", $t[2], $t[1], $t[0]);
}

# -----------------------------
# API base URLs and constants
# -----------------------------
my $UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search";
my $ENSEMBL_REST_URL   = "https://rest.ensembl.org";
my $NCBI_TAXONOMY_URL  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";

my $PRIMATES_TAXON_ID  = 9443;

# Delay between API calls (seconds)
my $API_CALL_DELAY_S   = 0.2;

my $UNIPROT_QUERY = join(" ",
    "organism_id:9606",
    "AND reviewed:true",
    "AND keyword:KW-0675",
    'AND locations:(location:"Cell membrane")'
);

my $UNIPROT_FIELDS    = "accession,gene_names,protein_name,xref_ensembl";
my $UNIPROT_FORMAT    = "tsv";
my $UNIPROT_PAGE_SIZE = 500;

my @OUTPUT_FIELDNAMES = qw(
  accession
  gene_names
  protein_name
  num_primate_orthologues
  primate_species
);

# -----------------------------
# Shared HTTP client
# -----------------------------
my $http = HTTP::Tiny->new(
    timeout  => 30,
    agent    => "extract_primate_receptors.pl/1.0",
    verify_SSL => 1,
);

sub _sleep_rate_limit {
    # HTTP::Tiny doesn't do fractional sleep; use usleep
    usleep(int($API_CALL_DELAY_S * 1_000_000));
}

sub _get {
    my (%args) = @_;
    my $url     = $args{url}     // die "_get requires url";
    my $headers = $args{headers} // {};

    _sleep_rate_limit();

    my $res = $http->get($url, { headers => $headers });

    if (!$res->{success}) {
        my $status = $res->{status} // 0;
        my $reason = $res->{reason} // "";
        die "HTTP GET failed ($status $reason) for $url\n";
    }
    return $res;
}

# Follow the 'next' URL in an HTTP Link header.
# $link_header may be a plain string or an array-ref (HTTP::Tiny joins
# repeated headers with ", " for most fields, but being defensive here
# protects against any version that returns an array-ref instead).
sub _next_link {
    my ($link_header) = @_;
    return undef unless defined($link_header);

    # Normalise to a flat list of comma-separated parts.
    my @parts;
    if (ref($link_header) eq 'ARRAY') {
        @parts = map { split(/\s*,\s*/, $_) } @$link_header;
    } else {
        return undef if $link_header eq "";
        @parts = split(/\s*,\s*/, $link_header);
    }

    for my $part (@parts) {
        if ($part =~ /<([^>]+)>\s*;\s*rel="next"/) {
            return $1;
        }
    }
    return undef;
}

# -----------------------------
# Step 1 – fetch UniProt receptors
# -----------------------------
sub fetch_uniprot_receptors {
    log_info("Querying UniProt for human membrane-bound receptors ...");

    my @results;
    my @col_headers;   # captured from the first page; reused for all subsequent pages
    my %seen_urls;     # guard against infinite-loop if the same next URL repeats
    my $page_num = 0;

    # First request uses query params. UniProt returns a full next-link for pagination.
    my $url = $UNIPROT_SEARCH_URL . "?"
        . "query="  . uri_escape_utf8($UNIPROT_QUERY)
        . "&format=" . uri_escape_utf8($UNIPROT_FORMAT)
        . "&fields=" . uri_escape_utf8($UNIPROT_FIELDS)
        . "&size="   . uri_escape_utf8($UNIPROT_PAGE_SIZE);

    while (defined $url) {
        if ($seen_urls{$url}++) {
            log_warn("Pagination loop detected – same URL seen twice; stopping.");
            last;
        }

        $page_num++;
        log_debug("Fetching UniProt page $page_num: $url");

        my $res = _get(url => $url);

        my $text = $res->{content} // "";
        $text =~ s/\s+\z//;

        my @lines = split(/\r?\n/, $text);
        last if !@lines;

        # UniProt repeats the TSV header row on every page.  Capture it from
        # the first page; on subsequent pages simply discard it.
        my $header_line = shift @lines;
        if (!@col_headers) {
            @col_headers = split(/\t/, $header_line);
        }

        for my $line (@lines) {
            next if $line !~ /\S/;
            my @cols = split(/\t/, $line, -1);

            # pad to headers length
            push @cols, ("") x (@col_headers - @cols) if @cols < @col_headers;

            my %row;
            @row{@col_headers} = @cols;
            push @results, \%row;
        }

        # HTTP::Tiny normalises all header names to lowercase, so only check
        # 'link'.  It may be a plain string or an array-ref – _next_link
        # handles both forms.
        $url = _next_link($res->{headers}{link});

        if (defined $url) {
            log_debug("Next page link found; continuing pagination.");
        } else {
            log_debug("No next page link; pagination complete after $page_num page(s).");
        }
    }

    log_info("Found " . scalar(@results) . " human membrane-bound receptor entries in UniProt.");
    return \@results;
}

sub parse_ensembl_ids {
    my ($cross_ref_field) = @_;
    return [] if !defined($cross_ref_field) || $cross_ref_field eq "";

    my @ids;
    for my $token (split(/;/, $cross_ref_field)) {
        $token =~ s/^\s+|\s+$//g;
        next if $token eq "";

        my ($ens_id) = split(/\s+/, $token);
        next if !defined($ens_id);

        if ($ens_id =~ /^(ENSG|ENST)/) {
            push @ids, $ens_id;
        }
    }

    # dedupe preserve order
    my %seen;
    my @deduped = grep { !$seen{$_}++ } @ids;
    return \@deduped;
}

# -----------------------------
# Step 2 – fetch orthologues from Ensembl
# -----------------------------
sub fetch_ensembl_orthologues {
    my ($ensembl_id) = @_;
    my $url = $ENSEMBL_REST_URL . "/homology/id/" . uri_escape_utf8($ensembl_id)
            . "?type=orthologues&content-type=application/json";

    my $res;
    eval {
        $res = _get(url => $url, headers => { "Content-Type" => "application/json" });
        1;
    } or do {
        my $err = $@ || "unknown error";
        # Ensembl returns 404 for missing IDs; HTTP::Tiny turns it into a die above.
        # We can't inspect status easily after die, so do a lightweight re-GET without dying:
        my $tmp = HTTP::Tiny->new(timeout => 30, agent => "extract_primate_receptors.pl/1.0", verify_SSL => 1);
        _sleep_rate_limit();
        my $r2 = $tmp->get($url, { headers => { "Content-Type" => "application/json" } });
        if (!$r2->{success} && ($r2->{status} // 0) == 404) {
            log_debug("Ensembl ID $ensembl_id not found (404).");
            return [];
        }
        die $err;
    };

    my $data = decode_json($res->{content} // "{}");
    my $entries = $data->{data} // [];
    return [] if ref($entries) ne 'ARRAY' || !@$entries;

    my $homologies = $entries->[0]{homologies} // [];
    return (ref($homologies) eq 'ARRAY') ? $homologies : [];
}

# -----------------------------
# Step 3 – NCBI taxonomy primate check (cached)
# -----------------------------
my %IS_PRIMATE_CACHE;

sub _extract_lineage_ids {
    my ($xml_text) = @_;
    my @ids = ($xml_text =~ m/<TaxId>(\d+)<\/TaxId>/g);
    return [ map { int($_) } @ids ];
}

sub is_primate_taxon {
    my ($taxon_id) = @_;
    $taxon_id = int($taxon_id);

    if (exists $IS_PRIMATE_CACHE{$taxon_id}) {
        return $IS_PRIMATE_CACHE{$taxon_id} ? 1 : 0;
    }

    my $url = $NCBI_TAXONOMY_URL
        . "?db=taxonomy&id=" . uri_escape_utf8($taxon_id)
        . "&retmode=xml";

    my $result = 0;
    eval {
        my $res = _get(url => $url);
        my $ids = _extract_lineage_ids($res->{content} // "");
        for my $id (@$ids) {
            if ($id == $PRIMATES_TAXON_ID) {
                $result = 1;
                last;
            }
        }
        1;
    } or do {
        my $err = $@ || "unknown error";
        chomp $err;
        log_warn("Taxonomy lookup failed for taxon $taxon_id: $err");
        $result = 0;
    };

    $IS_PRIMATE_CACHE{$taxon_id} = $result ? 1 : 0;
    return $result;
}

# -----------------------------
# Step 4 – apply primate-only filter
# -----------------------------
sub is_primate_only {
    my ($orthologues) = @_;
    return 0 if !defined($orthologues) || ref($orthologues) ne 'ARRAY' || !@$orthologues;

    for my $homology (@$orthologues) {
        my $target  = (ref($homology) eq 'HASH') ? ($homology->{target} // {}) : {};
        my $taxon_id = $target->{taxon_id};
        my $species  = $target->{species} // "unknown";

        if (!defined $taxon_id) {
            log_debug("Missing taxon_id for species '$species'; treating as non-primate.");
            return 0;
        }
        if (!is_primate_taxon($taxon_id)) {
            return 0;
        }
    }
    return 1;
}

# -----------------------------
# output formatting
# -----------------------------
sub build_output_row {
    my ($entry, $orthologues) = @_;

    my %species_set;
    for my $h (@$orthologues) {
        next if ref($h) ne 'HASH';
        my $sp = $h->{target}{species} // "";
        next if $sp eq "";
        $species_set{$sp} = 1;
    }
    my @primate_species = sort keys %species_set;

    my $gene_names    = ($entry->{"Gene Names"}    // "");
    my $protein_names = ($entry->{"Protein names"} // "");
    $gene_names    =~ s/\t/ /g;
    $protein_names =~ s/\t/ /g;

    return {
        accession               => ($entry->{"Entry"} // ""),
        gene_names              => $gene_names,
        protein_name            => $protein_names,
        num_primate_orthologues => scalar(@$orthologues),
        primate_species         => join("; ", @primate_species),
    };
}

sub write_tsv {
    my ($path, $rows) = @_;
    open(my $fh, ">:encoding(UTF-8)", $path) or die "Cannot open '$path' for writing: $!\n";

    print $fh join("\t", @OUTPUT_FIELDNAMES) . "\n";
    for my $row (@$rows) {
        my @vals = map { defined($row->{$_}) ? $row->{$_} : "" } @OUTPUT_FIELDNAMES;
        # remove newlines just in case
        for (@vals) { s/\R/ /g; }
        print $fh join("\t", @vals) . "\n";
    }
    close $fh;
}

sub print_summary {
    my ($rows, $output_path) = @_;

    print "\n" . ("=" x 80) . "\n";
    print "Human membrane-bound receptors present only in primates: " . scalar(@$rows) . "\n";
    print ("=" x 80) . "\n";

    if (@$rows) {
        printf "% -12s %-20s %s\n", "Accession", "Gene", "Protein name";
        print "-" x 80 . "\n";

        for my $row (@$rows) {
            my $gene = "—";
            if (defined $row->{gene_names} && $row->{gene_names} =~ /(\S+)/) {
                $gene = $1;
            }

            my $protein = $row->{protein_name} // "";
            $protein = substr($protein, 0, 45) . "..." if length($protein) > 45;

            printf "% -12s %-20s %s\n",
                ($row->{accession} // ""),
                $gene,
                $protein;
        }
    }

    print ("=" x 80) . "\n\n";
    print "Full results saved to: $output_path\n";
}

# -----------------------------
# Main pipeline
# -----------------------------
sub run_pipeline {
    my ($output_path) = @_;

    my $uniprot_entries = fetch_uniprot_receptors();
    my $total = scalar(@$uniprot_entries);

    my @primate_only_receptors;

    for (my $i = 0; $i < $total; $i++) {
        my $idx = $i + 1;
        my $entry = $uniprot_entries->[$i];

        my $accession = $entry->{"Entry"} // "?";

        my $cross_ref = $entry->{"Ensembl transcript"} // $entry->{"Cross-reference (Ensembl)"} // "";
        my $ensembl_ids = parse_ensembl_ids($cross_ref);

        if (!@$ensembl_ids) {
            log_debug("[$idx/$total] $accession - no Ensembl ID, skipping.");
            next;
        }

        log_info("[$idx/$total] $accession - checking " . scalar(@$ensembl_ids) . " Ensembl ID(s) ...");

        my %seen_target_ids;
        my @all_orthologues;

        for my $eid (@$ensembl_ids) {
            my $orthologs = fetch_ensembl_orthologues($eid);
            for my $ortho (@$orthologs) {
                next if ref($ortho) ne 'HASH';
                my $target_id = $ortho->{target}{id} // "";
                next if $target_id eq "";
                next if $seen_target_ids{$target_id}++;
                push @all_orthologues, $ortho;
            }
        }

        if (is_primate_only(\@all_orthologues)) {
            my $row = build_output_row($entry, \@all_orthologues);
            push @primate_only_receptors, $row;

            log_info("  ✓  $accession ($row->{gene_names}) is primate-only ($row->{num_primate_orthologues} orthologues).\n");
        }
    }

    write_tsv($output_path, \@primate_only_receptors);

    log_info("Done. " . scalar(@primate_only_receptors) . " primate-only membrane receptors written to '$output_path'.");

    print_summary(\@primate_only_receptors, $output_path);
}

# -----------------------------
# CLI
# -----------------------------
my $output = "primate_membrane_receptors.tsv";

GetOptions(
    "debug" => \$DEBUG,
) or die "Usage: perl extract_primate_receptors.pl [OUTPUT_FILE] [--debug]\n";

# Remaining arg (optional) is output file
if (@ARGV) {
    $output = shift @ARGV;
}

run_pipeline($output);