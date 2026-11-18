#!/usr/bin/env Rscript
# Validate schema, primary keys, foreign keys, and controlled vocabularies across catalogue/*.tsv.
# Exits with code 1 if any check fails (used in GitHub Actions PR gates).

cat_dir <- "catalogue"
errors <- character(0)

fail <- function(...) {
  errors <<- c(errors, paste0(...))
}

read_cat <- function(name) {
  path <- file.path(cat_dir, name)
  if (!file.exists(path)) {
    fail("Missing required file: ", path)
    return(NULL)
  }
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, quote = "")
}

studies   <- read_cat("studies.tsv")
strata    <- read_cat("strata.tsv")
cohorts   <- read_cat("cohort_membership.tsv")
proteins  <- read_cat("proteins_index.tsv")
runs      <- read_cat("search_runs.tsv")
queue     <- read_cat("candidate_queue.tsv")
releases  <- read_cat("releases.tsv")

check_cols <- function(df, name, req) {
  if (is.null(df)) return()
  missing_cols <- setdiff(req, names(df))
  if (length(missing_cols) > 0) {
    fail(name, " missing columns: ", paste(missing_cols, collapse = ", "))
  }
}

check_cols(studies, "studies.tsv", c(
  "study_id", "first_author", "year", "status", "superseded_by",
  "pmid", "doi", "journal", "title", "repository", "accession_or_url",
  "access_level", "added_release", "verified_full_text"
))

check_cols(strata, "strata.tsv", c(
  "stratum_id", "study_id", "cohorts_display", "ancestry_group", "ancestry_detail",
  "sex", "n_samples", "sample_matrix", "platform_category", "platform_version",
  "n_assays", "n_unique_proteins", "genome_build", "imputation_panel", "coverage",
  "cis_window_kb", "transformation", "effect_scale", "covariates", "gwas_software",
  "file_format", "has_beta_se", "has_eaf", "indel_encoding", "gwasprocessing_yaml",
  "harmonisation_notes"
))

check_cols(cohorts, "cohort_membership.tsv", c("study_id", "cohort", "n_cohort"))
check_cols(proteins, "proteins_index.tsv", c(
  "uniprot", "gene_symbol", "protein_name", "assay_id",
  "stratum_id", "study_id", "cis_pqtl_reported", "pav_epitope_flag"
))
check_cols(runs, "search_runs.tsv", c(
  "run_id", "run_date", "trigger", "raw_records", "deduplicated_new",
  "candidates_for_review", "studies_added", "status"
))
check_cols(releases, "releases.tsv", c("version", "release_date", "n_studies", "n_strata", "summary"))

if (!is.null(studies) && !is.null(strata)) {
  # Uniqueness
  if (anyDuplicated(studies$study_id)) {
    fail("Duplicate study_id in studies.tsv: ", paste(unique(studies$study_id[duplicated(studies$study_id)]), collapse = ", "))
  }
  if (anyDuplicated(strata$stratum_id)) {
    fail("Duplicate stratum_id in strata.tsv: ", paste(unique(strata$stratum_id[duplicated(strata$stratum_id)]), collapse = ", "))
  }

  # Foreign keys
  bad_fk <- setdiff(strata$study_id, studies$study_id)
  if (length(bad_fk) > 0) {
    fail("strata.tsv references unknown study_id: ", paste(bad_fk, collapse = ", "))
  }

  sup <- studies$superseded_by[nzchar(studies$superseded_by)]
  bad_sup <- setdiff(sup, studies$study_id)
  if (length(bad_sup) > 0) {
    fail("studies.tsv superseded_by references unknown study_id: ", paste(bad_sup, collapse = ", "))
  }

  if (!is.null(cohorts)) {
    bad_coh <- setdiff(cohorts$study_id, studies$study_id)
    if (length(bad_coh) > 0) {
      fail("cohort_membership.tsv references unknown study_id: ", paste(bad_coh, collapse = ", "))
    }
  }

  if (!is.null(proteins)) {
    bad_prot_st <- setdiff(proteins$stratum_id, strata$stratum_id)
    if (length(bad_prot_st) > 0) {
      fail("proteins_index.tsv references unknown stratum_id: ", paste(bad_prot_st, collapse = ", "))
    }
  }

  # Controlled vocabularies
  check_vocab <- function(vals, allowed, label) {
    bad <- setdiff(unique(vals), allowed)
    if (length(bad) > 0) {
      fail("Invalid values in ", label, ": ", paste(bad, collapse = ", "),
           " (allowed: ", paste(allowed, collapse = " | "), ")")
    }
  }

  check_vocab(studies$status, c("Published", "Preprint", "Superseded Preprint"), "studies$status")
  check_vocab(studies$access_level, c("Open download", "Portal registration", "Controlled access", "Upon request", "Not available"), "studies$access_level")
  check_vocab(strata$ancestry_group, c("EUR", "AFR", "EAS", "CSA", "AMR", "MID", "Multi-ancestry"), "strata$ancestry_group")
  check_vocab(strata$sex, c("Combined", "Female", "Male"), "strata$sex")
  check_vocab(strata$sample_matrix, c("EDTA plasma", "Citrate plasma", "Serum", "CSF", "Brain tissue", "Urine", "Saliva", "Other tissue"), "strata$sample_matrix")
  check_vocab(strata$platform_category, c("SomaScan", "Olink", "Mass Spectrometry", "NULISA", "Luminex", "Other Immunoassay"), "strata$platform_category")
  check_vocab(strata$genome_build, c("GRCh36", "GRCh37", "GRCh38", "T2T-CHM13", "Unknown"), "strata$genome_build")
  check_vocab(strata$coverage, c("genome_wide", "cis_only", "top_hits_only", "none"), "strata$coverage")
  check_vocab(strata$effect_scale, c("SD_IRNT", "SD_Z", "NPX_log2", "log10_RFU", "log2_intensity", "raw_units", "Z_score_only"), "strata$effect_scale")

  # Sample size sanity (positive integers, explicitly NO minimum threshold)
  if (any(is.na(strata$n_samples) | strata$n_samples < 1L)) {
    fail("strata$n_samples must be >= 1 (no sample-size floor, but must be positive)")
  }

  # Cis-only strata should document cis_window_kb
  cis_missing <- strata$stratum_id[strata$coverage == "cis_only" & is.na(strata$cis_window_kb)]
  if (length(cis_missing) > 0) {
    fail("cis_only strata missing cis_window_kb: ", paste(cis_missing, collapse = ", "))
  }
}

if (length(errors) > 0) {
  message("Catalogue validation FAILED with ", length(errors), " issue(s):")
  for (e in errors) message("  - ", e)
  quit(status = 1L)
} else {
  message(sprintf(
    "Catalogue validation OK: %d studies, %d strata, %d cohort links, %d protein-assay records.",
    nrow(studies), nrow(strata), nrow(cohorts), nrow(proteins)
  ))
}
