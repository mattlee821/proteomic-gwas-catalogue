#!/usr/bin/env Rscript
# Validate catalogue/*.tsv and compile site/data/catalogue.json + downloadable TSV copies
# for the Quarto + Observable JS website.

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required to build site/data/catalogue.json", call. = FALSE)
}

val_status <- system2("Rscript", "src/validate_catalogue.R")
if (val_status != 0L) stop("Aborting site data build because catalogue validation failed.", call. = FALSE)

cat_dir <- "catalogue"
site_data_dir <- file.path("site", "data")
dir.create(site_data_dir, showWarnings = FALSE, recursive = TRUE)

read_cat <- function(name) {
  utils::read.delim(file.path(cat_dir, name), stringsAsFactors = FALSE, check.names = FALSE, quote = "")
}

copy_if_changed <- function(src, dst) {
  if (!file.exists(dst) || unname(tools::md5sum(src)) != unname(tools::md5sum(dst))) {
    file.copy(src, dst, overwrite = TRUE)
  }
}

studies  <- read_cat("studies.tsv")
strata   <- read_cat("strata.tsv")
cohorts  <- read_cat("cohort_membership.tsv")
proteins <- read_cat("proteins_index.tsv")
runs     <- read_cat("search_runs.tsv")
queue    <- read_cat("candidate_queue.tsv")
releases <- read_cat("releases.tsv")

for (f in c("studies.tsv", "strata.tsv", "cohort_membership.tsv", "proteins_index.tsv", "search_runs.tsv", "releases.tsv")) {
  copy_if_changed(file.path(cat_dir, f), file.path(site_data_dir, f))
}

# Sync docs/methods.md -> site/_methods_content.md (if local docs/methods.md exists) for site/search-strategy.qmd
methods_src <- file.path("docs", "methods.md")
methods_dst <- file.path("site", "_methods_content.md")
if (file.exists(methods_src)) {
  m_lines <- readLines(methods_src, warn = FALSE, encoding = "UTF-8")
  if (length(m_lines) > 0L && grepl("^#\\s+", m_lines[1L])) {
    m_lines <- m_lines[-1L]
  }
  tmp_m <- tempfile(fileext = ".md")
  writeLines(m_lines, tmp_m, useBytes = TRUE)
  copy_if_changed(tmp_m, methods_dst)
  unlink(tmp_m)
}

# Compute pairwise study overlap from cohort_membership.tsv (excluding superseded preprints)
active_study_ids <- studies$study_id[studies$status != "Superseded Preprint"]
active_cohorts   <- cohorts[cohorts$study_id %in% active_study_ids, , drop = FALSE]

overlap_info <- lapply(studies$study_id, function(sid) {
  my_cohorts <- unique(cohorts$cohort[cohorts$study_id == sid])
  if (length(my_cohorts) == 0L) {
    return(list(count = 0L, compact = "None", ids = "None"))
  }
  shared <- active_cohorts[active_cohorts$cohort %in% my_cohorts & active_cohorts$study_id != sid, , drop = FALSE]
  if (nrow(shared) == 0L) {
    return(list(count = 0L, compact = "None", ids = "None"))
  }
  other_ids <- unique(shared$study_id)
  n_ov <- length(other_ids)
  compact_lbl <- if (n_ov == 1L) "1 study" else sprintf("%d studies", n_ov)
  list(
    count = n_ov,
    compact = compact_lbl,
    ids = paste(other_ids, collapse = ", ")
  )
})
names(overlap_info) <- studies$study_id

studies$overlap_count   <- vapply(overlap_info, `[[`, integer(1L), "count")
studies$overlap_compact <- vapply(overlap_info, `[[`, character(1L), "compact")
studies$overlap_ids     <- vapply(overlap_info, `[[`, character(1L), "ids")

# Enrich strata with study-level fields
strata_enriched <- merge(
  strata,
  studies,
  by = "study_id",
  all.x = TRUE,
  sort = FALSE
)
strata_enriched <- strata_enriched[match(strata$stratum_id, strata_enriched$stratum_id), ]

# Super-population abbreviation for POPULATION column (EUR, AFR, EAS, CSA, AMR, MID, MULTI)
strata_enriched$pop_super <- ifelse(
  strata_enriched$ancestry_group == "Multi-ancestry",
  "MULTI",
  strata_enriched$ancestry_group
)

# 4-tier DATA availability label: Genome-wide | Gene-region | Top-hits | Unavailable
strata_enriched$data_label <- ifelse(
  strata_enriched$coverage == "none" | strata_enriched$access_level == "Not available",
  "Unavailable",
  ifelse(
    strata_enriched$coverage == "genome_wide",
    "Genome-wide",
    ifelse(
      strata_enriched$coverage == "cis_only",
      "Gene-region",
      "Top-hits"
    )
  )
)

# 3-tier ACCESS label for expanded OTHER column: open | request | FALSE
strata_enriched$access_short <- ifelse(
  strata_enriched$coverage == "none" | strata_enriched$access_level == "Not available",
  "FALSE",
  ifelse(
    strata_enriched$access_level == "Open download",
    "open",
    "request"
  )
)

# Write full flattened catalogue export TSV for "Download Full Table"
export_cols <- c(
  "stratum_id", "study_id", "first_author", "year", "pmid", "doi", "status",
  "cohorts_display", "overlap_compact", "overlap_ids", "pop_super", "ancestry_detail",
  "sex", "n_samples", "sample_matrix", "platform_category", "platform_version",
  "n_assays", "n_unique_proteins", "data_label", "coverage", "cis_window_kb",
  "genome_build", "strand", "imputation_panel", "transformation", "effect_scale",
  "gwas_software", "covariates", "repository", "accession_or_url", "access_short",
  "curation_tier"
)
tmp_full <- tempfile(fileext = ".tsv")
utils::write.table(
  strata_enriched[, export_cols, drop = FALSE],
  file = tmp_full,
  sep = "\t", row.names = FALSE, quote = FALSE, na = ""
)
copy_if_changed(tmp_full, file.path(site_data_dir, "proteomic_gwas_catalogue_full.tsv"))
unlink(tmp_full)

# Enrich protein index with full stratum/study context for Lookup
proteins_enriched <- merge(
  proteins,
  strata_enriched[, setdiff(names(strata_enriched), "study_id"), drop = FALSE],
  by = "stratum_id",
  all.x = TRUE,
  sort = FALSE
)

if (nrow(proteins_enriched) > 0L) {
  prot_split <- split(proteins_enriched, proteins_enriched$uniprot)
  protein_summary <- do.call(rbind, lapply(prot_split, function(df) {
    active_df <- df[df$status != "Superseded Preprint", , drop = FALSE]
    data.frame(
      uniprot = df$uniprot[1],
      gene_symbol = df$gene_symbol[1],
      protein_name = df$protein_name[1],
      n_studies = length(unique(active_df$study_id)),
      n_strata = length(unique(active_df$stratum_id)),
      n_assays_total = nrow(active_df),
      n_genome_wide = sum(active_df$data_label == "Genome-wide"),
      max_n = if (nrow(active_df) > 0) max(active_df$n_samples) else 0L,
      platforms = paste(unique(active_df$platform_category), collapse = ", "),
      matrices = paste(unique(active_df$sample_matrix), collapse = ", "),
      ancestries = paste(unique(active_df$pop_super), collapse = ", "),
      any_cis_pqtl = any(active_df$cis_pqtl_reported),
      any_pav_flag = any(active_df$pav_epitope_flag),
      stringsAsFactors = FALSE
    )
  }))
  protein_summary <- protein_summary[order(protein_summary$gene_symbol), ]
  rownames(protein_summary) <- NULL
} else {
  protein_summary <- data.frame(
    uniprot = character(0), gene_symbol = character(0), protein_name = character(0),
    n_studies = integer(0), n_strata = integer(0), n_assays_total = integer(0),
    n_genome_wide = integer(0), max_n = integer(0), platforms = character(0),
    matrices = character(0), ancestries = character(0), any_cis_pqtl = logical(0),
    any_pav_flag = logical(0), stringsAsFactors = FALSE
  )
}

# Ensure optional search_runs detail columns exist even on minimal 13-col automated search logs
if (!("date_from" %in% names(runs))) runs$date_from <- runs$run_date
if (!("date_to" %in% names(runs)))   runs$date_to   <- runs$run_date
if (!("studies_added_ids" %in% names(runs))) runs$studies_added_ids <- ""
if (!("preprints_replaced_details" %in% names(runs))) runs$preprints_replaced_details <- ""

runs$excluded_records <- pmax(0L, as.integer(runs$deduplicated_new) - as.integer(runs$studies_added))
runs$manual_review_required <- as.integer(runs$candidates_for_review)
runs_sorted <- runs[order(runs$run_date, runs$run_id, decreasing = TRUE), , drop = FALSE]
rownames(runs_sorted) <- NULL

# Explicit catalogue-wide Preprint -> Published replacement records
sup_rows <- studies[studies$status == "Superseded Preprint" & !is.na(studies$superseded_by) & nzchar(studies$superseded_by), , drop = FALSE]
preprints_replaced_list <- if (nrow(sup_rows) > 0L) {
  do.call(rbind, lapply(seq_len(nrow(sup_rows)), function(i) {
    pr <- sup_rows[i, ]
    pb <- studies[studies$study_id == pr$superseded_by, , drop = FALSE]
    data.frame(
      preprint_study_id = pr$study_id,
      preprint_citation = sprintf("%s (%s)", pr$first_author, pr$year),
      preprint_doi = pr$doi,
      published_study_id = pr$superseded_by,
      published_citation = if (nrow(pb) > 0L) sprintf("%s (%s, %s)", pb$first_author[1], pb$year[1], pb$journal[1]) else pr$superseded_by,
      published_pmid = if (nrow(pb) > 0L) pb$pmid[1] else "",
      published_doi = if (nrow(pb) > 0L) pb$doi[1] else "",
      stringsAsFactors = FALSE
    )
  }))
} else {
  data.frame()
}

latest_run_tag <- if (nrow(runs_sorted) > 0 && "release_tag" %in% names(runs_sorted) && any(!is.na(runs_sorted$release_tag) & nzchar(runs_sorted$release_tag))) {
  runs_sorted$release_tag[!is.na(runs_sorted$release_tag) & nzchar(runs_sorted$release_tag)][1]
} else {
  "v0.1.0-pre"
}
cur_ver  <- if (nrow(releases) > 0 && !is.na(releases$version[1]) && nzchar(releases$version[1])) releases$version[1] else latest_run_tag
last_upd <- if (nrow(releases) > 0 && !is.na(releases$release_date[1]) && nzchar(releases$release_date[1])) releases$release_date[1] else if (nrow(runs_sorted) > 0) runs_sorted$run_date[1] else format(Sys.Date(), "%Y-%m-%d")
min_from <- if (nrow(runs) > 0 && any(!is.na(runs$date_from) & nzchar(runs$date_from))) min(runs$date_from[!is.na(runs$date_from) & nzchar(runs$date_from)]) else "—"
max_to   <- if (nrow(runs) > 0 && any(!is.na(runs$date_to) & nzchar(runs$date_to)))   max(runs$date_to[!is.na(runs$date_to) & nzchar(runs$date_to)])     else "—"
tot_screened <- sum(as.integer(runs$deduplicated_new), na.rm = TRUE)
tot_included <- sum(studies$status != "Superseded Preprint")

payload <- list(
  metadata = list(
    data_mode = "production",
    current_version = cur_ver,
    last_updated = last_upd,
    date_from = min_from,
    date_to = max_to,
    total_screened = tot_screened,
    total_included = tot_included
  ),
  studies = studies,
  strata = strata_enriched,
  cohort_membership = cohorts,
  protein_summary = protein_summary,
  proteins_index = proteins_enriched,
  search_runs = runs_sorted,
  preprints_replaced = preprints_replaced_list,
  candidate_queue = queue,
  releases = releases
)

json_path <- file.path(site_data_dir, "catalogue.json")
tmp_json  <- tempfile(fileext = ".json")
jsonlite::write_json(payload, path = tmp_json, auto_unbox = TRUE, pretty = TRUE, na = "null")
copy_if_changed(tmp_json, json_path)
unlink(tmp_json)
message("Built ", json_path, " (", round(file.info(json_path)$size / 1024, 1), " KB)")
