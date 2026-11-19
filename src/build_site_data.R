#!/usr/bin/env Rscript
# Validate catalogue/*.tsv and compile site/data/catalogue.json + downloadable TSV copies
# for the Quarto + Observable JS website.

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required to build site/data/catalogue.json", call. = FALSE)
}

# Run validation first
val_status <- system2("Rscript", "src/validate_catalogue.R")
if (val_status != 0L) stop("Aborting site data build because catalogue validation failed.", call. = FALSE)

cat_dir <- "catalogue"
site_data_dir <- file.path("site", "data")
dir.create(site_data_dir, showWarnings = FALSE, recursive = TRUE)

read_cat <- function(name) {
  utils::read.delim(file.path(cat_dir, name), stringsAsFactors = FALSE, check.names = FALSE, quote = "")
}

studies  <- read_cat("studies.tsv")
strata   <- read_cat("strata.tsv")
cohorts  <- read_cat("cohort_membership.tsv")
proteins <- read_cat("proteins_index.tsv")
runs     <- read_cat("search_runs.tsv")
queue    <- read_cat("candidate_queue.tsv")
releases <- read_cat("releases.tsv")

# Copy raw TSVs into site/data/ for static downloads
for (f in c("studies.tsv", "strata.tsv", "cohort_membership.tsv", "proteins_index.tsv", "search_runs.tsv", "releases.tsv")) {
  file.copy(file.path(cat_dir, f), file.path(site_data_dir, f), overwrite = TRUE)
}

# Compute pairwise study overlap from cohort_membership.tsv (excluding superseded preprints)
active_study_ids <- studies$study_id[studies$status != "Superseded Preprint"]
active_cohorts   <- cohorts[cohorts$study_id %in% active_study_ids, , drop = FALSE]

overlap_info <- lapply(studies$study_id, function(sid) {
  my_cohorts <- unique(cohorts$cohort[cohorts$study_id == sid])
  if (length(my_cohorts) == 0L) {
    return(list(count = 0L, compact = "None", details = "None", ids = ""))
  }
  shared <- active_cohorts[active_cohorts$cohort %in% my_cohorts & active_cohorts$study_id != sid, , drop = FALSE]
  if (nrow(shared) == 0L) {
    return(list(count = 0L, compact = "None", details = "None", ids = ""))
  }
  by_other <- split(shared$cohort, shared$study_id)
  other_ids <- names(by_other)
  n_ov <- length(other_ids)
  details_vec <- vapply(other_ids, function(oid) {
    st_row <- studies[studies$study_id == oid, , drop = FALSE]
    lbl <- if (nrow(st_row) > 0L) sprintf("%s (%s)", st_row$first_author[1], st_row$year[1]) else oid
    sprintf("%s [%s: %s]", oid, lbl, paste(unique(by_other[[oid]]), collapse = ", "))
  }, character(1L))
  compact_lbl <- if (n_ov == 1L) "1 study" else sprintf("%d studies", n_ov)
  list(
    count = n_ov,
    compact = compact_lbl,
    details = paste(details_vec, collapse = "; "),
    ids = paste(other_ids, collapse = ", ")
  )
})
names(overlap_info) <- studies$study_id

studies$overlap_count   <- vapply(overlap_info, `[[`, integer(1L), "count")
studies$overlap_compact <- vapply(overlap_info, `[[`, character(1L), "compact")
studies$overlap_details <- vapply(overlap_info, `[[`, character(1L), "details")
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

# Write full flattened catalogue export TSV for the one-click "Download Full Table" button
export_cols <- c(
  "stratum_id", "study_id", "first_author", "year", "pmid", "doi", "status",
  "cohorts_display", "overlap_compact", "overlap_ids", "ancestry_detail", "ancestry_group",
  "sex", "n_samples", "sample_matrix", "platform_category", "platform_version",
  "n_assays", "n_unique_proteins", "data_label", "coverage", "cis_window_kb",
  "genome_build", "strand", "imputation_panel", "transformation", "effect_scale",
  "gwas_software", "covariates", "repository", "accession_or_url", "access_level",
  "curation_tier", "gwasprocessing_yaml", "harmonisation_notes"
)
utils::write.table(
  strata_enriched[, export_cols, drop = FALSE],
  file = file.path(site_data_dir, "proteomic_gwas_catalogue_full.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = ""
)

# Enrich protein index with full stratum/study context so Lookup reuses the exact Catalogue table + expandable row
proteins_enriched <- merge(
  proteins,
  strata_enriched[, setdiff(names(strata_enriched), "study_id"), drop = FALSE],
  by = "stratum_id",
  all.x = TRUE,
  sort = FALSE
)

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
    ancestries = paste(unique(active_df$ancestry_group), collapse = ", "),
    any_cis_pqtl = any(active_df$cis_pqtl_reported),
    any_pav_flag = any(active_df$pav_epitope_flag),
    stringsAsFactors = FALSE
  )
}))
protein_summary <- protein_summary[order(protein_summary$gene_symbol), ]
rownames(protein_summary) <- NULL

# Enrich search_runs with explicit excluded & manual-review counts for the Log tab (newest at top)
runs$excluded_records <- pmax(0L, as.integer(runs$deduplicated_new) - as.integer(runs$studies_added))
runs$manual_review_required <- as.integer(runs$candidates_for_review)
runs_sorted <- runs[order(runs$run_date, runs$run_id, decreasing = TRUE), , drop = FALSE]
rownames(runs_sorted) <- NULL

cur_ver  <- if (nrow(releases) > 0 && nzchar(releases$version[1])) releases$version[1] else "unreleased (dev)"
last_upd <- if (nrow(releases) > 0 && nzchar(releases$release_date[1])) releases$release_date[1] else format(Sys.Date(), "%Y-%m-%d")
last_src <- if (nrow(runs_sorted) > 0) runs_sorted$run_date[1] else format(Sys.Date(), "%Y-%m-%d")

payload <- list(
  metadata = list(
    data_mode = "simulated_prototype",
    current_version = cur_ver,
    last_updated = last_upd,
    last_search_date = last_src,
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  ),
  studies = studies,
  strata = strata_enriched,
  cohort_membership = cohorts,
  protein_summary = protein_summary,
  proteins_index = proteins_enriched,
  search_runs = runs_sorted,
  candidate_queue = queue,
  releases = releases
)

json_path <- file.path(site_data_dir, "catalogue.json")
jsonlite::write_json(payload, path = json_path, auto_unbox = TRUE, pretty = TRUE, na = "null")
message("Built ", json_path, " (", round(file.info(json_path)$size / 1024, 1), " KB)")
