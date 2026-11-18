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

# Copy raw TSVs into site/data/ so static GitHub Pages visitors can download them directly
for (f in c("studies.tsv", "strata.tsv", "cohort_membership.tsv", "proteins_index.tsv", "search_runs.tsv", "releases.tsv")) {
  file.copy(file.path(cat_dir, f), file.path(site_data_dir, f), overwrite = TRUE)
}

# Enrich strata with study-level fields so the main explorer has one flat fast array
strata_enriched <- merge(
  strata,
  studies,
  by = "study_id",
  all.x = TRUE,
  sort = FALSE
)
# Preserve original strata order (chronological / study order)
strata_enriched <- strata_enriched[match(strata$stratum_id, strata_enriched$stratum_id), ]

# Compute MR & Coloc readiness tier per stratum
# - Full Ready: genome_wide or cis_only + has_beta_se + has_eaf + Open/Portal access
# - Cis Ready: cis_only with beta/se + eaf
# - Caution (Missing EAF / Indel codes): missing EAF or D/I indels
# - Limited: top_hits_only or Controlled/Upon request
strata_enriched$mr_coloc_readiness <- ifelse(
  strata_enriched$coverage == "top_hits_only" | strata_enriched$access_level %in% c("Controlled access", "Upon request"),
  "Limited (Access / Top-hits)",
  ifelse(
    !strata_enriched$has_eaf,
    "Caution (No signed EAF)",
    ifelse(
      strata_enriched$coverage == "genome_wide",
      "Full (Cis + Trans Ready)",
      "Cis-Only Ready"
    )
  )
)

# Attach cohort list per study onto studies
cohort_map <- split(cohorts, cohorts$study_id)
studies$cohorts_list <- lapply(studies$study_id, function(sid) {
  sub <- cohort_map[[sid]]
  if (is.null(sub)) list() else sub[, c("cohort", "n_cohort")]
})

# Attach max sample size, platforms, ancestries, and matrices summary to each study
strata_by_study <- split(strata, strata$study_id)
studies$max_n <- vapply(studies$study_id, function(sid) max(strata_by_study[[sid]]$n_samples), integer(1))
studies$n_strata <- vapply(studies$study_id, function(sid) nrow(strata_by_study[[sid]]), integer(1))
studies$platforms_summary <- vapply(studies$study_id, function(sid) {
  paste(unique(strata_by_study[[sid]]$platform_category), collapse = ", ")
}, character(1))
studies$ancestries_summary <- vapply(studies$study_id, function(sid) {
  paste(unique(strata_by_study[[sid]]$ancestry_group), collapse = ", ")
}, character(1))
studies$matrices_summary <- vapply(studies$study_id, function(sid) {
  paste(unique(strata_by_study[[sid]]$sample_matrix), collapse = ", ")
}, character(1))

# Enrich protein index with compact stratum/study context for fast client-side rendering
proteins_enriched <- merge(
  proteins,
  strata_enriched[, c(
    "stratum_id", "first_author", "year", "status", "superseded_by",
    "ancestry_group", "sex", "n_samples", "sample_matrix",
    "platform_category", "platform_version", "genome_build", "coverage",
    "access_level", "accession_or_url", "mr_coloc_readiness", "gwasprocessing_yaml"
  )],
  by = "stratum_id",
  all.x = TRUE,
  sort = FALSE
)

# Protein directory summary (1 row per unique protein)
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

cur_ver  <- if (nrow(releases) > 0 && nzchar(releases$version[1])) releases$version[1] else "unreleased (dev)"
last_upd <- if (nrow(releases) > 0 && nzchar(releases$release_date[1])) releases$release_date[1] else format(Sys.Date(), "%Y-%m-%d")
last_src <- if (nrow(runs) > 0) runs$run_date[nrow(runs)] else format(Sys.Date(), "%Y-%m-%d")

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
  search_runs = runs,
  candidate_queue = queue,
  releases = releases
)

json_path <- file.path(site_data_dir, "catalogue.json")
jsonlite::write_json(payload, path = json_path, auto_unbox = TRUE, pretty = TRUE, na = "null")
message("Built ", json_path, " (", round(file.info(json_path)$size / 1024, 1), " KB)")
