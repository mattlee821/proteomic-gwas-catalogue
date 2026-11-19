# Tier 4: Strictly Atomic Field Provenance, Conflict/Attrition Reconciliation,
# Preprint Supersede Patching, and NEEDS_MANUAL_REVIEW Escalation ([C-2], [C-5], [M-1], [M-6], [M-7]).

make_prov_row <- function(study_id, target_table, record_id, field, required_tier,
                          value, status, source_tier, location, evidence_snippet, conflict_notes = "") {
  data.frame(
    study_id = clean_str(study_id),
    target_table = clean_str(target_table),
    record_id = clean_str(record_id),
    field = clean_str(field),
    required_tier = clean_str(required_tier),
    value = clean_str(value),
    status = clean_str(status),
    source_tier = clean_str(source_tier),
    location = substr(clean_str(location), 1L, 180L),
    evidence_snippet = substr(clean_str(evidence_snippet), 1L, 240L),
    conflict_notes = substr(clean_str(conflict_notes), 1L, 240L),
    stringsAsFactors = FALSE
  )
}

pivot_provenance_to_table <- function(prov_df, tbl_name, col_order) {
  sub_df <- prov_df[prov_df$target_table == tbl_name, , drop = FALSE]
  if (nrow(sub_df) == 0L) {
    empty_l <- stats::setNames(replicate(length(col_order), character(0), simplify = FALSE), col_order)
    return(as.data.frame(empty_l, stringsAsFactors = FALSE, check.names = FALSE))
  }
  rec_ids <- unique(sub_df$record_id)
  rows <- lapply(rec_ids, function(rid) {
    rsub <- sub_df[sub_df$record_id == rid, , drop = FALSE]
    vals <- vapply(col_order, function(cn) {
      idx <- which(rsub$field == cn)
      if (length(idx) > 0L) rsub$value[[idx[[1L]]]] else ""
    }, character(1L))
    as.data.frame(as.list(vals), stringsAsFactors = FALSE, check.names = FALSE)
  })
  do.call(rbind, rows)
}

run_tier4_provenance_and_escalate <- function(t1, t2, t3_supp, t3_repo, out_dir,
                                              existing_studies = NULL,
                                              manifests_dir = NULL,
                                              added_release = "v0.3.0-staging") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  sid <- t1$study_id
  prov <- list()
  add_p <- function(...) prov[[length(prov) + 1L]] <<- make_prov_row(...)

  # --- 1. Resolve Preprint Lineage & `superseded_by` Foreign Key ([M-7]) ---
  superseded_by_val <- ""
  superseded_by_status <- "API_EXTRACTED"
  superseded_by_note <- ""
  preprint_patch_df <- NULL

  if (identical(t1$status_label, "Published")) {
    # Published rows ALWAYS have superseded_by = ""!
    # If this published study supersedes a tracked preprint, emit a patch for the *preprint* row.
    if (nzchar(t1$supersedes_preprint_doi) || nzchar(t1$supersedes_epmc_id)) {
      superseded_by_note <- sprintf(
        "Published study supersedes preprint (%s %s); patch emitted for preprint row.",
        t1$supersedes_preprint_doi, t1$supersedes_epmc_id
      )
      matched_pre_id <- ""
      if (!is.null(existing_studies) && nrow(existing_studies) > 0L && nzchar(t1$supersedes_preprint_doi)) {
        idx <- which(tolower(trimws(existing_studies$doi)) == t1$supersedes_preprint_doi)
        if (length(idx) > 0L) matched_pre_id <- existing_studies$study_id[[idx[[1L]]]]
      }
      preprint_patch_df <- data.frame(
        action = "UPDATE_SUPERSEDED_PREPRINT",
        preprint_study_id = ifelse(nzchar(matched_pre_id), matched_pre_id, paste0("UNTRACKED:", t1$supersedes_preprint_doi)),
        preprint_doi = t1$supersedes_preprint_doi,
        preprint_epmc_id = t1$supersedes_epmc_id,
        new_status = "Superseded Preprint",
        new_superseded_by = sid,
        stringsAsFactors = FALSE
      )
    }
  } else if (identical(t1$status_label, "Superseded Preprint")) {
    # Candidate itself is a preprint that already has a published version!
    matched_pub_id <- ""
    if (!is.null(existing_studies) && nrow(existing_studies) > 0L) {
      if (nzchar(t1$published_upgrade_pmid)) {
        idx <- which(as.character(existing_studies$pmid) == t1$published_upgrade_pmid)
        if (length(idx) > 0L) matched_pub_id <- existing_studies$study_id[[idx[[1L]]]]
      }
      if (!nzchar(matched_pub_id) && nzchar(t1$published_upgrade_doi)) {
        idx <- which(tolower(trimws(existing_studies$doi)) == t1$published_upgrade_doi)
        if (length(idx) > 0L) matched_pub_id <- existing_studies$study_id[[idx[[1L]]]]
      }
    }
    if (nzchar(matched_pub_id)) {
      superseded_by_val <- matched_pub_id
      superseded_by_status <- "API_EXTRACTED"
      superseded_by_note <- sprintf("Linked to catalogued published study %s (DOI: %s)", matched_pub_id, t1$published_upgrade_doi)
    } else {
      superseded_by_val <- ""
      superseded_by_status <- "AMBIGUOUS"
      superseded_by_note <- sprintf(
        "Published version detected (DOI: %s %s) — curate published study first so foreign key exists!",
        t1$published_upgrade_doi, t1$published_upgrade_pmid
      )
    }
  }

  # --- 2. Emit all 14 atomic rows for `studies.tsv` ---
  add_p(sid, "studies.tsv", sid, "study_id", "REQUIRED", sid, "API_EXTRACTED", "Tier 1", "Canonical study_id generator", sid)
  add_p(sid, "studies.tsv", sid, "first_author", "REQUIRED",
        t1$first_author, ifelse(nzchar(t1$first_author), "API_EXTRACTED", "MISSING"),
        "Tier 1", "EuropePMC / Crossref authorList[1]", t1$first_author)
  add_p(sid, "studies.tsv", sid, "year", "REQUIRED",
        ifelse(is.na(t1$year), "", as.character(t1$year)),
        ifelse(!is.na(t1$year), "API_EXTRACTED", "MISSING"),
        "Tier 1", "EuropePMC / bioRxiv pubYear", as.character(t1$year %||% ""))
  add_p(sid, "studies.tsv", sid, "status", "REQUIRED",
        t1$status_label, "API_EXTRACTED", "Tier 1",
        "EuropePMC source + Crossref/bioRxiv relation",
        sprintf("status=%s (published_upgrade_doi=%s)", t1$status_label, t1$published_upgrade_doi),
        superseded_by_note)
  add_p(sid, "studies.tsv", sid, "superseded_by",
        ifelse(identical(t1$status_label, "Superseded Preprint"), "CONDITIONALLY_REQUIRED", "OPTIONAL"),
        superseded_by_val, superseded_by_status, "Tier 1",
        "Crossref is-preprint-of / bioRxiv published", superseded_by_val, superseded_by_note)
  add_p(sid, "studies.tsv", sid, "pmid",
        ifelse(identical(t1$status_label, "Published"), "CONDITIONALLY_REQUIRED", "OPTIONAL"),
        t1$pmid, ifelse(nzchar(t1$pmid) || t1$is_preprint, "API_EXTRACTED", "MISSING"),
        "Tier 1", "EuropePMC / PubMed PMID", t1$pmid)
  add_p(sid, "studies.tsv", sid, "doi", "REQUIRED",
        t1$doi, ifelse(nzchar(t1$doi), "API_EXTRACTED", "MISSING"),
        "Tier 1", "EuropePMC / Crossref DOI", t1$doi)
  add_p(sid, "studies.tsv", sid, "journal", "REQUIRED",
        t1$journal, ifelse(nzchar(t1$journal), "API_EXTRACTED", "MISSING"),
        "Tier 1", "EuropePMC / bioRxiv journal", t1$journal)
  add_p(sid, "studies.tsv", sid, "title", "REQUIRED",
        t1$title, ifelse(nzchar(t1$title), "API_EXTRACTED", "MISSING"),
        "Tier 1", "EuropePMC / Crossref title", t1$title)

  repo_req <- "REQUIRED"
  add_p(sid, "studies.tsv", sid, "repository", repo_req,
        ifelse(nzchar(t3_repo$repository), t3_repo$repository, "Not specified"),
        t3_repo$access_status, t3_repo$access_tier, t3_repo$location, t3_repo$snippet, t3_repo$access_note)
  add_p(sid, "studies.tsv", sid, "accession_or_url", repo_req,
        t3_repo$accession_or_url, t3_repo$access_status, t3_repo$access_tier, t3_repo$location, t3_repo$snippet)
  add_p(sid, "studies.tsv", sid, "access_level", repo_req,
        ifelse(nzchar(t3_repo$access_level), t3_repo$access_level, "Not available"),
        t3_repo$access_status, t3_repo$access_tier, t3_repo$location, t3_repo$snippet, t3_repo$access_note)
  add_p(sid, "studies.tsv", sid, "added_release", "OPTIONAL",
        added_release, "API_EXTRACTED", "Tier 1", "Pipeline staging tag", added_release)
  add_p(sid, "studies.tsv", sid, "verified_full_text", "REQUIRED",
        ifelse(isTRUE(t2$verified_full_text), "TRUE", "FALSE"),
        ifelse(isTRUE(t2$verified_full_text), "EXTRACTED", "MISSING"),
        ifelse(isTRUE(t2$verified_full_text), "Tier 2", "Unresolved"),
        t2$ft_source,
        ifelse(isTRUE(t2$verified_full_text), paste0("Verified legal full text via ", t2$ft_source), t2$ft_gate_note),
        t2$ft_gate_note)

  # --- 3. Reconcile study-wide stratum attributes across Tier 1, Tier 2, and Tier 3 ---
  cohorts_vec <- t2$cohorts
  cohorts_disp <- if (length(cohorts_vec) > 0L) paste(cohorts_vec, collapse = "; ") else "Unspecified Cohort"
  cohorts_status <- if (length(cohorts_vec) > 0L) (if (t2$verified_full_text) "EXTRACTED" else "API_EXTRACTED") else "MISSING"

  # Sample matrix ([C-2] 12-value enum)
  mat_val <- "Plasma (unspecified)"
  mat_status <- "MISSING"
  mat_tier <- "Unresolved"
  mat_loc <- ""
  mat_snip <- ""
  mat_note <- ""

  if (nzchar(t3_supp$sample_matrix$value) && !isTRUE(t3_supp$sample_matrix$ambiguous)) {
    mat_val <- t3_supp$sample_matrix$value
    mat_status <- "SUPPLEMENT_EXTRACTED"
    mat_tier <- "Tier 3"
    mat_loc <- t3_supp$sample_matrix$location
    mat_snip <- t3_supp$sample_matrix$snippet
  } else if (nzchar(t2$sample_matrix$value)) {
    mat_val <- t2$sample_matrix$value
    mat_loc <- t2$sample_matrix$location
    mat_snip <- t2$sample_matrix$snippet
    mat_tier <- if (grepl("^Tier 2", mat_loc)) "Tier 2" else "Tier 1"
    if (isTRUE(t2$sample_matrix$ambiguous)) {
      mat_status <- "AMBIGUOUS"
      mat_note <- "Only generic 'plasma' detected; verify tube anticoagulant (EDTA vs Citrate vs Heparin) in Methods."
    } else {
      mat_status <- if (identical(mat_tier, "Tier 2")) "EXTRACTED" else "API_EXTRACTED"
    }
  } else {
    mat_val <- "Other tissue"
    mat_status <- "MISSING"
  }

  # Platform category & version
  plat_cat <- if (nzchar(t2$platform$category)) t2$platform$category else "Other Immunoassay"
  plat_cat_status <- if (nzchar(t2$platform$category)) (if (grepl("^Tier 2", t2$platform$location)) "EXTRACTED" else "API_EXTRACTED") else "MISSING"
  plat_ver <- t2$platform$version
  plat_ver_tier <- if (grepl("^Tier 2", t2$platform$location)) "Tier 2" else "Tier 1"
  plat_ver_status <- plat_cat_status
  plat_ver_loc <- t2$platform$location
  plat_ver_snip <- t2$platform$snippet

  if (!is.null(t3_supp$xlsx_parsed) && nzchar(t3_supp$xlsx_parsed$platform_version) &&
      (grepl("unspecified", plat_ver, ignore.case = TRUE) || !nzchar(plat_ver))) {
    plat_ver <- t3_supp$xlsx_parsed$platform_version
    plat_ver_tier <- "Tier 3"
    plat_ver_status <- "SUPPLEMENT_EXTRACTED"
    plat_ver_loc <- t3_supp$xlsx_parsed$assay_location
    plat_ver_snip <- t3_supp$xlsx_parsed$assay_snippet
  }
  if (grepl("unspecified", plat_ver, ignore.case = TRUE)) plat_ver_status <- "AMBIGUOUS"

  # n_assays & n_unique_proteins
  xl <- t3_supp$xlsx_parsed
  n_assays_val <- ""
  n_assays_status <- "MISSING"
  n_assays_tier <- "Unresolved"
  n_assays_loc <- ""
  n_assays_snip <- ""
  n_assays_note <- ""

  if (!is.null(xl) && !is.na(xl$n_assays)) {
    n_assays_val <- as.character(xl$n_assays)
    n_assays_status <- "SUPPLEMENT_EXTRACTED"
    n_assays_tier <- "Tier 3"
    n_assays_loc <- xl$assay_location
    n_assays_snip <- xl$assay_snippet
  } else if (!is.na(t2$methods$n_assays)) {
    n_assays_val <- as.character(t2$methods$n_assays)
    n_assays_status <- if (t2$verified_full_text) "EXTRACTED" else "API_EXTRACTED"
    n_assays_tier <- if (t2$verified_full_text) "Tier 2" else "Tier 1"
    n_assays_loc <- if (t2$verified_full_text) "JATS PROTEOMIC_QC / Overview" else "Tier 1 Abstract"
    n_assays_snip <- t2$methods$prot_snippet
  }

  n_uprot_val <- ""
  n_uprot_status <- "MISSING"
  n_uprot_tier <- "Unresolved"
  n_uprot_loc <- ""
  n_uprot_snip <- ""
  n_uprot_note <- ""

  if (t2$verified_full_text && !is.na(t2$methods$n_unique_proteins)) {
    n_uprot_val <- as.character(t2$methods$n_unique_proteins)
    n_uprot_status <- "EXTRACTED"
    n_uprot_tier <- "Tier 2"
    n_uprot_loc <- "JATS Overview / PROTEOMIC_QC"
    n_uprot_snip <- t2$methods$prot_snippet
    if (!is.null(xl) && !is.na(xl$n_unique_proteins) && xl$n_unique_proteins != t2$methods$n_unique_proteins) {
      n_uprot_note <- sprintf("Supplement table lists %d unique IDs vs %d in main text.", xl$n_unique_proteins, t2$methods$n_unique_proteins)
    }
  } else if (!is.null(xl) && !is.na(xl$n_unique_proteins)) {
    n_uprot_val <- as.character(xl$n_unique_proteins)
    n_uprot_tier <- "Tier 3"
    n_uprot_loc <- xl$assay_location
    n_uprot_snip <- xl$assay_snippet
    # When full text is unverified and n_assays > n_unique_proteins (e.g. SomaScan 4,907 aptamers -> 4,667 non-NA UniProts vs 4,719 prose), flag AMBIGUOUS for curator confirmation
    if (!t2$verified_full_text && !is.na(xl$n_assays) && xl$n_assays != xl$n_unique_proteins) {
      n_uprot_status <- "AMBIGUOUS"
      n_uprot_note <- sprintf("Counted %d non-NA unique UniProt IDs in supplement; verify author-reported unique protein count in Methods.", xl$n_unique_proteins)
    } else {
      n_uprot_status <- "SUPPLEMENT_EXTRACTED"
    }
  } else if (!is.na(t2$methods$n_unique_proteins)) {
    n_uprot_val <- as.character(t2$methods$n_unique_proteins)
    n_uprot_status <- "AMBIGUOUS"
    n_uprot_tier <- "Tier 1"
    n_uprot_loc <- "Tier 1 Abstract"
    n_uprot_snip <- t2$methods$prot_snippet
    n_uprot_note <- "Verify assay count vs unique protein count in full text."
  }

  # Genome build (Tier 2 prose + Tier 3 explicit header + Tier 3e Coordinate Build Resolver)
  gb_val <- "Unknown"
  gb_status <- "MISSING"
  gb_tier <- "Unresolved"
  gb_loc <- ""
  gb_snip <- ""
  gb_note <- ""

  if (nzchar(t2$methods$genome_build)) {
    gb_val <- t2$methods$genome_build
    gb_status <- "EXTRACTED"
    gb_tier <- "Tier 2"
    gb_loc <- "JATS GENOMIC_IMPUTATION"
    gb_snip <- t2$methods$build_snippet
    if (isTRUE(t3_supp$build_resolved$resolved) && !identical(t3_supp$build_resolved$genome_build, gb_val)) {
      gb_status <- "CONFLICT"
      gb_note <- sprintf("Prose states %s but supplement rsID coordinates match %s!", gb_val, t3_supp$build_resolved$genome_build)
    }
  } else if (isTRUE(t3_supp$build_resolved$resolved)) {
    gb_val <- t3_supp$build_resolved$genome_build
    gb_status <- "SUPPLEMENT_EXTRACTED"
    gb_tier <- "Tier 3"
    gb_loc <- t3_supp$build_resolved$location
    gb_snip <- t3_supp$build_resolved$evidence_snippet
    gb_note <- "Resolved deterministically via sentinel rsID coordinate lookup (genome build omitted from prose)."
  } else if (!is.null(xl) && isTRUE(xl$hg38_header_explicit)) {
    gb_val <- "GRCh38"
    gb_status <- "SUPPLEMENT_EXTRACTED"
    gb_tier <- "Tier 3"
    gb_loc <- xl$pqtl_cols_location
    gb_snip <- xl$pqtl_cols_snippet
  }

  # Coverage & cis_window_kb
  cov_val <- if (nzchar(t2$methods$coverage)) t2$methods$coverage else if (nzchar(t3_repo$repository)) "genome_wide" else "none"
  cov_status <- if (nzchar(t2$methods$coverage) || nzchar(t3_repo$repository)) {
    if (t2$verified_full_text) "EXTRACTED" else if (nzchar(t3_repo$repository)) "SUPPLEMENT_EXTRACTED" else "API_EXTRACTED"
  } else {
    "AMBIGUOUS"
  }
  cov_tier <- if (t2$verified_full_text) "Tier 2" else if (nzchar(t3_repo$repository)) "Tier 3" else "Tier 1"

  # Booleans has_beta_se & has_eaf ([C-2]: NEVER emit NA!)
  has_bs_bool <- (!is.null(xl) && isTRUE(xl$has_beta_se))
  has_eaf_bool <- (!is.null(xl) && isTRUE(xl$has_eaf)) || grepl("assocvariants\\.annotated", t3_repo$helper_files_note)
  bs_status <- if (has_bs_bool) "SUPPLEMENT_EXTRACTED" else "AMBIGUOUS"
  eaf_status <- if (has_eaf_bool) "SUPPLEMENT_EXTRACTED" else "AMBIGUOUS"
  bs_req <- if (identical(t3_repo$access_level, "Open download")) "CONDITIONALLY_REQUIRED" else "OPTIONAL"

  indel_enc <- if (!is.null(xl) && isTRUE(xl$multiallelic_bang)) {
    "sequence_plus_multiallelic"
  } else if (has_bs_bool) {
    "sequence"
  } else {
    ""
  }

  yaml_name <- ""
  if (!is.null(manifests_dir) && dir.exists(manifests_dir)) {
    cand_yaml <- paste0(sid, ".yaml")
    if (file.exists(file.path(manifests_dir, cand_yaml))) yaml_name <- cand_yaml
  }

  # --- 4. Emit all 26 atomic rows per stratum in `strata.tsv` ---
  strata_items <- t2$strata_info$strata
  if (length(strata_items) == 0L) {
    strata_items <- list(
      UNK = list(stratum_suffix = "UNK", ancestry_group = "EUR", ancestry_detail = "Unknown (requires manual check)",
                 sex = "Combined", n_samples = 1L, location = "Unresolved", snippet = "", n_note = "Placeholder N=1; missing sample size")
    )
  }

  au_yr_prefix <- sub("_[0-9]{6,9}$", "", sid)
  for (st_it in strata_items) {
    stid <- paste0(au_yr_prefix, "_", st_it$stratum_suffix)
    n_status <- if (identical(st_it$stratum_suffix, "UNK")) {
      "MISSING"
    } else if (!t2$verified_full_text && t1$is_preprint) {
      "AMBIGUOUS"
    } else if (grepl("^Tier 1", st_it$location)) {
      "API_EXTRACTED"
    } else {
      "EXTRACTED"
    }
    n_tier <- if (grepl("^Tier 1", st_it$location)) "Tier 1" else "Tier 2"

    harm_parts <- c(
      st_it$n_note,
      if (!is.null(xl)) xl$mixed_id_build_note else "",
      t3_repo$helper_files_note,
      t3_repo$access_note
    )
    harm_parts <- harm_parts[nzchar(harm_parts)]
    harm_txt <- paste(unique(harm_parts), collapse = " ")

    add_p(sid, "strata.tsv", stid, "stratum_id", "REQUIRED", stid, "API_EXTRACTED", "Tier 1", "Stratum ID builder", stid)
    add_p(sid, "strata.tsv", stid, "study_id", "REQUIRED", sid, "API_EXTRACTED", "Tier 1", "Foreign key to studies.tsv", sid)
    add_p(sid, "strata.tsv", stid, "cohorts_display", "REQUIRED", cohorts_disp, cohorts_status,
          ifelse(t2$verified_full_text, "Tier 2", "Tier 1"), "Cohort entity detector", cohorts_disp)
    add_p(sid, "strata.tsv", stid, "ancestry_group", "REQUIRED", st_it$ancestry_group, n_status, n_tier, st_it$location, st_it$snippet)
    add_p(sid, "strata.tsv", stid, "ancestry_detail", "REQUIRED", st_it$ancestry_detail, n_status, n_tier, st_it$location, st_it$snippet)
    add_p(sid, "strata.tsv", stid, "sex", "REQUIRED", st_it$sex,
          ifelse(t2$verified_full_text, "EXTRACTED", "API_EXTRACTED"), n_tier, st_it$location, "Sex-combined analysis")
    add_p(sid, "strata.tsv", stid, "n_samples", "REQUIRED", as.character(st_it$n_samples), n_status, n_tier, st_it$location, st_it$snippet, st_it$n_note)
    add_p(sid, "strata.tsv", stid, "sample_matrix", "REQUIRED", mat_val, mat_status, mat_tier, mat_loc, mat_snip, mat_note)
    add_p(sid, "strata.tsv", stid, "platform_category", "REQUIRED", plat_cat, plat_cat_status, plat_ver_tier, plat_ver_loc, plat_ver_snip)
    add_p(sid, "strata.tsv", stid, "platform_version", "REQUIRED", plat_ver, plat_ver_status, plat_ver_tier, plat_ver_loc, plat_ver_snip)
    add_p(sid, "strata.tsv", stid, "n_assays", "REQUIRED", n_assays_val, n_assays_status, n_assays_tier, n_assays_loc, n_assays_snip, n_assays_note)
    add_p(sid, "strata.tsv", stid, "n_unique_proteins", "REQUIRED", n_uprot_val, n_uprot_status, n_uprot_tier, n_uprot_loc, n_uprot_snip, n_uprot_note)
    add_p(sid, "strata.tsv", stid, "genome_build", "REQUIRED", gb_val, gb_status, gb_tier, gb_loc, gb_snip, gb_note)
    add_p(sid, "strata.tsv", stid, "imputation_panel", "OPTIONAL",
          t2$methods$imputation_panel,
          ifelse(nzchar(t2$methods$imputation_panel), "EXTRACTED", "MISSING"),
          ifelse(nzchar(t2$methods$imputation_panel), "Tier 2", "Unresolved"),
          "JATS GENOMIC_IMPUTATION", t2$methods$imputation_panel,
          if (!nzchar(t2$methods$imputation_panel)) "Optional: often delegated to prior cohort genotyping reference." else "")
    add_p(sid, "strata.tsv", stid, "coverage", "REQUIRED", cov_val, cov_status, cov_tier, "Coverage classifier", cov_val)
    add_p(sid, "strata.tsv", stid, "cis_window_kb",
          ifelse(identical(cov_val, "cis_only"), "CONDITIONALLY_REQUIRED", "OPTIONAL"),
          ifelse(is.na(t2$methods$cis_window_kb), "", as.character(t2$methods$cis_window_kb)),
          ifelse(!identical(cov_val, "cis_only") || !is.na(t2$methods$cis_window_kb), "EXTRACTED", "MISSING"),
          "Tier 2", "GWAS cis-window parser", "")
    add_p(sid, "strata.tsv", stid, "transformation", "REQUIRED",
          ifelse(nzchar(t2$methods$transformation), t2$methods$transformation, "Unspecified"),
          ifelse(nzchar(t2$methods$transformation), "EXTRACTED", "MISSING"),
          ifelse(nzchar(t2$methods$transformation), "Tier 2", "Unresolved"),
          "JATS GWAS_ASSOCIATION", t2$methods$tr_snippet)
    add_p(sid, "strata.tsv", stid, "effect_scale", "REQUIRED",
          ifelse(nzchar(t2$methods$effect_scale), t2$methods$effect_scale, "SD_IRNT"),
          ifelse(nzchar(t2$methods$effect_scale), "EXTRACTED", "MISSING"),
          ifelse(nzchar(t2$methods$effect_scale), "Tier 2", "Unresolved"),
          "JATS GWAS_ASSOCIATION", t2$methods$tr_snippet,
          if (!nzchar(t2$methods$effect_scale)) "Defaulted placeholder SD_IRNT; MISSING status forces curator verification in Methods." else "")
    add_p(sid, "strata.tsv", stid, "covariates", "CONDITIONALLY_REQUIRED",
          t2$methods$covariates,
          ifelse(nzchar(t2$methods$covariates), "EXTRACTED", "MISSING"),
          ifelse(nzchar(t2$methods$covariates), "Tier 2", "Unresolved"),
          "JATS GWAS_ASSOCIATION", t2$methods$cov_snippet)
    add_p(sid, "strata.tsv", stid, "gwas_software", "CONDITIONALLY_REQUIRED",
          t2$methods$gwas_software,
          ifelse(nzchar(t2$methods$gwas_software), "EXTRACTED", "MISSING"),
          ifelse(nzchar(t2$methods$gwas_software), "Tier 2", "Unresolved"),
          "JATS GWAS_ASSOCIATION", t2$methods$sw_snippet)
    add_p(sid, "strata.tsv", stid, "file_format", "OPTIONAL",
          t3_repo$detected_file_format,
          ifelse(nzchar(t3_repo$detected_file_format), "SUPPLEMENT_EXTRACTED", "MISSING"),
          ifelse(nzchar(t3_repo$detected_file_format), "Tier 3", "Unresolved"),
          t3_repo$location, t3_repo$detected_file_format)
    add_p(sid, "strata.tsv", stid, "has_beta_se", bs_req,
          ifelse(has_bs_bool, "TRUE", "FALSE"), bs_status,
          ifelse(has_bs_bool, "Tier 3", "Unresolved"),
          if (!is.null(xl)) xl$pqtl_cols_location else "",
          if (!is.null(xl)) xl$pqtl_cols_snippet else "Boolean non-NA enforced ([C-2])")
    add_p(sid, "strata.tsv", stid, "has_eaf", bs_req,
          ifelse(has_eaf_bool, "TRUE", "FALSE"), eaf_status,
          ifelse(has_eaf_bool, "Tier 3", "Unresolved"),
          if (!is.null(xl)) xl$pqtl_cols_location else t3_repo$location,
          if (!is.null(xl)) xl$pqtl_cols_snippet else t3_repo$helper_files_note)
    add_p(sid, "strata.tsv", stid, "indel_encoding", "OPTIONAL",
          indel_enc,
          ifelse(nzchar(indel_enc), "SUPPLEMENT_EXTRACTED", "MISSING"),
          ifelse(nzchar(indel_enc), "Tier 3", "Unresolved"),
          if (!is.null(xl)) xl$pqtl_cols_location else "", indel_enc)
    add_p(sid, "strata.tsv", stid, "gwasprocessing_yaml", "OPTIONAL",
          yaml_name, ifelse(nzchar(yaml_name), "EXTRACTED", "MISSING"),
          "Tier 3", "GWASprocessing manifests check", yaml_name)
    add_p(sid, "strata.tsv", stid, "harmonisation_notes", "OPTIONAL",
          harm_txt, ifelse(nzchar(harm_txt), "SUPPLEMENT_EXTRACTED", "EXTRACTED"),
          "Tier 3", "Automated harmonisation summary", harm_txt)
  }

  # --- 5. Emit all 3 atomic rows per cohort in `cohort_membership.tsv` ([C-5], [M-7]) ---
  max_gwas_n <- max(vapply(strata_items, function(x) as.integer(x$n_samples), integer(1L)))
  coh_list <- if (length(cohorts_vec) > 0L) cohorts_vec else "Unspecified Cohort"
  for (ci in seq_along(coh_list)) {
    cnm <- coh_list[[ci]]
    cid_rec <- paste0(sid, ":", cnm)
    # Single-cohort conservation rule ([C-5]): n_cohort equals max post-QC GWAS N
    coh_n <- if (length(coh_list) == 1L) max_gwas_n else NA_integer_
    coh_n_status <- if (!is.na(coh_n) && !identical(cohorts_status, "MISSING")) {
      if (t2$verified_full_text) "EXTRACTED" else "API_EXTRACTED"
    } else {
      "MISSING"
    }
    add_p(sid, "cohort_membership.tsv", cid_rec, "study_id", "REQUIRED", sid, "API_EXTRACTED", "Tier 1", "Foreign key", sid)
    add_p(sid, "cohort_membership.tsv", cid_rec, "cohort", "REQUIRED", cnm, cohorts_status,
          ifelse(t2$verified_full_text, "Tier 2", "Tier 1"), "Cohort detector", cnm)
    add_p(sid, "cohort_membership.tsv", cid_rec, "n_cohort", "REQUIRED",
          ifelse(is.na(coh_n), as.character(max_gwas_n), as.character(coh_n)),
          coh_n_status,
          ifelse(t2$verified_full_text, "Tier 2", "Tier 1"),
          "Post-QC GWAS sample conservation ([C-5])",
          sprintf("n_cohort=%s (max stratum GWAS N=%d)", ifelse(is.na(coh_n), "verify per-cohort split", as.character(coh_n)), max_gwas_n))
  }

  prov_df <- do.call(rbind, prov)

  # --- 6. Evaluate Escalation Rubric ---
  blocking_mask <- (prov_df$required_tier %in% c("REQUIRED", "CONDITIONALLY_REQUIRED")) &
    (prov_df$status %in% c("MISSING", "AMBIGUOUS", "CONFLICT"))
  blocking_df <- prov_df[blocking_mask, , drop = FALSE]

  needs_review <- (!isTRUE(t2$verified_full_text)) || (nrow(blocking_df) > 0L)
  decision <- if (needs_review) "NEEDS_MANUAL_REVIEW" else "READY_FOR_CURATOR_SIGNOFF"
  queue_curation_status <- if (needs_review) {
    "Pending Curator Review (Needs Manual Check)"
  } else {
    "Pending Curator Review (Auto-Extracted)"
  }

  # --- 7. Pivot atomic provenance into exact staging TSVs ---
  studies_cols <- c(
    "study_id", "first_author", "year", "status", "superseded_by",
    "pmid", "doi", "journal", "title", "repository", "accession_or_url",
    "access_level", "added_release", "verified_full_text"
  )
  strata_cols <- c(
    "stratum_id", "study_id", "cohorts_display", "ancestry_group", "ancestry_detail",
    "sex", "n_samples", "sample_matrix", "platform_category", "platform_version",
    "n_assays", "n_unique_proteins", "genome_build", "imputation_panel", "coverage",
    "cis_window_kb", "transformation", "effect_scale", "covariates", "gwas_software",
    "file_format", "has_beta_se", "has_eaf", "indel_encoding", "gwasprocessing_yaml",
    "harmonisation_notes"
  )
  cohorts_cols <- c("study_id", "cohort", "n_cohort")

  draft_studies <- pivot_provenance_to_table(prov_df, "studies.tsv", studies_cols)
  draft_strata  <- pivot_provenance_to_table(prov_df, "strata.tsv", strata_cols)
  draft_cohorts <- pivot_provenance_to_table(prov_df, "cohort_membership.tsv", cohorts_cols)

  prov_path    <- file.path(out_dir, paste0(sid, "_provenance.tsv"))
  st_path      <- file.path(out_dir, paste0(sid, "_draft_studies.tsv"))
  str_path     <- file.path(out_dir, paste0(sid, "_draft_strata.tsv"))
  coh_path     <- file.path(out_dir, paste0(sid, "_draft_cohort_membership.tsv"))
  chk_path     <- file.path(out_dir, paste0(sid, "_review_checklist.md"))

  utils::write.table(prov_df, prov_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  utils::write.table(draft_studies, st_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  utils::write.table(draft_strata, str_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  utils::write.table(draft_cohorts, coh_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")

  patch_path <- ""
  if (!is.null(preprint_patch_df)) {
    patch_path <- file.path(out_dir, paste0(sid, "_preprint_supersede_patch.tsv"))
    utils::write.table(preprint_patch_df, patch_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  }

  # Deduplicate blocking items across strata when study-wide fields share identical issue
  chk_lines <- c(
    sprintf("## Extraction Summary: `%s` (`%s`)", sid, decision),
    "",
    sprintf("- **Recommended `candidate_queue$curation_status`:** `%s`", queue_curation_status),
    sprintf("- **Title:** %s", t1$title),
    sprintf("- **Identifiers:** `PMID: %s` | `DOI: %s` | `PMCID: %s` | `EPMC: %s`",
            ifelse(nzchar(t1$pmid), t1$pmid, "None"),
            ifelse(nzchar(t1$doi), t1$doi, "None"),
            ifelse(nzchar(t1$pmcid), t1$pmcid, "None"),
            ifelse(nzchar(t1$epmc_id), t1$epmc_id, "None")),
    sprintf("- **Verified Legal Full Text:** `%s` (%s)", draft_studies$verified_full_text[[1L]],
            ifelse(isTRUE(t2$verified_full_text), t2$ft_source, t2$ft_gate_note)),
    sprintf("- **Discovered Strata (%d):** %s", nrow(draft_strata),
            paste(sprintf("`%s` (%s, N=%s)", draft_strata$stratum_id, draft_strata$ancestry_group, draft_strata$n_samples), collapse = ", ")),
    sprintf("- **Repository & Access Level:** `%s` — `%s` (`%s`)",
            draft_studies$repository[[1L]], draft_studies$access_level[[1L]], draft_studies$accession_or_url[[1L]]),
    ""
  )

  if (!is.null(preprint_patch_df)) {
    chk_lines <- c(
      chk_lines,
      "### Preprint Lineage Patch Detected",
      sprintf("- [ ] Update `%s` (`%s`) in `catalogue/studies.tsv`: set `status = \"Superseded Preprint\"` and `superseded_by = \"%s\"`.",
              preprint_patch_df$preprint_study_id[[1L]], preprint_patch_df$preprint_doi[[1L]], sid),
      ""
    )
  }

  if (nrow(blocking_df) > 0L) {
    # Collapse identical (target_table, field, status, conflict_notes) across multiple strata for crisp checklist
    blocking_df$dedup_key <- paste(blocking_df$target_table, blocking_df$field, blocking_df$status, blocking_df$value, sep = "||")
    uniq_blk <- blocking_df[!duplicated(blocking_df$dedup_key), , drop = FALSE]
    chk_lines <- c(chk_lines, sprintf("### Actionable Curator Checklist (%d Unresolved Required Item(s))", nrow(uniq_blk)))
    for (i in seq_len(nrow(uniq_blk))) {
      b <- uniq_blk[i, ]
      note_str <- if (nzchar(b$conflict_notes)) b$conflict_notes else b$evidence_snippet
      chk_lines <- c(chk_lines, sprintf(
        "- [ ] **`%s$%s` (`%s`, current=`\"%s\"`):** %s",
        b$target_table, b$field, b$status, b$value, ifelse(nzchar(note_str), note_str, "Verify in full text / supplement.")
      ))
    }
  } else {
    chk_lines <- c(chk_lines, "### Actionable Curator Checklist", "- [x] All required & conditionally required fields resolved cleanly! Ready for curator sign-off.")
  }

  writeLines(chk_lines, chk_path)

  list(
    study_id = sid,
    decision = decision,
    recommended_curation_status = queue_curation_status,
    n_provenance_rows = nrow(prov_df),
    n_strata = nrow(draft_strata),
    n_blocking_fields = nrow(blocking_df),
    files = list(
      provenance = prov_path,
      draft_studies = st_path,
      draft_strata = str_path,
      draft_cohorts = coh_path,
      checklist = chk_path,
      preprint_patch = patch_path
    )
  )
}
