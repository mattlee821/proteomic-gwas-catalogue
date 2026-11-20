#!/usr/bin/env Rscript
# Standalone Post-Hoc Blinded Benchmark Evaluator for proteomic-gwas-catalogue.
#
# STRICT ISOLATION CONTRACT:
#   - Never invoked by src/run_literature_search.R or CI discovery workflows.
#   - Defaults to --mode blinded-summary, which emits ONLY aggregate recall,
#     per-source unique lift, conditional triage survival rates, stratified
#     breakdowns, and 8-class failure root-cause counts — ZERO study IDs,
#     DOIs, authors, or titles are printed or written to disk!
#   - Row-level unblinded inspection requires BOTH --mode unblinded-audit AND
#     the explicit --unblind-audit confirmation flag (logged to disk).
#
# Usage:
#   Rscript src/evaluate_positive_controls.R \
#     [--controls data/positive_controls.tsv] \
#     [--dump data/search_dump_LATEST.tsv] \
#     [--studies catalogue/studies.tsv] \
#     [--queue catalogue/candidate_queue.tsv] \
#     [--mode blinded-summary|unblinded-audit] [--unblind-audit] \
#     [--out-dir output]

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Required R package 'jsonlite' is not installed.", call. = FALSE)
}

resolve_lexicons_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0L) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
    cand <- file.path(script_dir, "search", "lexicons.R")
    if (file.exists(cand)) return(cand)
  }
  file.path("src", "search", "lexicons.R")
}
source(resolve_lexicons_path())

args <- commandArgs(trailingOnly = TRUE)
get_opt <- function(flag, default = NA_character_) {
  idx <- match(flag, args)
  if (is.na(idx) || idx >= length(args)) default else args[idx + 1L]
}

controls_path <- get_opt("--controls", file.path("data", "positive_controls.tsv"))
dump_path     <- get_opt("--dump", "")
studies_path  <- get_opt("--studies", file.path("catalogue", "studies.tsv"))
queue_path    <- get_opt("--queue", file.path("catalogue", "candidate_queue.tsv"))
eval_mode     <- tolower(get_opt("--mode", "blinded-summary"))
unblind_ok    <- "--unblind-audit" %in% args
out_dir       <- get_opt("--out-dir", "output")
year_from_arg <- suppressWarnings(as.integer(get_opt("--year-from", NA_character_)))
year_to_arg   <- suppressWarnings(as.integer(get_opt("--year-to", NA_character_)))

if (!(eval_mode %in% c("blinded-summary", "unblinded-audit"))) {
  stop("--mode must be 'blinded-summary' or 'unblinded-audit'.", call. = FALSE)
}
if (identical(eval_mode, "unblinded-audit") && !isTRUE(unblind_ok)) {
  stop(
    "Refusing to output unblinded row-level identifiers without explicit '--unblind-audit' confirmation flag. ",
    "Use default '--mode blinded-summary' during search/triage tuning to prevent benchmark leakage.",
    call. = FALSE
  )
}

if (!file.exists(controls_path)) {
  stop("Benchmark file not found at: ", controls_path, call. = FALSE)
}

# Auto-discover latest search_dump_*.tsv in data/ if --dump omitted
if (!nzchar(dump_path) && dir.exists("data")) {
  dumps <- sort(Sys.glob(file.path("data", "search_dump_*.tsv")), decreasing = TRUE)
  if (length(dumps) > 0L) dump_path <- dumps[1L]
}

read_table_flex <- function(path) {
  if (!nzchar(path) || identical(tolower(path), "none") || !file.exists(path)) return(NULL)
  sep_char <- if (grepl("\\.csv$", tolower(path))) "," else "\t"
  utils::read.delim(path, sep = sep_char, stringsAsFactors = FALSE, check.names = FALSE,
                    quote = "\"", colClasses = "character", comment.char = "")
}

ctrl_raw   <- read_table_flex(controls_path)
studies_df <- read_table_flex(studies_path) %||% data.frame()
queue_df   <- read_table_flex(queue_path)   %||% data.frame()
dump_df    <- read_table_flex(dump_path)    %||% data.frame()

if (is.null(ctrl_raw) || nrow(ctrl_raw) == 0L) {
  stop("Benchmark file is empty: ", controls_path, call. = FALSE)
}

# ---------------------------------------------------------------------------
# Flexible Case-Insensitive Header Mapping (M-5)
# ---------------------------------------------------------------------------
norm_hdr <- gsub("[^a-z0-9]+", "_", tolower(trimws(names(ctrl_raw))))
pick_col <- function(aliases, default_val = "") {
  hit <- which(norm_hdr %in% aliases)
  if (length(hit) > 0L) {
    v <- as.character(ctrl_raw[[hit[1L]]])
    v[is.na(v)] <- ""
    trimws(v)
  } else {
    rep(default_val, nrow(ctrl_raw))
  }
}

ctrl <- data.frame(
  row_idx        = seq_len(nrow(ctrl_raw)),
  control_id     = pick_col(c("control_id", "id", "study_id"), ""),
  pmid_raw       = pick_col(c("pmid", "pmids", "pubmed", "pubmed_id"), ""),
  doi_raw        = pick_col(c("doi", "dois", "preprint_doi", "identifier"), ""),
  first_author   = pick_col(c("first_author", "author", "authors", "lead_author"), ""),
  year           = pick_col(c("year", "pub_year", "publication_year"), ""),
  title          = pick_col(c("title", "study_title", "paper_title"), ""),
  platform       = pick_col(c("platform", "platform_category", "technology", "assay"), "Unspecified"),
  sample_matrix  = pick_col(c("sample_matrix", "matrix", "tissue", "biofluid", "fluid"), "Unspecified"),
  n_samples      = pick_col(c("n_samples", "sample_size", "n", "max_n"), ""),
  expected_scope = pick_col(c("expected_scope", "scope", "in_scope"), "Primary_pGWAS"),
  stringsAsFactors = FALSE
)
ctrl$control_id <- ifelse(nzchar(ctrl$control_id), ctrl$control_id, sprintf("ctrl_%03d", ctrl$row_idx))
ctrl$expected_scope <- ifelse(nzchar(ctrl$expected_scope), ctrl$expected_scope, "Primary_pGWAS")

# Filter to in-scope primary pGWAS controls (tracking out-of-scope count separately)
is_in_scope <- !(tolower(ctrl$expected_scope) %in% c("out_of_scope", "excluded", "false", "0", "non_pgwas"))
n_out_of_scope <- sum(!is_in_scope)
ctrl <- ctrl[is_in_scope, , drop = FALSE]

# Optional publication-year window filter (--year-from / --year-to) for chronological branch-by-branch checks
n_outside_year_window <- 0L
if (!is.na(year_from_arg) || !is.na(year_to_arg)) {
  ctrl_yr_num <- suppressWarnings(as.integer(gsub("[^0-9]", "", ctrl$year)))
  if (all(is.na(ctrl_yr_num))) {
    stop("--year-from / --year-to specified, but no numeric years found in the 'year' column of ", controls_path, call. = FALSE)
  }
  keep_yr <- !is.na(ctrl_yr_num) &
    (is.na(year_from_arg) | ctrl_yr_num >= year_from_arg) &
    (is.na(year_to_arg)   | ctrl_yr_num <= year_to_arg)
  n_outside_year_window <- sum(!keep_yr)
  ctrl <- ctrl[keep_yr, , drop = FALSE]
}

if (nrow(ctrl) == 0L) stop("0 in-scope benchmark controls in requested window after filtering.", call. = FALSE)

is_consortium_author <- function(a) {
  a_low <- tolower(ascii_fold(a %||% ""))
  !nzchar(trimws(a_low)) || grepl("\\b(consortium|project|biobank|group|initiative|network|collaboration|ukb-ppp|scallop|decode|gtex|aric|interval)\\b", a_low, perl = TRUE)
}

# ---------------------------------------------------------------------------
# Build Unified Indexing & Preprint <-> Published Lineage Graph
# ---------------------------------------------------------------------------
# Index 1: Curated studies.tsv
cat_pmid_to_sid <- character(0)
cat_doi_to_sid  <- character(0)
lineage_map_sid <- list() # study_id -> character vector of equivalent DOIs & PMIDs

if (nrow(studies_df) > 0L) {
  for (i in seq_len(nrow(studies_df))) {
    sid <- trimws(studies_df$study_id[i] %||% "")
    p   <- normalize_pmid_one(studies_df$pmid[i] %||% "")
    d   <- normalize_doi_one(studies_df$doi[i] %||% "")
    if (nzchar(p)) cat_pmid_to_sid[p] <- sid
    if (nzchar(d)) cat_doi_to_sid[d]  <- sid
    lineage_map_sid[[sid]] <- unique(c(
      lineage_map_sid[[sid]],
      if (nzchar(p)) paste0("pmid:", p) else character(0),
      if (nzchar(d)) paste0("doi:", d)  else character(0)
    ))
  }
  # Connect superseded_by pairs bidirectionally
  for (i in seq_len(nrow(studies_df))) {
    sid <- trimws(studies_df$study_id[i] %||% "")
    sup <- trimws(studies_df$superseded_by[i] %||% "")
    if (nzchar(sid) && nzchar(sup) && !is.null(lineage_map_sid[[sup]])) {
      comb <- unique(c(lineage_map_sid[[sid]], lineage_map_sid[[sup]]))
      lineage_map_sid[[sid]] <- comb
      lineage_map_sid[[sup]] <- comb
    }
  }
}

# Index 2: Search Dump (preferred uncapped ledger) + Candidate Queue fallback
build_cand_records <- function(dump, queue) {
  recs <- list()
  if (!is.null(dump) && nrow(dump) > 0L) {
    for (i in seq_len(nrow(dump))) {
      p <- normalize_pmid_one(dump$pmid[i] %||% dump$identifier[i] %||% "")
      d <- normalize_doi_one(dump$doi[i] %||% dump$identifier[i] %||% "")
      lin <- dump$preprint_lineage_note[i] %||% ""
      lin_dois  <- extract_dois_from_cell(lin)
      lin_pmids <- extract_pmids_from_cell(lin)
      recs[[length(recs) + 1L]] <- list(
        pmid                = p,
        doi                 = d,
        first_author        = dump$first_author[i] %||% "",
        year                = dump$pub_year[i] %||% "",
        title               = dump$title[i] %||% "",
        matched_sources     = dump$matched_sources[i] %||% "",
        auto_classification = dump$auto_classification[i] %||% "",
        p_score             = suppressWarnings(as.integer(dump$p_score[i] %||% NA)),
        s_score             = suppressWarnings(as.integer(dump$s_score[i] %||% NA)),
        triage_rules        = dump$triage_rules[i] %||% "",
        lineage_dois        = lin_dois,
        lineage_pmids       = lin_pmids
      )
    }
  } else if (!is.null(queue) && nrow(queue) > 0L) {
    for (i in seq_len(nrow(queue))) {
      p <- normalize_pmid_one((if ("pmid" %in% names(queue)) queue$pmid[i] else "") %||% queue$identifier[i] %||% "")
      d <- normalize_doi_one((if ("doi" %in% names(queue)) queue$doi[i] else "") %||% queue$identifier[i] %||% "")
      lin <- if ("preprint_lineage_note" %in% names(queue)) queue$preprint_lineage_note[i] %||% "" else ""
      recs[[length(recs) + 1L]] <- list(
        pmid                = p,
        doi                 = d,
        first_author        = "",
        year                = "",
        title               = queue$title[i] %||% "",
        matched_sources     = if ("matched_sources" %in% names(queue)) queue$matched_sources[i] %||% "" else queue$source[i] %||% "",
        auto_classification = queue$auto_classification[i] %||% "",
        p_score             = if ("p_score" %in% names(queue)) suppressWarnings(as.integer(queue$p_score[i])) else NA_integer_,
        s_score             = if ("s_score" %in% names(queue)) suppressWarnings(as.integer(queue$s_score[i])) else NA_integer_,
        triage_rules        = "",
        lineage_dois        = extract_dois_from_cell(lin),
        lineage_pmids       = extract_pmids_from_cell(lin)
      )
    }
  }
  recs
}

cand_recs <- build_cand_records(dump_df, queue_df)

match_fuzzy_title_author <- function(c_au, c_yr, c_ttl, ref_au, ref_yr, ref_ttl) {
  if (!nzchar(c_ttl) || !nzchar(ref_ttl)) return(FALSE)
  y1 <- suppressWarnings(as.integer(gsub("[^0-9]", "", c_yr)))
  y2 <- suppressWarnings(as.integer(gsub("[^0-9]", "", ref_yr)))
  yr_ok <- is.na(y1) || is.na(y2) || abs(y1 - y2) <= 1L
  if (!yr_ok) return(FALSE)

  jac <- title_jaccard(c_ttl, ref_ttl)
  if (is_consortium_author(c_au) || is_consortium_author(ref_au)) {
    return(jac >= 0.65)
  }
  s1 <- normalize_author_surname(c_au)
  s2 <- normalize_author_surname(ref_au)
  nzchar(s1) && identical(s1, s2) && jac >= 0.40
}

classify_cand_outcome <- function(cls, via_lineage = FALSE) {
  if (grepl("^Auto-Deprioritised:\\s*Secondary", cls)) return("FAIL_TRIAGE_FALSE_2SMR")
  if (grepl("^Auto-Deprioritised:\\s*Low-Signal", cls)) return("FAIL_TRIAGE_LOW_SIGNAL")
  if (grepl("^Auto-Deprioritised", cls)) return("FAIL_TRIAGE_STAGE1_GATE")
  if (isTRUE(via_lineage)) return("PASS_LINEAGE_PREPRINT_PUBLISHED")
  if (grepl("^Needs Human Check", cls)) return("PASS_ACTIONABLE_HUMAN_CHECK")
  "PASS_ACTIONABLE_PRIMARY"
}

# ---------------------------------------------------------------------------
# Evaluate Every Benchmark Row (Supporting 1-to-Many IDs, Lineage & Fuzzy ASCII)
# ---------------------------------------------------------------------------
eval_rows <- vector("list", nrow(ctrl))
demoted_rule_tally <- character(0)

for (i in seq_len(nrow(ctrl))) {
  r <- ctrl[i, ]
  # M-5: Multi-PMID & Multi-DOI extraction (also checking pmid_raw for embedded DOIs & vice versa)
  pmids_i <- unique(c(extract_pmids_from_cell(r$pmid_raw, strip_dois_first = TRUE),
                      extract_pmids_from_cell(r$doi_raw, strip_dois_first = TRUE)))
  dois_i  <- unique(c(extract_dois_from_cell(r$doi_raw),
                      extract_dois_from_cell(r$pmid_raw)))
  tok_keys <- c(if (length(pmids_i)) paste0("pmid:", pmids_i) else character(0),
                if (length(dois_i))  paste0("doi:", dois_i)   else character(0))

  # 1. Check search_dump / candidate_queue FIRST when dump is provided (tests active pipeline!),
  #    and check studies.tsv alongside it
  matched_cand_idxs <- integer(0)
  matched_via_lineage <- FALSE

  for (ci in seq_along(cand_recs)) {
    cr <- cand_recs[[ci]]
    direct_hit <- (nzchar(cr$pmid) && (cr$pmid %in% pmids_i)) ||
                  (nzchar(cr$doi)  && (cr$doi  %in% dois_i))
    lin_hit <- (!direct_hit) && (
      any(cr$lineage_pmids %in% pmids_i) || any(cr$lineage_dois %in% dois_i)
    )
    fuzzy_hit <- (!direct_hit && !lin_hit) &&
      match_fuzzy_title_author(r$first_author, r$year, r$title, cr$first_author, cr$year, cr$title)

    if (direct_hit || lin_hit || fuzzy_hit) {
      matched_cand_idxs <- c(matched_cand_idxs, ci)
      if (lin_hit) matched_via_lineage <- TRUE
    }
  }

  # 2. Check studies.tsv (exact, lineage-equivalent via superseded_by, or ASCII-folded fuzzy)
  in_studies_exact   <- any(pmids_i %in% names(cat_pmid_to_sid)) || any(dois_i %in% names(cat_doi_to_sid))
  in_studies_lineage <- FALSE
  if (!in_studies_exact && length(tok_keys) > 0L && length(lineage_map_sid) > 0L) {
    for (eq_vec in lineage_map_sid) {
      if (any(tok_keys %in% eq_vec)) { in_studies_lineage <- TRUE; break }
    }
  }
  in_studies_fuzzy <- FALSE
  if (!in_studies_exact && !in_studies_lineage && nrow(studies_df) > 0L && nzchar(r$title)) {
    for (si in seq_len(nrow(studies_df))) {
      if (match_fuzzy_title_author(r$first_author, r$year, r$title,
                                   studies_df$first_author[si], studies_df$year[si], studies_df$title[si])) {
        in_studies_fuzzy <- TRUE; break
      }
    }
  }

  # Also check if control gave preprint DOI -> studies.tsv maps it to published PMID -> matched in cand_recs!
  if (length(matched_cand_idxs) == 0L && in_studies_lineage && length(cand_recs) > 0L) {
    eq_all <- unique(unlist(Filter(function(v) any(tok_keys %in% v), lineage_map_sid)))
    for (ci in seq_along(cand_recs)) {
      cr <- cand_recs[[ci]]
      cr_keys <- c(if (nzchar(cr$pmid)) paste0("pmid:", cr$pmid) else character(0),
                   if (nzchar(cr$doi))  paste0("doi:", cr$doi)   else character(0))
      if (any(cr_keys %in% eq_all)) {
        matched_cand_idxs <- c(matched_cand_idxs, ci)
        matched_via_lineage <- TRUE
      }
    }
  }

  # Determine outcome class & matched sources
  srcs_hit <- character(0)
  outcome  <- "FAIL_NOT_RETRIEVED"
  best_p   <- NA_integer_
  best_s   <- NA_integer_

  if (length(matched_cand_idxs) > 0L) {
    sub_cands <- cand_recs[matched_cand_idxs]
    srcs_raw  <- unlist(strsplit(paste(vapply(sub_cands, `[[`, "", "matched_sources"), collapse = ";"), ";"))
    srcs_hit  <- unique(trimws(srcs_raw[nzchar(trimws(srcs_raw))]))

    outcomes_vec <- vapply(sub_cands, function(x) classify_cand_outcome(x$auto_classification, matched_via_lineage), character(1))
    # Prefer PASS outcomes if multiple companion records matched
    pass_priority <- c(
      "PASS_ACTIONABLE_PRIMARY", "PASS_LINEAGE_PREPRINT_PUBLISHED", "PASS_ACTIONABLE_HUMAN_CHECK",
      "FAIL_TRIAGE_FALSE_2SMR", "FAIL_TRIAGE_STAGE1_GATE", "FAIL_TRIAGE_LOW_SIGNAL"
    )
    outcome <- pass_priority[pass_priority %in% outcomes_vec][1L]
    best_p  <- max(vapply(sub_cands, function(x) x$p_score %||% NA_integer_, integer(1)), na.rm = TRUE)
    best_s  <- min(vapply(sub_cands, function(x) x$s_score %||% NA_integer_, integer(1)), na.rm = TRUE)

    if (grepl("^FAIL_TRIAGE", outcome)) {
      rls <- unlist(strsplit(paste(vapply(sub_cands, `[[`, "", "triage_rules"), collapse = ","), ","))
      demoted_rule_tally <- c(demoted_rule_tally, trimws(rls[nzchar(trimws(rls))]))
    }
  } else if (in_studies_exact || in_studies_fuzzy) {
    outcome <- "PASS_IN_STUDIES_TSV"
    srcs_hit <- "Catalogued_studies_tsv"
  } else if (in_studies_lineage) {
    outcome <- "PASS_LINEAGE_PREPRINT_PUBLISHED"
    srcs_hit <- "Catalogued_studies_tsv(Lineage)"
  }

  # Multi-ID companion completeness tracking
  n_sub_ids <- max(1L, length(unique(c(pmids_i, dois_i))))
  yr_num <- suppressWarnings(as.integer(gsub("[^0-9]", "", r$year)))
  era <- if (is.na(yr_num)) "Unspecified" else if (yr_num <= 2016L) "<=2016" else if (yr_num <= 2020L) "2017-2020" else if (yr_num <= 2023L) "2021-2023" else "2024+"

  n_num <- suppressWarnings(as.numeric(gsub("[^0-9]", "", r$n_samples)))
  n_bin <- if (is.na(n_num) || n_num <= 0) "Unspecified" else if (n_num < 200) "<200" else if (n_num <= 2000) "200-2,000" else ">2,000"

  eval_rows[[i]] <- data.frame(
    control_id      = r$control_id,
    pmid_norm       = paste(pmids_i, collapse = ";"),
    doi_norm        = paste(dois_i, collapse = ";"),
    n_identifiers   = n_sub_ids,
    first_author    = r$first_author,
    year            = r$year,
    era             = era,
    platform        = ifelse(nzchar(r$platform), r$platform, "Unspecified"),
    sample_matrix   = ifelse(nzchar(r$sample_matrix), r$sample_matrix, "Unspecified"),
    n_bin           = n_bin,
    title           = r$title,
    retrieved       = outcome != "FAIL_NOT_RETRIEVED",
    actionable_pass = grepl("^PASS_", outcome),
    outcome_class   = outcome,
    matched_sources = paste(srcs_hit, collapse = ";"),
    p_score         = if (is.finite(best_p)) as.integer(best_p) else NA_integer_,
    s_score         = if (is.finite(best_s)) as.integer(best_s) else NA_integer_,
    stringsAsFactors = FALSE
  )
}

res_df <- do.call(rbind, eval_rows)

# ---------------------------------------------------------------------------
# Compute Headline, Stratified & Source-Lift Metrics
# ---------------------------------------------------------------------------
n_total <- nrow(res_df)
n_disc  <- sum(res_df$retrieved)
n_act   <- sum(res_df$actionable_pass)

r_disc   <- 100 * n_disc / n_total
r_e2e    <- 100 * n_act  / n_total
s_triage <- if (n_disc > 0L) 100 * n_act / n_disc else NA_real_

n_dump_total   <- nrow(dump_df)
n_dump_deprior <- if (n_dump_total > 0L) sum(grepl("^Auto-Deprioritised", dump_df$auto_classification)) else NA_integer_
compression_pct <- if (!is.na(n_dump_deprior) && n_dump_total > 0L) 100 * n_dump_deprior / n_dump_total else NA_real_

ALL_OUTCOME_CLASSES <- c(
  "PASS_IN_STUDIES_TSV",
  "PASS_LINEAGE_PREPRINT_PUBLISHED",
  "PASS_ACTIONABLE_PRIMARY",
  "PASS_ACTIONABLE_HUMAN_CHECK",
  "FAIL_TRIAGE_FALSE_2SMR",
  "FAIL_TRIAGE_STAGE1_GATE",
  "FAIL_TRIAGE_LOW_SIGNAL",
  "FAIL_NOT_RETRIEVED"
)
outcome_counts <- setNames(integer(length(ALL_OUTCOME_CLASSES)), ALL_OUTCOME_CLASSES)
tab_out <- table(res_df$outcome_class)
outcome_counts[names(tab_out)] <- as.integer(tab_out)

# Per-Source Recall & Marginal Unique Lift
all_src_tokens <- unique(unlist(strsplit(res_df$matched_sources[nzchar(res_df$matched_sources)], ";")))
source_lift_df <- if (length(all_src_tokens) > 0L) {
  do.call(rbind, lapply(sort(all_src_tokens), function(st) {
    has_st <- vapply(strsplit(res_df$matched_sources, ";"), function(v) st %in% v, logical(1))
    only_st <- vapply(strsplit(res_df$matched_sources, ";"), function(v) identical(v, st), logical(1))
    data.frame(
      source          = st,
      recalled_n      = sum(has_st),
      recall_pct      = round(100 * sum(has_st) / n_total, 1),
      unique_lift_n   = sum(only_st),
      unique_lift_pct = round(100 * sum(only_st) / n_total, 1),
      stringsAsFactors = FALSE
    )
  }))
} else {
  data.frame(source = character(0), recalled_n = integer(0), recall_pct = numeric(0),
             unique_lift_n = integer(0), unique_lift_pct = numeric(0), stringsAsFactors = FALSE)
}

stratify_table <- function(df, col_name) {
  sp <- split(df, df[[col_name]])
  do.call(rbind, lapply(names(sp), function(k) {
    sub <- sp[[k]]
    nt <- nrow(sub); nd <- sum(sub$retrieved); na <- sum(sub$actionable_pass)
    data.frame(
      stratum          = k,
      n_controls       = nt,
      discovered_n     = nd,
      discovery_pct    = round(100 * nd / nt, 1),
      actionable_n     = na,
      actionable_pct   = round(100 * na / nt, 1),
      triage_surv_pct  = if (nd > 0L) round(100 * na / nd, 1) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}

strat_era      <- stratify_table(res_df, "era")
strat_platform <- stratify_table(res_df, "platform")
strat_matrix   <- stratify_table(res_df, "sample_matrix")
strat_nbin     <- stratify_table(res_df, "n_bin")

# ---------------------------------------------------------------------------
# Render Strictly Blinded Summary Report (Zero IDs/Titles/Authors!)
# ---------------------------------------------------------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S", tz = "UTC")

fmt_md_table <- function(df) {
  if (nrow(df) == 0L) return("_None_")
  hdr <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  rows <- apply(df, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(hdr, sep, rows), collapse = "\n")
}

outcome_df <- data.frame(
  outcome_class = names(outcome_counts),
  count         = as.integer(outcome_counts),
  pct_of_total  = round(100 * as.integer(outcome_counts) / n_total, 1),
  stringsAsFactors = FALSE
)

md_lines <- c(
  sprintf("# Positive-Control Benchmark Evaluation (%s)", toupper(eval_mode)),
  "",
  sprintf("- **Timestamp (UTC):** `%s`", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
  sprintf("- **Blinding Mode:** `%s` %s", eval_mode,
          if (identical(eval_mode, "blinded-summary")) "(STRICT: study identifiers, authors & titles suppressed)" else "(UNBLINDED POST-FREEZE AUDIT)"),
  sprintf("- **In-Scope Controls Evaluated:** `%d` *(Year filter: `%s..%s`, outside-window skipped: `%d`; Out-of-scope excluded: `%d`; Multi-ID rows: `%d`)*",
          n_total,
          if (!is.na(year_from_arg)) as.character(year_from_arg) else "min",
          if (!is.na(year_to_arg)) as.character(year_to_arg) else "max",
          n_outside_year_window, n_out_of_scope, sum(res_df$n_identifiers > 1L)),
  sprintf("- **Search Dump Evaluated:** `%s` (`%d` total deduplicated candidates)",
          if (nzchar(dump_path)) basename(dump_path) else "studies.tsv + candidate_queue.tsv only", n_dump_total),
  "",
  "## 1. Headline Recall, Triage Survival & Workload Compression",
  "",
  sprintf("- **Discovery Recall (`R_disc`):** `%d / %d` (**`%.2f%%`**)", n_disc, n_total, r_disc),
  sprintf("- **End-to-End Actionable Recall (`R_e2e`):** `%d / %d` (**`%.2f%%`**)", n_act, n_total, r_e2e),
  sprintf("- **Conditional Triage Survival Rate (`S_triage`):** %s",
          if (!is.na(s_triage)) sprintf("`%d / %d` (**`%.2f%%`**)", n_act, n_disc, s_triage) else "`N/A`"),
  sprintf("- **Triage Noise Compression on Dump:** %s",
          if (!is.na(compression_pct)) sprintf("`%d / %d` (**`%.1f%%`** auto-deprioritised)", n_dump_deprior, n_dump_total, compression_pct) else "`N/A`"),
  "",
  "## 2. Aggregate 8-Class Outcome & Failure Root-Cause Distribution",
  "",
  fmt_md_table(outcome_df),
  "",
  if (length(demoted_rule_tally) > 0L) {
    rt <- as.data.frame(table(rule = demoted_rule_tally), stringsAsFactors = FALSE)
    rt <- rt[order(-rt$Freq), , drop = FALSE]
    paste0("\n### Aggregate Rules Triggered on Falsely Demoted Controls (Identity-Free)\n\n", fmt_md_table(rt), "\n")
  } else "",
  "## 3. Per-Source Recall & Marginal Unique Lift",
  "",
  fmt_md_table(source_lift_df),
  "",
  "## 4. Stratified Recall Breakdowns",
  "",
  "### By Publication Era",
  fmt_md_table(strat_era),
  "",
  "### By Proteomic Platform",
  fmt_md_table(strat_platform),
  "",
  "### By Biological Matrix",
  fmt_md_table(strat_matrix),
  "",
  "### By Sample-Size Bin",
  fmt_md_table(strat_nbin)
)

summary_md_path   <- file.path(out_dir, sprintf("positive_control_%s_%s.md", gsub("-", "_", eval_mode), stamp))
summary_json_path <- file.path(out_dir, sprintf("positive_control_%s_%s.json", gsub("-", "_", eval_mode), stamp))

writeLines(md_lines, con = summary_md_path)
jsonlite::write_json(
  list(
    mode                   = eval_mode,
    generated_at           = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    n_in_scope             = n_total,
    n_out_of_scope         = n_out_of_scope,
    discovery_recall_pct   = round(r_disc, 3),
    actionable_recall_pct  = round(r_e2e, 3),
    triage_survival_pct    = round(s_triage, 3),
    dump_compression_pct   = round(compression_pct, 2),
    outcome_counts         = as.list(outcome_counts),
    source_lift            = source_lift_df,
    stratified_era         = strat_era,
    stratified_platform    = strat_platform,
    stratified_matrix      = strat_matrix,
    stratified_sample_size = strat_nbin
  ),
  path = summary_json_path,
  auto_unbox = TRUE,
  pretty = TRUE
)

cat(paste(md_lines, collapse = "\n"), "\n")
message(sprintf("\nWrote %s report to %s", eval_mode, summary_md_path))

# Only write row-level table when explicitly authorized via --mode unblinded-audit --unblind-audit
if (identical(eval_mode, "unblinded-audit") && isTRUE(unblind_ok)) {
  audit_log <- file.path(out_dir, "unblind_audit_log.txt")
  cat(sprintf("%s\tmode=unblinded-audit\tcontrols=%s\tdump=%s\tn=%d\n",
              format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
              controls_path, dump_path, n_total),
      file = audit_log, append = TRUE)

  rows_path <- file.path(out_dir, sprintf("positive_control_unblinded_rows_%s.tsv", stamp))
  utils::write.table(res_df, file = rows_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  message(sprintf("UNBLINDED AUDIT LOGGED (%s); row-level details written to %s", audit_log, rows_path))
}
