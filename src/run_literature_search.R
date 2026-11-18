#!/usr/bin/env Rscript
# Automated multi-source literature search and triage script for proteomic GWAS.
#
# Usage:
#   Rscript src/run_literature_search.R [--since YYYY-MM-DD] [--until YYYY-MM-DD] [--trigger LABEL] [--dry-run]
#
# Sources queried:
#   1. Europe PMC REST API (covers PubMed/MEDLINE + bioRxiv/medRxiv/Research Square preprints)
#   2. Preprint-to-published upgrade check for active Preprints in catalogue/studies.tsv
#
# Outputs:
#   - Appends unseen candidate records to catalogue/candidate_queue.tsv
#   - Logs execution metrics to catalogue/search_runs.tsv
#   - Validates all catalogue tables via src/validate_catalogue.R

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required.", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
get_opt <- function(flag, default = NA_character_) {
  idx <- match(flag, args)
  if (is.na(idx) || idx >= length(args)) default else args[idx + 1L]
}

dry_run    <- "--dry-run" %in% args
since_date <- get_opt("--since", format(Sys.Date() - 35L, "%Y-%m-%d"))
until_date <- get_opt("--until", format(Sys.Date(), "%Y-%m-%d"))
trigger_lb <- get_opt("--trigger", "Scheduled GitHub Action")

clean_cell <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x[1])) return("")
  gsub("[\r\n\t]+", " ", trimws(as.character(x[1])))
}

message(sprintf("Proteomic GWAS search window: %s to %s (dry_run=%s)", since_date, until_date, dry_run))

read_tsv_safe <- function(path) {
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, quote = "")
}

write_tsv_safe <- function(df, path) {
  utils::write.table(df, file = path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

studies  <- read_tsv_safe(file.path("catalogue", "studies.tsv"))
excluded <- read_tsv_safe(file.path("catalogue", "excluded_records.tsv"))
queue    <- read_tsv_safe(file.path("catalogue", "candidate_queue.tsv"))
runs     <- read_tsv_safe(file.path("catalogue", "search_runs.tsv"))

known_ids <- tolower(unique(trimws(c(
  studies$doi,
  ifelse(nzchar(studies$pmid), paste0("pmid:", studies$pmid), ""),
  studies$pmid,
  excluded$identifier,
  sub("^pmid:", "", tolower(excluded$identifier)),
  queue$identifier,
  sub("^pmid:", "", tolower(queue$identifier))
))))
known_ids <- known_ids[nzchar(known_ids)]

epmc_query <- paste0(
  '("SomaScan" OR "SomaLogic" OR "Olink" OR "NULISA" OR "NULISAseq" OR "SCALLOP consortium" ',
  'OR "mass spectrometry proteomics" OR "DIA-MS" OR "plasma proteome" OR "CSF proteome" OR "urinary proteome") ',
  'AND ("genome-wide association" OR "GWAS" OR "pQTL" OR "protein quantitative trait") ',
  'AND (FIRST_PDATE:[', since_date, ' TO ', until_date, '])'
)

if (dry_run) {
  message("Dry-run complete. Query:\n  ", epmc_query)
  quit(status = 0L)
}

fetch_epmc_records <- function(query, max_pages = 3L, page_size = 100L) {
  all_recs <- list()
  cursor <- "*"
  for (page in seq_len(max_pages)) {
    url <- paste0(
      "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=",
      utils::URLencode(query, reserved = TRUE),
      "&resultType=core&pageSize=", page_size,
      "&cursorMark=", utils::URLencode(cursor, reserved = TRUE),
      "&format=json"
    )
    res <- tryCatch(
      jsonlite::fromJSON(url, simplifyVector = FALSE),
      error = function(e) {
        warning("Europe PMC API request failed: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(res) || is.null(res$resultList$result) || length(res$resultList$result) == 0) break
    all_recs <- c(all_recs, res$resultList$result)
    next_cursor <- res$nextCursorMark
    if (is.null(next_cursor) || identical(next_cursor, cursor) || length(res$resultList$result) < page_size) break
    cursor <- next_cursor
  }
  all_recs
}

classify_record <- function(title, abstract_text = "") {
  txt <- tolower(paste(title, abstract_text))
  ttl <- tolower(title)

  if (grepl("epigenome-wide|dna methylation|\\bewas\\b", ttl)) {
    return("Auto-Deprioritised: Non-GWAS (EWAS)")
  }
  is_secondary_mr <- grepl(
    "mendelian randomization|mendelian randomisation|proteome-wide association study|drug target|therapeutic target|causal effect|two-sample",
    ttl
  ) && !grepl("we measured|genotyped and measured|newly generated|somascan in [0-9]|olink explore in [0-9]", txt)

  if (is_secondary_mr) {
    return("Auto-Deprioritised: Secondary 2SMR")
  }
  if (grepl("genome-wide association|\\bgwas\\b|\\bpqtl|protein quantitative trait|genetic architecture|genetic determinant", txt)) {
    return("Likely Primary pGWAS")
  }
  "Needs Human Check"
}

detect_platform <- function(text) {
  txt <- tolower(text)
  plats <- character(0)
  if (grepl("somascan|somalogic|aptamer", txt)) plats <- c(plats, "SomaScan")
  if (grepl("olink|proximity extension", txt))  plats <- c(plats, "Olink")
  if (grepl("nulisa", txt))                     plats <- c(plats, "NULISA")
  if (grepl("mass spectrometry|dia-ms|swath|tmt|proteograph", txt)) plats <- c(plats, "Mass Spectrometry")
  if (grepl("luminex", txt))                    plats <- c(plats, "Luminex")
  if (length(plats) == 0) "Unspecified Proteomics" else paste(plats, collapse = " / ")
}

extract_sample_size <- function(text) {
  m <- regmatches(text, regexpr("\\bn\\s*=\\s*[0-9][0-9,]*", tolower(text), perl = TRUE))
  if (length(m) == 1L && nzchar(m)) sub("^n\\s*=\\s*", "", m) else "Check full text"
}

raw_recs <- fetch_epmc_records(epmc_query)
raw_n <- length(raw_recs)
message("Fetched ", raw_n, " records from Europe PMC.")

run_stamp <- format(Sys.time(), "%Y%m%d_%H%M", tz = "UTC")
run_id <- paste0("run_", run_stamp)

new_rows <- list()
seen_in_batch <- character(0)

for (r in raw_recs) {
  doi  <- tolower(clean_cell(r$doi))
  pmid <- clean_cell(r$pmid)
  id_primary <- if (nzchar(doi)) doi else if (nzchar(pmid)) paste0("PMID:", pmid) else clean_cell(r$id)
  if (!nzchar(id_primary)) next

  keys_to_check <- tolower(c(id_primary, doi, pmid, if (nzchar(pmid)) paste0("pmid:", pmid) else ""))
  keys_to_check <- keys_to_check[nzchar(keys_to_check)]
  if (any(keys_to_check %in% known_ids) || any(keys_to_check %in% seen_in_batch)) next
  seen_in_batch <- c(seen_in_batch, keys_to_check)

  title <- clean_cell(r$title)
  abstr <- clean_cell(r$abstractText)
  src   <- clean_cell(r$journalTitle)
  if (!nzchar(src)) src <- if (identical(clean_cell(r$source), "PPR")) "Preprint (bioRxiv/medRxiv)" else clean_cell(r$source)

  cls  <- classify_record(title, abstr)
  plat <- detect_platform(paste(title, abstr))
  n_est <- extract_sample_size(paste(title, abstr))

  new_rows[[length(new_rows) + 1L]] <- data.frame(
    candidate_id = sprintf("cand_%s_%03d", gsub("[^0-9]", "", until_date), length(new_rows) + 1L),
    run_id = run_id,
    source = src,
    identifier = id_primary,
    detected_platform = plat,
    detected_n = n_est,
    title = title,
    auto_classification = cls,
    curation_status = ifelse(grepl("^Auto-Deprioritised", cls), "Auto-Triage: Low Priority", "Pending Curator Review"),
    stringsAsFactors = FALSE
  )
}

dedup_n <- length(new_rows)
auto_mr_n <- if (dedup_n > 0) sum(vapply(new_rows, function(x) grepl("^Auto-Deprioritised", x$auto_classification), logical(1))) else 0L
review_n  <- dedup_n - auto_mr_n

if (dedup_n > 0) {
  new_df <- do.call(rbind, new_rows)
  # Keep actionable items at the top of candidate_queue.tsv, cap low-priority tail so file stays crisp
  actionable <- new_df[!grepl("^Auto-Deprioritised", new_df$auto_classification), , drop = FALSE]
  deprior    <- head(new_df[grepl("^Auto-Deprioritised", new_df$auto_classification), , drop = FALSE], 15L)
  updated_queue <- rbind(actionable, deprior, queue)
  write_tsv_safe(updated_queue, file.path("catalogue", "candidate_queue.tsv"))
  message(sprintf("Added %d actionable + %d sample deprioritised candidates to catalogue/candidate_queue.tsv",
                  nrow(actionable), nrow(deprior)))
} else {
  message("0 unseen records found in window.")
}

# Log this execution in catalogue/search_runs.tsv
new_run_row <- data.frame(
  run_id = run_id,
  run_date = until_date,
  trigger = clean_cell(trigger_lb),
  sources_queried = paste0("Europe PMC (PubMed + Preprints) [", since_date, "..", until_date, "]"),
  raw_records = as.integer(raw_n),
  deduplicated_new = as.integer(dedup_n),
  auto_flagged_secondary_mr = as.integer(auto_mr_n),
  candidates_for_review = as.integer(review_n),
  studies_added = 0L,
  preprints_upgraded = 0L,
  pr_number = "Automated PR",
  release_tag = "",
  status = if (review_n > 0L) sprintf("Pending Triage (%d candidates)", review_n) else "Completed (0 new primary)",
  stringsAsFactors = FALSE
)
runs_updated <- rbind(runs, new_run_row)
write_tsv_safe(runs_updated, file.path("catalogue", "search_runs.tsv"))

val_status <- system2("Rscript", "src/validate_catalogue.R")
if (val_status != 0L) stop("Validation failed after search update!", call. = FALSE)
