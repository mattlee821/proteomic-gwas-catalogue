#!/usr/bin/env Rscript
# Detect and apply preprint-to-published replacements in catalogue/studies.tsv,
# and record explicit study ID replacements ("preprint_id -> published_id") in catalogue/search_runs.tsv.
#
# Usage:
#   # 1. Apply an explicit replacement pair:
#   Rscript src/apply_preprint_upgrades.R --preprint-id sun_2022_preprint --published-id sun_2023_37794186 [--run-id run_2026_07_01]
#
#   # 2. Apply all *_preprint_supersede_patch.tsv files emitted by src/extract_candidate_study.R in staging/:
#   Rscript src/apply_preprint_upgrades.R --staging-dir staging [--run-id <run_id>]
#
#   # 3. Check active Preprints in catalogue/studies.tsv against bioRxiv/EuropePMC/Crossref APIs (--check-live):
#   Rscript src/apply_preprint_upgrades.R --check-live

args <- commandArgs(trailingOnly = TRUE)
get_opt <- function(flag, default = "") {
  idx <- match(flag, args)
  if (is.na(idx) || idx >= length(args)) default else args[[idx + 1L]]
}
has_flag <- function(flag) flag %in% args
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

studies_path <- file.path("catalogue", "studies.tsv")
strata_path  <- file.path("catalogue", "strata.tsv")
runs_path    <- file.path("catalogue", "search_runs.tsv")

read_tsv <- function(p) utils::read.delim(p, stringsAsFactors = FALSE, check.names = FALSE, quote = "")
write_tsv <- function(df, p) utils::write.table(df, file = p, sep = "\t", row.names = FALSE, quote = FALSE, na = "")

studies <- read_tsv(studies_path)
runs    <- read_tsv(runs_path)

if (!("preprints_replaced_details" %in% names(runs))) {
  runs$preprints_replaced_details <- ""
}

pre_arg <- get_opt("--preprint-id", "")
pub_arg <- get_opt("--published-id", "")
stg_dir <- get_opt("--staging-dir", "")
run_arg <- get_opt("--run-id", if (nrow(runs) > 0L) runs$run_id[nrow(runs)] else "")

pairs <- list()

if (nzchar(pre_arg) && nzchar(pub_arg)) {
  pairs[[length(pairs) + 1L]] <- list(pre = pre_arg, pub = pub_arg)
}

if (nzchar(stg_dir) && dir.exists(stg_dir)) {
  patch_files <- list.files(stg_dir, pattern = "_preprint_supersede_patch\\.tsv$", full.names = TRUE)
  for (pf in patch_files) {
    pdf <- read_tsv(pf)
    for (i in seq_len(nrow(pdf))) {
      pre_id <- pdf$preprint_study_id[i]
      pub_id <- pdf$new_superseded_by[i]
      if (nzchar(pre_id) && !grepl("^UNTRACKED:", pre_id) && nzchar(pub_id)) {
        pairs[[length(pairs) + 1L]] <- list(pre = pre_id, pub = pub_id)
      }
    }
  }
}

if (has_flag("--check-live") && requireNamespace("httr2", quietly = TRUE) && requireNamespace("jsonlite", quietly = TRUE)) {
  active_pre <- studies[studies$status == "Preprint" & nzchar(studies$doi), , drop = FALSE]
  message(sprintf("Checking %d active preprint(s) in catalogue/studies.tsv for published journal versions...", nrow(active_pre)))
  for (i in seq_len(nrow(active_pre))) {
    p_doi <- tolower(trimws(active_pre$doi[i]))
    p_id  <- active_pre$study_id[i]
    pub_doi <- ""
    if (grepl("^10\\.(1101|64898)/", p_doi)) {
      for (srv in c("biorxiv", "medrxiv")) {
        u <- sprintf("https://api.biorxiv.org/pubs/%s/%s/na/json", srv, p_doi)
        resp <- tryCatch(httr2::request(u) |> httr2::req_timeout(15) |> httr2::req_error(is_error = function(r) FALSE) |> httr2::req_perform(), error = function(e) NULL)
        if (!is.null(resp) && httr2::resp_status(resp) == 200L) {
          js <- tryCatch(httr2::resp_body_json(resp, simplifyVector = FALSE), error = function(e) NULL)
          if (!is.null(js$collection) && length(js$collection) > 0L && nzchar(js$collection[[1]]$published_doi %||% "")) {
            pub_doi <- tolower(trimws(js$collection[[1]]$published_doi))
            break
          }
        }
      }
    }
    if (nzchar(pub_doi)) {
      hit_idx <- which(tolower(trimws(studies$doi)) == pub_doi)
      if (length(hit_idx) > 0L) {
        pairs[[length(pairs) + 1L]] <- list(pre = p_id, pub = studies$study_id[hit_idx[1L]])
      } else {
        message(sprintf("  [UPGRADE DETECTED] %s (%s) -> published DOI %s (run extract_candidate_study.R --doi %s to add published record)", p_id, p_doi, pub_doi, pub_doi))
      }
    }
  }
}

applied <- character(0)
for (pr in pairs) {
  i_pre <- which(studies$study_id == pr$pre)
  i_pub <- which(studies$study_id == pr$pub)
  if (length(i_pre) == 1L && length(i_pub) == 1L) {
    studies$status[i_pre] <- "Superseded Preprint"
    studies$superseded_by[i_pre] <- pr$pub
    applied <- c(applied, sprintf("%s -> %s", pr$pre, pr$pub))
  }
}

if (length(applied) > 0L) {
  applied <- unique(applied)
  write_tsv(studies, studies_path)
  if (nzchar(run_arg) && any(runs$run_id == run_arg)) {
    r_idx <- which(runs$run_id == run_arg)[1L]
    existing_det <- trimws(runs$preprints_replaced_details[r_idx])
    merged_det <- unique(c(if (nzchar(existing_det)) trimws(strsplit(existing_det, ";")[[1L]]) else character(0), applied))
    runs$preprints_replaced_details[r_idx] <- paste(merged_det, collapse = "; ")
    runs$preprints_upgraded[r_idx] <- length(merged_det)
    write_tsv(runs, runs_path)
  }
  message("Applied preprint->published replacements: ", paste(applied, collapse = "; "))
} else {
  message("No pending preprint->published study replacements applied.")
}
