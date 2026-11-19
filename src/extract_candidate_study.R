#!/usr/bin/env Rscript
# 4-Tier Cascading Paper, Supplement, Repository API & Provenance Extraction CLI
#
# Usage examples:
#   Rscript src/extract_candidate_study.R --pmid 37794186 --cache-dir <dir> --out-dir <dir>
#   Rscript src/extract_candidate_study.R --doi 10.1101/2022.06.17.496443 --epmc-id PPR508031 --out-dir <dir>
#   Rscript src/extract_candidate_study.R --queue-file catalogue/candidate_queue.tsv --candidate-id cand_2026_11_01
#   Rscript src/extract_candidate_study.R --queue-file catalogue/candidate_queue.tsv --batch-pending --max-extract 10

req_pkgs <- c("httr2", "jsonlite", "xml2", "readxl")
missing_pkgs <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}

# Locate script directory so submodules load cleanly regardless of caller CWD
cmd_args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args_all, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
} else {
  file.path(getwd(), "src")
}
repo_root_default <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

source(file.path(script_dir, "extract", "http_client.R"))
source(file.path(script_dir, "extract", "tier1_apis.R"))
source(file.path(script_dir, "extract", "tier2_fulltext.R"))
source(file.path(script_dir, "extract", "tier3_build_resolver.R"))
source(file.path(script_dir, "extract", "tier3_supplements.R"))
source(file.path(script_dir, "extract", "tier3_repositories.R"))
source(file.path(script_dir, "extract", "tier4_provenance_escalate.R"))

args <- commandArgs(trailingOnly = TRUE)
get_opt <- function(flag, default = "") {
  idx <- match(flag, args)
  if (is.na(idx) || idx >= length(args)) default else args[[idx + 1L]]
}
has_flag <- function(flag) flag %in% args

if (has_flag("--help") || has_flag("-h") || length(args) == 0L) {
  cat(paste(
    "Usage: Rscript src/extract_candidate_study.R [options]",
    "  --pmid <PMID>                 PubMed ID",
    "  --doi <DOI>                   DOI (published or preprint)",
    "  --pmcid <PMCID>               PubMed Central ID (e.g. PMC10567551)",
    "  --epmc-id <PPR_ID>            Europe PMC preprint ID (e.g. PPR508031)",
    "  --identifier <ID>             Generic DOI / PMID / PMCID identifier",
    "  --queue-file <TSV>            Path to candidate_queue.tsv (9-col or 17-col)",
    "  --candidate-id <CAND_ID>      Extract a specific candidate_id from --queue-file",
    "  --batch-pending               Extract pending actionable rows from --queue-file",
    "  --max-extract <INT>           Cap batch extraction count (default: 10)",
    "  --cache-dir <DIR>             Raw HTTP download cache outside Git/OneDrive",
    "  --out-dir <DIR>               Output directory for draft TSVs, provenance & checklist",
    "  --local-pdf <PDF>             Optional local curator PDF override",
    "  --local-xml <XML>             Optional local curator JATS XML override",
    "  --refresh-cache               Bypass existing cached HTTP responses",
    sep = "\n"
  ), "\n")
  quit(status = 0L)
}

xdg_cache <- Sys.getenv("XDG_CACHE_HOME", file.path(Sys.getenv("HOME", tempdir()), ".cache"))
cache_dir <- get_opt("--cache-dir", file.path(xdg_cache, "proteomic-gwas-catalogue", "extractions"))
out_dir   <- get_opt("--out-dir", file.path(repo_root_default, "staging"))
repo_root <- get_opt("--repo-root", repo_root_default)
manifests_dir <- get_opt("--manifests-dir", normalizePath(file.path(repo_root, "..", "GWASprocessing", "manifests"), mustWork = FALSE))
local_pdf <- get_opt("--local-pdf", "")
local_xml <- get_opt("--local-xml", "")
refresh_cache <- has_flag("--refresh-cache")
max_extract <- suppressWarnings(as.integer(get_opt("--max-extract", "10")))
if (is.na(max_extract) || max_extract < 1L) max_extract <- 10L

studies_path <- file.path(repo_root, "catalogue", "studies.tsv")
existing_studies <- if (file.exists(studies_path)) {
  utils::read.delim(studies_path, stringsAsFactors = FALSE, check.names = FALSE, quote = "")
} else {
  NULL
}

# Build list of candidate records to process ([C-1] supports 9-col and 17-col candidate_queue.tsv)
candidates <- list()
queue_file <- get_opt("--queue-file", "")
cand_id_opt <- get_opt("--candidate-id", "")

if (nzchar(queue_file)) {
  if (!file.exists(queue_file)) stop("Queue file not found: ", queue_file, call. = FALSE)
  qdf <- utils::read.delim(queue_file, stringsAsFactors = FALSE, check.names = FALSE, quote = "")
  row_to_cand <- function(r) {
    list(
      candidate_id = clean_str(r$candidate_id %||% ""),
      identifier = clean_str(r$identifier %||% ""),
      pmid = clean_str(if ("pmid" %in% names(r)) r$pmid else ""),
      doi = clean_str(if ("doi" %in% names(r)) r$doi else ""),
      pmcid = clean_str(if ("pmcid" %in% names(r)) r$pmcid else ""),
      epmc_id = clean_str(if ("epmc_id" %in% names(r)) r$epmc_id else ""),
      preprint_lineage_note = clean_str(if ("preprint_lineage_note" %in% names(r)) r$preprint_lineage_note else ""),
      title = clean_str(r$title %||% "")
    )
  }

  if (nzchar(cand_id_opt)) {
    hit <- which(qdf$candidate_id == cand_id_opt)
    if (length(hit) == 0L) stop("candidate_id not found in queue: ", cand_id_opt, call. = FALSE)
    candidates[[1L]] <- row_to_cand(qdf[hit[[1L]], , drop = FALSE])
  } else if (has_flag("--batch-pending")) {
    keep <- grepl("^Pending", qdf$curation_status %||% "") & !grepl("^Auto-Deprioritised", qdf$auto_classification %||% "")
    sub_q <- qdf[keep, , drop = FALSE]
    if ("p_score" %in% names(sub_q)) {
      ord <- order(
        !grepl("Preprint -> Published", sub_q$auto_classification %||% ""),
        -suppressWarnings(as.numeric(sub_q$p_score)),
        na.last = TRUE
      )
      sub_q <- sub_q[ord, , drop = FALSE]
    }
    sub_q <- utils::head(sub_q, max_extract)
    for (i in seq_len(nrow(sub_q))) candidates[[i]] <- row_to_cand(sub_q[i, , drop = FALSE])
  } else {
    stop("When passing --queue-file, specify either --candidate-id <ID> or --batch-pending.", call. = FALSE)
  }
} else {
  candidates[[1L]] <- list(
    candidate_id = "",
    identifier = get_opt("--identifier", ""),
    pmid = get_opt("--pmid", ""),
    doi = get_opt("--doi", ""),
    pmcid = get_opt("--pmcid", ""),
    epmc_id = get_opt("--epmc-id", ""),
    preprint_lineage_note = get_opt("--lineage-note", ""),
    title = get_opt("--title", "")
  )
}

session <- new_http_session(cache_dir = cache_dir, refresh_cache = refresh_cache, strict_oa_only = TRUE)

for (cand in candidates) {
  t0 <- proc.time()[["elapsed"]]
  t1 <- run_tier1_apis(cand, session, existing_studies = existing_studies)
  t2 <- run_tier2_fulltext(t1, session, local_pdf = local_pdf, local_xml = local_xml)
  t3_supp <- run_tier3_supplements(t1, t2, session, repo_root = repo_root)
  t3_repo <- run_tier3_repositories(t1, t2, session)
  t4 <- run_tier4_provenance_and_escalate(
    t1, t2, t3_supp, t3_repo,
    out_dir = out_dir,
    existing_studies = existing_studies,
    manifests_dir = manifests_dir
  )
  http_log_path <- file.path(out_dir, paste0(t4$study_id, "_http_log.tsv"))
  export_http_log(session, http_log_path)
  elapsed <- proc.time()[["elapsed"]] - t0

  message(sprintf(
    "[%s] decision=%s | queue_status=\"%s\" | strata=%d | prov_rows=%d | unresolved_req=%d | elapsed=%.2fs\n  -> %s",
    t4$study_id, t4$decision, t4$recommended_curation_status,
    t4$n_strata, t4$n_provenance_rows, t4$n_blocking_fields, elapsed,
    t4$files$checklist
  ))
}
