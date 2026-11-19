#!/usr/bin/env Rscript
# Multi-resource literature discovery, connected-components deduplication,
# 4-tier preprint-to-published upgrade detector, and Two-Stage Dual-Axis Triage
# engine (v2) for proteomic-gwas-catalogue.
#
# Usage:
#   Rscript src/run_literature_search.R [--since YYYY-MM-DD] [--until YYYY-MM-DD]
#     [--mode incremental|retrospective] [--trigger LABEL] [--dry-run]
#     [--verify-pmids 37794186,34857953] [--out-dir DIR] [--write-catalogue]
#     [--include-gwas-bulk-tsv]

for (pkg in c("jsonlite", "xml2", "httr2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required R package '%s' is not installed.", pkg), call. = FALSE)
  }
}

# Locate and source shared canonical lexicons & triage scorer (src/search/lexicons.R)
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

# ---------------------------------------------------------------------------
# CLI Argument Parsing
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_opt <- function(flag, default = NA_character_) {
  idx <- match(flag, args)
  if (is.na(idx) || idx >= length(args)) default else args[idx + 1L]
}

dry_run       <- "--dry-run" %in% args
mode          <- tolower(get_opt("--mode", "incremental"))
if (!(mode %in% c("incremental", "retrospective"))) {
  stop("--mode must be 'incremental' or 'retrospective'.", call. = FALSE)
}

default_since <- if (identical(mode, "retrospective")) "2005-01-01" else format(Sys.Date() - 35L, "%Y-%m-%d")
since_date    <- get_opt("--since", default_since)
until_date    <- get_opt("--until", format(Sys.Date(), "%Y-%m-%d"))
trigger_lb    <- get_opt("--trigger", "Scheduled GitHub Action")
verify_raw    <- get_opt("--verify-pmids", "")
out_dir       <- get_opt("--out-dir", "")
use_gwas_bulk <- "--include-gwas-bulk-tsv" %in% args

verify_pmids <- if (nzchar(verify_raw)) extract_pmids_from_cell(verify_raw, strip_dois_first = FALSE) else character(0)
is_verify    <- length(verify_pmids) > 0L

write_catalogue <- if ("--write-catalogue" %in% args) {
  TRUE
} else if (is_verify || nzchar(out_dir) || dry_run) {
  FALSE
} else {
  TRUE
}

clean_cell <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[1L])) return("")
  gsub("[\r\n\t]+", " ", trimws(as.character(x[1L])))
}

read_tsv_safe <- function(path) {
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, quote = "", colClasses = "character")
}

write_tsv_safe <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(df, file = path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

UA_STRING <- "proteomic-gwas-catalogue/2.0 (https://github.com/matthewlee/proteomic-gwas-catalogue)"

# Compile queries from single source of truth (src/search/lexicons.R)
pm_date_field <- if (identical(mode, "retrospective")) "dp" else "edat"
q_pubmed      <- build_pubmed_query(since_date, until_date, date_field = pm_date_field,
                                    verify_pmids = if (is_verify) verify_pmids else NULL)
q_epmc_2a     <- build_epmc_prong2a_query(since_date, until_date,
                                          include_idate = identical(mode, "incremental"),
                                          verify_pmids = if (is_verify) verify_pmids else NULL)
q_epmc_2b     <- build_epmc_prong2b_query(since_date, until_date,
                                          include_idate = identical(mode, "incremental"),
                                          verify_pmids = if (is_verify) verify_pmids else NULL)
q_epmc_snow   <- build_epmc_snowball_query(LANDMARK_PGWAS_SEEDS$pmid, since_date, until_date,
                                           include_idate = identical(mode, "incremental"),
                                           verify_pmids = if (is_verify) verify_pmids else NULL)

message(sprintf(
  "Proteomic GWAS Search v2 | mode=%s | window=%s..%s | verify_pmids=%s | dry_run=%s | write_catalogue=%s",
  mode, since_date, until_date,
  if (is_verify) paste(verify_pmids, collapse = ",") else "none",
  dry_run, write_catalogue
))

if (dry_run) {
  cat("\n--- [Source 1: PubMed E-utilities POST Query] ---\n", q_pubmed, "\n", sep = "")
  cat("\n--- [Source 2A: Europe PMC Metadata + Preprint POST Query] ---\n", q_epmc_2a, "\n", sep = "")
  cat("\n--- [Source 2B: Europe PMC Full-Text METHODS Rescue Query] ---\n", q_epmc_2b, "\n", sep = "")
  cat("\n--- [Source 4: Europe PMC Server-Side Snowball Query] ---\n", q_epmc_snow, "\n", sep = "")
  cat("\n--- [Source 3: NHGRI-EBI GWAS Catalog REST API v2] ---\n",
      "GET https://www.ebi.ac.uk/gwas/rest/api/v2/studies?efo_id=EFO_0004747&show_child_traits=true&page=0&size=200\n", sep = "")
  quit(status = 0L)
}

# ---------------------------------------------------------------------------
# Load Existing Catalogue State (Read-Only unless write_catalogue == TRUE)
# ---------------------------------------------------------------------------
studies  <- read_tsv_safe(file.path("catalogue", "studies.tsv"))
excluded <- read_tsv_safe(file.path("catalogue", "excluded_records.tsv"))
queue    <- read_tsv_safe(file.path("catalogue", "candidate_queue.tsv"))
runs     <- read_tsv_safe(file.path("catalogue", "search_runs.tsv"))

# Standardized intermediate candidate record constructor
make_raw_rec <- function(source_tag,
                         pmid = "",
                         doi = "",
                         pmcid = "",
                         epmc_id = "",
                         journal = "",
                         pub_year = "",
                         first_author = "",
                         last_author = "",
                         title = "",
                         abstract = "",
                         pub_types = character(0),
                         mesh_headings = character(0),
                         preprint_links = character(0),
                         erratum_targets = character(0),
                         gwas_study_count = 0L) {
  list(
    source_tag       = clean_cell(source_tag),
    pmid             = normalize_pmid_one(pmid),
    doi              = normalize_doi_one(doi),
    pmcid            = normalize_pmcid_one(pmcid),
    epmc_id          = clean_cell(epmc_id),
    journal          = clean_cell(journal),
    pub_year         = clean_cell(pub_year),
    first_author     = clean_cell(first_author),
    last_author      = clean_cell(last_author),
    title            = clean_cell(title),
    abstract         = clean_cell(abstract),
    pub_types        = unique(trimws(as.character(pub_types[nzchar(pub_types)]))),
    mesh_headings    = unique(trimws(as.character(mesh_headings[nzchar(mesh_headings)]))),
    preprint_links   = unique(trimws(as.character(preprint_links[nzchar(preprint_links)]))),
    erratum_targets  = unique(trimws(as.character(erratum_targets[nzchar(erratum_targets)]))),
    gwas_study_count = as.integer(gwas_study_count %||% 0L)
  )
}

# ---------------------------------------------------------------------------
# Source 1: PubMed via NCBI E-utilities (esearch POST + efetch XML POST)
# ---------------------------------------------------------------------------
fetch_pubmed_efetch_xml <- function(pmids, source_tag = "PubMed") {
  pmids <- unique(pmids[nzchar(pmids)])
  if (length(pmids) == 0L) return(list())
  ncbi_key <- Sys.getenv("NCBI_API_KEY", unset = "")
  sleep_s  <- if (nzchar(ncbi_key)) 0.15 else 0.40
  out <- list()

  chunks <- split(pmids, ceiling(seq_along(pmids) / 200L))
  for (ch in chunks) {
    Sys.sleep(sleep_s)
    form_params <- list(
      db      = "pubmed",
      id      = paste(ch, collapse = ","),
      retmode = "xml",
      tool    = "proteomic_gwas_catalogue"
    )
    if (nzchar(ncbi_key)) form_params$api_key <- ncbi_key

    resp <- tryCatch(
      httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi") |>
        httr2::req_user_agent(UA_STRING) |>
        httr2::req_timeout(45) |>
        httr2::req_body_form(!!!form_params) |>
        httr2::req_error(is_error = function(r) FALSE) |>
        httr2::req_perform(),
      error = function(e) {
        warning("PubMed efetch failed: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(resp) || httr2::resp_status(resp) != 200L) next
    doc <- tryCatch(xml2::read_xml(httr2::resp_body_string(resp)), error = function(e) NULL)
    if (is.null(doc)) next

    arts <- xml2::xml_find_all(doc, ".//PubmedArticle")
    for (art in arts) {
      pmid_v <- xml2::xml_text(xml2::xml_find_first(art, ".//MedlineCitation/PMID"))
      ttl_v  <- xml2::xml_text(xml2::xml_find_first(art, ".//Article/ArticleTitle"))
      jrn_v  <- xml2::xml_text(xml2::xml_find_first(art, ".//Article/Journal/Title"))
      yr_v   <- xml2::xml_text(xml2::xml_find_first(art, ".//Article/Journal/JournalIssue/PubDate/Year"))
      if (is.na(yr_v) || !nzchar(yr_v)) {
        md <- xml2::xml_text(xml2::xml_find_first(art, ".//Article/Journal/JournalIssue/PubDate/MedlineDate"))
        yr_v <- regmatches(md %||% "", regexpr("[0-9]{4}", md %||% ""))
      }

      abs_nodes <- xml2::xml_find_all(art, ".//Article/Abstract/AbstractText")
      abs_parts <- vapply(abs_nodes, function(n) {
        lbl <- xml2::xml_attr(n, "Label")
        tx  <- xml2::xml_text(n)
        if (!is.na(lbl) && nzchar(lbl)) paste0(lbl, ": ", tx) else tx
      }, character(1))
      abs_v <- paste(abs_parts, collapse = " ")

      au_surnames <- xml2::xml_text(xml2::xml_find_all(art, ".//AuthorList/Author/LastName"))
      au_collect  <- xml2::xml_text(xml2::xml_find_all(art, ".//AuthorList/Author/CollectiveName"))
      first_au <- if (length(au_surnames) > 0L) au_surnames[1L] else if (length(au_collect) > 0L) au_collect[1L] else ""
      last_au  <- if (length(au_surnames) > 0L) au_surnames[length(au_surnames)] else first_au

      id_nodes <- xml2::xml_find_all(art, ".//PubmedData/ArticleIdList/ArticleId")
      doi_v <- ""; pmc_v <- ""
      for (idn in id_nodes) {
        idt <- xml2::xml_attr(idn, "IdType") %||% ""
        idv <- xml2::xml_text(idn)
        if (identical(idt, "doi") && !nzchar(doi_v)) doi_v <- idv
        if (identical(idt, "pmc") && !nzchar(pmc_v)) pmc_v <- idv
      }

      pts_v  <- xml2::xml_text(xml2::xml_find_all(art, ".//PublicationTypeList/PublicationType"))
      mesh_v <- xml2::xml_text(xml2::xml_find_all(art, ".//MeshHeadingList/MeshHeading/DescriptorName"))

      cc_nodes <- xml2::xml_find_all(art, ".//CommentsCorrectionsList/CommentsCorrections")
      pre_links <- character(0); err_links <- character(0)
      for (cc in cc_nodes) {
        rt <- xml2::xml_attr(cc, "RefType") %||% ""
        rs <- xml2::xml_text(xml2::xml_find_first(cc, "./RefSource")) %||% ""
        rp <- xml2::xml_text(xml2::xml_find_first(cc, "./PMID")) %||% ""
        if (rt %in% c("UpdateOf", "PreprintIn", "PreprintOf")) {
          pre_links <- c(pre_links, extract_dois_from_cell(rs), if (!is.na(rp) && nzchar(rp)) paste0("PMID:", rp) else "")
        }
        if (rt %in% c("ErratumFor", "CorrectionFor")) {
          err_links <- c(err_links, extract_dois_from_cell(rs), if (!is.na(rp) && nzchar(rp)) rp else "")
        }
      }

      out[[length(out) + 1L]] <- make_raw_rec(
        source_tag = source_tag, pmid = pmid_v, doi = doi_v, pmcid = pmc_v,
        journal = jrn_v, pub_year = yr_v, first_author = first_au, last_author = last_au,
        title = ttl_v, abstract = abs_v, pub_types = pts_v, mesh_headings = mesh_v,
        preprint_links = pre_links, erratum_targets = err_links
      )
    }
  }
  out
}

query_pubmed_v2 <- function(query_str) {
  ncbi_key <- Sys.getenv("NCBI_API_KEY", unset = "")
  sleep_s  <- if (nzchar(ncbi_key)) 0.15 else 0.40
  all_pmids <- character(0)
  retstart <- 0L
  retmax   <- 500L

  repeat {
    Sys.sleep(sleep_s)
    form_params <- list(
      db       = "pubmed",
      term     = query_str,
      retmode  = "json",
      retstart = as.character(retstart),
      retmax   = as.character(retmax),
      tool     = "proteomic_gwas_catalogue"
    )
    if (nzchar(ncbi_key)) form_params$api_key <- ncbi_key

    resp <- tryCatch(
      httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi") |>
        httr2::req_user_agent(UA_STRING) |>
        httr2::req_timeout(45) |>
        httr2::req_body_form(!!!form_params) |>
        httr2::req_error(is_error = function(r) FALSE) |>
        httr2::req_perform(),
      error = function(e) {
        warning("PubMed esearch failed: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(resp) || httr2::resp_status(resp) != 200L) break
    js <- tryCatch(httr2::resp_body_json(resp, simplifyVector = TRUE), error = function(e) NULL)
    if (is.null(js) || is.null(js$esearchresult)) break

    ids <- as.character(js$esearchresult$idlist %||% character(0))
    tot <- suppressWarnings(as.integer(js$esearchresult$count %||% "0"))
    all_pmids <- unique(c(all_pmids, ids))
    retstart  <- retstart + retmax
    if (length(ids) == 0L || retstart >= tot) break
  }

  message(sprintf("[Source 1: PubMed] esearch matched %d PMID(s).", length(all_pmids)))
  fetch_pubmed_efetch_xml(all_pmids, source_tag = "PubMed")
}

# ---------------------------------------------------------------------------
# Source 2: Europe PMC REST API (searchPOST across MED + PMC + PPR)
# ---------------------------------------------------------------------------
parse_epmc_result_item <- function(r, default_tag = "EuropePMC") {
  src_db <- clean_cell(r$source)
  tag <- if (identical(default_tag, "EuropePMC")) {
    if (identical(src_db, "PPR")) "EuropePMC_PPR" else "EuropePMC_MED"
  } else {
    default_tag
  }

  epmc_id <- if (identical(src_db, "PPR")) clean_cell(r$id) else ""
  pre_links <- character(0)
  err_targets <- character(0)

  cclist <- r$commentCorrectionList$commentCorrection
  if (!is.null(cclist) && length(cclist) > 0L) {
    for (cc in cclist) {
      cct <- clean_cell(cc$type)
      cci <- clean_cell(cc$id)
      ccs <- clean_cell(cc$source)
      ccr <- clean_cell(cc$reference)
      if (cct %in% c("Preprint in", "Preprint of")) {
        if (identical(ccs, "PPR") && nzchar(cci)) {
          if (!nzchar(epmc_id)) epmc_id <- cci
          pre_links <- c(pre_links, cci)
        } else if (identical(ccs, "MED") && nzchar(cci)) {
          pre_links <- c(pre_links, paste0("PMID:", cci))
        }
        pre_links <- c(pre_links, extract_dois_from_cell(ccr))
      }
      if (cct %in% c("Erratum for", "Correction for")) {
        if (nzchar(cci)) err_targets <- c(err_targets, cci)
        err_targets <- c(err_targets, extract_dois_from_cell(ccr))
      }
    }
  }

  pts <- character(0)
  if (!is.null(r$pubTypeList$pubType)) {
    pts <- unlist(r$pubTypeList$pubType)
  }
  meshes <- character(0)
  if (!is.null(r$meshHeadingList$meshHeading)) {
    meshes <- vapply(r$meshHeadingList$meshHeading, function(m) clean_cell(m$descriptorName), character(1))
  }

  au_str <- clean_cell(r$authorString)
  au_parts <- trimws(unlist(strsplit(au_str, ",")))
  first_au <- if (length(au_parts) > 0L) au_parts[1L] else ""
  last_au  <- if (length(au_parts) > 0L) au_parts[length(au_parts)] else first_au

  jrn <- clean_cell(r$journalTitle)
  if (!nzchar(jrn)) {
    jrn <- if (identical(src_db, "PPR")) clean_cell(r$bookOrReportDetails$publisher %||% "Preprint") else src_db
  }

  make_raw_rec(
    source_tag      = tag,
    pmid            = r$pmid %||% "",
    doi             = r$doi %||% "",
    pmcid           = r$pmcid %||% "",
    epmc_id         = epmc_id,
    journal         = jrn,
    pub_year        = r$pubYear %||% "",
    first_author    = first_au,
    last_author     = last_au,
    title           = r$title %||% "",
    abstract        = r$abstractText %||% "",
    pub_types       = pts,
    mesh_headings   = meshes,
    preprint_links  = pre_links,
    erratum_targets = err_targets
  )
}

fetch_epmc_post <- function(query_str, source_tag = "EuropePMC", page_size = 200L, max_pages = 25L) {
  recs <- list()
  cursor <- "*"
  for (pg in seq_len(max_pages)) {
    Sys.sleep(0.35)
    resp <- tryCatch(
      httr2::request("https://www.ebi.ac.uk/europepmc/webservices/rest/searchPOST") |>
        httr2::req_user_agent(UA_STRING) |>
        httr2::req_timeout(45) |>
        httr2::req_headers(`Content-Type` = "application/x-www-form-urlencoded") |>
        httr2::req_body_form(
          query      = query_str,
          resultType = "core",
          pageSize   = as.character(page_size),
          cursorMark = cursor,
          format     = "json"
        ) |>
        httr2::req_error(is_error = function(r) FALSE) |>
        httr2::req_perform(),
      error = function(e) {
        warning("Europe PMC searchPOST failed: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(resp) || httr2::resp_status(resp) != 200L) break
    js <- tryCatch(httr2::resp_body_json(resp, simplifyVector = FALSE), error = function(e) NULL)
    items <- js$resultList$result
    if (is.null(items) || length(items) == 0L) break

    parsed <- lapply(items, parse_epmc_result_item, default_tag = source_tag)
    recs <- c(recs, parsed)

    next_cursor <- js$nextCursorMark
    if (is.null(next_cursor) || identical(next_cursor, cursor) || length(items) < page_size) break
    cursor <- next_cursor
  }
  recs
}

run_epmc_windowed <- function(builder_fn, source_tag, since_d, until_d, is_retro, verify_ids = NULL) {
  if (!is.null(verify_ids) && length(verify_ids) > 0L) {
    q <- builder_fn(since_d, until_d, include_idate = FALSE, verify_pmids = verify_ids)
    return(fetch_epmc_post(q, source_tag = source_tag, page_size = 100L, max_pages = 2L))
  }
  # M-4 Fix: Year-chunk retrospective sweeps with pageSize = 1000 to prevent deep-cursor expiration
  if (isTRUE(is_retro)) {
    y_start <- as.integer(substr(since_d, 1L, 4L))
    y_end   <- as.integer(substr(until_d, 1L, 4L))
    all_out <- list()
    for (yr in seq(y_start, y_end)) {
      w_from <- if (yr == y_start) since_d else sprintf("%04d-01-01", yr)
      w_to   <- if (yr == y_end)   until_d else sprintf("%04d-12-31", yr)
      q_yr   <- builder_fn(w_from, w_to, include_idate = FALSE, verify_pmids = NULL)
      yr_res <- fetch_epmc_post(q_yr, source_tag = source_tag, page_size = 1000L, max_pages = 15L)
      all_out <- c(all_out, yr_res)
    }
    return(all_out)
  }
  q_inc <- builder_fn(since_d, until_d, include_idate = TRUE, verify_pmids = NULL)
  fetch_epmc_post(q_inc, source_tag = source_tag, page_size = 500L, max_pages = 10L)
}

# Source 2b: Narrow 14/21-day bioRxiv/medRxiv category catch-up scan (ONLY when window <= 21 days)
query_biorxiv_narrow_v2 <- function(since_d, until_d) {
  win_days <- as.numeric(as.Date(until_d) - as.Date(since_d))
  if (is.na(win_days) || win_days > 21 || win_days < 0) return(list())

  cats <- list(
    biorxiv = c("genetics", "genomics", "bioinformatics", "systems_biology"),
    medrxiv = c("genetic_and_genomic_medicine", "epidemiology")
  )
  out <- list()
  for (srv in names(cats)) {
    for (cat_name in cats[[srv]]) {
      cursor <- 0L
      for (pg in seq_len(5L)) {
        Sys.sleep(0.40)
        url <- sprintf("https://api.biorxiv.org/details/%s/%s/%s/%d?category=%s",
                       srv, since_d, until_d, cursor, cat_name)
        resp <- tryCatch(
          httr2::request(url) |>
            httr2::req_user_agent(UA_STRING) |>
            httr2::req_timeout(25) |>
            httr2::req_error(is_error = function(r) FALSE) |>
            httr2::req_perform(),
          error = function(e) NULL
        )
        if (is.null(resp) || httr2::resp_status(resp) != 200L) break
        js <- tryCatch(httr2::resp_body_json(resp, simplifyVector = FALSE), error = function(e) NULL)
        coll <- js$collection
        if (is.null(coll) || length(coll) == 0L) break

        for (item in coll) {
          ttl <- clean_cell(item$title)
          abs <- clean_cell(item$abstract)
          comb <- tolower(paste(ttl, abs))
          if (grepl(BROAD_PROTEIN_REGEX, comb, perl = TRUE) && grepl(BROAD_GENETICS_REGEX, comb, perl = TRUE)) {
            au_list <- trimws(unlist(strsplit(clean_cell(item$authors), ";")))
            pub_doi <- normalize_doi_one(item$published %||% "")
            out[[length(out) + 1L]] <- make_raw_rec(
              source_tag     = paste0("bioRxivAPI_", srv),
              doi            = item$doi %||% "",
              journal        = srv,
              pub_year       = substr(clean_cell(item$date), 1L, 4L),
              first_author   = if (length(au_list) > 0L) au_list[1L] else "",
              last_author    = if (length(au_list) > 0L) au_list[length(au_list)] else "",
              title          = ttl,
              abstract       = abs,
              preprint_links = if (nzchar(pub_doi)) pub_doi else character(0)
            )
          }
        }
        if (length(coll) < 100L) break
        cursor <- cursor + 100L
      }
    }
  }
  message(sprintf("[Source 2b: bioRxiv/medRxiv Narrow API] matched %d recent preprint(s).", length(out)))
  out
}

# ---------------------------------------------------------------------------
# Source 3: NHGRI-EBI GWAS Catalog REST API v2 (+ Optional Bulk TSV)
# ---------------------------------------------------------------------------
query_gwas_catalog_v2 <- function(since_d, until_d, is_retro = FALSE,
                                  verify_ids = NULL, include_bulk = FALSE) {
  pmid_counts <- integer(0)

  if (!is.null(verify_ids) && length(verify_ids) > 0L) {
    for (pid in verify_ids) {
      Sys.sleep(0.25)
      u_pub <- sprintf("https://www.ebi.ac.uk/gwas/rest/api/v2/publications/%s", pid)
      r_pub <- tryCatch(
        httr2::request(u_pub) |> httr2::req_user_agent(UA_STRING) |>
          httr2::req_timeout(20) |> httr2::req_error(is_error = function(r) FALSE) |> httr2::req_perform(),
        error = function(e) NULL
      )
      st_pub <- if (!is.null(r_pub)) httr2::resp_status(r_pub) else NA_integer_

      Sys.sleep(0.25)
      u_std <- sprintf("https://www.ebi.ac.uk/gwas/rest/api/v2/studies?pubmed_id=%s&page=0&size=20", pid)
      r_std <- tryCatch(
        httr2::request(u_std) |> httr2::req_user_agent(UA_STRING) |>
          httr2::req_timeout(20) |> httr2::req_error(is_error = function(r) FALSE) |> httr2::req_perform(),
        error = function(e) NULL
      )
      n_el <- 0L
      if (!is.null(r_std) && httr2::resp_status(r_std) == 200L) {
        js <- tryCatch(httr2::resp_body_json(r_std, simplifyVector = FALSE), error = function(e) NULL)
        n_el <- as.integer(js$page$totalElements %||% 0L)
      }
      message(sprintf("[Source 3: GWAS Catalog v2] verify PMID %s -> /v2/publications HTTP %s, /v2/studies totalElements=%d",
                      pid, as.character(st_pub), n_el))
      if (n_el > 0L) pmid_counts[pid] <- n_el
    }
  } else {
    max_pg <- if (isTRUE(is_retro)) 40L else 10L
    for (pg in seq(0L, max_pg - 1L)) {
      Sys.sleep(0.25)
      u <- sprintf(
        "https://www.ebi.ac.uk/gwas/rest/api/v2/studies?efo_id=EFO_0004747&show_child_traits=true&page=%d&size=200",
        pg
      )
      resp <- tryCatch(
        httr2::request(u) |> httr2::req_user_agent(UA_STRING) |>
          httr2::req_timeout(30) |> httr2::req_error(is_error = function(r) FALSE) |> httr2::req_perform(),
        error = function(e) NULL
      )
      if (is.null(resp) || httr2::resp_status(resp) != 200L) break
      js <- tryCatch(httr2::resp_body_json(resp, simplifyVector = FALSE), error = function(e) NULL)
      stds <- js$`_embedded`$studies %||% js$studies %||% list()
      if (length(stds) == 0L) break

      for (s in stds) {
        pid <- normalize_pmid_one(s$pubmed_id %||% s$publicationInfo$pubmedId %||% "")
        pdt <- clean_cell(s$publication_date %||% s$publicationInfo$publicationDate %||% "")
        if (!nzchar(pid)) next
        if (nzchar(pdt) && (pdt < since_d || pdt > until_d)) next
        pmid_counts[pid] <- (pmid_counts[pid] %||% 0L) + 1L
      }
      tot_pg <- as.integer(js$page$totalPages %||% 1L)
      if (pg + 1L >= tot_pg) break
    }

    # Optional M-4 streaming bulk studies check for retrospective sweeps
    if (isTRUE(is_retro) && isTRUE(include_bulk) && requireNamespace("data.table", quietly = TRUE)) {
      bulk_url <- "https://www.ebi.ac.uk/gwas/api/search/downloads/studies_alternative"
      r_b <- tryCatch(
        httr2::request(bulk_url) |> httr2::req_user_agent(UA_STRING) |>
          httr2::req_timeout(60) |> httr2::req_error(is_error = function(r) FALSE) |> httr2::req_perform(),
        error = function(e) NULL
      )
      if (!is.null(r_b) && httr2::resp_status(r_b) == 200L &&
          !grepl("html", httr2::resp_content_type(r_b), ignore.case = TRUE)) {
        tf <- tempfile(fileext = ".tsv")
        writeBin(httr2::resp_body_raw(r_b), tf)
        dt <- tryCatch(
          data.table::fread(tf, sep = "\t", quote = "", select = c("PUBMEDID", "DATE", "DISEASE/TRAIT", "MAPPED_TRAIT")),
          error = function(e) NULL
        )
        unlink(tf)
        if (!is.null(dt) && nrow(dt) > 0L) {
          keep <- grepl("protein|proteom|cytokin|chemokin|peptid|aptamer|olink|somascan",
                        tolower(paste(dt[["DISEASE/TRAIT"]], dt[["MAPPED_TRAIT"]])), perl = TRUE)
          sub_dt <- dt[keep & dt$DATE >= since_d & dt$DATE <= until_d]
          if (nrow(sub_dt) > 0L) {
            tab <- table(as.character(sub_dt$PUBMEDID))
            for (nm in names(tab)) {
              pmid_counts[nm] <- max(pmid_counts[nm] %||% 0L, as.integer(tab[[nm]]))
            }
          }
        }
      }
    }
  }

  # Retain publications with >= 5 distinct protein GCST studies (multi-protein panel/proteome signal)
  qual_pmids <- names(pmid_counts[pmid_counts >= 5L])
  if (length(qual_pmids) == 0L) {
    message("[Source 3: GWAS Catalog v2] 0 multi-protein PMIDs matched in window.")
    return(list())
  }
  message(sprintf("[Source 3: GWAS Catalog v2] %d multi-protein PMID(s) (>=5 GCSTs); hydrating via efetch.", length(qual_pmids)))
  recs <- fetch_pubmed_efetch_xml(qual_pmids, source_tag = "GWASCatalog_v2")
  for (i in seq_along(recs)) {
    pid <- recs[[i]]$pmid
    recs[[i]]$gwas_study_count <- as.integer(pmid_counts[pid] %||% 0L)
  }
  recs
}

# ---------------------------------------------------------------------------
# Source 4: Blinded Citation Graph Snowballing (6 Landmark Seeds)
# ---------------------------------------------------------------------------
query_snowball_backward_epmc <- function(seed_pmids) {
  cited_med <- character(0)
  for (sp in seed_pmids) {
    Sys.sleep(0.35)
    u <- sprintf("https://www.ebi.ac.uk/europepmc/webservices/rest/MED/%s/references?page=1&pageSize=1000&format=json", sp)
    resp <- tryCatch(
      httr2::request(u) |> httr2::req_user_agent(UA_STRING) |>
        httr2::req_timeout(30) |> httr2::req_error(is_error = function(r) FALSE) |> httr2::req_perform(),
      error = function(e) NULL
    )
    if (is.null(resp) || httr2::resp_status(resp) != 200L) next
    js <- tryCatch(httr2::resp_body_json(resp, simplifyVector = FALSE), error = function(e) NULL)
    refs <- js$referenceList$reference %||% list()
    for (rf in refs) {
      if (identical(clean_cell(rf$source), "MED") && nzchar(clean_cell(rf$id))) {
        cited_med <- c(cited_med, clean_cell(rf$id))
      }
    }
  }
  cited_med <- unique(cited_med)
  if (length(cited_med) == 0L) return(list())
  raw_refs <- fetch_pubmed_efetch_xml(cited_med, source_tag = "Snowball_Backward")
  # C-3 Fix: Strictly require BOTH broad protein context AND broad genetics context
  Filter(function(r) {
    comb <- tolower(paste(r$title, r$abstract, paste(r$mesh_headings, collapse = " ")))
    grepl(BROAD_PROTEIN_REGEX, comb, perl = TRUE) && grepl(BROAD_GENETICS_REGEX, comb, perl = TRUE)
  }, raw_refs)
}

query_snowball_v2 <- function(since_d, until_d, is_retro = FALSE, verify_ids = NULL) {
  fwd <- run_epmc_windowed(
    builder_fn = function(s, u, include_idate, verify_pmids) {
      build_epmc_snowball_query(LANDMARK_PGWAS_SEEDS$pmid, s, u,
                                include_idate = include_idate, verify_pmids = verify_pmids)
    },
    source_tag = "Snowball_Forward",
    since_d    = since_d,
    until_d    = until_d,
    is_retro   = is_retro,
    verify_ids = verify_ids
  )
  message(sprintf("[Source 4: Snowball Forward CITES:] matched %d record(s).", length(fwd)))

  bwd <- if (isTRUE(is_retro) && is.null(verify_ids)) {
    res_b <- query_snowball_backward_epmc(LANDMARK_PGWAS_SEEDS$pmid)
    message(sprintf("[Source 4: Snowball Backward /references] matched %d proteomic-genetic record(s).", length(res_b)))
    res_b
  } else {
    list()
  }
  c(fwd, bwd)
}

# ---------------------------------------------------------------------------
# Execute Multi-Resource Discovery
# ---------------------------------------------------------------------------
recs_s1  <- query_pubmed_v2(q_pubmed)
recs_s2a <- run_epmc_windowed(build_epmc_prong2a_query, "EuropePMC", since_date, until_date,
                              is_retro = identical(mode, "retrospective"),
                              verify_ids = if (is_verify) verify_pmids else NULL)
message(sprintf("[Source 2A: Europe PMC MED+PMC+PPR] matched %d record(s).", length(recs_s2a)))

recs_s2b <- run_epmc_windowed(build_epmc_prong2b_query, "EuropePMC_MethodsRescue", since_date, until_date,
                              is_retro = identical(mode, "retrospective"),
                              verify_ids = if (is_verify) verify_pmids else NULL)
message(sprintf("[Source 2B: Europe PMC METHODS Rescue] matched %d record(s).", length(recs_s2b)))

recs_s2_bio <- if (!is_verify && identical(mode, "incremental")) {
  query_biorxiv_narrow_v2(since_date, until_date)
} else {
  list()
}

recs_s3 <- query_gwas_catalog_v2(since_date, until_date,
                                 is_retro = identical(mode, "retrospective"),
                                 verify_ids = if (is_verify) verify_pmids else NULL,
                                 include_bulk = use_gwas_bulk)

recs_s4 <- query_snowball_v2(since_date, until_date,
                             is_retro = identical(mode, "retrospective"),
                             verify_ids = if (is_verify) verify_pmids else NULL)

all_raw_recs <- c(recs_s1, recs_s2a, recs_s2b, recs_s2_bio, recs_s3, recs_s4)
raw_n <- length(all_raw_recs)
message(sprintf("Total raw multi-source hits prior to graph deduplication: %d", raw_n))

# ---------------------------------------------------------------------------
# Stage B: Connected-Components Deduplication & Provenance Merger
# ---------------------------------------------------------------------------
merge_raw_records <- function(recs) {
  n <- length(recs)
  if (n == 0L) return(list())

  parent <- seq_len(n)
  find_root <- function(i) {
    while (parent[i] != i) {
      parent[i] <<- parent[parent[i]]
      i <- parent[i]
    }
    i
  }
  union_nodes <- function(i, j) {
    ri <- find_root(i); rj <- find_root(j)
    if (ri != rj) parent[rj] <<- ri
  }

  key_map <- new.env(hash = TRUE, parent = emptyenv())
  for (i in seq_len(n)) {
    r <- recs[[i]]
    tfp <- build_title_fingerprint(r$first_author, r$pub_year, r$title)
    keys <- c(
      if (nzchar(r$doi))     paste0("doi:", r$doi) else "",
      if (nzchar(r$pmid))    paste0("pmid:", r$pmid) else "",
      if (nzchar(r$pmcid))   paste0("pmcid:", r$pmcid) else "",
      if (nzchar(r$epmc_id)) paste0("epmc:", r$epmc_id) else "",
      if (nzchar(tfp))       paste0("tfp:", tfp) else ""
    )
    keys <- keys[nzchar(keys)]
    for (k in keys) {
      if (exists(k, envir = key_map, inherits = FALSE)) {
        union_nodes(get(k, envir = key_map, inherits = FALSE), i)
      } else {
        assign(k, i, envir = key_map)
      }
    }
  }

  roots <- vapply(seq_len(n), find_root, integer(1))
  groups <- split(seq_len(n), roots)

  pick_nonempty <- function(vals) {
    vals <- vals[nzchar(vals)]
    if (length(vals) == 0L) "" else vals[1L]
  }
  pick_longest <- function(vals) {
    vals <- vals[nzchar(vals)]
    if (length(vals) == 0L) "" else vals[which.max(nchar(vals))]
  }

  lapply(groups, function(idxs) {
    sub <- recs[idxs]
    list(
      matched_sources  = paste(unique(vapply(sub, `[[`, "", "source_tag")), collapse = ";"),
      pmid             = pick_nonempty(vapply(sub, `[[`, "", "pmid")),
      doi              = pick_nonempty(vapply(sub, `[[`, "", "doi")),
      pmcid            = pick_nonempty(vapply(sub, `[[`, "", "pmcid")),
      epmc_id          = pick_nonempty(vapply(sub, `[[`, "", "epmc_id")),
      journal          = pick_nonempty(vapply(sub, `[[`, "", "journal")),
      pub_year         = pick_nonempty(vapply(sub, `[[`, "", "pub_year")),
      first_author     = pick_nonempty(vapply(sub, `[[`, "", "first_author")),
      last_author      = pick_nonempty(vapply(sub, `[[`, "", "last_author")),
      title            = pick_longest(vapply(sub, `[[`, "", "title")),
      abstract         = pick_longest(vapply(sub, `[[`, "", "abstract")),
      pub_types        = unique(unlist(lapply(sub, `[[`, "pub_types"))),
      mesh_headings    = unique(unlist(lapply(sub, `[[`, "mesh_headings"))),
      preprint_links   = unique(unlist(lapply(sub, `[[`, "preprint_links"))),
      erratum_targets  = unique(unlist(lapply(sub, `[[`, "erratum_targets"))),
      gwas_study_count = max(vapply(sub, `[[`, 0L, "gwas_study_count"))
    )
  })
}

merged_recs <- merge_raw_records(all_raw_recs)
message(sprintf("Unique works after cross-source connected-components merge: %d", length(merged_recs)))

# ---------------------------------------------------------------------------
# 4-Tier Preprint-to-Published Upgrade Detector
# ---------------------------------------------------------------------------
resolve_ppr_to_doi <- function(ppr_id) {
  if (!grepl("^PPR[0-9]+$", ppr_id)) return("")
  Sys.sleep(0.30)
  u <- sprintf("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:%s%%20AND%%20SRC:PPR&resultType=lite&format=json", ppr_id)
  resp <- tryCatch(
    httr2::request(u) |> httr2::req_user_agent(UA_STRING) |>
      httr2::req_timeout(20) |> httr2::req_error(is_error = function(r) FALSE) |> httr2::req_perform(),
    error = function(e) NULL
  )
  if (is.null(resp) || httr2::resp_status(resp) != 200L) return("")
  js <- tryCatch(httr2::resp_body_json(resp, simplifyVector = FALSE), error = function(e) NULL)
  res <- js$resultList$result
  if (is.null(res) || length(res) == 0L) return("")
  normalize_doi_one(res[[1L]]$doi %||% "")
}

check_active_preprints_direct <- function(studies_df) {
  active_idx <- which(studies_df$status == "Preprint")
  upgrades <- list()
  if (length(active_idx) == 0L) return(upgrades)

  for (idx in active_idx) {
    sid <- clean_cell(studies_df$study_id[idx])
    p_doi <- normalize_doi_one(studies_df$doi[idx])
    if (!nzchar(p_doi)) next
    pub_doi <- ""; pub_pmid <- ""; via_tier <- ""

    # Tier 1: bioRxiv / medRxiv /pubs/ endpoint
    if (grepl("^10\\.1101/", p_doi)) {
      for (srv in c("biorxiv", "medrxiv")) {
        Sys.sleep(0.30)
        u <- sprintf("https://api.biorxiv.org/pubs/%s/%s/na/json", srv, p_doi)
        r <- tryCatch(
          httr2::request(u) |> httr2::req_user_agent(UA_STRING) |>
            httr2::req_timeout(20) |> httr2::req_error(is_error = function(x) FALSE) |> httr2::req_perform(),
          error = function(e) NULL
        )
        if (!is.null(r) && httr2::resp_status(r) == 200L) {
          js <- tryCatch(httr2::resp_body_json(r, simplifyVector = FALSE), error = function(e) NULL)
          if (!is.null(js$collection) && length(js$collection) > 0L) {
            cand_d <- normalize_doi_one(js$collection[[1L]]$published_doi %||% "")
            if (nzchar(cand_d)) { pub_doi <- cand_d; via_tier <- paste0("Tier1_biorxiv_pubs(", srv, ")"); break }
          }
        }
      }
    }

    # Tier 2: Europe PMC commentCorrectionList ("Preprint of")
    if (!nzchar(pub_doi) && !nzchar(pub_pmid)) {
      Sys.sleep(0.30)
      u2 <- paste0("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:",
                   utils::URLencode(p_doi, reserved = TRUE), "&resultType=core&format=json")
      r2 <- tryCatch(
        httr2::request(u2) |> httr2::req_user_agent(UA_STRING) |>
          httr2::req_timeout(20) |> httr2::req_error(is_error = function(x) FALSE) |> httr2::req_perform(),
        error = function(e) NULL
      )
      if (!is.null(r2) && httr2::resp_status(r2) == 200L) {
        js2 <- tryCatch(httr2::resp_body_json(r2, simplifyVector = FALSE), error = function(e) NULL)
        items <- js2$resultList$result %||% list()
        for (it in items) {
          for (cc in (it$commentCorrectionList$commentCorrection %||% list())) {
            if (identical(clean_cell(cc$type), "Preprint of") && nzchar(clean_cell(cc$id))) {
              pub_pmid <- normalize_pmid_one(cc$id)
              via_tier <- "Tier2_EuropePMC_PreprintOf"
            }
          }
        }
      }
    }

    # Tier 3: Crossref is-preprint-of
    if (!nzchar(pub_doi) && !nzchar(pub_pmid)) {
      Sys.sleep(0.30)
      u3 <- paste0("https://api.crossref.org/works/", utils::URLencode(p_doi, reserved = TRUE))
      r3 <- tryCatch(
        httr2::request(u3) |> httr2::req_user_agent(UA_STRING) |>
          httr2::req_timeout(20) |> httr2::req_error(is_error = function(x) FALSE) |> httr2::req_perform(),
        error = function(e) NULL
      )
      if (!is.null(r3) && httr2::resp_status(r3) == 200L) {
        js3 <- tryCatch(httr2::resp_body_json(r3, simplifyVector = FALSE), error = function(e) NULL)
        rel <- js3$message$relation[["is-preprint-of"]]
        if (!is.null(rel) && length(rel) > 0L) {
          pub_doi <- normalize_doi_one(rel[[1L]]$id %||% "")
          if (nzchar(pub_doi)) via_tier <- "Tier3_Crossref_isPreprintOf"
        }
      }
    }

    if (nzchar(pub_doi) || nzchar(pub_pmid)) {
      upgrades[[sid]] <- list(
        preprint_study_id = sid, preprint_doi = p_doi,
        published_doi = pub_doi, published_pmid = pub_pmid, via = via_tier
      )
    }
  }
  upgrades
}

direct_preprint_upgrades <- if (!is_verify) check_active_preprints_direct(studies) else list()

# Build lookup tables of catalogued studies & preprints for lineage + erratum detection
cat_pmids <- unique( vapply(studies$pmid, normalize_pmid_one, character(1)) )
cat_pmids <- cat_pmids[nzchar(cat_pmids)]
cat_dois  <- unique( vapply(studies$doi, normalize_doi_one, character(1)) )
cat_dois  <- cat_dois[nzchar(cat_dois)]

preprint_rows_df <- studies[studies$status %in% c("Preprint", if (is_verify) "Superseded Preprint" else character(0)), , drop = FALSE]

known_ids <- tolower(unique(trimws(c(
  vapply(studies$doi, normalize_doi_one, character(1)),
  vapply(studies$pmid, normalize_pmid_one, character(1)),
  ifelse(nzchar(studies$pmid), paste0("pmid:", vapply(studies$pmid, normalize_pmid_one, character(1))), ""),
  vapply(excluded$identifier, normalize_doi_one, character(1)),
  extract_pmids_from_cell(paste(excluded$identifier, collapse = " ")),
  tolower(excluded$identifier),
  vapply(queue$identifier, normalize_doi_one, character(1)),
  extract_pmids_from_cell(paste(queue$identifier, collapse = " ")),
  tolower(queue$identifier)
))))
known_ids <- known_ids[nzchar(known_ids)]

# ---------------------------------------------------------------------------
# Stage C: Triage Scoring, Preprint Lineage Attachment & Queue Formatting
# ---------------------------------------------------------------------------
run_stamp <- format(Sys.time(), "%Y%m%d_%H%M", tz = "UTC")
run_id    <- paste0("run_", run_stamp)

dump_rows  <- list()
queue_rows <- list()

for (m in merged_recs) {
  doi  <- m$doi
  pmid <- m$pmid
  pmcid <- m$pmcid
  epmc_id <- m$epmc_id

  id_primary <- if (nzchar(doi)) doi else if (nzchar(pmid)) paste0("PMID:", pmid) else epmc_id
  if (!nzchar(id_primary)) next

  # Resolve PPR* preprint link to DOI when present so lineage note is maximally informative
  lineage_parts <- m$preprint_links
  ppr_tokens <- grep("^PPR[0-9]+$", lineage_parts, value = TRUE)
  for (ptok in ppr_tokens) {
    resolved_d <- resolve_ppr_to_doi(ptok)
    if (nzchar(resolved_d)) lineage_parts <- unique(c(lineage_parts, resolved_d))
  }

  # Match against catalogued Preprint rows (Tiers 1-3 exact DOI + Tier 4 author/N/Jaccard fallback)
  matched_preprint_sid <- ""
  if (nrow(preprint_rows_df) > 0L && (nzchar(pmid) || !grepl("^10\\.(1101|21203)/", doi))) {
    for (pi in seq_len(nrow(preprint_rows_df))) {
      p_sid <- clean_cell(preprint_rows_df$study_id[pi])
      p_doi <- normalize_doi_one(preprint_rows_df$doi[pi])
      p_st  <- clean_cell(preprint_rows_df$status[pi])

      # Exact link via Tiers 1-3
      if (nzchar(p_doi) && (p_doi %in% lineage_parts)) {
        lineage_parts <- unique(c(lineage_parts, p_sid))
        if (identical(p_st, "Preprint")) matched_preprint_sid <- p_sid
        next
      }
      if (p_sid %in% names(direct_preprint_upgrades)) {
        up <- direct_preprint_upgrades[[p_sid]]
        if ((nzchar(up$published_doi) && identical(up$published_doi, doi)) ||
            (nzchar(up$published_pmid) && identical(up$published_pmid, pmid))) {
          lineage_parts <- unique(c(lineage_parts, p_doi, p_sid, up$via))
          if (identical(p_st, "Preprint")) matched_preprint_sid <- p_sid
          next
        }
      }

      # Tier 4 Fuzzy Author + (Title Jaccard >= 0.35 OR Platform + Sample Size within +/-5%)
      if (identical(p_st, "Preprint")) {
        p_au  <- normalize_author_surname(preprint_rows_df$first_author[pi])
        m_au1 <- normalize_author_surname(m$first_author)
        m_au2 <- normalize_author_surname(m$last_author)
        if (nzchar(p_au) && (identical(p_au, m_au1) || identical(p_au, m_au2))) {
          jac <- title_jaccard(preprint_rows_df$title[pi], m$title)
          if (jac >= 0.35) {
            matched_preprint_sid <- p_sid
            lineage_parts <- unique(c(lineage_parts, p_doi, paste0("Tier4_TitleJaccard=", round(jac, 2))))
          }
        }
      }
    }
  }

  err_on_cat <- any(m$erratum_targets %in% c(cat_pmids, cat_dois))

  sc <- score_triage_v2(
    title                      = m$title,
    abstract_text              = m$abstract,
    pub_types                  = m$pub_types,
    mesh_headings              = m$mesh_headings,
    gwas_catalog_study_count   = m$gwas_study_count,
    erratum_targets_catalogued = err_on_cat
  )

  auto_cls <- if (nzchar(matched_preprint_sid) && sc$P_score >= 5L) {
    sprintf("Preprint -> Published Match (%s)", matched_preprint_sid)
  } else {
    sc$classification
  }

  plat_det <- detect_platform_v2(paste(m$title, m$abstract))
  n_det    <- extract_sample_size_v2(paste(m$title, m$abstract))
  lineage_note <- paste(unique(lineage_parts[nzchar(lineage_parts)]), collapse = "|")
  rules_str    <- paste(c(sc$p_rules, sc$s_rules), collapse = ",")

  # M-1 Fix: Keep curation_status starting with "Pending" for all actionable items
  cur_status <- if (grepl("^Auto-Deprioritised", auto_cls)) {
    "Auto-Triage: Low Priority"
  } else {
    "Pending Curator Review"
  }

  # Record every merged candidate in the uncapped search dump (for post-hoc evaluation)
  dump_rows[[length(dump_rows) + 1L]] <- data.frame(
    run_id                = run_id,
    identifier            = id_primary,
    pmid                  = pmid,
    doi                   = doi,
    pmcid                 = pmcid,
    epmc_id               = epmc_id,
    first_author          = m$first_author,
    pub_year              = m$pub_year,
    journal               = m$journal,
    title                 = m$title,
    detected_platform     = plat_det,
    detected_n            = n_det,
    matched_sources       = m$matched_sources,
    p_score               = sc$P_score,
    s_score               = sc$S_score,
    triage_rules          = rules_str,
    auto_classification   = auto_cls,
    curation_status       = cur_status,
    preprint_lineage_note = lineage_note,
    stringsAsFactors      = FALSE
  )

  # Unseen check for candidate_queue.tsv (bypassed in --verify-pmids mode, or when active preprint upgraded)
  keys_check <- tolower(c(id_primary, doi, pmid, if (nzchar(pmid)) paste0("pmid:", pmid) else ""))
  keys_check <- keys_check[nzchar(keys_check)]
  if (!is_verify && !nzchar(matched_preprint_sid) && any(keys_check %in% known_ids)) next

  src_display <- if (nzchar(m$journal)) m$journal else strsplit(m$matched_sources, ";")[[1L]][1L]

  queue_rows[[length(queue_rows) + 1L]] <- data.frame(
    candidate_id          = sprintf("cand_%s_%03d", gsub("[^0-9]", "", until_date), length(queue_rows) + 1L),
    run_id                = run_id,
    source                = src_display,
    identifier            = id_primary,
    detected_platform     = plat_det,
    detected_n            = n_det,
    title                 = m$title,
    auto_classification   = auto_cls,
    curation_status       = cur_status,
    pmid                  = pmid,
    doi                   = doi,
    pmcid                 = pmcid,
    epmc_id               = epmc_id,
    matched_sources       = m$matched_sources,
    p_score               = sc$P_score,
    s_score               = sc$S_score,
    preprint_lineage_note = lineage_note,
    stringsAsFactors      = FALSE
  )
}

empty_queue_df <- function() {
  df <- as.data.frame(setNames(replicate(length(CANONICAL_CANDIDATE_QUEUE_COLS), character(0), simplify = FALSE),
                               CANONICAL_CANDIDATE_QUEUE_COLS), stringsAsFactors = FALSE)
  df$p_score <- integer(0); df$s_score <- integer(0)
  df
}

new_queue_df <- if (length(queue_rows) > 0L) do.call(rbind, queue_rows) else empty_queue_df()
dump_df      <- if (length(dump_rows) > 0L)  do.call(rbind, dump_rows)  else data.frame()

dedup_n   <- nrow(new_queue_df)
auto_mr_n <- if (dedup_n > 0L) sum(grepl("^Auto-Deprioritised", new_queue_df$auto_classification)) else 0L
review_n  <- dedup_n - auto_mr_n
upg_n     <- if (dedup_n > 0L) sum(grepl("^Preprint -> Published Match", new_queue_df$auto_classification)) else 0L

# Write complete uncapped audit dump to --out-dir or git-ignored data/
target_dump_dir <- if (nzchar(out_dir)) out_dir else "data"
dump_path <- file.path(target_dump_dir, sprintf("search_dump_%s.tsv", run_id))
if (nrow(dump_df) > 0L) {
  write_tsv_safe(dump_df, dump_path)
  message(sprintf("Wrote uncapped search ledger (%d row(s)) to %s", nrow(dump_df), dump_path))
}

sources_desc <- sprintf(
  "PubMed(%s);EuropePMC(MED+PMC+PPR);GWASCatalog_v2;Snowball [%s..%s]",
  pm_date_field, since_date, until_date
)

new_run_row <- data.frame(
  run_id                    = run_id,
  run_date                  = until_date,
  trigger                   = clean_cell(trigger_lb),
  sources_queried           = sources_desc,
  raw_records               = as.integer(raw_n),
  deduplicated_new          = as.integer(dedup_n),
  auto_flagged_secondary_mr = as.integer(auto_mr_n),
  candidates_for_review     = as.integer(review_n),
  studies_added             = 0L,
  preprints_upgraded        = as.integer(upg_n),
  pr_number                 = if (is_verify) "Verify-Only" else "Automated PR",
  release_tag               = "",
  status                    = if (review_n > 0L) sprintf("Pending Triage (%d candidates)", review_n) else "Completed (0 new primary)",
  stringsAsFactors          = FALSE
)

if (nzchar(out_dir)) {
  write_tsv_safe(new_queue_df, file.path(out_dir, "candidate_queue_preview.tsv"))
  write_tsv_safe(new_run_row,  file.path(out_dir, "search_runs_preview.tsv"))
  message(sprintf("Wrote preview tables (%d candidate(s)) to %s", nrow(new_queue_df), out_dir))
}

if (isTRUE(write_catalogue)) {
  if (dedup_n > 0L) {
    actionable <- new_queue_df[!grepl("^Auto-Deprioritised", new_queue_df$auto_classification), , drop = FALSE]
    deprior    <- utils::head(new_queue_df[grepl("^Auto-Deprioritised", new_queue_df$auto_classification), , drop = FALSE], 15L)

    # Backfill any missing v2 handoff columns on existing queue before rbind (C-1 safe upgrade)
    for (col in CANONICAL_CANDIDATE_QUEUE_COLS) {
      if (!(col %in% names(queue))) queue[[col]] <- ""
    }
    queue_aligned <- queue[, CANONICAL_CANDIDATE_QUEUE_COLS, drop = FALSE]
    updated_queue <- rbind(actionable[, CANONICAL_CANDIDATE_QUEUE_COLS, drop = FALSE],
                           deprior[, CANONICAL_CANDIDATE_QUEUE_COLS, drop = FALSE],
                           queue_aligned)
    write_tsv_safe(updated_queue, file.path("catalogue", "candidate_queue.tsv"))
    message(sprintf("Updated catalogue/candidate_queue.tsv (+%d actionable, +%d sample deprioritised).",
                    nrow(actionable), nrow(deprior)))
  } else {
    message("0 unseen records found in window.")
  }

  runs_updated <- rbind(runs, new_run_row)
  write_tsv_safe(runs_updated, file.path("catalogue", "search_runs.tsv"))

  val_status <- system2("Rscript", "src/validate_catalogue.R")
  if (val_status != 0L) stop("Catalogue validation failed after search update!", call. = FALSE)
}
