#!/usr/bin/env Rscript
# Multi-resource literature discovery, connected-components + preprint-lineage
# deduplication, 4-tier preprint-to-published upgrade detector, strict canonical
# first-publication-date windowing, selective Stage B.5 PMC full-text multi-omic
# modality verification, and Unified Two-Stage Dual-Axis Triage engine (v4.11)
# for proteomic-gwas-catalogue.
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

parse_month_num <- function(m) {
  if (is.null(m) || length(m) == 0L || is.na(m[1L])) return("")
  s <- tolower(trimws(as.character(m[1L])))
  if (!nzchar(s)) return("")
  if (grepl("^[0-9]{1,2}$", s)) {
    v <- as.integer(s)
    if (!is.na(v) && v >= 1L && v <= 12L) return(sprintf("%02d", v))
  }
  abbr <- substr(s, 1L, 3L)
  idx <- match(abbr, c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"))
  if (!is.na(idx)) sprintf("%02d", idx) else ""
}

format_ymd_safe <- function(yr, mo = "", dy = "") {
  y <- regmatches(as.character(yr %||% ""), regexpr("[0-9]{4}", as.character(yr %||% "")))
  if (length(y) == 0L || !nzchar(y)) return("")
  m <- parse_month_num(mo)
  if (!nzchar(m)) m <- "01"
  d_raw <- trimws(as.character(dy %||% ""))
  d <- if (grepl("^[0-9]{1,2}$", d_raw)) {
    dv <- as.integer(d_raw)
    if (!is.na(dv) && dv >= 1L && dv <= 31L) sprintf("%02d", dv) else "01"
  } else {
    "01"
  }
  sprintf("%s-%s-%s", y, m, d)
}

UA_STRING <- "proteomic-gwas-catalogue/4.0 (https://github.com/matthewlee/proteomic-gwas-catalogue)"

# Polite HTTP engine with jittered exponential backoff, Retry-After support, and
# `expect_json` validation against transient EBI/proxy HTTP 200 text/html error pages.
perform_polite_http <- function(req, sleep_s = 0.55, max_tries = 5L, label = "HTTP", expect_json = FALSE) {
  req_safe <- req |> httr2::req_error(is_error = function(r) FALSE)
  for (attempt in seq_len(max_tries)) {
    Sys.sleep(sleep_s)
    resp <- tryCatch(
      httr2::req_perform(req_safe),
      error = function(e) {
        if (attempt == max_tries) {
          warning(sprintf("%s failed after %d tries: %s", label, max_tries, conditionMessage(e)))
        }
        NULL
      }
    )
    if (!is.null(resp)) {
      st <- httr2::resp_status(resp)
      ct <- httr2::resp_content_type(resp) %||% ""
      bad_json <- isTRUE(expect_json) && st == 200L && !grepl("json", ct, ignore.case = TRUE)
      if (!(st %in% c(408L, 429L, 500L, 502L, 503L, 504L)) && !bad_json) {
        return(resp)
      }
      ra <- suppressWarnings(as.numeric(httr2::resp_header(resp, "Retry-After") %||% NA_real_))
      wait_s <- if (!is.na(ra) && ra > 0 && ra <= 60) {
        ra
      } else {
        min(45, 1.5 * (2 ^ (attempt - 1L)) + stats::runif(1L, 0.1, 0.7))
      }
      if (attempt < max_tries) {
        message(sprintf("[%s] Transient HTTP %d (ct=%s) on attempt %d/%d; backing off %.1fs...",
                        label, st, ct, attempt, max_tries, wait_s))
        Sys.sleep(wait_s)
      } else {
        warning(sprintf("%s returned HTTP %d (ct=%s) after %d tries.", label, st, ct, max_tries))
        return(if (bad_json) NULL else resp)
      }
    } else if (attempt < max_tries) {
      wait_s <- min(45, 1.5 * (2 ^ (attempt - 1L)) + stats::runif(1L, 0.1, 0.7))
      Sys.sleep(wait_s)
    }
  }
  NULL
}

# Compile queries from single source of truth (src/search/lexicons.R)
pm_date_field <- if (identical(mode, "retrospective")) "dp_or_edat" else "edat"
q_pubmed      <- build_pubmed_query(since_date, until_date, date_field = pm_date_field,
                                    verify_pmids = if (is_verify) verify_pmids else NULL)
q_epmc_2a_p1  <- build_epmc_prong2a_platforms_query(since_date, until_date,
                                                    include_idate = TRUE,
                                                    verify_pmids = if (is_verify) verify_pmids else NULL)
q_epmc_2a_p2  <- build_epmc_prong2a_proteome_query(since_date, until_date,
                                                   include_idate = TRUE,
                                                   verify_pmids = if (is_verify) verify_pmids else NULL)
q_epmc_2b     <- build_epmc_prong2b_query(since_date, until_date,
                                          include_idate = TRUE,
                                          verify_pmids = if (is_verify) verify_pmids else NULL)
q_epmc_snow   <- build_epmc_snowball_query(LANDMARK_PGWAS_SEEDS$pmid, since_date, until_date,
                                           include_idate = TRUE,
                                           verify_pmids = if (is_verify) verify_pmids else NULL)

message(sprintf(
  "Proteomic GWAS Search v4.11 | mode=%s | window=%s..%s | verify_pmids=%s | dry_run=%s | write_catalogue=%s",
  mode, since_date, until_date,
  if (is_verify) paste(verify_pmids, collapse = ",") else "none",
  dry_run, write_catalogue
))

if (dry_run) {
  cat("\n--- [Source 1: PubMed E-utilities POST Query] ---\n", q_pubmed, "\n", sep = "")
  cat("\n--- [Source 2A-1: Europe PMC Platforms + pQTL Anchors Shard] ---\n", q_epmc_2a_p1, "\n", sep = "")
  cat("\n--- [Source 2A-2: Europe PMC Proteome + Multi-Omic Anchors Shard] ---\n", q_epmc_2a_p2, "\n", sep = "")
  cat("\n--- [Source 2B: Europe PMC Full-Text METHODS Rescue Query] ---\n", q_epmc_2b, "\n", sep = "")
  cat("\n--- [Source 4: Europe PMC 10-Seed Server-Side Snowball Query] ---\n", q_epmc_snow, "\n", sep = "")
  cat("\n--- [Source 3: NHGRI-EBI GWAS Catalog REST API v2] ---\n",
      "GET https://www.ebi.ac.uk/gwas/rest/api/v2/publications?page=0&size=500&sort=publication_date&direction=desc -> /v2/studies?pubmed_id=... (is_gwas_catalog_protein_trait + EFO_0004747 OR-fallback)\n", sep = "")
  quit(status = 0L)
}

# ---------------------------------------------------------------------------
# Load Existing Catalogue State (Read-Only unless write_catalogue == TRUE)
# ---------------------------------------------------------------------------
studies  <- read_tsv_safe(file.path("catalogue", "studies.tsv"))
strata_cat <- if (file.exists(file.path("catalogue", "strata.tsv"))) read_tsv_safe(file.path("catalogue", "strata.tsv")) else data.frame()
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
                         first_pub_date = "",
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
    first_pub_date   = clean_cell(first_pub_date),
    first_author     = clean_cell(first_author),
    last_author      = clean_cell(last_author),
    title            = clean_cell(strip_html_tags(title)),
    abstract         = clean_cell(strip_html_tags(abstract)),
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
parse_pubmed_xml_doc <- function(doc, source_tag = "PubMed") {
  out <- list()
  if (is.null(doc)) return(out)
  # Hydrate both <PubmedArticle> and <PubmedBookArticle> so zero PubMed hits are dropped.
  # Safeguard 2: Strictly scope ALL XPath queries to ./MedlineCitation, ./PubmedData,
  # ./BookDocument, or ./PubmedBookData so cited <ReferenceList><Reference><ArticleIdList>
  # PMCIDs/DOIs/PMIDs NEVER leak into article metadata!
  arts <- xml2::xml_find_all(doc, ".//PubmedArticle | .//PubmedBookArticle")
  for (art in arts) {
    is_book <- identical(xml2::xml_name(art), "PubmedBookArticle")
    pmid_v <- xml2::xml_text(xml2::xml_find_first(
      art,
      if (is_book) "./BookDocument/PMID" else "./MedlineCitation/PMID"
    ))
    ttl_v  <- xml2::xml_text(xml2::xml_find_first(
      art,
      if (is_book) "./BookDocument/ArticleTitle | ./BookDocument/Book/BookTitle" else "./MedlineCitation/Article/ArticleTitle"
    ))
    jrn_v  <- xml2::xml_text(xml2::xml_find_first(
      art,
      if (is_book) "./BookDocument/Book/BookTitle" else "./MedlineCitation/Article/Journal/Title"
    ))

    # Parse electronic ArticleDate (epub), PubmedData History dates, and JournalIssue/Book PubDate (ppub)
    ad_node <- xml2::xml_find_first(
      art,
      if (is_book) "./BookDocument/ArticleDate" else "./MedlineCitation/Article/ArticleDate[@DateType='Electronic'] | ./MedlineCitation/Article/ArticleDate"
    )
    art_dt  <- if (!inherits(ad_node, "xml_missing")) {
      format_ymd_safe(
        xml2::xml_text(xml2::xml_find_first(ad_node, "./Year")),
        xml2::xml_text(xml2::xml_find_first(ad_node, "./Month")),
        xml2::xml_text(xml2::xml_find_first(ad_node, "./Day"))
      )
    } else {
      ""
    }

    pd_y <- xml2::xml_text(xml2::xml_find_first(art, if (is_book) "./BookDocument/Book/PubDate/Year" else "./MedlineCitation/Article/Journal/JournalIssue/PubDate/Year"))
    pd_m <- xml2::xml_text(xml2::xml_find_first(art, if (is_book) "./BookDocument/Book/PubDate/Month" else "./MedlineCitation/Article/Journal/JournalIssue/PubDate/Month"))
    pd_d <- xml2::xml_text(xml2::xml_find_first(art, if (is_book) "./BookDocument/Book/PubDate/Day" else "./MedlineCitation/Article/Journal/JournalIssue/PubDate/Day"))
    if (is.na(pd_y) || !nzchar(pd_y)) {
      md <- xml2::xml_text(xml2::xml_find_first(art, if (is_book) "./BookDocument/Book/PubDate/MedlineDate" else "./MedlineCitation/Article/Journal/JournalIssue/PubDate/MedlineDate")) %||% ""
      pd_y <- regmatches(md, regexpr("[0-9]{4}", md))
    }
    pub_dt <- format_ymd_safe(pd_y, pd_m, pd_d)

    id_nodes <- xml2::xml_find_all(
      art,
      if (is_book) "./PubmedBookData/ArticleIdList/ArticleId" else "./PubmedData/ArticleIdList/ArticleId"
    )
    doi_v <- ""; pmc_v <- ""
    for (idn in id_nodes) {
      idt <- xml2::xml_attr(idn, "IdType") %||% ""
      idv <- xml2::xml_text(idn)
      if (identical(idt, "doi") && !nzchar(doi_v)) doi_v <- idv
      if (identical(idt, "pmc") && !nzchar(pmc_v)) pmc_v <- idv
    }

    doi_norm  <- normalize_doi_one(doi_v)
    prep_dt   <- extract_preprint_doi_date(doi_norm)
    is_prep_d <- any(is_preprint_doi(doi_norm))

    # Safeguard 1: When publishers (e.g. OUP) omit <ArticleDate> upon issue assignment,
    # collect online History dates (aheadofprint, epublish, pubmed, entrez) into pm_hist_dt
    # and apply strictly monotonic min(pub_dt, pm_hist_dt) so older print dates are NEVER
    # shifted forward by retroactive PMC/NLM history entries!
    pm_hist_dt <- ""
    if (!nzchar(art_dt) && !is_prep_d) {
      h_xpath <- if (is_book) {
        "./PubmedBookData/History/PubMedPubDate[@PubStatus='aheadofprint' or @PubStatus='epublish' or @PubStatus='pubmed' or @PubStatus='entrez']"
      } else {
        "./PubmedData/History/PubMedPubDate[@PubStatus='aheadofprint' or @PubStatus='epublish' or @PubStatus='pubmed' or @PubStatus='entrez']"
      }
      h_nodes <- xml2::xml_find_all(art, h_xpath)
      if (length(h_nodes) > 0L) {
        h_dts <- vapply(h_nodes, function(hn) {
          format_ymd_safe(
            xml2::xml_text(xml2::xml_find_first(hn, "./Year")),
            xml2::xml_text(xml2::xml_find_first(hn, "./Month")),
            xml2::xml_text(xml2::xml_find_first(hn, "./Day"))
          )
        }, character(1))
        h_dts <- h_dts[nzchar(h_dts)]
        if (length(h_dts) > 0L) pm_hist_dt <- min(h_dts)
      }
    }

    # Canonical raw-record first_pub_date precedence:
    # preprint_doi_date > epub ArticleDate > min(PubDate, pm_hist_dt)
    first_dt <- if (is_prep_d && nzchar(prep_dt)) {
      prep_dt
    } else if (nzchar(art_dt)) {
      art_dt
    } else if (nzchar(pub_dt) && nzchar(pm_hist_dt)) {
      min(pub_dt, pm_hist_dt)
    } else if (nzchar(pub_dt)) {
      pub_dt
    } else {
      pm_hist_dt
    }
    yr_v <- if (nzchar(first_dt)) substr(first_dt, 1L, 4L) else (pd_y %||% "")

    abs_nodes <- xml2::xml_find_all(
      art,
      if (is_book) "./BookDocument/Abstract/AbstractText" else "./MedlineCitation/Article/Abstract/AbstractText"
    )
    abs_parts <- vapply(abs_nodes, function(n) {
      lbl <- xml2::xml_attr(n, "Label")
      tx  <- xml2::xml_text(n)
      if (!is.na(lbl) && nzchar(lbl)) paste0(lbl, ": ", tx) else tx
    }, character(1))
    abs_v <- paste(abs_parts, collapse = " ")

    au_surnames <- xml2::xml_text(xml2::xml_find_all(
      art,
      if (is_book) "./BookDocument/AuthorList/Author/LastName" else "./MedlineCitation/Article/AuthorList/Author/LastName"
    ))
    au_collect  <- xml2::xml_text(xml2::xml_find_all(
      art,
      if (is_book) "./BookDocument/AuthorList/Author/CollectiveName" else "./MedlineCitation/Article/AuthorList/Author/CollectiveName"
    ))
    first_au <- if (length(au_surnames) > 0L) au_surnames[1L] else if (length(au_collect) > 0L) au_collect[1L] else ""
    last_au  <- if (length(au_surnames) > 0L) au_surnames[length(au_surnames)] else first_au

    pts_v <- xml2::xml_text(xml2::xml_find_all(
      art,
      if (is_book) "./BookDocument/PublicationType" else "./MedlineCitation/Article/PublicationTypeList/PublicationType"
    ))
    if (is_book) pts_v <- unique(c(pts_v, "Book Chapter"))
    mesh_v <- xml2::xml_text(xml2::xml_find_all(art, "./MedlineCitation/MeshHeadingList/MeshHeading/DescriptorName"))

    cc_nodes <- xml2::xml_find_all(art, "./MedlineCitation/CommentsCorrectionsList/CommentsCorrections")
    pre_links <- character(0); err_links <- character(0)
    for (cc in cc_nodes) {
      rt <- xml2::xml_attr(cc, "RefType") %||% ""
      rs <- xml2::xml_text(xml2::xml_find_first(cc, "./RefSource")) %||% ""
      rp <- xml2::xml_text(xml2::xml_find_first(cc, "./PMID")) %||% ""
      if (rt %in% c("UpdateOf", "UpdateIn", "PreprintIn", "PreprintOf")) {
        pre_links <- c(pre_links, extract_dois_from_cell(rs), if (!is.na(rp) && nzchar(rp)) paste0("PMID:", rp) else "")
      }
      if (rt %in% c("ErratumFor", "CorrectionFor")) {
        err_links <- c(err_links, extract_dois_from_cell(rs), if (!is.na(rp) && nzchar(rp)) rp else "")
      }
    }

    out[[length(out) + 1L]] <- make_raw_rec(
      source_tag = source_tag, pmid = pmid_v, doi = doi_norm, pmcid = pmc_v,
      journal = jrn_v, pub_year = yr_v, first_pub_date = first_dt,
      first_author = first_au, last_author = last_au,
      title = ttl_v, abstract = abs_v, pub_types = pts_v, mesh_headings = mesh_v,
      preprint_links = pre_links, erratum_targets = err_links
    )
  }
  out
}

fetch_pubmed_efetch_xml <- function(pmids, source_tag = "PubMed") {
  pmids <- unique(pmids[nzchar(pmids)])
  if (length(pmids) == 0L) return(list())
  ncbi_key <- Sys.getenv("NCBI_API_KEY", unset = "")
  sleep_s  <- 0.55
  out <- list()

  chunks <- split(pmids, ceiling(seq_along(pmids) / 200L))
  for (ch in chunks) {
    form_params <- list(
      db      = "pubmed",
      id      = paste(ch, collapse = ","),
      retmode = "xml",
      tool    = "proteomic_gwas_catalogue"
    )
    if (nzchar(ncbi_key)) form_params$api_key <- ncbi_key

    resp <- perform_polite_http(
      httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi") |>
        httr2::req_user_agent(UA_STRING) |>
        httr2::req_timeout(60) |>
        httr2::req_body_form(!!!form_params),
      sleep_s = sleep_s,
      label   = "PubMed efetch"
    )
    if (is.null(resp) || httr2::resp_status(resp) != 200L) next
    doc <- tryCatch(xml2::read_xml(httr2::resp_body_string(resp)), error = function(e) NULL)
    if (is.null(doc)) next
    out <- c(out, parse_pubmed_xml_doc(doc, source_tag = source_tag))
  }
  out
}

query_pubmed_v2 <- function(query_str) {
  ncbi_key <- Sys.getenv("NCBI_API_KEY", unset = "")
  sleep_s  <- 0.55
  all_pmids <- character(0)
  retstart <- 0L
  retmax   <- 500L

  repeat {
    form_params <- list(
      db       = "pubmed",
      term     = query_str,
      retmode  = "json",
      retstart = as.character(retstart),
      retmax   = as.character(retmax),
      tool     = "proteomic_gwas_catalogue"
    )
    if (nzchar(ncbi_key)) form_params$api_key <- ncbi_key

    resp <- perform_polite_http(
      httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi") |>
        httr2::req_user_agent(UA_STRING) |>
        httr2::req_timeout(60) |>
        httr2::req_body_form(!!!form_params),
      sleep_s     = sleep_s,
      label       = "PubMed esearch",
      expect_json = TRUE
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
# Source 2: Europe PMC REST API (Sharded TITLE_ABS searchPOST across MED+PMC+PPR)
# ---------------------------------------------------------------------------
parse_epmc_result_item <- function(r, default_tag = "EuropePMC") {
  src_db <- clean_cell(r$source)
  tag <- if (identical(default_tag, "EuropePMC")) {
    if (identical(src_db, "PPR")) "EuropePMC_PPR" else "EuropePMC_MED"
  } else {
    default_tag
  }

  # Preserve id on both PPR and PMC-only records so zero SRC:PMC records are silently dropped
  epmc_id <- if (src_db %in% c("PPR", "PMC")) clean_cell(r$id) else ""
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

  jrn <- clean_cell(r$journalTitle %||% r$journalInfo$journal$title %||% "")
  if (!nzchar(jrn)) {
    jrn <- if (identical(src_db, "PPR")) clean_cell(r$bookOrReportDetails$publisher %||% "Preprint") else src_db
  }

  doi_norm  <- normalize_doi_one(r$doi %||% "")
  fpd       <- clean_cell(r$firstPublicationDate %||% "")
  epd       <- clean_cell(r$electronicPublicationDate %||% "")
  fid       <- clean_cell(r$firstIndexDate %||% "")
  prep_dt   <- extract_preprint_doi_date(doi_norm)
  is_prep_d <- any(is_preprint_doi(doi_norm)) || identical(src_db, "PPR")
  valid_pub_dts <- c(
    if (grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", fpd)) fpd else character(0),
    if (grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", epd)) epd else character(0)
  )
  first_dt <- if (is_prep_d && nzchar(prep_dt)) {
    prep_dt
  } else if (!is_prep_d && length(valid_pub_dts) > 0L && !nzchar(epd) &&
             grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", fid) && fid < min(valid_pub_dts)) {
    # Advance-Online publication indexed before print issue without electronicPublicationDate
    fid
  } else if (length(valid_pub_dts) > 0L) {
    min(valid_pub_dts)
  } else {
    ""
  }
  yr_v <- if (nzchar(first_dt)) substr(first_dt, 1L, 4L) else clean_cell(r$pubYear %||% "")

  make_raw_rec(
    source_tag      = tag,
    pmid            = r$pmid %||% "",
    doi             = doi_norm,
    pmcid           = r$pmcid %||% (if (identical(src_db, "PMC")) r$id else ""),
    epmc_id         = epmc_id,
    journal         = jrn,
    pub_year        = yr_v,
    first_pub_date  = first_dt,
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
    resp <- perform_polite_http(
      httr2::request("https://www.ebi.ac.uk/europepmc/webservices/rest/searchPOST") |>
        httr2::req_user_agent(UA_STRING) |>
        httr2::req_timeout(60) |>
        httr2::req_headers(`Content-Type` = "application/x-www-form-urlencoded") |>
        httr2::req_body_form(
          query      = query_str,
          resultType = "core",
          pageSize   = as.character(page_size),
          cursorMark = cursor,
          format     = "json"
        ),
      sleep_s     = 0.50,
      label       = paste0("EuropePMC searchPOST (", source_tag, ")"),
      expect_json = TRUE
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
  if (isTRUE(is_retro)) {
    y_start <- as.integer(substr(since_d, 1L, 4L))
    y_end   <- as.integer(substr(until_d, 1L, 4L))
    all_out <- list()
    for (yr in seq(y_start, y_end)) {
      w_from <- if (yr == y_start) since_d else sprintf("%04d-01-01", yr)
      w_to   <- if (yr == y_end)   until_d else sprintf("%04d-12-31", yr)
      # Safeguard 1: Include FIRST_IDATE OR-union so Advance-Online Epubs whose print issue moved
      # to the next calendar year are captured; post-merge canonical first_pub_date drops prior-year leaks.
      q_yr   <- builder_fn(w_from, w_to, include_idate = TRUE, verify_pmids = NULL)
      yr_res <- fetch_epmc_post(q_yr, source_tag = source_tag, page_size = 100L, max_pages = 35L)
      all_out <- c(all_out, yr_res)
    }
    return(all_out)
  }
  q_inc <- builder_fn(since_d, until_d, include_idate = TRUE, verify_pmids = NULL)
  fetch_epmc_post(q_inc, source_tag = source_tag, page_size = 100L, max_pages = 30L)
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
        url <- sprintf("https://api.biorxiv.org/details/%s/%s/%s/%d?category=%s",
                       srv, since_d, until_d, cursor, cat_name)
        resp <- perform_polite_http(
          httr2::request(url) |> httr2::req_user_agent(UA_STRING) |> httr2::req_timeout(35),
          sleep_s     = 0.50,
          label       = "bioRxiv details API",
          expect_json = TRUE
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
            dt_v    <- clean_cell(item$date)
            out[[length(out) + 1L]] <- make_raw_rec(
              source_tag     = paste0("bioRxivAPI_", srv),
              doi            = item$doi %||% "",
              journal        = srv,
              pub_year       = substr(dt_v, 1L, 4L),
              first_pub_date = dt_v,
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
# Source 3: NHGRI-EBI GWAS Catalog REST API v2
#   Date-sorted /v2/publications + EFO-independent is_gwas_catalog_protein_trait()
#   verification + EFO_0004747 OR-fallback for bare gene-symbol protein traits.
# ---------------------------------------------------------------------------
NON_PROTEOMIC_BIOMARKER_REGEX <- paste0(
  "\\b(lipoprotein|apolipoprotein a|apolipoprotein b|c-reactive protein|serum albumin|",
  "total protein|urinary albumin|albumin-to-creatinine|glycoprotein acetyls|n-glycan|",
  "n-glycosyl|glycated hemoglobin|hemoglobin a1c)\\b"
)

count_verified_gwas_protein_traits <- function(pid, sleep_s = 0.50) {
  u_st <- sprintf("https://www.ebi.ac.uk/gwas/rest/api/v2/studies?pubmed_id=%s&page=0&size=50", pid)
  resp_s <- perform_polite_http(
    httr2::request(u_st) |> httr2::req_user_agent(UA_STRING) |> httr2::req_timeout(25),
    sleep_s     = sleep_s,
    label       = "GWASCatalog pubmed studies",
    expect_json = TRUE
  )
  if (is.null(resp_s) || httr2::resp_status(resp_s) != 200L) return(0L)
  js_s <- tryCatch(httr2::resp_body_json(resp_s, simplifyVector = FALSE), error = function(e) NULL)
  stds <- js_s$`_embedded`$studies %||% js_s$studies %||% list()
  if (length(stds) == 0L) return(0L)

  trs <- vapply(stds, function(s) clean_cell(s$disease_trait %||% ""), character(1))
  # Pre-2017 GWAS Catalog curation bundled entire multi-protein pGWAS studies into a single
  # pooled GCST study with plural umbrella traits ("Protein quantitative trait loci",
  # "Cerebrospinal fluid levels of Alzheimer's disease-related proteins", "Blood protein levels").
  trs_low <- tolower(ascii_fold(trs))
  if (any(grepl("^protein quantitative trait loci$|\\blevels of [a-z0-9' -]+ proteins\\b|^blood protein levels$|^plasma protein levels$|^serum protein levels$|^cerebrospinal fluid protein levels$", trs_low, perl = TRUE))) {
    return(max(5L, length(stds)))
  }
  if (length(stds) < 5L) return(0L)

  # Prong 3A: EFO-independent `disease_trait` protein verification (rescues uncurated deposits)
  keep_prot <- is_gwas_catalog_protein_trait(trs)
  n_dist_prot <- length(unique(tolower(trs[keep_prot])))
  tot_el <- as.integer(js_s$page$totalElements %||% length(stds))
  if (n_dist_prot >= 5L) {
    est_n <- if (sum(keep_prot) >= 0.8 * length(stds) && tot_el > length(stds)) {
      as.integer(round(tot_el * (sum(keep_prot) / length(stds))))
    } else {
      as.integer(n_dist_prot)
    }
    return(est_n)
  }

  # Prong 3B: EFO_0004747 OR-fallback when traits use bare protein/gene symbols without "protein"
  u_efo <- sprintf(
    "https://www.ebi.ac.uk/gwas/rest/api/v2/studies?pubmed_id=%s&efo_id=EFO_0004747&show_child_traits=true&page=0&size=50",
    pid
  )
  resp_e <- perform_polite_http(
    httr2::request(u_efo) |> httr2::req_user_agent(UA_STRING) |> httr2::req_timeout(25),
    sleep_s     = sleep_s,
    label       = "GWASCatalog EFO_0004747 fallback",
    expect_json = TRUE
  )
  if (is.null(resp_e) || httr2::resp_status(resp_e) != 200L) return(0L)
  js_e <- tryCatch(httr2::resp_body_json(resp_e, simplifyVector = FALSE), error = function(e) NULL)
  stds_e <- js_e$`_embedded`$studies %||% js_e$studies %||% list()
  if (length(stds_e) < 5L) return(0L)

  trs_e <- tolower(ascii_fold(vapply(stds_e, function(s) clean_cell(s$disease_trait %||% ""), character(1))))
  trs_e_keep <- trs_e[nzchar(trs_e) & !grepl(NON_PROTEOMIC_BIOMARKER_REGEX, trs_e, perl = TRUE)]
  n_dist_efo <- length(unique(trs_e_keep))
  if (n_dist_efo >= 5L) {
    tot_efo <- as.integer(js_e$page$totalElements %||% length(stds_e))
    return(if (length(trs_e_keep) >= 0.8 * length(stds_e) && tot_efo > length(stds_e)) tot_efo else n_dist_efo)
  }
  0L
}

query_gwas_catalog_v2 <- function(since_d, until_d, is_retro = FALSE,
                                  verify_ids = NULL, include_bulk = FALSE) {
  pmid_counts <- integer(0)
  window_pubs <- list()

  if (!is.null(verify_ids) && length(verify_ids) > 0L) {
    for (pid in verify_ids) {
      u_pub <- sprintf("https://www.ebi.ac.uk/gwas/rest/api/v2/publications?pubmed_id=%s", pid)
      r_pub <- perform_polite_http(
        httr2::request(u_pub) |> httr2::req_user_agent(UA_STRING) |> httr2::req_timeout(20),
        sleep_s     = 0.50,
        label       = "GWASCatalog verify pub",
        expect_json = TRUE
      )
      st_pub <- if (!is.null(r_pub)) httr2::resp_status(r_pub) else NA_integer_
      n_prot <- count_verified_gwas_protein_traits(pid, sleep_s = 0.50)
      message(sprintf("[Source 3: GWAS Catalog v2] verify PMID %s -> /v2/publications?pubmed_id HTTP %s, verified_protein_traits=%d",
                      pid, as.character(st_pub), n_prot))
      if (n_prot >= 5L) pmid_counts[pid] <- n_prot
    }
  } else {
    max_pub_pages <- if (isTRUE(is_retro)) 40L else 10L
    sort_dir      <- if (suppressWarnings(as.integer(substr(until_d, 1L, 4L))) <= 2018L) "asc" else "desc"
    for (pg in seq(0L, max_pub_pages - 1L)) {
      u_pubs <- sprintf(
        "https://www.ebi.ac.uk/gwas/rest/api/v2/publications?page=%d&size=500&sort=publication_date&direction=%s",
        pg, sort_dir
      )
      resp_p <- perform_polite_http(
        httr2::request(u_pubs) |> httr2::req_user_agent(UA_STRING) |> httr2::req_timeout(45),
        sleep_s     = 0.50,
        label       = "GWASCatalog publications",
        expect_json = TRUE
      )
      if (is.null(resp_p) || httr2::resp_status(resp_p) != 200L) break
      js_p <- tryCatch(httr2::resp_body_json(resp_p, simplifyVector = FALSE), error = function(e) NULL)
      pubs <- js_p$`_embedded`$publications %||% js_p$publications %||% list()
      if (length(pubs) == 0L) break

      page_dts <- character(0)
      for (p in pubs) {
        pid <- normalize_pmid_one(p$pubmed_id %||% "")
        pdt <- clean_cell(p$publication_date %||% "")
        ttl <- clean_cell(p$title %||% "")
        if (nzchar(pdt)) page_dts <- c(page_dts, pdt)
        if (!nzchar(pid) || !nzchar(pdt)) next
        if (pdt >= since_d && pdt <= until_d) {
          window_pubs[[pid]] <- list(pmid = pid, pub_date = pdt, title = ttl)
        }
      }

      tot_pg <- as.integer(js_p$page$totalPages %||% 1L)
      if (pg + 1L >= tot_pg) break
      if (identical(sort_dir, "desc") && length(page_dts) > 0L && max(page_dts) < since_d) break
      if (identical(sort_dir, "asc")  && length(page_dts) > 0L && min(page_dts) > until_d) break
    }

    cand_pat <- paste0(
      BROAD_PROTEIN_REGEX, "|",
      "\\b(multi-?omics?|multi-?trait|multitrait|phenome-wide|biobank|atlas|quantitative traits?|",
      "mixed-model|variational inference|whole-genome sequencing|exome sequencing|pan-uk|scalable|pleiotrop|",
      "heritability|adipokin[a-z]*|adipocytokin[a-z]*|placental|pregnancy|maternal|",
      "blood|plasma|serum|cerebrospinal|csf|urine|urinary)\\b"
    )
    cand_pids <- names(Filter(function(x) {
      grepl(cand_pat, tolower(ascii_fold(x$title)), perl = TRUE)
    }, window_pubs))

    message(sprintf(
      "[Source 3: GWAS Catalog v2] %d publication(s) in [%s..%s]; verifying %d candidate PMID(s) for >=5 distinct protein traits.",
      length(window_pubs), since_d, until_date, length(cand_pids)
    ))

    for (pid in cand_pids) {
      n_prot <- count_verified_gwas_protein_traits(pid, sleep_s = 0.50)
      if (n_prot >= 5L) pmid_counts[pid] <- n_prot
    }

    # Optional streaming bulk TSV fallback for multi-decade retrospective runs
    if (isTRUE(is_retro) && isTRUE(include_bulk) && requireNamespace("data.table", quietly = TRUE)) {
      bulk_url <- "https://www.ebi.ac.uk/gwas/api/search/downloads/studies_alternative"
      r_b <- perform_polite_http(
        httr2::request(bulk_url) |> httr2::req_user_agent(UA_STRING) |> httr2::req_timeout(60),
        sleep_s = 0.50,
        label   = "GWASCatalog bulk TSV"
      )
      if (!is.null(r_b) && httr2::resp_status(r_b) == 200L &&
          !grepl("html", httr2::resp_content_type(r_b) %||% "", ignore.case = TRUE)) {
        tf <- tempfile(fileext = ".tsv")
        writeBin(httr2::resp_body_raw(r_b), tf)
        dt <- tryCatch(
          data.table::fread(tf, sep = "\t", quote = "", select = c("PUBMEDID", "DATE", "DISEASE/TRAIT", "MAPPED_TRAIT")),
          error = function(e) NULL
        )
        unlink(tf)
        if (!is.null(dt) && nrow(dt) > 0L) {
          keep <- is_gwas_catalog_protein_trait(paste(dt[["DISEASE/TRAIT"]], dt[["MAPPED_TRAIT"]]))
          sub_dt <- dt[keep & dt$DATE >= since_d & dt$DATE <= until_d]
          if (nrow(sub_dt) > 0L) {
            by_pid <- split(tolower(sub_dt[["DISEASE/TRAIT"]]), as.character(sub_dt$PUBMEDID))
            for (nm in names(by_pid)) {
              nd <- length(unique(by_pid[[nm]]))
              if (nd >= 5L) pmid_counts[nm] <- max(pmid_counts[nm] %||% 0L, as.integer(nd))
            }
          }
        }
      }
    }
  }

  qual_pmids <- names(pmid_counts[pmid_counts >= 5L])
  if (length(qual_pmids) == 0L) {
    message("[Source 3: GWAS Catalog v2] 0 multi-protein PMIDs matched in window.")
    return(list())
  }
  message(sprintf("[Source 3: GWAS Catalog v2] %d multi-protein PMID(s) (>=5 distinct protein traits); hydrating via efetch.", length(qual_pmids)))
  recs <- fetch_pubmed_efetch_xml(qual_pmids, source_tag = "GWASCatalog_v2")
  for (i in seq_along(recs)) {
    pid <- recs[[i]]$pmid
    recs[[i]]$gwas_study_count <- as.integer(pmid_counts[pid] %||% 0L)
    if (!nzchar(recs[[i]]$first_pub_date) && !is.null(window_pubs[[pid]]$pub_date)) {
      recs[[i]]$first_pub_date <- window_pubs[[pid]]$pub_date
      recs[[i]]$pub_year       <- substr(window_pubs[[pid]]$pub_date, 1L, 4L)
    }
  }
  recs
}

# ---------------------------------------------------------------------------
# Source 4: Blinded Citation Graph Snowballing (10 Landmark Seeds + Date Window)
# ---------------------------------------------------------------------------
query_snowball_backward_epmc <- function(seed_pmids, since_d = NULL, until_d = NULL) {
  if (!is.null(since_d) && nzchar(since_d)) {
    y_from <- suppressWarnings(as.integer(substr(since_d, 1L, 4L)))
    max_seed_yr <- max(LANDMARK_PGWAS_SEEDS$pub_year, na.rm = TRUE)
    if (!is.na(y_from) && y_from > max_seed_yr) {
      message(sprintf("[Source 4: Snowball Backward] skipped (since_year=%d > max_seed_year=%d).", y_from, max_seed_yr))
      return(list())
    }
  }

  y_min <- if (!is.null(since_d) && nzchar(since_d)) as.integer(substr(since_d, 1L, 4L)) - 1L else -Inf
  y_max <- if (!is.null(until_d) && nzchar(until_d)) as.integer(substr(until_d, 1L, 4L)) + 1L else Inf

  cited_med <- character(0)
  for (sp in seed_pmids) {
    u <- sprintf("https://www.ebi.ac.uk/europepmc/webservices/rest/MED/%s/references?page=1&pageSize=1000&format=json", sp)
    resp <- perform_polite_http(
      httr2::request(u) |> httr2::req_user_agent(UA_STRING) |> httr2::req_timeout(45),
      sleep_s     = 0.50,
      label       = "EuropePMC references",
      expect_json = TRUE
    )
    if (is.null(resp) || httr2::resp_status(resp) != 200L) next
    js <- tryCatch(httr2::resp_body_json(resp, simplifyVector = FALSE), error = function(e) NULL)
    refs <- js$referenceList$reference %||% list()
    for (rf in refs) {
      if (identical(clean_cell(rf$source), "MED") && nzchar(clean_cell(rf$id))) {
        ry <- suppressWarnings(as.integer(clean_cell(rf$pubYear %||% "")))
        if (!is.na(ry) && (ry < y_min || ry > y_max)) next
        cited_med <- c(cited_med, clean_cell(rf$id))
      }
    }
  }
  cited_med <- unique(cited_med)
  if (length(cited_med) == 0L) return(list())
  raw_refs <- fetch_pubmed_efetch_xml(cited_med, source_tag = "Snowball_Backward")
  Filter(function(r) {
    in_win <- TRUE
    if (!is.null(since_d) && nzchar(since_d) && !is.null(until_d) && nzchar(until_d)) {
      if (nzchar(r$first_pub_date)) {
        in_win <- (r$first_pub_date >= since_d && r$first_pub_date <= until_d)
      } else if (nzchar(r$pub_year)) {
        in_win <- (r$pub_year >= substr(since_d, 1L, 4L) && r$pub_year <= substr(until_d, 1L, 4L))
      }
    }
    comb <- tolower(paste(r$title, r$abstract, paste(r$mesh_headings, collapse = " ")))
    in_win && grepl(BROAD_PROTEIN_REGEX, comb, perl = TRUE) && grepl(BROAD_GENETICS_REGEX, comb, perl = TRUE)
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
  message(sprintf("[Source 4: Snowball Forward CITES: (%d seeds)] matched %d record(s).", nrow(LANDMARK_PGWAS_SEEDS), length(fwd)))

  bwd <- if (isTRUE(is_retro) && is.null(verify_ids)) {
    res_b <- query_snowball_backward_epmc(LANDMARK_PGWAS_SEEDS$pmid, since_d = since_d, until_d = until_d)
    message(sprintf("[Source 4: Snowball Backward /references] matched %d in-window record(s).", length(res_b)))
    res_b
  } else {
    list()
  }
  c(fwd, bwd)
}

# ---------------------------------------------------------------------------
# Stage B: Connected-Components Deduplication, Preprint-to-Published &
#          Cross-Server Preprint Consolidation (Safeguards 2, 3, 4)
# ---------------------------------------------------------------------------
resolve_biorxiv_published_doi <- function(p_doi) {
  if (!grepl("^10\\.(1101|64898)/", p_doi)) return("")
  for (srv in c("biorxiv", "medrxiv")) {
    u <- sprintf("https://api.biorxiv.org/pubs/%s/%s/na/json", srv, p_doi)
    r <- perform_polite_http(
      httr2::request(u) |> httr2::req_user_agent(UA_STRING) |> httr2::req_timeout(25),
      sleep_s     = 0.50,
      max_tries   = 3L,
      label       = "bioRxiv pubs batch",
      expect_json = TRUE
    )
    if (!is.null(r) && httr2::resp_status(r) == 200L) {
      js <- tryCatch(httr2::resp_body_json(r, simplifyVector = FALSE), error = function(e) NULL)
      if (!is.null(js$collection) && length(js$collection) > 0L) {
        cand_d <- normalize_doi_one(js$collection[[1L]]$published_doi %||% "")
        if (nzchar(cand_d)) return(cand_d)
      }
    }
  }
  ""
}

is_preprint_record <- function(r) {
  any(is_preprint_doi(r$doi)) ||
    r$source_tag %in% c("EuropePMC_PPR", "bioRxivAPI_biorxiv", "bioRxivAPI_medrxiv") ||
    (grepl("^PPR[0-9]+$", r$epmc_id %||% "") && !nzchar(r$pmid %||% "") && !nzchar(r$doi %||% "")) ||
    any(tolower(r$pub_types) == "preprint")
}

merge_raw_records <- function(recs, resolve_remote_pubs = TRUE) {
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

  # Pass 1: Exact IDs (using `doi_dedup_key` so Preprints.org .v1/.v2 merge while preserving .v1 in `doi`)
  # + title fingerprint + explicit `preprint_links` graph edges
  key_map <- new.env(hash = TRUE, parent = emptyenv())
  bind_key <- function(k, idx) {
    if (!nzchar(k)) return()
    if (exists(k, envir = key_map, inherits = FALSE)) {
      union_nodes(get(k, envir = key_map, inherits = FALSE), idx)
    } else {
      assign(k, idx, envir = key_map)
    }
  }

  for (i in seq_len(n)) {
    r <- recs[[i]]
    dk  <- doi_dedup_key(r$doi)
    tfp <- build_title_fingerprint(r$first_author, r$pub_year, r$title)
    keys <- c(
      if (nzchar(dk))        paste0("doi:", dk)        else "",
      if (nzchar(r$pmid))    paste0("pmid:", r$pmid)   else "",
      if (nzchar(r$pmcid))   paste0("pmcid:", r$pmcid) else "",
      if (nzchar(r$epmc_id)) paste0("epmc:", r$epmc_id) else "",
      if (nzchar(tfp))       paste0("tfp:", tfp)       else ""
    )
    for (pl in r$preprint_links) {
      if (grepl("^10\\.", pl)) {
        d_pl <- doi_dedup_key(pl)
        if (nzchar(d_pl)) keys <- c(keys, paste0("doi:", d_pl))
      } else if (grepl("^PMID:[0-9]+$", pl, ignore.case = TRUE)) {
        p_pl <- normalize_pmid_one(sub("^PMID:", "", pl, ignore.case = TRUE))
        if (nzchar(p_pl)) keys <- c(keys, paste0("pmid:", p_pl))
      } else if (grepl("^PPR[0-9]+$", pl)) {
        keys <- c(keys, paste0("epmc:", pl))
      }
    }
    for (k in unique(keys[nzchar(keys)])) bind_key(k, i)
  }

  # Pass 2A: Resolve bioRxiv/medRxiv /pubs/ links on standalone actionable preprints
  if (isTRUE(resolve_remote_pubs)) {
    roots_p1 <- vapply(seq_len(n), find_root, integer(1))
    groups_p1 <- split(seq_len(n), roots_p1)
    for (idxs in groups_p1) {
      sub <- recs[idxs]
      if (!all(vapply(sub, is_preprint_record, logical(1)))) next
      dois <- unique(vapply(sub, `[[`, "", "doi"))
      dois <- dois[grepl("^10\\.(1101|64898)/", dois)]
      if (length(dois) == 0L) next
      sc_tmp <- score_triage_v2(sub[[1L]]$title, sub[[1L]]$abstract, sub[[1L]]$pub_types, sub[[1L]]$mesh_headings, 0L, doi = dois[1L])
      if (!grepl("^Auto-Deprioritised", sc_tmp$classification)) {
        pub_d <- resolve_biorxiv_published_doi(dois[1L])
        if (!nzchar(pub_d)) {
          u_cr <- paste0("https://api.crossref.org/works/", utils::URLencode(dois[1L], reserved = TRUE))
          r_cr <- perform_polite_http(
            httr2::request(u_cr) |> httr2::req_user_agent(UA_STRING) |> httr2::req_timeout(20),
            sleep_s     = 0.45,
            max_tries   = 3L,
            label       = "Crossref batch preprint->pub",
            expect_json = TRUE
          )
          if (!is.null(r_cr) && httr2::resp_status(r_cr) == 200L) {
            js_cr <- tryCatch(httr2::resp_body_json(r_cr, simplifyVector = FALSE), error = function(e) NULL)
            rel_cr <- js_cr$message$relation[["is-preprint-of"]]
            if (!is.null(rel_cr) && length(rel_cr) > 0L) {
              pub_d <- normalize_doi_one(rel_cr[[1L]]$id %||% "")
            }
          }
        }
        if (nzchar(pub_d)) {
          recs[[idxs[1L]]]$preprint_links <- unique(c(recs[[idxs[1L]]]$preprint_links, pub_d, dois[1L]))
          k_pub <- paste0("doi:", doi_dedup_key(pub_d))
          if (exists(k_pub, envir = key_map, inherits = FALSE)) {
            union_nodes(get(k_pub, envir = key_map, inherits = FALSE), idxs[1L])
          }
        }
      }
    }
  }

  # Pass 2B: Fuzzy cross-server preprint & preprint->published consolidation
  # Merge Safeguard 4: MR_BOILERPLATE_STOPWORDS stripped in `title_content_tokens()`,
  # `has_compatible_phenotype()` enforced, and stricter author guard applied!
  roots_p2 <- vapply(seq_len(n), find_root, integer(1))
  groups_p2 <- split(seq_len(n), roots_p2)
  g_len <- length(groups_p2)
  if (g_len > 1L) {
    g_meta <- lapply(groups_p2, function(idxs) {
      sub <- recs[idxs]
      is_p <- all(vapply(sub, is_preprint_record, logical(1)))
      t_best <- sub[[which.max(nchar(vapply(sub, `[[`, "", "title")))]]$title
      a_best <- tolower(ascii_fold(sub[[which.max(nchar(vapply(sub, `[[`, "", "abstract")))]]$abstract %||% ""))
      p_phr  <- regmatches(a_best, regexpr("\\b[1-9][0-9,]+\\s+[a-z-]+\\s+proteins\\b", a_best, perl = TRUE))
      au1 <- unique(vapply(sub, function(x) normalize_author_surname(x$first_author), character(1)))
      au2 <- unique(vapply(sub, function(x) normalize_author_surname(x$last_author), character(1)))
      list(
        idx = idxs[1L],
        is_prep = is_p,
        title = t_best,
        panel = if (length(p_phr) == 1L) p_phr else "",
        au1 = au1[nzchar(au1)],
        au2 = au2[nzchar(au2)],
        toks = title_content_tokens(t_best, strip_mr_boilerplate = TRUE)
      )
    })
    prep_g <- which(vapply(g_meta, `[[`, logical(1), "is_prep"))
    for (pi in prep_g) {
      m1 <- g_meta[[pi]]
      if (length(m1$toks) < 3L) next
      for (qj in seq_len(g_len)) {
        if (pi == qj) next
        m2 <- g_meta[[qj]]
        if (length(m2$toks) < 3L) next
        if (m2$is_prep && qj < pi) next
        if (find_root(m1$idx) == find_root(m2$idx)) next
        inter_n <- length(intersect(m1$toks, m2$toks))
        if (inter_n < 3L) next
        jac <- inter_n / length(union(m1$toks, m2$toks))
        same_first <- length(intersect(m1$au1, m2$au1)) > 0L
        same_last  <- length(intersect(m1$au2, m2$au2)) > 0L
        any_au     <- length(intersect(c(m1$au1, m1$au2), c(m2$au1, m2$au2))) > 0L
        same_panel <- (m1$is_prep != m2$is_prep) && nzchar(m1$panel) && identical(m1$panel, m2$panel)
        if (((m1$is_prep && m2$is_prep && jac >= 0.85) ||
             (any_au && jac >= 0.65) ||
             (same_first && same_last && jac >= 0.40) ||
             (same_first && same_panel && jac >= 0.28)) &&
            has_compatible_phenotype(m1$title, m2$title)) {
          union_nodes(m1$idx, m2$idx)
        }
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
    is_prep_vec <- vapply(sub, is_preprint_record, logical(1))
    pub_sub  <- sub[!is_prep_vec]
    prep_sub <- sub[is_prep_vec]
    # Safeguards 2 & 5: Prefer peer-reviewed journal records (`pub_sub`) when present;
    # when consolidating cross-server preprints (`pub_sub` empty), prioritize permanent
    # openRxiv/bioRxiv/medRxiv DOIs (`10.1101/`, `10.64898/`) and earliest posting date!
    pref_sub <- if (length(pub_sub) > 0L) {
      pub_sub
    } else if (length(sub) > 1L) {
      s_dois <- vapply(sub, `[[`, "", "doi")
      s_dts  <- vapply(sub, `[[`, "", "first_pub_date")
      s_pmid <- vapply(sub, `[[`, "", "pmid")
      ord_p  <- order(
        !nzchar(s_dois),
        !grepl("^10\\.(1101|64898)/", s_dois),
        ifelse(nzchar(s_dts), s_dts, "9999-99-99"),
        !nzchar(s_pmid)
      )
      sub[ord_p]
    } else {
      sub
    }

    chosen_pmid <- pick_nonempty(vapply(pref_sub, `[[`, "", "pmid"))
    chosen_doi  <- pick_nonempty(vapply(pref_sub, `[[`, "", "doi"))

    merged_prep_links <- unique(unlist(lapply(sub, `[[`, "preprint_links")))
    if (length(pub_sub) > 0L && length(prep_sub) > 0L) {
      p_dois <- unique(vapply(prep_sub, `[[`, "", "doi"))
      p_epmc <- unique(vapply(prep_sub, `[[`, "", "epmc_id"))
      p_pmid <- unique(vapply(prep_sub, `[[`, "", "pmid"))
      merged_prep_links <- unique(c(
        merged_prep_links,
        p_dois[nzchar(p_dois)],
        p_epmc[grepl("^PPR", p_epmc)],
        ifelse(nzchar(p_pmid), paste0("PMID:", p_pmid), "")
      ))
    } else if (length(prep_sub) > 1L) {
      other_dois  <- setdiff(unique(vapply(prep_sub, `[[`, "", "doi")), chosen_doi)
      other_pmids <- setdiff(unique(vapply(prep_sub, `[[`, "", "pmid")), chosen_pmid)
      merged_prep_links <- unique(c(
        merged_prep_links,
        other_dois[nzchar(other_dois)],
        ifelse(nzchar(other_pmids), paste0("PMID:", other_pmids), "")
      ))
    }
    # Safeguard 5: Strip self-referential own-DOI/own-PMID tokens from `merged_prep_links`
    self_refs <- c(
      chosen_doi,
      doi_dedup_key(chosen_doi),
      if (nzchar(chosen_pmid)) paste0("PMID:", chosen_pmid) else ""
    )
    self_refs <- self_refs[nzchar(self_refs)]
    merged_prep_links <- setdiff(merged_prep_links, self_refs)

    # Canonical `first_pub_date` from `pref_sub`: peer-reviewed journal date (`epub` > `ppub`)
    # takes precedence over superseded preprint dates, while standalone preprints use `extract_preprint_doi_date()`!
    pref_dates <- vapply(pref_sub, `[[`, "", "first_pub_date")
    pref_dates <- pref_dates[nzchar(pref_dates)]
    canon_dt   <- if (length(pref_dates) > 0L) min(pref_dates) else ""
    canon_yr   <- if (nzchar(canon_dt)) substr(canon_dt, 1L, 4L) else pick_nonempty(vapply(pref_sub, `[[`, "", "pub_year"))

    jrns <- vapply(pref_sub, `[[`, "", "journal")
    jrns_clean <- jrns[nzchar(jrns) & !(jrns %in% c("MED", "PMC", "Preprint"))]
    best_jrn <- if (length(jrns_clean) > 0L) jrns_clean[1L] else pick_nonempty(c(jrns, vapply(sub, `[[`, "", "journal")))

    all_pts <- unique(unlist(lapply(pref_sub, `[[`, "pub_types")))
    if (length(pub_sub) > 0L) {
      all_pts <- setdiff(all_pts, c("Preprint", "preprint"))
    }

    best_abs <- if (any(nzchar(vapply(pref_sub, `[[`, "", "abstract")))) {
      pick_longest(vapply(pref_sub, `[[`, "", "abstract"))
    } else {
      pick_longest(vapply(sub, `[[`, "", "abstract"))
    }
    if (isTRUE(resolve_remote_pubs) && !nzchar(best_abs) && length(pub_sub) == 0L && grepl("^10\\.1101/", chosen_doi)) {
      for (srv_h in c("biorxiv", "medrxiv")) {
        u_bh <- sprintf("https://api.biorxiv.org/details/%s/%s/na/json", srv_h, chosen_doi)
        r_bh <- perform_polite_http(
          httr2::request(u_bh) |> httr2::req_user_agent(UA_STRING) |> httr2::req_timeout(15),
          sleep_s = 0.25, max_tries = 2L, label = "bioRxiv merge abstract hydrate", expect_json = TRUE
        )
        if (!is.null(r_bh) && httr2::resp_status(r_bh) == 200L) {
          jbh <- tryCatch(httr2::resp_body_json(r_bh, simplifyVector = FALSE), error = function(e) NULL)
          if (!is.null(jbh$collection) && length(jbh$collection) > 0L) {
            best_abs <- clean_cell(jbh$collection[[1L]]$abstract %||% "")
            if (nzchar(best_abs)) break
          }
        }
      }
    }

    list(
      matched_sources  = paste(unique(vapply(sub, `[[`, "", "source_tag")), collapse = ";"),
      pmid             = chosen_pmid,
      doi              = chosen_doi,
      pmcid            = pick_nonempty(c(vapply(pref_sub, `[[`, "", "pmcid"), vapply(sub, `[[`, "", "pmcid"))),
      epmc_id          = pick_nonempty(c(vapply(pref_sub, `[[`, "", "epmc_id"), vapply(sub, `[[`, "", "epmc_id"))),
      journal          = best_jrn,
      pub_year         = canon_yr,
      first_pub_date   = canon_dt,
      first_author     = pick_nonempty(vapply(pref_sub, `[[`, "", "first_author")),
      last_author      = pick_nonempty(vapply(pref_sub, `[[`, "", "last_author")),
      title            = pick_longest(vapply(pref_sub, `[[`, "", "title")),
      abstract         = best_abs,
      pub_types        = all_pts,
      mesh_headings    = unique(unlist(lapply(sub, `[[`, "mesh_headings"))),
      preprint_links   = merged_prep_links[nzchar(merged_prep_links)],
      erratum_targets  = unique(unlist(lapply(sub, `[[`, "erratum_targets"))),
      gwas_study_count = max(vapply(sub, `[[`, 0L, "gwas_study_count"))
    )
  })
}

# ---------------------------------------------------------------------------
# Execute Multi-Resource Discovery
# ---------------------------------------------------------------------------
recs_s1 <- query_pubmed_v2(q_pubmed)

recs_s2a_1 <- run_epmc_windowed(build_epmc_prong2a_platforms_query, "EuropePMC", since_date, until_date,
                                is_retro = identical(mode, "retrospective"),
                                verify_ids = if (is_verify) verify_pmids else NULL)
recs_s2a_2 <- run_epmc_windowed(build_epmc_prong2a_proteome_query, "EuropePMC", since_date, until_date,
                                is_retro = identical(mode, "retrospective"),
                                verify_ids = if (is_verify) verify_pmids else NULL)

# Rescue preprints first posted in `[since_d..until_d]` via Crossref `type:posted-content`
# (immune to Europe PMC overwriting `FIRST_PDATE` when `v2`/`v3` revisions are uploaded in `yr+1`)
query_crossref_posted_preprints <- function(since_d, until_d, is_retro, verify_ids = NULL) {
  if (!isTRUE(is_retro) || (!is.null(verify_ids) && length(verify_ids) > 0L)) return(list())
  if (is.null(since_d) || is.null(until_d) || !nzchar(since_d) || !nzchar(until_d)) return(list())
  out <- list()
  seen_d <- character(0)
  cr_queries <- c("pQTL proteome", "SomaScan GWAS", "Olink GWAS", "aptamer pQTL")
  for (cq in cr_queries) {
    u_cr <- paste0(
      "https://api.crossref.org/works?filter=type:posted-content,from-posted-date:", since_d,
      ",until-posted-date:", until_d, "&query=", utils::URLencode(cq, reserved = TRUE), "&rows=50"
    )
    r_cr <- perform_polite_http(
      httr2::request(u_cr) |> httr2::req_user_agent(UA_STRING) |> httr2::req_timeout(20),
      sleep_s     = 0.35,
      max_tries   = 3L,
      label       = "Crossref posted-content preprints",
      expect_json = TRUE
    )
    if (is.null(r_cr) || httr2::resp_status(r_cr) != 200L) next
    js_cr <- tryCatch(httr2::resp_body_json(r_cr, simplifyVector = FALSE), error = function(e) NULL)
    items <- js_cr$message$items %||% list()
    for (it in items) {
      d_raw <- normalize_doi_one(it$DOI %||% "")
      if (!nzchar(d_raw) || d_raw %in% seen_d || !is_preprint_doi(d_raw)) next
      seen_d <- c(seen_d, d_raw)
      ttl <- clean_cell(if (length(it$title) > 0L) it$title[[1L]] else "")
      abs <- strip_html_tags(it$abstract %||% "")
      # Hydrate missing abstract on 10.1101/ preprints from bioRxiv/medRxiv details API
      if (!nzchar(abs) && grepl("^10\\.1101/", d_raw)) {
        for (srv in c("biorxiv", "medrxiv")) {
          u_bx <- sprintf("https://api.biorxiv.org/details/%s/%s/na/json", srv, d_raw)
          r_bx <- perform_polite_http(
            httr2::request(u_bx) |> httr2::req_user_agent(UA_STRING) |> httr2::req_timeout(15),
            sleep_s = 0.25, max_tries = 2L, label = "bioRxiv details hydrate", expect_json = TRUE
          )
          if (!is.null(r_bx) && httr2::resp_status(r_bx) == 200L) {
            jb <- tryCatch(httr2::resp_body_json(r_bx, simplifyVector = FALSE), error = function(e) NULL)
            if (!is.null(jb$collection) && length(jb$collection) > 0L) {
              abs <- clean_cell(jb$collection[[1L]]$abstract %||% "")
              if (!nzchar(ttl)) ttl <- clean_cell(jb$collection[[1L]]$title %||% "")
              break
            }
          }
        }
      }
      comb_lc <- tolower(paste(ttl, abs))
      if (!grepl(BROAD_PROTEIN_REGEX, comb_lc, perl = TRUE) || !grepl(BROAD_GENETICS_REGEX, comb_lc, perl = TRUE)) next
      au_l <- it$author %||% list()
      f_au <- if (length(au_l) > 0L) clean_cell(au_l[[1L]]$family %||% "") else ""
      l_au <- if (length(au_l) > 0L) clean_cell(au_l[[length(au_l)]]$family %||% f_au) else ""
      dp_parts <- tryCatch(unlist(it$posted$`date-parts`[[1L]]), error = function(e) integer(0))
      dt_str <- if (length(dp_parts) >= 3L) {
        sprintf("%04d-%02d-%02d", as.integer(dp_parts[1]), as.integer(dp_parts[2]), as.integer(dp_parts[3]))
      } else if (nzchar(extract_preprint_doi_date(d_raw))) {
        extract_preprint_doi_date(d_raw)
      } else {
        sprintf("%s-01-01", substr(since_d, 1L, 4L))
      }
      rel_pub <- it$relation[["is-preprint-of"]]
      pub_lnk <- if (!is.null(rel_pub) && length(rel_pub) > 0L) normalize_doi_one(rel_pub[[1L]]$id %||% "") else ""
      out[[length(out) + 1L]] <- make_raw_rec(
        source_tag     = "EuropePMC_PPR",
        doi            = d_raw,
        journal        = clean_cell(it$`group-title` %||% it$institution[[1L]]$name %||% "Preprint"),
        pub_year       = substr(dt_str, 1L, 4L),
        first_pub_date = dt_str,
        first_author   = f_au,
        last_author    = l_au,
        title          = ttl,
        abstract       = abs,
        preprint_links = if (nzchar(pub_lnk)) pub_lnk else character(0)
      )
    }
  }
  out
}
recs_s2a_rev <- query_crossref_posted_preprints(since_date, until_date,
                                                is_retro = identical(mode, "retrospective"),
                                                verify_ids = if (is_verify) verify_pmids else NULL)
recs_s2a   <- c(recs_s2a_1, recs_s2a_2, recs_s2a_rev)
message(sprintf("[Source 2A: Europe PMC MED+PMC+PPR (Platforms=%d + Proteome=%d + RevisedPPR=%d)] matched %d raw record(s).",
                length(recs_s2a_1), length(recs_s2a_2), length(recs_s2a_rev), length(recs_s2a)))

recs_s2b <- run_epmc_windowed(build_epmc_prong2b_query, "EuropePMC_MethodsRescue", since_date, until_date,
                              is_retro = identical(mode, "retrospective"),
                              verify_ids = if (is_verify) verify_pmids else NULL)

# Fast NCBI PMC full-text rescue backup (resilient to transient Europe PMC full-text search timeouts;
# note `[pt]` does not exist in `db=pmc` and must NOT be appended here — PubMed efetch handles pub_types!)
query_pmc_methods_rescue_v2 <- function(since_d, until_d, verify_ids = NULL) {
  if (!is.null(verify_ids) && length(verify_ids) > 0L) return(list())
  if (is.null(since_d) || is.null(until_d) || !nzchar(since_d) || !nzchar(until_d)) return(list())
  d_from <- gsub("-", "/", since_d)
  y_to_plus1 <- suppressWarnings(as.integer(substr(until_d, 1L, 4L))) + 1L
  d_to_pmc   <- if (!is.na(y_to_plus1)) sprintf("%d/12/31", y_to_plus1) else gsub("-", "/", until_d)
  q_pmc <- paste0(
    '("', d_from, '"[PubDate] : "', d_to_pmc, '"[PubDate]) AND ',
    '("SomaScan"[Body - All Words] OR "SOMAscan"[Body - All Words] OR "SomaLogic"[Body - All Words] OR "SOMAmer"[Body - All Words] OR ',
    '"aptamer-based proteomics"[Body - All Words] OR "aptamer-based proteomic"[Body - All Words] OR "Olink"[Body - All Words] OR "proximity extension assay"[Body - All Words] OR ',
    '"NULISA"[Body - All Words] OR "Proteograph"[Body - All Words] OR "SWATH-MS"[Body - All Words] OR "DIA-MS"[Body - All Words] OR "OpenSWATH"[Body - All Words]) AND ',
    '("pQTL"[Body - All Words] OR "pQTLs"[Body - All Words] OR "cis-pQTL"[Body - All Words] OR "cis-pQTLs"[Body - All Words] OR ',
    '"protein QTL"[Body - All Words] OR "protein QTLs"[Body - All Words] OR "protein quantitative trait loci"[Body - All Words] OR ',
    '"GWAF"[Body - All Words] OR "seqMeta"[Body - All Words] OR "MatrixEQTL"[Body - All Words] OR "tensorQTL"[Body - All Words] OR ',
    '"meta-analysis of genome-wide association"[Body - All Words])'
  )
  ncbi_key <- Sys.getenv("NCBI_API_KEY", unset = "")
  req_es <- httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi") |>
    httr2::req_user_agent(UA_STRING) |>
    httr2::req_timeout(30) |>
    httr2::req_body_form(
      db      = "pmc",
      term    = q_pmc,
      retmax  = "2000",
      retmode = "json",
      api_key = if (nzchar(ncbi_key)) ncbi_key else NULL
    )
  r_es <- perform_polite_http(req_es, sleep_s = 0.45, max_tries = 3L, label = "NCBI PMC MethodsRescue esearch", expect_json = TRUE)
  if (is.null(r_es) || httr2::resp_status(r_es) != 200L) return(list())
  js_es <- tryCatch(httr2::resp_body_json(r_es, simplifyVector = FALSE), error = function(e) NULL)
  pmc_ids <- as.character(unlist(js_es$esearchresult$idlist %||% list()))
  if (length(pmc_ids) == 0L) return(list())

  pmids_out <- character(0)
  for (chunk in split(pmc_ids, ceiling(seq_along(pmc_ids) / 150L))) {
    req_sm <- httr2::request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi") |>
      httr2::req_user_agent(UA_STRING) |>
      httr2::req_timeout(30) |>
      httr2::req_body_form(
        db      = "pmc",
        id      = paste(chunk, collapse = ","),
        retmode = "json",
        api_key = if (nzchar(ncbi_key)) ncbi_key else NULL
      )
    r_sm <- perform_polite_http(req_sm, sleep_s = 0.45, max_tries = 3L, label = "NCBI PMC MethodsRescue esummary", expect_json = TRUE)
    if (is.null(r_sm) || httr2::resp_status(r_sm) != 200L) next
    js_sm <- tryCatch(httr2::resp_body_json(r_sm, simplifyVector = FALSE)$result, error = function(e) NULL)
    if (is.null(js_sm)) next
    for (cid in chunk) {
      it <- js_sm[[cid]]
      if (is.null(it)) next
      for (aid in it$articleids %||% list()) {
        if (identical(aid$idtype, "pmid") && nzchar(aid$value %||% "")) {
          pmids_out <- c(pmids_out, as.character(aid$value))
        }
      }
    }
  }
  pmids_out <- unique(pmids_out[nzchar(pmids_out)])
  if (length(pmids_out) == 0L) return(list())
  recs_pmc_ef <- fetch_pubmed_efetch_xml(pmids_out, source_tag = "EuropePMC_MethodsRescue")
  Filter(function(r) !nzchar(r$first_pub_date) || r$first_pub_date <= until_d, recs_pmc_ef)
}

recs_s2b_pmc <- query_pmc_methods_rescue_v2(since_date, until_date, verify_ids = if (is_verify) verify_pmids else NULL)
recs_s2b     <- c(recs_s2b, recs_s2b_pmc)
message(sprintf("[Source 2B: Europe PMC + NCBI PMC METHODS Rescue] matched %d record(s).", length(recs_s2b)))

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

# ---------------------------------------------------------------------------
# 4-Tier Preprint-to-Published Upgrade Detector & Pre-Merge Hydration (v4.9)
#   Resolves active catalogued preprints via bioRxiv/medRxiv /pubs/, EuropePMC
#   commentCorrectionList, and Crossref is-preprint-of BEFORE merge_raw_records(),
#   hydrating any newly published journal counterparts via PubMed efetch so
#   generalized journal titles ("blood molecular phenotypes", Vinuela 2023) or
#   transient EuropePMC gateway errors cannot drop a published preprint upgrade.
# ---------------------------------------------------------------------------
resolve_ppr_to_doi <- function(ppr_id) {
  if (!grepl("^PPR[0-9]+$", ppr_id)) return("")
  u <- sprintf("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:%s%%20AND%%20SRC:PPR&resultType=lite&format=json", ppr_id)
  resp <- perform_polite_http(
    httr2::request(u) |> httr2::req_user_agent(UA_STRING) |> httr2::req_timeout(25),
    sleep_s     = 0.50,
    max_tries   = 3L,
    label       = "EuropePMC resolve PPR",
    expect_json = TRUE
  )
  if (is.null(resp) || httr2::resp_status(resp) != 200L) return("")
  js <- tryCatch(httr2::resp_body_json(resp, simplifyVector = FALSE), error = function(e) NULL)
  res <- js$resultList$result
  if (is.null(res) || length(res) == 0L) return("")
  normalize_doi_one(res[[1L]]$doi %||% "")
}

resolve_pmid_to_doi <- function(pmid_tok) {
  pid <- sub("^PMID:", "", trimws(pmid_tok), ignore.case = TRUE)
  if (!grepl("^[0-9]+$", pid)) return("")
  ncbi_key <- Sys.getenv("NCBI_API_KEY", unset = "")
  u <- sprintf(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=%s&retmode=json%s",
    pid, if (nzchar(ncbi_key)) paste0("&api_key=", ncbi_key) else ""
  )
  resp <- perform_polite_http(
    httr2::request(u) |> httr2::req_user_agent(UA_STRING) |> httr2::req_timeout(20),
    sleep_s     = 0.40,
    max_tries   = 3L,
    label       = "NCBI esummary resolve PMID->DOI",
    expect_json = TRUE
  )
  if (is.null(resp) || httr2::resp_status(resp) != 200L) return("")
  js <- tryCatch(httr2::resp_body_json(resp, simplifyVector = FALSE), error = function(e) NULL)
  aids <- js$result[[pid]]$articleids %||% list()
  for (aid in aids) {
    if (identical(tolower(aid$idtype %||% ""), "doi") && nzchar(aid$value %||% "")) {
      return(normalize_doi_one(aid$value))
    }
  }
  ""
}

resolve_doi_to_pmid <- function(doi_str) {
  d_clean <- normalize_doi_one(doi_str)
  if (!nzchar(d_clean)) return("")
  ncbi_key <- Sys.getenv("NCBI_API_KEY", unset = "")
  u <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=",
    utils::URLencode(paste0('"', d_clean, '"[DOI]'), reserved = TRUE),
    "&retmode=json", if (nzchar(ncbi_key)) paste0("&api_key=", ncbi_key) else ""
  )
  resp <- perform_polite_http(
    httr2::request(u) |> httr2::req_user_agent(UA_STRING) |> httr2::req_timeout(20),
    sleep_s     = 0.35,
    max_tries   = 3L,
    label       = "NCBI esearch resolve DOI->PMID",
    expect_json = TRUE
  )
  if (is.null(resp) || httr2::resp_status(resp) != 200L) return("")
  js <- tryCatch(httr2::resp_body_json(resp, simplifyVector = FALSE), error = function(e) NULL)
  ids <- as.character(unlist(js$esearchresult$idlist %||% list()))
  if (length(ids) > 0L) normalize_pmid_one(ids[1L]) else ""
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
    if (grepl("^10\\.(1101|64898)/", p_doi)) {
      for (srv in c("biorxiv", "medrxiv")) {
        u <- sprintf("https://api.biorxiv.org/pubs/%s/%s/na/json", srv, p_doi)
        r <- perform_polite_http(
          httr2::request(u) |> httr2::req_user_agent(UA_STRING) |> httr2::req_timeout(25),
          sleep_s     = 0.40,
          max_tries   = 4L,
          label       = "bioRxiv pubs API",
          expect_json = TRUE
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
      u2 <- paste0("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:",
                   utils::URLencode(p_doi, reserved = TRUE), "&resultType=core&format=json")
      r2 <- perform_polite_http(
        httr2::request(u2) |> httr2::req_user_agent(UA_STRING) |> httr2::req_timeout(25),
        sleep_s     = 0.40,
        max_tries   = 3L,
        label       = "EuropePMC preprint upgrade",
        expect_json = TRUE
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
      u3 <- paste0("https://api.crossref.org/works/", utils::URLencode(p_doi, reserved = TRUE))
      r3 <- perform_polite_http(
        httr2::request(u3) |> httr2::req_user_agent(UA_STRING) |> httr2::req_timeout(25),
        sleep_s     = 0.40,
        max_tries   = 3L,
        label       = "Crossref preprint upgrade",
        expect_json = TRUE
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

    if (nzchar(pub_doi) && !nzchar(pub_pmid)) {
      pub_pmid <- resolve_doi_to_pmid(pub_doi)
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
if (!is_verify && length(direct_preprint_upgrades) > 0L) {
  existing_raw_pmids <- unique(vapply(all_raw_recs, function(r) normalize_pmid_one(r$pmid %||% ""), character(1)))
  up_pmids_to_fetch  <- character(0)
  up_pmid_to_pdoi    <- list()
  for (up in direct_preprint_upgrades) {
    if (nzchar(up$published_pmid)) {
      up_pmid_to_pdoi[[up$published_pmid]] <- up$preprint_doi
      if (!(up$published_pmid %in% existing_raw_pmids)) {
        up_pmids_to_fetch <- unique(c(up_pmids_to_fetch, up$published_pmid))
      }
    }
  }
  if (length(up_pmids_to_fetch) > 0L) {
    up_recs <- fetch_pubmed_efetch_xml(up_pmids_to_fetch, source_tag = "PreprintUpgrade")
    for (ui in seq_along(up_recs)) {
      up_pid <- normalize_pmid_one(up_recs[[ui]]$pmid %||% "")
      if (nzchar(up_pid) && !is.null(up_pmid_to_pdoi[[up_pid]])) {
        up_recs[[ui]]$preprint_links <- unique(c(up_recs[[ui]]$preprint_links, up_pmid_to_pdoi[[up_pid]]))
      }
    }
    all_raw_recs <- c(all_raw_recs, up_recs)
  }
}

if (identical(mode, "retrospective") && !is_verify && nzchar(until_date)) {
  # Drop future-window (`> until_date`) journal articles before `merge_raw_records()` so `yr+1`
  # peer-reviewed publications cannot absorb in-window (`<= until_date`) preprints prior to date filtering.
  all_raw_recs <- Filter(function(r) {
    is_preprint_record(r) || !nzchar(r$first_pub_date) || r$first_pub_date <= until_date
  }, all_raw_recs)
}
raw_n <- length(all_raw_recs)
message(sprintf("Total raw multi-source hits prior to graph deduplication: %d", raw_n))

merged_recs_all <- merge_raw_records(all_raw_recs, resolve_remote_pubs = !is_verify)
message(sprintf("Unique works after cross-source & preprint-lineage merge: %d", length(merged_recs_all)))

# Strict retrospective date-window enforcement on canonical post-merge `first_pub_date`
# (preserving unseen Q4-prior-year Epub articles whose journal print `[dp]` / MEDLINE `[edat]`
# fell inside `[since_date..until_date]` so late-year Epub -> Jan-print papers never fall between periods).
prior_seen_ids <- tolower(unique(trimws(c(
  vapply(studies$doi, normalize_doi_one, character(1)),
  vapply(studies$doi, doi_dedup_key, character(1)),
  vapply(studies$pmid, normalize_pmid_one, character(1)),
  vapply(excluded$identifier, normalize_doi_one, character(1)),
  extract_pmids_from_cell(paste(excluded$identifier, collapse = " ")),
  vapply(queue$identifier, normalize_doi_one, character(1)),
  extract_pmids_from_cell(paste(queue$identifier, collapse = " "))
))))
prior_seen_ids <- prior_seen_ids[nzchar(prior_seen_ids)]

if (identical(mode, "retrospective") && !is_verify && nzchar(since_date) && nzchar(until_date)) {
  y_s <- substr(since_date, 1L, 4L)
  y_u <- substr(until_date, 1L, 4L)
  q4_floor <- sprintf("%04d-10-01", as.integer(y_s) - 1L)
  merged_recs <- Filter(function(m) {
    if (nzchar(m$first_pub_date)) {
      if (m$first_pub_date >= since_date && m$first_pub_date <= until_date) return(TRUE)
      if (!is_preprint_doi(m$doi) && m$first_pub_date >= q4_floor && m$first_pub_date < since_date) {
        m_ids <- tolower(c(m$pmid, normalize_doi_one(m$doi), doi_dedup_key(m$doi)))
        m_ids <- m_ids[nzchar(m_ids)]
        if (length(m_ids) > 0L && !any(m_ids %in% prior_seen_ids)) return(TRUE)
      }
      FALSE
    } else if (nzchar(m$pub_year)) {
      m$pub_year >= y_s && m$pub_year <= y_u
    } else {
      TRUE
    }
  }, merged_recs_all)
  n_dropped_dt <- length(merged_recs_all) - length(merged_recs)
  message(sprintf("Strict window filter [%s..%s] retained %d work(s) (removed %d out-of-window Epub/preprint leak(s)).",
                  since_date, until_date, length(merged_recs), n_dropped_dt))
} else {
  merged_recs <- merged_recs_all
}

cat_pmids <- unique(vapply(studies$pmid, normalize_pmid_one, character(1)))
cat_pmids <- cat_pmids[nzchar(cat_pmids)]
cat_dois  <- unique(c(
  vapply(studies$doi, normalize_doi_one, character(1)),
  vapply(studies$doi, doi_dedup_key, character(1))
))
cat_dois  <- cat_dois[nzchar(cat_dois)]

preprint_rows_df <- studies[studies$status %in% c("Preprint", if (is_verify) "Superseded Preprint" else character(0)), , drop = FALSE]

known_ids <- tolower(unique(trimws(c(
  vapply(studies$doi, normalize_doi_one, character(1)),
  vapply(studies$doi, doi_dedup_key, character(1)),
  vapply(studies$pmid, normalize_pmid_one, character(1)),
  ifelse(nzchar(studies$pmid), paste0("pmid:", vapply(studies$pmid, normalize_pmid_one, character(1))), ""),
  vapply(excluded$identifier, normalize_doi_one, character(1)),
  vapply(excluded$identifier, doi_dedup_key, character(1)),
  extract_pmids_from_cell(paste(excluded$identifier, collapse = " ")),
  tolower(excluded$identifier),
  vapply(queue$identifier, normalize_doi_one, character(1)),
  vapply(queue$identifier, doi_dedup_key, character(1)),
  extract_pmids_from_cell(paste(queue$identifier, collapse = " ")),
  tolower(queue$identifier)
))))
known_ids <- known_ids[nzchar(known_ids)]

# ---------------------------------------------------------------------------
# Stage B.5: Hardened Selective PMC Full-Text Modality Verifier (Safeguard 3)
#   Resolves umbrella-abstract xQTL / molQTL / TR-xQTL / qQTL / SV-QTL studies
#   that omit explicit protein terms from title+abstract while mapping primary
#   human pQTLs in full-text Methods/Results. Strips <ref-list>, <back>, and
#   Introduction/Background/Discussion/Perspectives/Limitations sections,
#   enforces same-paragraph primary co-occurrence, skips records in `known_ids`,
#   and caps checks at `max_ft_checks_per_run = 60L`.
# ---------------------------------------------------------------------------
prune_non_evidence_secs <- function(parent_node) {
  child_secs <- xml2::xml_find_all(parent_node, "./sec")
  for (sc_node in child_secs) {
    st_attr  <- tolower(trimws(xml2::xml_attr(sc_node, "sec-type") %||% ""))
    st_title <- tolower(trimws(ascii_fold(xml2::xml_text(xml2::xml_find_first(sc_node, "./title")) %||% "")))
    is_bad_type <- nzchar(st_attr) &&
      grepl("\\b(intro|introduction|background|discussion|future|perspectives?|limitations?)\\b", st_attr, perl = TRUE)
    is_bad_title <- nzchar(st_title) &&
      grepl("^(?:[0-9.ivx]+\\s+)?(introduction|background(?:\\s+and\\s+[a-z ]+)?|discussion(?:\\s+and\\s+[a-z ]+)?|future\\s+[a-z ]+|perspectives?(?:\\s+and\\s+[a-z ]+)?|(?:study\\s+)?limitations?(?:\\s+and\\s+[a-z ]+)?)\\s*$",
            st_title, perl = TRUE)
    if (is_bad_type || is_bad_title) {
      xml2::xml_remove(sc_node)
    } else {
      prune_non_evidence_secs(sc_node)
    }
  }
}

verify_pmc_fulltext_pqtl <- function(pmcid, xml_override = NULL) {
  xml_str <- ""
  if (!is.null(xml_override) && length(xml_override) > 0L && nzchar(xml_override[1L])) {
    xml_str <- as.character(xml_override[1L])
  } else {
    if (!nzchar(pmcid %||% "")) return(list(verified = FALSE, platform = ""))
    ncbi_key <- Sys.getenv("NCBI_API_KEY", unset = "")
    u_pmc <- sprintf(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=%s&retmode=xml%s",
      pmcid, if (nzchar(ncbi_key)) paste0("&api_key=", ncbi_key) else ""
    )
    r1 <- perform_polite_http(
      httr2::request(u_pmc) |> httr2::req_user_agent(UA_STRING) |> httr2::req_timeout(35),
      sleep_s   = 0.55,
      max_tries = 3L,
      label     = paste0("NCBI PMC efetch (", pmcid, ")")
    )
    if (!is.null(r1) && httr2::resp_status(r1) == 200L) {
      cand_xml <- httr2::resp_body_string(r1)
      if (grepl("<article", cand_xml, fixed = TRUE)) xml_str <- cand_xml
    }
    if (!nzchar(xml_str)) {
      u_epmc <- sprintf("https://www.ebi.ac.uk/europepmc/webservices/rest/%s/fullTextXML", pmcid)
      r2 <- perform_polite_http(
        httr2::request(u_epmc) |> httr2::req_user_agent(UA_STRING) |> httr2::req_timeout(35),
        sleep_s   = 0.55,
        max_tries = 2L,
        label     = paste0("EuropePMC fullTextXML (", pmcid, ")")
      )
      if (!is.null(r2) && httr2::resp_status(r2) == 200L) {
        cand_xml <- httr2::resp_body_string(r2)
        if (grepl("<article", cand_xml, fixed = TRUE)) xml_str <- cand_xml
      }
    }
  }
  if (!nzchar(xml_str)) return(list(verified = FALSE, platform = ""))

  xdoc <- tryCatch(xml2::read_xml(xml_str), error = function(e) NULL)
  if (is.null(xdoc)) return(list(verified = FALSE, platform = ""))

  # Safeguard 3: Strip bibliography, back-matter, AND Intro/Background/Discussion/Limitations sections
  xml2::xml_remove(xml2::xml_find_all(xdoc, ".//ref-list | .//back"))
  body_node <- xml2::xml_find_first(xdoc, ".//body")
  if (inherits(body_node, "xml_missing")) return(list(verified = FALSE, platform = ""))

  has_body_secs <- length(xml2::xml_find_all(body_node, "./sec")) > 0L
  if (has_body_secs) {
    # Strip only the first untitled Nature-style leading Introduction <p> directly under <body>,
    # preserving subsequent main-text Results <p> nodes in Brief Communications (e.g. Nat Genet)
    direct_ps <- xml2::xml_find_all(body_node, "./p")
    if (length(direct_ps) >= 3L) {
      xml2::xml_remove(direct_ps[1L])
    } else if (length(direct_ps) > 0L && length(direct_ps) <= 2L) {
      xml2::xml_remove(direct_ps)
    }
    prune_non_evidence_secs(body_node)
  }

  cand_nodes  <- xml2::xml_find_all(body_node, ".//p | .//caption | .//table-wrap-foot")
  floats_caps <- xml2::xml_find_all(xdoc, ".//floats-group//caption | .//floats-group//p")
  para_txts   <- tolower(ascii_fold(gsub("\\s+", " ", trimws(c(xml2::xml_text(cand_nodes), xml2::xml_text(floats_caps))))))
  para_txts   <- para_txts[nzchar(para_txts)]
  if (length(para_txts) == 0L) {
    fb_txt <- tolower(ascii_fold(gsub("\\s+", " ", trimws(xml2::xml_text(body_node)))))
    if (nzchar(fb_txt)) para_txts <- fb_txt
  }
  if (length(para_txts) == 0L) return(list(verified = FALSE, platform = ""))

  pqtl_pat <- "\\b(tr-?pqtls?|sv-?pqtls?|cis-?pqtls?|trans-?pqtls?|sd-?pqtls?|pqtls?|protein[- ]qtls?|protein quantitative trait loc[a-z]*)\\b"
  prot_assay_pat <- "\\b(proteom[a-z]*|protein (?:abundance|expression|levels|quantification|profiles?)|olink|somascan|somalogic|somamers?|aptamer[- ]based|nulisa[a-z]*|mass spectrometry|openswath[0-9]*|dia-ms|dia-pasef|swath(?:-ms)?|tandem mass tags?|lc-ms/ms|proteograph|orbitrap|myriad[- ]rbm|luminex)\\b"
  primary_evid_pat <- paste0(
    "\\b(we (?:aimed to |next |then |also |systematically |comprehensively )?(?:performed|applied|mapped|computed|tested|conducted|generated|identified|discovered) [^.;]{0,110}\\b(?:tr-?pqtls?|sv-?pqtls?|cis-?pqtls?|trans-?pqtls?|sd-?pqtls?|pqtls?|protein[- ]qtls?|protein quantitative trait loc[a-z]*)|",
    "we used plink [0-9 ]*to estimate (?:protein quantitative trait loc[a-z]*|pqtls?)|",
    "to test for (?:tr-?|sv-?|cis-?|trans-?)?pqtls?|",
    "was also used for cis[- ]?pqtl analysis|",
    "conducted genome-wide association studies \\(gwas\\) of these proteins|",
    "meta-analysis of genome-wide association analyses was performed to identify variants associated with circulating protein levels|",
    "gwas analyses were performed using the r gwaf package|",
    "(?:of the [1-9][0-9,]* pqtls? in cis|n\\s*=\\s*[1-9][0-9,]*\\s+(?:[a-z-]+\\s+){0,2}proteom[a-z]*|[1-9][0-9,]*\\s+(?:[a-z-]+\\s+){0,2}proteomics? samples)(?=[\\s\\S]{0,220}\\b(?:tr-?pqtls?|sv-?pqtls?|tensorqtl|matrixeqtl|qtltools|fastgwa|regenie|qrank|quail|susie))|",
    "(?:tr-?pqtls?|sv-?pqtls?|cis-?pqtls?|pqtls?|protein[- ]qtls?) (?:were|was) (?:performed|mapped|computed|tested|generated)|",
    "tensorqtl|matrixeqtl|qtltools|fastgwa|regenie|qrank|quail|fastlmm|susie)\\b"
  )
  ext_only_pat <- paste0(
    "\\b(protein qtls? from the [a-z0-9 -]{0,35}generated by|curated (?:tumor )?database integrating|",
    "(?:previous|prior|published|existing|publicly available|publically available) (?:eqtl and )?pqtl (?:stud(?:y|ies)|analys(?:is|es)|datasets?)|",
    "pre-computed associations for (?:cis/?trans-?|cis-?|trans-?)?pqtls?|",
    "obtained locus-specific[^.;]{0,45}summary statistics from our genome-wide association|",
    "pqtls? acquired from phenoscanner|",
    "previously identified (?:genetic loci \\([^)]*p/?m?qtls?|p/?m?qtls?)|had a p/?m?qtls? previously identified|",
    "selected from (?:a )?pqtls? analys(?:is|es) conducted in|",
    "identified by relating genetic variants to protein levels in [a-z0-9-]+ as previously described|",
    "pqtls? come from interval|extracted the significant snps? reported in interval|",
    "also influenced by pqtls? in the|",
    "protein quantitative trait locus analysis for (?:two|three|four|[1-9]) (?:previously reported )?snps|",
    "determine the genetic drivers of the [a-z0-9-]+ and [a-z0-9-]+ protein levels|",
    "using the data from sun and colleagues|",
    "builds upon previously published datasets,\\s*including pqtl|",
    "leveraged summary statistics for protein quantitative trait loci|",
    "used pqtls? from the somascan|",
    "using two large somascan cis-?pqtl studies|",
    "gene expression at [0-9]+ h was compared to protein abundance|",
    "analysed top snps from [1-9][0-9]? [a-z0-9 -]+ gwas loci|",
    "performed gwass? for the metabolites|",
    "pan[- ]genomic pqtl analysis using [a-z0-9-]+ protein levels|",
    "[a-z0-9-]+ protein qtl analysis was performed on log10 transformed [a-z0-9-]+ protein levels)\\b"
  )

  if (any(grepl("\\b(the gwas was conducted to determine the genetic drivers of the [a-z0-9-]+ and [a-z0-9-]+ protein levels|gene expression at [0-9]+ h was compared to protein abundance|analysed top snps from [1-9][0-9]? [a-z0-9 -]+ gwas loci|builds upon previously published datasets,\\s*including pqtl|leveraged summary statistics for protein quantitative trait loci|used pqtls? from the somascan|using two large somascan cis-?pqtl studies)\\b", para_txts, perl = TRUE))) {
    return(list(verified = FALSE, platform = ""))
  }

  has_cis_window_protein_regression <- any(
    grepl("\\b(?:snps?|single nucleotide polymorphisms?(?:\\s*\\(snps?\\))?|variants?)\\s+within\\s+a\\s+cis\\s+window\\s+of\\s+[0-9]+\\s*(?:kb|mb)[^.;]{0,120}were\\s+tested\\s+for\\s+(?:an\\s+)?association\\s+with\\s+protein\\s+levels\\b",
          para_txts, perl = TRUE) &
      !grepl(ext_only_pat, para_txts, perl = TRUE)
  ) && any(grepl(prot_assay_pat, para_txts, perl = TRUE))

  has_in_cohort_protein_gwas <- any(
    grepl(paste0(
      "\\b(meta-analysis of genome-wide association analyses was performed to identify variants associated with circulating protein levels|",
      "gwas analyses were performed using the r gwaf package|",
      "association between each protein and genetic variant was performed using a linear regression model[^.;]{0,120}\\bplink\\b|",
      "we used fusion [0-9 ]*to estimate protein weights in the discovery|",
      "significantly associated with at least 1 of [1-9][0-9,]* genetic variants \\(protein quantitative trait loci)\\b"
    ), para_txts, perl = TRUE) &
      !grepl(ext_only_pat, para_txts, perl = TRUE)
  ) && any(grepl(prot_assay_pat, para_txts, perl = TRUE))

  para_ok <- (grepl(pqtl_pat, para_txts, perl = TRUE) &
    grepl(prot_assay_pat, para_txts, perl = TRUE) &
    grepl(primary_evid_pat, para_txts, perl = TRUE) &
    !grepl(ext_only_pat, para_txts, perl = TRUE)) |
    (any(grepl("\\bwe conducted genome-wide association studies \\(gwas\\) of these proteins in [1-9][0-9,]+", para_txts, perl = TRUE)) &
       grepl(pqtl_pat, para_txts, perl = TRUE) &
       !grepl(ext_only_pat, para_txts, perl = TRUE)) |
    has_cis_window_protein_regression |
    has_in_cohort_protein_gwas

  if (any(para_ok)) {
    kept_body_txt <- paste(para_txts, collapse = " ")
    return(list(verified = TRUE, platform = detect_platform_v2(kept_body_txt)))
  }
  list(verified = FALSE, platform = "")
}

needs_ft_pqtl_check <- function(m) {
  if (!nzchar(m$pmcid %||% "")) return(FALSE)
  pts_lc <- tolower(trimws(m$pub_types %||% character(0)))
  if (any(pts_lc %in% c("review", "systematic review", "published erratum", "erratum", "editorial", "comment"))) {
    return(FALSE)
  }
  ttl_lc  <- tolower(ascii_fold(strip_html_tags(m$title %||% "")))
  comb_lc <- tolower(ascii_fold(strip_html_tags(paste(m$title %||% "", m$abstract %||% ""))))
  if (grepl(NONHUMAN_ORGANISM_REGEX, ttl_lc, perl = TRUE) &&
      !grepl("\\b(human|humans|patients?|participants?|donors?|biobank|cohort)\\b", ttl_lc, perl = TRUE)) {
    return(FALSE)
  }
  if (grepl(BROAD_PROTEIN_REGEX, comb_lc, perl = TRUE)) {
    sc_pre <- score_triage_v2(
      title                    = m$title,
      abstract_text            = m$abstract,
      pub_types                = m$pub_types,
      mesh_headings            = m$mesh_headings,
      gwas_catalog_study_count = m$gwas_study_count,
      doi                      = m$doi
    )
    if (!grepl("^Auto-Deprioritised", sc_pre$classification)) return(FALSE)
    has_human_cohort_kw <- grepl("\\b(cohorts?|population-based|participants?|individuals?|subjects?|patients?|donors?|adults?|men|women|offspring|icelanders|twins?|human|humans)\\b",
                                 comb_lc, perl = TRUE)
    # Case 2B-1: Source 2B full-text Methods/Results rescue on human cohort proteomic papers without external sumstats reuse (S1)
    is_methods_rescue_low_s <- grepl("EuropePMC_MethodsRescue", m$matched_sources %||% "", fixed = TRUE) &&
      !("S1_PublicSumStatsReuse(+7)" %in% sc_pre$s_rules) &&
      sc_pre$S_score <= 6L && has_human_cohort_kw &&
      !grepl("\\b(two-sample mendelian|2-sample mendelian)\\b", comb_lc, perl = TRUE)
    # Case 2B-2: Primary cohort proteomics (P1 + P4, no external sumstats download S1) paired with bidirectional MR
    is_cohort_bidir_mr <- all(c("P1_DirectAssay(+5)", "P4_AssayScaleAndN(+4)") %in% sc_pre$p_rules) &&
      !("S1_PublicSumStatsReuse(+7)" %in% sc_pre$s_rules) &&
      sc_pre$S_score <= 6L && has_human_cohort_kw &&
      grepl("\\bbidirectional\\s+(?:two-sample\\s+)?mendelian\\s+randomi[sz]ation\\b", comb_lc, perl = TRUE)
    # Case 2B-3: Named proteomics platform + genomic architecture (P5, S=0) in a human cohort (e.g. Epub boundary papers missed by EPMC FIRST_PDATE)
    is_named_plat_gen_cohort <- sc_pre$P_score >= 3L && sc_pre$S_score == 0L && has_human_cohort_kw &&
      ("P5_GenomicArchitecture(+3)" %in% sc_pre$p_rules) &&
      !identical(detect_platform_v2(comb_lc), "Unspecified Proteomics")
    return(is_methods_rescue_low_s || is_cohort_bidir_mr || is_named_plat_gen_cohort)
  }
  has_multiomic_qtl_anchor <- grepl("\\b(molqtls?|xqtls?|tr-?xqtls?|sv-?qtls?|qqtls?|(?:molecular|multi-?omic|quantile) (?:quantitative trait loc[a-z]*|qtls?))\\b",
                                    comb_lc, perl = TRUE)
  has_primary_qtl_mapping  <- grepl("\\b(we (?:comprehensively |systematically )?(?:mapped|identified|discovered|quantified|performed|generated|profiled) [^.;]{0,140}\\b(?:qtls?|xqtls?|molqtls?|tr-?xqtls?|sv-?qtls?|qqtls?|quantitative trait loc[a-z]*)|(?:identified|discovered|mapped) [>~]?[0-9,.]+(?:\\s*million)? (?:[a-z-]+ ){0,3}(?:xqtls?|molqtls?|tr-?xqtls?|sv-?qtls?|qqtls?|quantitative trait loc[a-z]*))\\b",
                                    comb_lc, perl = TRUE)
  has_multiomic_qtl_anchor && has_primary_qtl_mapping &&
    grepl("\\b(donors?|participants?|individuals?|subjects?|patients?|multi-?omics? samples)\\b", comb_lc, perl = TRUE)
}

max_ft_checks_per_run <- 250L
n_ft_candidates       <- 0L
n_ft_checked          <- 0L
n_ft_skipped_known    <- 0L
n_ft_verified         <- 0L
for (idx_m in seq_along(merged_recs)) {
  merged_recs[[idx_m]]$ft_verified_pqtl <- FALSE
  merged_recs[[idx_m]]$ft_platform      <- ""
  m_cur <- merged_recs[[idx_m]]
  if (!needs_ft_pqtl_check(m_cur)) next

  m_idp <- if (nzchar(m_cur$doi)) m_cur$doi else if (nzchar(m_cur$pmid)) paste0("PMID:", m_cur$pmid) else if (nzchar(m_cur$pmcid)) m_cur$pmcid else m_cur$epmc_id
  m_keys <- tolower(c(
    m_idp, m_cur$doi, doi_dedup_key(m_cur$doi), m_cur$pmid,
    if (nzchar(m_cur$pmid)) paste0("pmid:", m_cur$pmid) else ""
  ))
  m_keys <- m_keys[nzchar(m_keys)]
  if (!is_verify && any(m_keys %in% known_ids)) {
    n_ft_skipped_known <- n_ft_skipped_known + 1L
    next
  }

  n_ft_candidates <- n_ft_candidates + 1L
  if (n_ft_checked >= max_ft_checks_per_run) next
  n_ft_checked <- n_ft_checked + 1L
  ft_res <- verify_pmc_fulltext_pqtl(m_cur$pmcid)
  if (isTRUE(ft_res$verified)) {
    merged_recs[[idx_m]]$ft_verified_pqtl <- TRUE
    merged_recs[[idx_m]]$ft_platform      <- ft_res$platform
    n_ft_verified <- n_ft_verified + 1L
  }
}
if (n_ft_candidates > max_ft_checks_per_run) {
  warning(sprintf(
    "[Stage B.5: Multi-Omic PMC Full-Text Verifier] %d candidate(s) matched needs_ft_pqtl_check(), exceeding max_ft_checks_per_run=%d; %d excess candidate(s) skipped.",
    n_ft_candidates, max_ft_checks_per_run, n_ft_candidates - max_ft_checks_per_run
  ))
}
if (n_ft_checked > 0L || n_ft_skipped_known > 0L) {
  message(sprintf(
    "[Stage B.5: Multi-Omic PMC Full-Text Verifier] checked %d umbrella xQTL/molQTL record(s) (skipped %d already in known_ids); verified %d with primary pQTL mapping.",
    n_ft_checked, n_ft_skipped_known, n_ft_verified
  ))
}

# ---------------------------------------------------------------------------
# Stage C: Triage Scoring, Preprint Lineage Attachment & Queue Formatting
# ---------------------------------------------------------------------------
run_stamp <- format(Sys.time(), "%Y%m%d_%H%M", tz = "UTC")
run_id    <- paste0("run_", run_stamp)

dump_rows  <- list()
queue_rows <- list()

for (m in merged_recs) {
  doi     <- m$doi
  pmid    <- m$pmid
  pmcid   <- m$pmcid
  epmc_id <- m$epmc_id

  id_primary <- if (nzchar(doi)) doi else if (nzchar(pmid)) paste0("PMID:", pmid) else if (nzchar(pmcid)) pmcid else epmc_id
  if (!nzchar(id_primary)) next

  err_targets_norm <- unique(c(m$erratum_targets, vapply(m$erratum_targets, doi_dedup_key, character(1))))
  err_on_cat <- any(err_targets_norm %in% c(cat_pmids, cat_dois))

  sc <- score_triage_v2(
    title                      = m$title,
    abstract_text              = m$abstract,
    pub_types                  = m$pub_types,
    mesh_headings              = m$mesh_headings,
    gwas_catalog_study_count   = m$gwas_study_count,
    erratum_targets_catalogued = err_on_cat,
    ft_verified_pqtl           = isTRUE(m$ft_verified_pqtl),
    doi                        = doi
  )

  # Selective PPR* -> DOI and PMID* -> DOI resolution on actionable / primary-signal records
  lineage_parts <- m$preprint_links
  has_pre_doi_already <- any(is_preprint_doi(lineage_parts))
  if (!has_pre_doi_already && !grepl("^Auto-Deprioritised", sc$classification)) {
    ppr_tokens <- grep("^PPR[0-9]+$", lineage_parts, value = TRUE)
    for (ptok in ppr_tokens) {
      resolved_d <- resolve_ppr_to_doi(ptok)
      if (nzchar(resolved_d)) lineage_parts <- unique(c(lineage_parts, resolved_d))
    }
  }
  lin_dois_now <- extract_dois_from_cell(paste(lineage_parts, collapse = "|"))
  has_pub_doi_already <- any(!is_preprint_doi(lin_dois_now))
  if (is_preprint_doi(doi) && !has_pub_doi_already && (!grepl("^Auto-Deprioritised", sc$classification) || sc$P_score >= 5L)) {
    pmid_tokens <- grep("^PMID:[0-9]+$", lineage_parts, value = TRUE, ignore.case = TRUE)
    for (pmtok in pmid_tokens) {
      resolved_pub_d <- resolve_pmid_to_doi(pmtok)
      if (nzchar(resolved_pub_d)) lineage_parts <- unique(c(lineage_parts, resolved_pub_d))
    }
  }

  lineage_dedup_keys <- unique(vapply(lineage_parts, doi_dedup_key, character(1)))
  lineage_dedup_keys <- lineage_dedup_keys[nzchar(lineage_dedup_keys)]

  matched_preprint_sid <- ""
  if (nrow(preprint_rows_df) > 0L && (nzchar(pmid) || !any(is_preprint_doi(doi)))) {
    for (pi in seq_len(nrow(preprint_rows_df))) {
      p_sid <- clean_cell(preprint_rows_df$study_id[pi])
      p_doi <- normalize_doi_one(preprint_rows_df$doi[pi])
      p_dk  <- doi_dedup_key(p_doi)
      p_st  <- clean_cell(preprint_rows_df$status[pi])

      if ((nzchar(p_doi) && (p_doi %in% lineage_parts)) || (nzchar(p_dk) && (p_dk %in% lineage_dedup_keys))) {
        lineage_parts <- unique(c(lineage_parts, p_sid))
        if (identical(p_st, "Preprint")) matched_preprint_sid <- p_sid
        next
      }
      if (p_sid %in% names(direct_preprint_upgrades)) {
        up <- direct_preprint_upgrades[[p_sid]]
        if ((nzchar(up$published_doi) && identical(doi_dedup_key(up$published_doi), doi_dedup_key(doi))) ||
            (nzchar(up$published_pmid) && identical(up$published_pmid, pmid))) {
          lineage_parts <- unique(c(lineage_parts, p_doi, p_sid, up$via))
          if (identical(p_st, "Preprint")) matched_preprint_sid <- p_sid
          next
        }
      }

      if (identical(p_st, "Preprint")) {
        p_au  <- normalize_author_surname(preprint_rows_df$first_author[pi])
        m_au1 <- normalize_author_surname(m$first_author)
        m_au2 <- normalize_author_surname(m$last_author)
        if (nzchar(p_au) && (identical(p_au, m_au1) || identical(p_au, m_au2))) {
          jac <- title_jaccard(preprint_rows_df$title[pi], m$title)
          if ((jac >= 0.35 || (identical(p_au, m_au1) && jac >= 0.30 && sc$P_score >= 12L && sc$S_score <= 4L)) &&
              has_compatible_phenotype(preprint_rows_df$title[pi], m$title)) {
            matched_preprint_sid <- p_sid
            lineage_parts <- unique(c(lineage_parts, p_doi, paste0("Tier4_TitleJaccard=", round(jac, 2))))
          } else if (identical(p_au, m_au1) && sc$P_score >= 4L && sc$S_score <= 6L && nrow(strata_cat) > 0L) {
            p_yr <- suppressWarnings(as.integer(preprint_rows_df$year[pi]))
            m_yr <- suppressWarnings(as.integer(m$pub_year))
            if (!is.na(p_yr) && !is.na(m_yr) && m_yr >= p_yr && (m_yr - p_yr) <= 3L) {
              st_rows <- strata_cat[strata_cat$study_id == p_sid, , drop = FALSE]
              p_assays <- suppressWarnings(as.integer(gsub("[^0-9]", "", c(st_rows$n_assays, st_rows$n_unique_proteins))))
              p_assays <- unique(p_assays[!is.na(p_assays) & (p_assays >= 100L | (p_assays >= 40L & sc$P_score >= 4L))])
              abs_digits <- gsub(",", "", normalize_si_thousands(m$abstract %||% ""))
              matched_panel <- Filter(function(na_v) {
                grepl(sprintf("\\b%d\\s+(?:[a-z-]+\\s+){0,2}(?:proteins?|protein\\s+measurements?|aptamers?|somamers?|circulating\\s+cytokines?|cytokines?)\\b", na_v),
                      tolower(abs_digits), perl = TRUE)
              }, p_assays)
              p_nsamp <- suppressWarnings(as.integer(gsub("[^0-9]", "", st_rows$n_samples)))
              p_nsamp <- unique(p_nsamp[!is.na(p_nsamp) & p_nsamp >= 500L])
              if (length(p_nsamp) >= 2L) p_nsamp <- unique(c(p_nsamp, sum(p_nsamp)))
              matched_nsamp <- Filter(function(ns_v) {
                grepl(sprintf("\\b%d\\b", ns_v), tolower(abs_digits), perl = TRUE)
              }, p_nsamp)
              pre_txt_lc <- tolower(ascii_fold(paste(preprint_rows_df$title[pi], paste(st_rows$harmonisation_notes, collapse = " "))))
              pub_txt_lc <- tolower(ascii_fold(paste(m$title, m$abstract)))
              same_sex_proteome <- grepl("\\bsex(?:es|\\s+differences)\\b", pre_txt_lc) &&
                grepl("\\bsex(?:es|\\s+differences)\\b", pub_txt_lc) &&
                grepl("\\bproteom[a-z]*\\b", pub_txt_lc)
              if (length(matched_panel) > 0L) {
                matched_preprint_sid <- p_sid
                lineage_parts <- unique(c(lineage_parts, p_doi, p_sid, sprintf("Tier4_SameAuthorPanel=%d", matched_panel[1L])))
              } else if (length(matched_nsamp) > 0L && sc$P_score >= 7L && sc$S_score <= 4L) {
                matched_preprint_sid <- p_sid
                lineage_parts <- unique(c(lineage_parts, p_doi, p_sid, sprintf("Tier4_SameAuthorSampleN=%d", matched_nsamp[1L])))
              } else if (same_sex_proteome && sc$P_score >= 12L && sc$S_score <= 4L) {
                matched_preprint_sid <- p_sid
                lineage_parts <- unique(c(lineage_parts, p_doi, p_sid, "Tier4_SameAuthorSexProteome"))
              }
            }
          }
        }
      }
    }
  }

  auto_cls <- if (nzchar(matched_preprint_sid) && sc$P_score >= 4L) {
    sprintf("Preprint -> Published Match (%s)", matched_preprint_sid)
  } else {
    sc$classification
  }

  plat_det <- detect_platform_v2(paste(m$title, m$abstract))
  if (identical(plat_det, "Unspecified Proteomics") && nzchar(m$ft_platform %||% "") &&
      !identical(m$ft_platform, "Unspecified Proteomics")) {
    plat_det <- m$ft_platform
  }
  if (identical(plat_det, "Unspecified Proteomics") && m$gwas_study_count >= 5L) {
    plat_det <- sprintf("GWAS Catalog (%d protein traits)", m$gwas_study_count)
  }
  n_det        <- extract_sample_size_v2(paste(m$title, m$abstract))
  lineage_note <- paste(unique(lineage_parts[nzchar(lineage_parts)]), collapse = "|")
  rules_str    <- paste(c(sc$p_rules, sc$s_rules), collapse = ",")

  cur_status <- if (grepl("^Auto-Deprioritised", auto_cls)) {
    "Auto-Triage: Low Priority"
  } else {
    "Pending Curator Review"
  }

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

  keys_check <- tolower(c(
    id_primary, doi, doi_dedup_key(doi), pmid,
    if (nzchar(pmid)) paste0("pmid:", pmid) else ""
  ))
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

# Sort candidate queue so Likely Primary & Needs Human Check appear first by P_score desc
if (nrow(new_queue_df) > 0L) {
  cls_rank <- ifelse(new_queue_df$auto_classification == "Likely Primary pGWAS", 1L,
              ifelse(grepl("^Preprint -> Published", new_queue_df$auto_classification), 2L,
              ifelse(grepl("^Needs Human Check", new_queue_df$auto_classification), 3L, 4L)))
  ord <- order(cls_rank, -(new_queue_df$p_score - new_queue_df$s_score), -new_queue_df$p_score)
  new_queue_df <- new_queue_df[ord, , drop = FALSE]
  rownames(new_queue_df) <- NULL
}

dedup_n   <- nrow(new_queue_df)
auto_mr_n <- if (dedup_n > 0L) sum(grepl("^Auto-Deprioritised", new_queue_df$auto_classification)) else 0L
review_n  <- dedup_n - auto_mr_n
upg_n     <- if (dedup_n > 0L) sum(grepl("^Preprint -> Published Match", new_queue_df$auto_classification)) else 0L

target_dump_dir <- if (nzchar(out_dir)) out_dir else "data"
dump_path <- file.path(target_dump_dir, sprintf("search_dump_%s.tsv", run_id))
if (nrow(dump_df) > 0L) {
  write_tsv_safe(dump_df, dump_path)
  message(sprintf("Wrote uncapped search ledger (%d row(s)) to %s", nrow(dump_df), dump_path))
}

sources_desc <- sprintf(
  "PubMed(%s);EuropePMC(2A_Platforms+2A_Proteome+2B_Methods);GWASCatalog_v2;Snowball_%dSeeds [%s..%s]",
  pm_date_field, nrow(LANDMARK_PGWAS_SEEDS), since_date, until_date
)

upg_details <- if (dedup_n > 0L) {
  upg_idx <- which(grepl("^Preprint -> Published Match", new_queue_df$auto_classification))
  if (length(upg_idx) > 0L) {
    paste(vapply(upg_idx, function(ii) {
      pre_sid <- sub("^Preprint -> Published Match \\(([^)]+)\\)$", "\\1", new_queue_df$auto_classification[ii])
      pub_ref <- if (nzchar(new_queue_df$pmid[ii])) paste0("PMID:", new_queue_df$pmid[ii]) else new_queue_df$identifier[ii]
      sprintf("%s -> %s", pre_sid, pub_ref)
    }, character(1)), collapse = "; ")
  } else ""
} else ""

new_run_row <- data.frame(
  run_id                     = run_id,
  run_date                   = until_date,
  date_from                  = since_date,
  date_to                    = until_date,
  trigger                    = clean_cell(trigger_lb),
  sources_queried            = sources_desc,
  raw_records                = as.integer(raw_n),
  deduplicated_new           = as.integer(dedup_n),
  auto_flagged_secondary_mr  = as.integer(auto_mr_n),
  candidates_for_review      = as.integer(review_n),
  studies_added              = 0L,
  preprints_upgraded         = as.integer(upg_n),
  pr_number                  = if (is_verify) "Verify-Only" else "Automated PR",
  release_tag                = "",
  status                     = if (review_n > 0L) sprintf("Pending Triage (%d candidates)", review_n) else "Completed (0 new primary)",
  studies_added_ids          = "",
  preprints_replaced_details = upg_details,
  stringsAsFactors           = FALSE
)

if (nzchar(out_dir)) {
  write_tsv_safe(new_queue_df, file.path(out_dir, "candidate_queue_preview.tsv"))
  write_tsv_safe(new_run_row,  file.path(out_dir, "search_runs_preview.tsv"))
  message(sprintf("Wrote preview tables (%d candidate(s), %d actionable) to %s",
                  nrow(new_queue_df), review_n, out_dir))
}

if (isTRUE(write_catalogue)) {
  if (dedup_n > 0L) {
    actionable <- new_queue_df[!grepl("^Auto-Deprioritised", new_queue_df$auto_classification), , drop = FALSE]
    deprior    <- utils::head(new_queue_df[grepl("^Auto-Deprioritised", new_queue_df$auto_classification), , drop = FALSE], 15L)

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

  all_run_cols <- union(names(new_run_row), names(runs))
  for (rc in all_run_cols) {
    if (!(rc %in% names(runs))) runs[[rc]] <- ""
    if (!(rc %in% names(new_run_row))) new_run_row[[rc]] <- ""
  }
  runs_updated <- rbind(runs[, all_run_cols, drop = FALSE], new_run_row[, all_run_cols, drop = FALSE])
  write_tsv_safe(runs_updated, file.path("catalogue", "search_runs.tsv"))

  val_status <- system2("Rscript", "src/validate_catalogue.R")
  if (val_status != 0L) stop("Catalogue validation failed after search update!", call. = FALSE)
}
