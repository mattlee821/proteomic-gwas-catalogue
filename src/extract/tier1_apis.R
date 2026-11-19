# Tier 1: Structured Metadata, Preprint Lineage, OA Licensing & GWAS Catalog v2 APIs

normalize_doi <- function(x) {
  d <- tolower(clean_str(x))
  if (!nzchar(d)) return("")
  d <- utils::URLdecode(d)
  d <- sub("^https?://(dx\\.)?doi\\.org/", "", d)
  d <- sub("^doi:\\s*", "", d)
  d <- sub("(\\.full(\\.pdf)?|\\.abstract|\\.short|/v[0-9]+|v[0-9]+)$", "", d, perl = TRUE)
  d <- sub("[.,;)\\]]+$", "", d)
  if (grepl("^10\\.[0-9]{4,9}/", d)) d else ""
}

normalize_pmid <- function(x) {
  p <- clean_str(x)
  p <- sub("^pmid:\\s*", "", p, ignore.case = TRUE)
  m <- regmatches(p, regexpr("\\b[1-9][0-9]{5,8}\\b", p))
  if (length(m) == 1L && nzchar(m)) m else ""
}

normalize_pmcid <- function(x) {
  p <- toupper(clean_str(x))
  m <- regmatches(p, regexpr("\\bPMC[0-9]+\\b", p))
  if (length(m) == 1L && nzchar(m)) m else ""
}

normalize_epmc_id <- function(x) {
  p <- toupper(clean_str(x))
  m <- regmatches(p, regexpr("\\bPPR[0-9]+\\b", p))
  if (length(m) == 1L && nzchar(m)) m else ""
}

ascii_slug <- function(x) {
  s <- clean_str(x)
  s_tr <- tryCatch(iconv(s, from = "", to = "ASCII//TRANSLIT"), error = function(e) s)
  if (is.na(s_tr)) s_tr <- s
  s_tr <- tolower(gsub("[^A-Za-z0-9]+", "_", s_tr))
  s_tr <- gsub("^_+|_+$", "", s_tr)
  if (nzchar(s_tr)) s_tr else "study"
}

legacy_prefix_for_study <- function(pmid, doi) {
  if (identical(pmid, "37794186") || identical(doi, "10.1038/s41586-023-06592-6")) return("sun_2023")
  if (identical(pmid, "34857953") || identical(doi, "10.1038/s41588-021-00978-w")) return("ferkingstad_2021")
  ""
}

run_tier1_apis <- function(cand, session, existing_studies = NULL) {
  pmid <- normalize_pmid(cand$pmid %||% "")
  doi <- normalize_doi(cand$doi %||% "")
  pmcid <- normalize_pmcid(cand$pmcid %||% "")
  epmc_id <- normalize_epmc_id(cand$epmc_id %||% "")
  lineage_hint <- clean_str(cand$preprint_lineage_note %||% "")

  # Fallback parsing if caller only supplied a generic `identifier`
  raw_ident <- clean_str(cand$identifier %||% "")
  if (nzchar(raw_ident)) {
    if (!nzchar(doi)) doi <- normalize_doi(raw_ident)
    if (!nzchar(pmid) && !grepl("^10\\.", raw_ident)) pmid <- normalize_pmid(raw_ident)
    if (!nzchar(pmcid)) pmcid <- normalize_pmcid(raw_ident)
    if (!nzchar(epmc_id)) epmc_id <- normalize_epmc_id(raw_ident)
  }

  leg_pref <- legacy_prefix_for_study(pmid, doi)
  slug_base <- if (nzchar(pmid)) paste0("pmid_", pmid) else if (nzchar(epmc_id)) tolower(epmc_id) else ascii_slug(doi)

  out <- list(
    pmid = pmid,
    doi = doi,
    pmcid = pmcid,
    epmc_id = epmc_id,
    epmc_source = "",
    first_author = "",
    year = NA_integer_,
    journal = "",
    title = clean_str(cand$title %||% ""),
    abstract_text = "",
    is_preprint = grepl("^10\\.1101/|researchsquare", doi) || nzchar(epmc_id),
    is_open_access_epmc = "N",
    in_epmc = "N",
    in_pmc = "N",
    has_suppl = "N",
    is_oa = FALSE,
    oa_status = "closed",
    is_cc_license = FALSE,
    license_url = "",
    best_oa_url = "",
    full_text_urls = list(),
    supersedes_preprint_doi = "",
    supersedes_epmc_id = "",
    published_upgrade_doi = "",
    published_upgrade_pmid = "",
    biostudies_acc = "",
    datalinks_urls = character(0),
    databank_accessions = character(0),
    biorxiv_jats_url = "",
    biorxiv_server = "",
    gwas_catalog_count = 0L,
    gwas_studies = list(),
    lineage_hint = lineage_hint,
    study_id = "",
    status_label = ""
  )

  # 1a. Europe PMC core lookup (NEVER append OPEN_ACCESS:y here!)
  epmc_q <- if (nzchar(pmid)) {
    paste0("EXT_ID:", pmid, " AND SRC:MED")
  } else if (nzchar(epmc_id)) {
    paste0("EXT_ID:", epmc_id, " AND SRC:PPR")
  } else if (nzchar(doi)) {
    paste0("DOI:\"", doi, "\"")
  } else if (nzchar(pmcid)) {
    paste0("PMCID:", pmcid)
  } else {
    ""
  }

  if (nzchar(epmc_q)) {
    epmc_url <- paste0(
      "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=",
      utils::URLencode(epmc_q, reserved = TRUE),
      "&resultType=core&format=json"
    )
    leg_epmc <- if (nzchar(leg_pref)) paste0(leg_pref, "_epmc_core.json") else character(0)
    r_epmc <- perform_cached_request(
      session, "Tier 1", "EuropePMC_core", epmc_url,
      cache_key = paste0(slug_base, "_epmc_core.json"),
      accept = "application/json",
      legacy_candidates = leg_epmc
    )
    if (r_epmc$ok) {
      j <- tryCatch(jsonlite::fromJSON(r_epmc$path, simplifyVector = FALSE), error = function(e) NULL)
      recs <- j$resultList$result %||% list()
      if (length(recs) > 0L) {
        # Prefer MED over PPR when querying by DOI if both exist, unless explicitly looking up PPR
        rec <- recs[[1L]]
        if (!nzchar(epmc_id) && length(recs) > 1L) {
          for (rc in recs) {
            if (identical(rc$source, "MED")) { rec <- rc; break }
          }
        }
        out$epmc_source <- clean_str(rec$source)
        if (!nzchar(out$pmid)) out$pmid <- normalize_pmid(rec$pmid %||% "")
        if (!nzchar(out$pmcid)) out$pmcid <- normalize_pmcid(rec$pmcid %||% "")
        if (!nzchar(out$doi)) out$doi <- normalize_doi(rec$doi %||% "")
        if (identical(out$epmc_source, "PPR") && !nzchar(out$epmc_id)) {
          out$epmc_id <- normalize_epmc_id(rec$id %||% "")
        }
        if (identical(out$epmc_source, "PPR")) out$is_preprint <- TRUE
        if (identical(out$epmc_source, "MED")) out$is_preprint <- FALSE

        auths <- rec$authorList$author %||% list()
        if (length(auths) > 0L) {
          out$first_author <- clean_str(auths[[1L]]$lastName %||% auths[[1L]]$collectiveName %||% "")
        }
        if (!nzchar(out$first_author) && nzchar(clean_str(rec$authorString))) {
          out$first_author <- sub("\\s+.*$", "", clean_str(rec$authorString))
        }
        out$year <- suppressWarnings(as.integer(rec$pubYear %||% NA_integer_))
        out$journal <- clean_str(rec$journalInfo$journal$title %||% "")
        if (!nzchar(out$journal) && out$is_preprint) {
          out$journal <- clean_str(rec$bookOrReportDetails$publisher %||% "bioRxiv")
        }
        out$title <- sub("\\.$", "", clean_str(rec$title %||% out$title))
        out$abstract_text <- clean_str(gsub("<[^>]+>", " ", rec$abstractText %||% ""))
        out$is_open_access_epmc <- clean_str(rec$isOpenAccess %||% "N")
        out$in_epmc <- clean_str(rec$inEPMC %||% "N")
        out$in_pmc <- clean_str(rec$inPMC %||% "N")
        out$has_suppl <- clean_str(rec$hasSuppl %||% "N")
        out$full_text_urls <- rec$fullTextUrlList$fullTextUrl %||% list()

        ccs <- rec$commentCorrectionList$commentCorrection %||% list()
        for (cc in ccs) {
          cc_type <- tolower(clean_str(cc$type))
          cc_src <- clean_str(cc$source)
          cc_id <- clean_str(cc$id)
          if (grepl("preprint in", cc_type) && identical(cc_src, "PPR")) {
            out$supersedes_epmc_id <- cc_id
          }
          if (grepl("published in|peer-reviewed|update in", cc_type) && identical(cc_src, "MED")) {
            out$published_upgrade_pmid <- normalize_pmid(cc_id)
          }
        }
      }
    }
  }

  # 1b. Fallback PubMed esearch by DOI if PMID still missing on non-preprint DOI ([C-1])
  if (!nzchar(out$pmid) && nzchar(out$doi) && !out$is_preprint) {
    es_url <- paste0(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=",
      utils::URLencode(paste0(out$doi, "[doi]"), reserved = TRUE),
      "&retmode=json"
    )
    r_es <- perform_cached_request(
      session, "Tier 1", "PubMed_esearch_DOI", es_url,
      cache_key = paste0(slug_base, "_pubmed_esearch.json"),
      accept = "application/json"
    )
    if (r_es$ok) {
      jes <- tryCatch(jsonlite::fromJSON(r_es$path, simplifyVector = TRUE), error = function(e) NULL)
      id_vec <- jes$esearchresult$idlist %||% character(0)
      if (length(id_vec) > 0L) out$pmid <- normalize_pmid(id_vec[[1L]])
    }
  }

  # Refresh legacy prefix / slug if PMID resolved via DOI
  if (!nzchar(leg_pref)) leg_pref <- legacy_prefix_for_study(out$pmid, out$doi)

  # 1c. Europe PMC datalinks (works for MED/{pmid} AND PPR/{epmc_id} — [C-6])
  dl_src <- if (nzchar(out$pmid) && !out$is_preprint) "MED" else if (nzchar(out$epmc_id)) "PPR" else ""
  dl_id <- if (identical(dl_src, "MED")) out$pmid else out$epmc_id
  if (nzchar(dl_src) && nzchar(dl_id)) {
    dl_url <- sprintf("https://www.ebi.ac.uk/europepmc/webservices/rest/%s/%s/datalinks?format=json", dl_src, dl_id)
    leg_dl <- if (nzchar(leg_pref) && identical(dl_src, "MED")) paste0(leg_pref, "_epmc_datalinks.json") else character(0)
    r_dl <- perform_cached_request(
      session, "Tier 1", "EuropePMC_datalinks", dl_url,
      cache_key = paste0(slug_base, "_epmc_datalinks.json"),
      accept = "application/json",
      legacy_candidates = leg_dl
    )
    if (r_dl$ok) {
      jdl <- tryCatch(jsonlite::fromJSON(r_dl$path, simplifyVector = FALSE), error = function(e) NULL)
      cats <- jdl$dataLinkList$Category %||% list()
      dl_targets <- character(0)
      for (cat_it in cats) {
        secs <- cat_it$Section %||% list()
        for (sc in secs) {
          lnks <- sc$Linklist$Link %||% list()
          for (lk in lnks) {
            tid <- clean_str(lk$Target$Identifier$ID %||% "")
            turl <- clean_str(lk$Target$Identifier$IDURL %||% "")
            if (grepl("S-EPMC[0-9]+", tid)) {
              out$biostudies_acc <- regmatches(tid, regexpr("S-EPMC[0-9]+", tid))
            } else if (grepl("S-EPMC[0-9]+", turl)) {
              out$biostudies_acc <- regmatches(turl, regexpr("S-EPMC[0-9]+", turl))
            }
            if (nzchar(turl)) dl_targets <- c(dl_targets, turl)
            if (nzchar(tid)) dl_targets <- c(dl_targets, tid)
          }
        }
      }
      out$datalinks_urls <- unique(dl_targets)
    }
  }
  if (!nzchar(out$biostudies_acc) && nzchar(out$pmcid) && identical(out$has_suppl, "Y")) {
    out$biostudies_acc <- paste0("S-E", out$pmcid)
  }

  # 1d. PubMed efetch XML (when PMID available)
  if (nzchar(out$pmid)) {
    pm_url <- sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=%s&retmode=xml", out$pmid)
    leg_pm <- if (nzchar(leg_pref)) paste0(leg_pref, "_pubmed_efetch.xml") else character(0)
    r_pm <- perform_cached_request(
      session, "Tier 1", "PubMed_efetch_XML", pm_url,
      cache_key = paste0(slug_base, "_pubmed_efetch.xml"),
      accept = "application/xml",
      legacy_candidates = leg_pm
    )
    if (r_pm$ok) {
      doc <- tryCatch(xml2::read_xml(r_pm$path), error = function(e) NULL)
      if (!is.null(doc)) {
        pmc_nd <- xml2::xml_find_first(doc, "//PubmedData/ArticleIdList/ArticleId[@IdType='pmc']")
        if (!inherits(pmc_nd, "xml_missing") && !nzchar(out$pmcid)) {
          out$pmcid <- normalize_pmcid(xml2::xml_text(pmc_nd))
        }
        doi_nd <- xml2::xml_find_first(doc, "//PubmedData/ArticleIdList/ArticleId[@IdType='doi']")
        if (!inherits(doi_nd, "xml_missing") && !nzchar(out$doi)) {
          out$doi <- normalize_doi(xml2::xml_text(doi_nd))
        }
        if (!nzchar(out$abstract_text)) {
          out$abstract_text <- clean_str(paste(xml2::xml_text(xml2::xml_find_all(doc, "//Abstract/AbstractText")), collapse = " "))
        }
        dbs <- xml2::xml_find_all(doc, "//DataBankList/DataBank")
        accs <- character(0)
        for (db in dbs) {
          db_nm <- clean_str(xml2::xml_text(xml2::xml_find_first(db, "DataBankName")))
          ac_vec <- xml2::xml_text(xml2::xml_find_all(db, ".//AccessionNumber"))
          if (length(ac_vec) > 0L) accs <- c(accs, paste0(db_nm, ":", ac_vec))
        }
        out$databank_accessions <- unique(accs)
      }
    }
  }

  # 1e. Crossref + Unpaywall + OpenAlex + bioRxiv/medRxiv details API
  if (nzchar(out$doi)) {
    cr_url <- paste0("https://api.crossref.org/works/", utils::URLencode(out$doi, reserved = TRUE))
    leg_cr <- if (nzchar(leg_pref)) paste0(leg_pref, "_crossref.json") else character(0)
    r_cr <- perform_cached_request(
      session, "Tier 1", "Crossref_works", cr_url,
      cache_key = paste0(slug_base, "_crossref.json"),
      accept = "application/json",
      legacy_candidates = leg_cr
    )
    if (r_cr$ok) {
      jcr <- tryCatch(jsonlite::fromJSON(r_cr$path, simplifyVector = FALSE)$message, error = function(e) NULL)
      if (!is.null(jcr)) {
        lics <- jcr$license %||% list()
        for (lc in lics) {
          u_lc <- clean_str(lc$URL %||% "")
          if (nzchar(u_lc) && !nzchar(out$license_url)) out$license_url <- u_lc
          if (grepl("creativecommons\\.org|cc0", u_lc, ignore.case = TRUE)) out$is_cc_license <- TRUE
        }
        hp <- jcr$relation$`has-preprint` %||% list()
        if (length(hp) > 0L) {
          out$supersedes_preprint_doi <- normalize_doi(hp[[1L]]$id %||% "")
        }
        ipo <- jcr$relation$`is-preprint-of` %||% list()
        if (length(ipo) > 0L) {
          out$published_upgrade_doi <- normalize_doi(ipo[[1L]]$id %||% "")
        }
        if (!nzchar(out$title) && length(jcr$title) > 0L) {
          out$title <- sub("\\.$", "", clean_str(jcr$title[[1L]]))
        }
        if (!nzchar(out$first_author) && length(jcr$author) > 0L) {
          out$first_author <- clean_str(jcr$author[[1L]]$family %||% jcr$author[[1L]]$name %||% "")
        }
      }
    }

    # Unpaywall v2
    upw_url <- paste0("https://api.unpaywall.org/v2/", out$doi, "?email=", utils::URLencode(session$contact_email, reserved = TRUE))
    leg_upw <- if (nzchar(leg_pref)) paste0(leg_pref, "_unpaywall.json") else character(0)
    r_upw <- perform_cached_request(
      session, "Tier 1", "Unpaywall_v2", upw_url,
      cache_key = paste0(slug_base, "_unpaywall.json"),
      accept = "application/json",
      legacy_candidates = leg_upw
    )
    if (r_upw$ok) {
      jup <- tryCatch(jsonlite::fromJSON(r_upw$path, simplifyVector = FALSE), error = function(e) NULL)
      if (!is.null(jup)) {
        out$is_oa <- isTRUE(jup$is_oa)
        out$oa_status <- clean_str(jup$oa_status %||% "closed")
        out$best_oa_url <- clean_str(jup$best_oa_location$url %||% "")
        lic_u <- clean_str(jup$best_oa_location$license %||% "")
        if (grepl("^cc-|cc0|public-domain", lic_u, ignore.case = TRUE)) out$is_cc_license <- TRUE
      }
    }

    # OpenAlex fallback / corroboration
    oa_url <- paste0("https://api.openalex.org/works/https://doi.org/", out$doi)
    leg_oa <- if (nzchar(leg_pref)) paste0(leg_pref, "_openalex.json") else character(0)
    r_oa <- perform_cached_request(
      session, "Tier 1", "OpenAlex_work", oa_url,
      cache_key = paste0(slug_base, "_openalex.json"),
      accept = "application/json",
      legacy_candidates = leg_oa
    )
    if (r_oa$ok) {
      joa <- tryCatch(jsonlite::fromJSON(r_oa$path, simplifyVector = FALSE), error = function(e) NULL)
      if (!is.null(joa) && !out$is_oa && isTRUE(joa$open_access$is_oa)) {
        out$is_oa <- TRUE
        out$oa_status <- clean_str(joa$open_access$oa_status %||% out$oa_status)
        if (!nzchar(out$best_oa_url)) out$best_oa_url <- clean_str(joa$open_access$oa_url %||% "")
      }
    }

    # bioRxiv / medRxiv details API when DOI is 10.1101/* ([C-6])
    if (grepl("^10\\.1101/", out$doi)) {
      out$is_preprint <- TRUE
      for (srv in c("biorxiv", "medrxiv")) {
        bx_url <- sprintf("https://api.biorxiv.org/details/%s/%s/na/json", srv, out$doi)
        r_bx <- perform_cached_request(
          session, "Tier 1", paste0("bioRxiv_API_", srv), bx_url,
          cache_key = paste0(slug_base, "_biorxiv_", srv, ".json"),
          accept = "application/json"
        )
        if (r_bx$ok) {
          jbx <- tryCatch(jsonlite::fromJSON(r_bx$path, simplifyVector = FALSE), error = function(e) NULL)
          coll <- jbx$collection %||% list()
          if (length(coll) > 0L) {
            rec_bx <- coll[[length(coll)]] # latest version
            out$biorxiv_server <- srv
            out$biorxiv_jats_url <- clean_str(rec_bx$jatsxml %||% "")
            if (!nzchar(out$journal) || identical(tolower(out$journal), "biorxiv")) {
              out$journal <- if (identical(srv, "medrxiv")) "medRxiv" else "bioRxiv"
            }
            if (!nzchar(out$abstract_text)) out$abstract_text <- clean_str(rec_bx$abstract %||% "")
            if (!nzchar(out$title)) out$title <- sub("\\.$", "", clean_str(rec_bx$title %||% ""))
            if (is.na(out$year) && nzchar(clean_str(rec_bx$date))) {
              out$year <- suppressWarnings(as.integer(substr(clean_str(rec_bx$date), 1, 4)))
            }
            pub_d <- normalize_doi(rec_bx$published %||% "")
            if (nzchar(pub_d) && !nzchar(out$published_upgrade_doi)) {
              out$published_upgrade_doi <- pub_d
            }
            lic_bx <- clean_str(rec_bx$license %||% "")
            if (grepl("cc_by|cc0|cc-by", lic_bx, ignore.case = TRUE)) out$is_cc_license <- TRUE
            break
          }
        }
      }
    }
  }

  # 1f. EBI GWAS Catalog REST API v2 (strictly ?pubmed_id=, never legacy v1 or ?pmid=)
  if (nzchar(out$pmid)) {
    gw_url <- sprintf("https://www.ebi.ac.uk/gwas/rest/api/v2/studies?pubmed_id=%s&page=0&size=100", out$pmid)
    leg_gw <- if (nzchar(leg_pref)) paste0(leg_pref, "_gwas_v2_studies_page0.json") else character(0)
    r_gw <- perform_cached_request(
      session, "Tier 1", "GWASCatalog_v2_studies", gw_url,
      cache_key = paste0(slug_base, "_gwas_v2_studies.json"),
      accept = "application/json",
      legacy_candidates = leg_gw
    )
    if (r_gw$ok) {
      jgw <- tryCatch(jsonlite::fromJSON(r_gw$path, simplifyVector = FALSE), error = function(e) NULL)
      if (!is.null(jgw)) {
        out$gwas_catalog_count <- as.integer(jgw$page$totalElements %||% 0L)
        out$gwas_studies <- jgw$`_embedded`$studies %||% list()
      }
    }
  }

  # Resolve canonical status & study_id
  if (out$is_preprint) {
    out$status_label <- if (nzchar(out$published_upgrade_doi) || nzchar(out$published_upgrade_pmid)) {
      "Superseded Preprint"
    } else {
      "Preprint"
    }
  } else {
    out$status_label <- "Published"
  }

  # Match existing study_id in catalogue/studies.tsv if already tracked
  matched_existing_id <- ""
  if (!is.null(existing_studies) && nrow(existing_studies) > 0L) {
    if (nzchar(out$pmid)) {
      idx <- which(as.character(existing_studies$pmid) == out$pmid)
      if (length(idx) > 0L) matched_existing_id <- existing_studies$study_id[[idx[[1L]]]]
    }
    if (!nzchar(matched_existing_id) && nzchar(out$doi)) {
      idx <- which(tolower(trimws(existing_studies$doi)) == out$doi)
      if (length(idx) > 0L) matched_existing_id <- existing_studies$study_id[[idx[[1L]]]]
    }
  }

  if (nzchar(matched_existing_id)) {
    out$study_id <- matched_existing_id
  } else {
    yr_str <- if (!is.na(out$year)) as.character(out$year) else "0000"
    au_slug <- ascii_slug(out$first_author)
    out$study_id <- if (nzchar(out$pmid)) {
      paste0(au_slug, "_", yr_str, "_", out$pmid)
    } else if (nzchar(out$epmc_id)) {
      paste0(au_slug, "_", yr_str, "_", tolower(out$epmc_id))
    } else {
      paste0(au_slug, "_", yr_str, "_pre")
    }
  }

  out
}
