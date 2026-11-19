# Tier 3d/3e: 1-Hop Portal Resolution, Public Repository Metadata Preflights,
# and Objective HTTP-Testable `access_level` Classification ([M-6], [m-4]).
#
# Objective `access_level` policy:
#   - "Open download": Direct unauthenticated HTTP/HTTPS/FTP binary/text file download
#     with no account login and no interactive JS/HTML click-through form.
#   - "Portal registration": Self-service researcher access gated by either a free portal
#     account/token (e.g. Synapse anonymous HTTP 401/403 Read-Only gate, NIAGADS, CNCB OMIX)
#     OR an interactive browser name/email + license click-through form (e.g. deCODE /form/ SPA).
#   - "Controlled access": Formal DAC proposal / DUA approval required (dbGaP phs*, EGA EGAD*).
#   - "Upon request": Email / collaboration request only.
#   - "Not available": No external summary statistics deposited.

run_tier3_repositories <- function(t1, t2, session) {
  leg_pref <- legacy_prefix_for_study(t1$pmid, t1$doi)
  slug_base <- t1$study_id

  search_blob <- paste(
    t1$title, t1$abstract_text, t2$data_availability_text,
    paste(t1$datalinks_urls, collapse = " "),
    paste(t1$databank_accessions, collapse = " "),
    t1$lineage_hint
  )

  res <- list(
    repository = "",
    accession_or_url = "",
    access_level = "",
    access_tier = "Unresolved",
    access_status = "MISSING",
    location = "",
    snippet = "",
    access_note = "",
    synapse_folders = character(0),
    detected_file_format = "",
    helper_files_note = ""
  )

  # 1. 1-Hop Portal Resolution for UKB-PPP portal (http://ukb-ppp.gwas.eu -> syn51364943)
  if (grepl("ukb-ppp\\.gwas\\.eu", search_blob, ignore.case = TRUE)) {
    leg_p <- if (nzchar(leg_pref)) paste0(leg_pref, "_ukbppp_portal.html") else character(0)
    r_ukb <- perform_cached_request(
      session, "Tier 3", "Portal_1hop_ukbppp", "http://ukb-ppp.gwas.eu",
      cache_key = paste0(slug_base, "_ukbppp_portal.html"),
      accept = "text/html",
      timeout_sec = 15,
      allow_transient_retry = FALSE,
      legacy_candidates = leg_p
    )
    if (r_ukb$ok && file.exists(r_ukb$path)) {
      ptxt <- paste(readLines(r_ukb$path, warn = FALSE), collapse = "\n")
      search_blob <- paste(search_blob, ptxt)
    }
  }

  # 2. Check Synapse accession (syn[0-9]{5,})
  m_syn <- regmatches(search_blob, regexpr("\\bsyn[0-9]{6,10}\\b", search_blob, ignore.case = TRUE))
  if (length(m_syn) == 1L && nzchar(m_syn)) {
    syn_id <- tolower(m_syn)
    ent_url <- paste0("https://repo-prod.prod.sagebase.org/repo/v1/entity/", syn_id)
    leg_ent <- if (nzchar(leg_pref)) paste0(leg_pref, "_synapse_entity.json") else character(0)
    r_ent <- perform_cached_request(
      session, "Tier 3", "Synapse_entity", ent_url,
      cache_key = paste0(slug_base, "_synapse_", syn_id, ".json"),
      accept = "application/json",
      legacy_candidates = leg_ent
    )

    res$repository <- paste0("Synapse (", syn_id, ")")
    res$accession_or_url <- paste0("https://www.synapse.org/Synapse:", syn_id)
    res$access_level <- "Portal registration"
    res$access_tier <- "Tier 3"
    res$access_status <- "SUPPLEMENT_EXTRACTED"
    res$location <- paste0("Synapse REST /repo/v1/entity/", syn_id)
    res$snippet <- paste0("Synapse project ", syn_id)

    # Enumerate top-level children and summary-statistics subfolder if public READ allowed
    ch_url <- "https://repo-prod.prod.sagebase.org/repo/v1/entity/children"
    leg_ch <- if (nzchar(leg_pref)) paste0(leg_pref, "_synapse_children.json") else character(0)
    r_ch <- perform_cached_request(
      session, "Tier 3", "Synapse_children_root", ch_url,
      cache_key = paste0(slug_base, "_synapse_children_", syn_id, ".json"),
      method = "POST",
      accept = "application/json",
      json_body = list(parentId = syn_id, includeTypes = list("folder", "file", "link")),
      legacy_candidates = leg_ch
    )
    if (r_ch$ok) {
      jch <- tryCatch(jsonlite::fromJSON(r_ch$path, simplifyVector = FALSE), error = function(e) NULL)
      pg <- jch$page %||% list()
      sumstat_folder_id <- ""
      for (it in pg) {
        if (grepl("summary|pgwas|gwas", it$name %||% "", ignore.case = TRUE)) {
          sumstat_folder_id <- clean_str(it$id)
          break
        }
      }
      if (nzchar(sumstat_folder_id)) {
        leg_sub <- if (nzchar(leg_pref)) paste0(leg_pref, "_synapse_children_", sumstat_folder_id, ".json") else character(0)
        r_sub <- perform_cached_request(
          session, "Tier 3", "Synapse_children_strata", ch_url,
          cache_key = paste0(slug_base, "_synapse_children_", sumstat_folder_id, ".json"),
          method = "POST",
          accept = "application/json",
          json_body = list(parentId = sumstat_folder_id, includeTypes = list("folder", "file")),
          legacy_candidates = leg_sub
        )
        if (r_sub$ok) {
          jsub <- tryCatch(jsonlite::fromJSON(r_sub$path, simplifyVector = FALSE), error = function(e) NULL)
          spg <- jsub$page %||% list()
          res$synapse_folders <- vapply(spg, function(x) clean_str(x$name), character(1L))
          # Inspect first subfolder to verify file format & anonymous download gate (HTTP 401/403)
          if (length(spg) > 0L && grepl("Folder$", spg[[1L]]$type %||% "")) {
            f1_id <- clean_str(spg[[1L]]$id)
            r_files <- perform_cached_request(
              session, "Tier 3", "Synapse_children_files", ch_url,
              cache_key = paste0(slug_base, "_synapse_files_", f1_id, ".json"),
              method = "POST",
              accept = "application/json",
              json_body = list(parentId = f1_id, includeTypes = list("file"))
            )
            if (r_files$ok) {
              jf <- tryCatch(jsonlite::fromJSON(r_files$path, simplifyVector = FALSE), error = function(e) NULL)
              fpg <- jf$page %||% list()
              if (length(fpg) > 0L) {
                fn1 <- clean_str(fpg[[1L]]$name)
                fid1 <- clean_str(fpg[[1L]]$id)
                if (grepl("\\.tar$", fn1, ignore.case = TRUE)) res$detected_file_format <- "TAR of space-delimited"
                else if (grepl("\\.tsv\\.gz$", fn1, ignore.case = TRUE)) res$detected_file_format <- "TSV.gz"
                else if (grepl("\\.txt\\.gz$", fn1, ignore.case = TRUE)) res$detected_file_format <- "TXT.gz"

                dl_url <- paste0("https://repo-prod.prod.sagebase.org/repo/v1/entity/", fid1, "/file?redirect=false")
                r_gate <- perform_cached_request(
                  session, "Tier 3", "Synapse_file_gate_check", dl_url,
                  cache_key = paste0(slug_base, "_synapse_gate_", fid1, ".json"),
                  accept = "application/json",
                  allow_transient_retry = FALSE
                )
                if (r_gate$status %in% c(401L, 403L)) {
                  res$access_level <- "Portal registration"
                  res$snippet <- sprintf("HTTP %d on /entity/%s/file (Anonymous users have READ metadata access; Synapse account required to download %s)", r_gate$status, fid1, fn1)
                } else if (identical(r_gate$status, 200L)) {
                  res$access_level <- "Open download"
                }
              }
            }
          }
        }
      }
    }
    return(res)
  }

  # 3. Check deCODE Summary Data Portal (decode.com/summarydata) — [M-6]
  if (grepl("decode\\.com/summarydata", search_blob, ignore.case = TRUE)) {
    dec_url <- "https://www.decode.com/summarydata/"
    leg_dec <- if (nzchar(leg_pref)) paste0(leg_pref, "_decode_summarydata.html") else character(0)
    r_dec <- perform_cached_request(
      session, "Tier 3", "deCODE_summarydata_index", dec_url,
      cache_key = paste0(slug_base, "_decode_summarydata.html"),
      accept = "text/html",
      legacy_candidates = leg_dec
    )
    res$repository <- "deCODE Genetics Portal"
    res$accession_or_url <- dec_url
    res$access_tier <- "Tier 3"
    res$access_status <- "SUPPLEMENT_EXTRACTED"
    res$location <- "Public Data Availability + https://www.decode.com/summarydata/"

    if (r_dec$ok && file.exists(r_dec$path)) {
      dtxt <- paste(readLines(r_dec$path, warn = FALSE), collapse = "\n")
      if (grepl("assocvariants\\.annotated\\.txt\\.gz", dtxt, ignore.case = TRUE)) {
        res$helper_files_note <- "Requires joining deCODE assocvariants.annotated.txt.gz (EAF + '!' multiallelic marker) and excluding assocvariants.excluded.txt.gz."
        res$detected_file_format <- "TXT.gz + Variant Anno"
      }
      # Probe small 2KB proteomics_readme.txt to test direct CLI download vs JS form gate ([M-6])
      rdm_url <- "https://download.decode.is/form/2021/proteomics_readme.txt"
      leg_rdm <- if (nzchar(leg_pref)) paste0(leg_pref, "_decode_proteomics_readme_resp.txt") else character(0)
      r_rdm <- perform_cached_request(
        session, "Tier 3", "deCODE_file_gate_check", rdm_url,
        cache_key = paste0(slug_base, "_decode_readme_probe.bin"),
        expected_data_file = TRUE,
        allow_transient_retry = FALSE,
        legacy_candidates = leg_rdm
      )
      if (isTRUE(r_rdm$gates$html_gate)) {
        res$access_level <- "Portal registration"
        res$snippet <- "GET download.decode.is/form/2021/proteomics_readme.txt -> HTTP 200 text/html (605B React SPA license/identity form gate)"
        res$access_note <- "Self-service browser name/email click-through license agreement on decode.com/summarydata/ (blocks unattended curl/httr2; classified as Portal registration per [M-6])."
      } else if (r_rdm$ok) {
        res$access_level <- "Open download"
        res$snippet <- "Direct HTTP 200 file download on decode.com"
      } else {
        res$access_level <- "Portal registration"
      }
    } else {
      res$access_level <- "Portal registration"
    }
    return(res)
  }

  # 4. Check Zenodo (10.5281/zenodo.{id} or zenodo.org/record(s)/{id})
  m_zen <- regmatches(search_blob, regexpr("zenodo\\.(org/records?/)?([0-9]{5,10})", search_blob, ignore.case = TRUE))
  if (length(m_zen) == 1L && nzchar(m_zen)) {
    zid <- regmatches(m_zen, regexpr("[0-9]{5,10}$", m_zen))
    z_api <- paste0("https://zenodo.org/api/records/", zid)
    r_z <- perform_cached_request(session, "Tier 3", "Zenodo_API", z_api, cache_key = paste0(slug_base, "_zenodo_", zid, ".json"), accept = "application/json")
    res$repository <- "Zenodo"
    res$accession_or_url <- paste0("https://doi.org/10.5281/zenodo.", zid)
    res$access_tier <- "Tier 3"
    res$access_status <- "SUPPLEMENT_EXTRACTED"
    res$location <- z_api
    if (r_z$ok) {
      jz <- tryCatch(jsonlite::fromJSON(r_z$path, simplifyVector = FALSE), error = function(e) NULL)
      ar <- clean_str(jz$metadata$access_right %||% "open")
      res$access_level <- if (identical(ar, "open")) "Open download" else "Upon request"
      res$snippet <- paste0("Zenodo record ", zid, " access_right=", ar)
    } else {
      res$access_level <- "Open download"
    }
    return(res)
  }

  # 5. Check Figshare / OSF ([m-4]) / GitHub ([m-4]) / dbGaP / EGA
  if (grepl("osf\\.io/([a-z0-9]{5})\\b", search_blob, ignore.case = TRUE)) {
    m_osf <- regmatches(search_blob, regexpr("osf\\.io/[a-z0-9]{5}", search_blob, ignore.case = TRUE))
    nid <- sub("^osf\\.io/", "", m_osf)
    osf_api <- paste0("https://api.osf.io/v2/nodes/", nid, "/files/osfstorage/")
    r_osf <- perform_cached_request(session, "Tier 3", "OSF_API", osf_api, cache_key = paste0(slug_base, "_osf_", nid, ".json"), accept = "application/json")
    res$repository <- "OSF"
    res$accession_or_url <- paste0("https://osf.io/", nid, "/")
    res$access_level <- if (r_osf$ok) "Open download" else "Portal registration"
    res$access_tier <- "Tier 3"
    res$access_status <- "SUPPLEMENT_EXTRACTED"
    res$location <- osf_api
    res$snippet <- paste0("OSF node ", nid)
    return(res)
  }

  if (grepl("\\bphs[0-9]{6}\\b", search_blob, ignore.case = TRUE)) {
    phs <- regmatches(search_blob, regexpr("\\bphs[0-9]{6}\\b", search_blob, ignore.case = TRUE))
    res$repository <- "dbGaP"
    res$accession_or_url <- paste0("https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=", phs)
    res$access_level <- "Controlled access"
    res$access_tier <- "Tier 3"
    res$access_status <- "SUPPLEMENT_EXTRACTED"
    res$location <- "Data Availability / NCBI DataBank"
    res$snippet <- paste0("dbGaP controlled-access accession ", phs)
    return(res)
  }

  if (grepl("\\bEGA[SD][0-9]{8,11}\\b", search_blob)) {
    ega <- regmatches(search_blob, regexpr("\\bEGA[SD][0-9]{8,11}\\b", search_blob))
    res$repository <- "EGA"
    res$accession_or_url <- paste0("https://ega-archive.org/studies/", ega)
    res$access_level <- "Controlled access"
    res$access_tier <- "Tier 3"
    res$access_status <- "SUPPLEMENT_EXTRACTED"
    res$location <- "Data Availability"
    res$snippet <- paste0("European Genome-phenome Archive ", ega)
    return(res)
  }

  # 6. Check EBI GWAS Catalog summary statistics from Tier 1f
  if (t1$gwas_catalog_count > 0L) {
    any_full <- any(vapply(t1$gwas_studies, function(s) isTRUE(s$full_summary_stats_available), logical(1L)))
    res$repository <- "EBI GWAS Catalog"
    res$accession_or_url <- paste0("https://www.ebi.ac.uk/gwas/publications/", t1$pmid)
    res$access_level <- if (any_full) "Open download" else "Not available"
    res$access_tier <- "Tier 1"
    res$access_status <- "API_EXTRACTED"
    res$location <- "EBI GWAS Catalog REST v2"
    res$snippet <- sprintf("%d study record(s); full_summary_stats_available=%s", t1$gwas_catalog_count, any_full)
    return(res)
  }

  if (grepl("github\\.com/[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", search_blob)) {
    gh <- regmatches(search_blob, regexpr("https?://github\\.com/[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", search_blob))
    if (length(gh) == 1L && nzchar(gh)) {
      res$access_note <- paste0("Code/manifest repository detected at ", gh, "; verify summary-statistics deposit.")
    }
  }

  res
}
