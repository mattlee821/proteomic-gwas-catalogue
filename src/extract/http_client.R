# Polite HTTP client, disk cache, rate-limiter, and zero-bypass gate detector
# for the Proteomic GWAS Catalogue 4-tier extraction pipeline.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

clean_str <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[[1L]])) return("")
  s <- as.character(x[[1L]])
  # Normalize Unicode thin/non-breaking spaces (e.g. U+2009 in Nature JATS XML `n = 931`) to ASCII space
  s <- gsub("[\u00a0\u2000-\u200b\u202f\u205f\u3000\r\n\t]+", " ", s, perl = TRUE)
  s <- gsub("  +", " ", s)
  trimws(s)
}

pgwas_default_ua <- function() {
  paste0(
    "proteomic-gwas-catalogue/0.4.0 ",
    "(https://github.com/matthewlee/proteomic-gwas-catalogue; polite scientific curation client)"
  )
}

new_http_session <- function(cache_dir,
                             refresh_cache = FALSE,
                             strict_oa_only = TRUE,
                             min_interval_sec = 0.45,
                             contact_email = Sys.getenv("PGWAS_CONTACT_EMAIL", "curator@proteomic-gwas-catalogue.org")) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  env <- new.env(parent = emptyenv())
  env$cache_dir <- normalizePath(cache_dir, mustWork = TRUE)
  env$refresh_cache <- isTRUE(refresh_cache)
  env$strict_oa_only <- isTRUE(strict_oa_only)
  env$min_interval_sec <- as.numeric(min_interval_sec)
  env$biorxiv_min_interval_sec <- 3.0
  env$contact_email <- clean_str(contact_email)
  env$last_req_time <- 0
  env$last_biorxiv_time <- 0
  env$http_log <- list()
  env
}

append_http_log <- function(session, tier, source, method, url, status, content_type, bytes, cached, gate_flag, note = "") {
  session$http_log[[length(session$http_log) + 1L]] <- data.frame(
    timestamp_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    tier = clean_str(tier),
    source = clean_str(source),
    method = clean_str(method),
    status = ifelse(is.na(status), NA_integer_, as.integer(status)),
    content_type = clean_str(content_type),
    bytes = as.numeric(bytes %||% 0),
    cached = isTRUE(cached),
    gate_flag = clean_str(gate_flag),
    note = clean_str(note),
    url = clean_str(url),
    stringsAsFactors = FALSE
  )
}

export_http_log <- function(session, out_path) {
  if (length(session$http_log) == 0L) {
    df <- data.frame(
      timestamp_utc = character(0), tier = character(0), source = character(0),
      method = character(0), status = integer(0), content_type = character(0),
      bytes = numeric(0), cached = logical(0), gate_flag = character(0),
      note = character(0), url = character(0), stringsAsFactors = FALSE
    )
  } else {
    df <- do.call(rbind, session$http_log)
  }
  utils::write.table(df, file = out_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  invisible(df)
}

throttle_session <- function(session, url = "") {
  if (grepl("(biorxiv|medrxiv)\\.org", url, ignore.case = TRUE) && !grepl("api\\.biorxiv\\.org", url, ignore.case = TRUE)) {
    now_bx <- as.numeric(Sys.time())
    el_bx <- now_bx - (session$last_biorxiv_time %||% 0)
    min_bx <- session$biorxiv_min_interval_sec %||% 3.0
    if (el_bx < min_bx) {
      Sys.sleep(min_bx - el_bx)
    }
    session$last_biorxiv_time <- as.numeric(Sys.time())
  }
  now <- as.numeric(Sys.time())
  elapsed <- now - session$last_req_time
  if (elapsed < session$min_interval_sec) {
    Sys.sleep(session$min_interval_sec - elapsed)
  }
  session$last_req_time <- as.numeric(Sys.time())
}

inspect_payload_gates <- function(file_path, content_type = "", expected_data_file = FALSE) {
  if (!file.exists(file_path)) {
    return(list(antibot = FALSE, html_gate = FALSE, bioc_error = FALSE, epmc_error_xml = FALSE, head_text = ""))
  }
  sz <- file.info(file_path)$size
  if (is.na(sz) || sz == 0) {
    return(list(antibot = FALSE, html_gate = FALSE, bioc_error = FALSE, epmc_error_xml = FALSE, head_text = ""))
  }
  con <- file(file_path, open = "rb")
  on.exit(close(con), add = TRUE)
  raw_head <- readBin(con, what = "raw", n = min(4096L, as.integer(sz)))

  # Fast-path binary magic headers (ZIP/XLSX: PK\x03\x04, GZIP: \x1f\x8b, PDF: %PDF-)
  if (length(raw_head) >= 4L) {
    is_zip <- identical(raw_head[1:4], as.raw(c(0x50L, 0x4bL, 0x03L, 0x04L)))
    is_gz  <- identical(raw_head[1:2], as.raw(c(0x1fL, 0x8bL)))
    is_pdf <- identical(raw_head[1:4], charToRaw("%PDF"))
    if (is_zip || is_gz || is_pdf) {
      return(list(antibot = FALSE, html_gate = FALSE, bioc_error = FALSE, epmc_error_xml = FALSE, head_text = ""))
    }
  }

  # Strip NULs and sanitize invalid UTF-8 bytes before regex checks
  raw_no_nul <- raw_head[raw_head != as.raw(0L)]
  head_txt <- tryCatch(rawToChar(raw_no_nul), error = function(e) "")
  head_txt <- iconv(head_txt, from = "UTF-8", to = "UTF-8", sub = "") %||% ""

  antibot <- grepl(
    "Just a moment\\.\\.\\.|cf-mitigated|challenge-platform|Attention Required! \\| Cloudflare|g-recaptcha|hcaptcha|Verify you are human",
    head_txt, ignore.case = TRUE, perl = TRUE
  )
  is_html <- grepl("text/html", content_type, ignore.case = TRUE) ||
    grepl("^\\s*(<!doctype html|<html)", head_txt, ignore.case = TRUE, perl = TRUE)
  html_gate <- isTRUE(expected_data_file) && is_html
  bioc_err <- grepl("^\\s*\\[Error\\]\\s*:", head_txt, perl = TRUE)
  epmc_err <- grepl("<(ns[0-9]+:)?errorBean\\b|<errMsg>", head_txt, perl = TRUE)

  list(
    antibot = antibot,
    html_gate = html_gate,
    bioc_error = bioc_err,
    epmc_error_xml = epmc_err,
    head_text = head_txt
  )
}

resolve_cached_file <- function(session, cache_key, legacy_candidates = character(0)) {
  primary <- file.path(session$cache_dir, cache_key)
  meta_path <- paste0(primary, ".meta.json")
  if (!session$refresh_cache && file.exists(primary) && file.info(primary)$size > 0) {
    meta <- if (file.exists(meta_path)) {
      tryCatch(jsonlite::fromJSON(meta_path, simplifyVector = TRUE), error = function(e) list())
    } else {
      list()
    }
    return(list(hit = TRUE, path = primary, meta = meta))
  }
  if (!session$refresh_cache && length(legacy_candidates) > 0L) {
    for (cand in legacy_candidates) {
      p <- file.path(session$cache_dir, cand)
      if (file.exists(p) && file.info(p)$size > 0) {
        ext <- tolower(tools::file_ext(p))
        ct_guess <- switch(ext,
          json = "application/json",
          xml = "application/xml",
          html = "text/html",
          pdf = "application/pdf",
          xlsx = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
          txt = "text/plain",
          "application/octet-stream"
        )
        return(list(hit = TRUE, path = p, meta = list(status = 200L, content_type = ct_guess, url_final = "")))
      }
    }
  }
  list(hit = FALSE, path = primary, meta = list())
}

perform_cached_request <- function(session,
                                   tier,
                                   source,
                                   url,
                                   cache_key,
                                   method = "GET",
                                   accept = NULL,
                                   json_body = NULL,
                                   timeout_sec = 30,
                                   expected_data_file = FALSE,
                                   allow_transient_retry = TRUE,
                                   legacy_candidates = character(0)) {
  cache_chk <- resolve_cached_file(session, cache_key, legacy_candidates = legacy_candidates)
  if (isTRUE(cache_chk$hit)) {
    p <- cache_chk$path
    st <- as.integer(cache_chk$meta$status %||% 200L)
    ct <- clean_str(cache_chk$meta$content_type %||% "")
    sz <- as.numeric(file.info(p)$size)
    gates <- inspect_payload_gates(p, content_type = ct, expected_data_file = expected_data_file)
    gate_lbl <- if (gates$antibot) "ANTIBOT_BLOCK" else if (gates$html_gate) "HTML_FORM_GATE" else if (gates$bioc_error) "BIOC_NO_RESULT" else if (gates$epmc_error_xml) "EPMC_ERROR_XML" else ""
    append_http_log(session, tier, source, method, url, st, ct, sz, TRUE, gate_lbl, "cache_hit")
    return(list(
      ok = (st >= 200L && st < 300L && !gates$antibot && !gates$html_gate && !gates$bioc_error && !gates$epmc_error_xml),
      status = st,
      content_type = ct,
      url_final = clean_str(cache_chk$meta$url_final %||% url),
      bytes = sz,
      path = p,
      cached = TRUE,
      gates = gates,
      error = NULL
    ))
  }

  throttle_session(session, url)
  dest <- cache_chk$path
  meta_path <- paste0(dest, ".meta.json")

  req <- httr2::request(url) |>
    httr2::req_user_agent(pgwas_default_ua()) |>
    httr2::req_method(method) |>
    httr2::req_timeout(timeout_sec) |>
    httr2::req_error(is_error = function(resp) FALSE)

  if (!is.null(accept) && nzchar(accept)) {
    req <- req |> httr2::req_headers(Accept = accept)
  }
  if (!is.null(json_body)) {
    req <- req |> httr2::req_body_json(json_body, auto_unbox = TRUE)
  }
  if (isTRUE(allow_transient_retry)) {
    # Retry at most once ONLY on transient 502/503/504; NEVER retry 401/403/404/410/429/500
    req <- req |> httr2::req_retry(
      max_tries = 2L,
      is_transient = function(resp) httr2::resp_status(resp) %in% c(502L, 503L, 504L),
      backoff = function(i) 1.5
    )
  }

  resp <- tryCatch(httr2::req_perform(req), error = function(e) e)
  if (inherits(resp, "error")) {
    err_msg <- clean_str(conditionMessage(resp))
    append_http_log(session, tier, source, method, url, NA_integer_, "", 0, FALSE, "NETWORK_ERROR", err_msg)
    return(list(
      ok = FALSE, status = NA_integer_, content_type = "", url_final = url,
      bytes = 0, path = "", cached = FALSE,
      gates = list(antibot = FALSE, html_gate = FALSE, bioc_error = FALSE, epmc_error_xml = FALSE, head_text = ""),
      error = err_msg
    ))
  }

  st <- httr2::resp_status(resp)
  ct <- tryCatch(clean_str(httr2::resp_content_type(resp)), error = function(e) "")
  uf <- clean_str(resp$url %||% url)
  raw_b <- if (httr2::resp_has_body(resp)) {
    tryCatch(httr2::resp_body_raw(resp), error = function(e) raw(0))
  } else {
    raw(0)
  }
  sz <- as.numeric(length(raw_b))

  if (sz > 0) {
    writeBin(raw_b, dest)
    meta_obj <- list(
      status = st,
      content_type = ct,
      url_requested = url,
      url_final = uf,
      bytes = sz,
      fetched_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    )
    writeLines(jsonlite::toJSON(meta_obj, auto_unbox = TRUE, pretty = TRUE), meta_path)
  }

  gates <- if (sz > 0) {
    inspect_payload_gates(dest, content_type = ct, expected_data_file = expected_data_file)
  } else {
    list(antibot = FALSE, html_gate = FALSE, bioc_error = FALSE, epmc_error_xml = FALSE, head_text = "")
  }

  gate_lbl <- if (gates$antibot || st %in% c(401L, 403L, 429L)) {
    if (gates$antibot) "ANTIBOT_BLOCK" else paste0("HTTP_", st)
  } else if (gates$html_gate) {
    "HTML_FORM_GATE"
  } else if (gates$bioc_error) {
    "BIOC_NO_RESULT"
  } else if (gates$epmc_error_xml) {
    "EPMC_ERROR_XML"
  } else if (st == 410L) {
    "HTTP_410_DEPRECATED"
  } else {
    ""
  }

  append_http_log(session, tier, source, method, url, st, ct, sz, FALSE, gate_lbl, "")
  list(
    ok = (st >= 200L && st < 300L && sz > 0 && !gates$antibot && !gates$html_gate && !gates$bioc_error && !gates$epmc_error_xml),
    status = st,
    content_type = ct,
    url_final = uf,
    bytes = sz,
    path = if (sz > 0) dest else "",
    cached = FALSE,
    gates = gates,
    error = NULL
  )
}
