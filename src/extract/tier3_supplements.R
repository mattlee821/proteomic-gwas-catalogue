# Tier 3a-3c: Supplementary Material Discovery, Download & Zero-DOM Low-RSS Parsers ([C-4], [C-6])
#
# Key memory invariant ([C-4]):
#   - NEVER loop readxl::read_excel() across all sheets of a workbook.
#   - NEVER call readxl::read_excel() on any worksheet whose uncompressed XML exceeds 15 MB.
#   - Inspect ZIP directory + xl/workbook.xml first (< 5 ms), stream only the first 64 KB of each
#     sheet via unz(), and resolve sharedStrings via 1 MB chunked streaming (< 10 MB R heap).

MAX_READXL_SHEET_XML_BYTES <- 15e6

xml_unescape_basic <- function(s) {
  s <- gsub("&lt;", "<", s, fixed = TRUE)
  s <- gsub("&gt;", ">", s, fixed = TRUE)
  s <- gsub("&quot;", "\"", s, fixed = TRUE)
  s <- gsub("&apos;", "'", s, fixed = TRUE)
  s <- gsub("&amp;", "&", s, fixed = TRUE)
  s
}

col_letters_to_int <- function(cell_ref) {
  lets <- gsub("[^A-Z]", "", toupper(cell_ref))
  if (!nzchar(lets)) return(NA_integer_)
  chars <- utf8ToInt(lets) - 64L
  val <- 0L
  for (ch in chars) val <- val * 26L + ch
  val
}

inspect_xlsx_zip_manifest <- function(xlsx_path) {
  zinfo <- utils::unzip(xlsx_path, list = TRUE)
  ss_len <- zinfo$Length[zinfo$Name == "xl/sharedStrings.xml"]
  ss_bytes <- if (length(ss_len) > 0L) as.numeric(ss_len[[1L]]) else 0

  con_wb <- unz(xlsx_path, "xl/workbook.xml", open = "rb")
  wb_raw <- readBin(con_wb, "raw", n = as.integer(zinfo$Length[zinfo$Name == "xl/workbook.xml"][[1L]]))
  close(con_wb)
  wb_doc <- xml2::read_xml(wb_raw)

  con_rel <- unz(xlsx_path, "xl/_rels/workbook.xml.rels", open = "rb")
  rel_raw <- readBin(con_rel, "raw", n = as.integer(zinfo$Length[zinfo$Name == "xl/_rels/workbook.xml.rels"][[1L]]))
  close(con_rel)
  rel_doc <- xml2::read_xml(rel_raw)

  rel_nodes <- xml2::xml_find_all(rel_doc, "//*[local-name()='Relationship']")
  rel_map <- character(0)
  for (rn in rel_nodes) {
    rid <- xml2::xml_attr(rn, "Id") %||% ""
    tgt <- xml2::xml_attr(rn, "Target") %||% ""
    tgt_norm <- sub("^/?(xl/)?", "xl/", tgt)
    if (nzchar(rid) && nzchar(tgt_norm)) rel_map[[rid]] <- tgt_norm
  }

  sh_nodes <- xml2::xml_find_all(wb_doc, "//*[local-name()='sheet']")
  rows <- list()
  for (i in seq_along(sh_nodes)) {
    sn <- sh_nodes[[i]]
    nm <- xml2::xml_attr(sn, "name") %||% paste0("Sheet", i)
    attrs <- xml2::xml_attrs(sn)
    rid_idx <- grep("(^|:)id$", names(attrs))
    rid <- if (length(rid_idx) > 0L) unname(attrs[[rid_idx[[1L]]]]) else ""
    xml_entry <- if (nzchar(rid) && rid %in% names(rel_map)) rel_map[[rid]] else paste0("xl/worksheets/sheet", i, ".xml")
    ub <- zinfo$Length[zinfo$Name == xml_entry]
    rows[[i]] <- data.frame(
      sheet_idx = i,
      sheet_name = nm,
      xml_entry = xml_entry,
      uncompressed_bytes = if (length(ub) > 0L) as.numeric(ub[[1L]]) else 0,
      stringsAsFactors = FALSE
    )
  }
  list(sheets = do.call(rbind, rows), shared_strings_bytes = ss_bytes)
}

stream_sheet_head_raw <- function(xlsx_path, xml_entry, entry_bytes, max_bytes = 65536L, max_rows = 35L) {
  n_read <- as.integer(min(entry_bytes, max_bytes))
  if (is.na(n_read) || n_read <= 0L) return(list(rows = list(), ss_needed = integer(0)))
  con <- unz(xlsx_path, xml_entry, open = "rb")
  raw_b <- tryCatch(readBin(con, "raw", n = n_read), finally = close(con))
  raw_b <- raw_b[raw_b != as.raw(0L)]
  txt <- iconv(tryCatch(rawToChar(raw_b), error = function(e) ""), from = "UTF-8", to = "UTF-8", sub = "") %||% ""
  if (!nzchar(txt)) return(list(rows = list(), ss_needed = integer(0)))

  row_matches <- regmatches(txt, gregexpr("<row\\b[^>]*>[\\s\\S]*?</row>", txt, perl = TRUE))[[1L]]
  if (length(row_matches) == 0L) return(list(rows = list(), ss_needed = integer(0)))
  row_matches <- head(row_matches, max_rows)

  parsed_rows <- list()
  ss_needed <- integer(0)

  for (ri in seq_along(row_matches)) {
    rm_str <- row_matches[[ri]]
    r_num_m <- regmatches(rm_str, regexpr("(?<=\\br=\")[0-9]+", rm_str, perl = TRUE))
    r_num <- if (length(r_num_m) == 1L) as.integer(r_num_m) else ri

    c_matches <- regmatches(rm_str, gregexpr("<c\\b[^/>]*(?:/>|>[\\s\\S]*?</c>)", rm_str, perl = TRUE))[[1L]]
    cells <- list()
    for (ci in seq_along(c_matches)) {
      cm <- c_matches[[ci]]
      ref_m <- regmatches(cm, regexpr("(?<=\\br=\")[A-Z]+[0-9]+", cm, perl = TRUE))
      col_i <- if (length(ref_m) == 1L) col_letters_to_int(ref_m) else ci
      t_m <- regmatches(cm, regexpr("(?<=\\bt=\")[a-zA-Z]+", cm, perl = TRUE))
      c_type <- if (length(t_m) == 1L) t_m else "n"

      if (identical(c_type, "s")) {
        v_m <- regmatches(cm, regexpr("(?<=<v>)[0-9]+(?=</v>)", cm, perl = TRUE))
        if (length(v_m) == 1L) {
          idx1 <- as.integer(v_m) + 1L
          ss_needed <- c(ss_needed, idx1)
          cells[[length(cells) + 1L]] <- list(col = col_i, is_ss = TRUE, ss_idx = idx1, val = "")
        }
      } else if (identical(c_type, "inlineStr")) {
        t_vals <- regmatches(cm, gregexpr("(?<=<t>|<t xml:space=\"preserve\">)[\\s\\S]*?(?=</t>)", cm, perl = TRUE))[[1L]]
        cells[[length(cells) + 1L]] <- list(col = col_i, is_ss = FALSE, ss_idx = NA_integer_, val = xml_unescape_basic(paste(t_vals, collapse = "")))
      } else {
        v_m <- regmatches(cm, regexpr("(?<=<v>)[\\s\\S]*?(?=</v>)", cm, perl = TRUE))
        if (length(v_m) == 1L) {
          cells[[length(cells) + 1L]] <- list(col = col_i, is_ss = FALSE, ss_idx = NA_integer_, val = xml_unescape_basic(v_m))
        }
      }
    }
    parsed_rows[[length(parsed_rows) + 1L]] <- list(row_num = r_num, cells = cells)
  }

  list(rows = parsed_rows, ss_needed = unique(ss_needed))
}

stream_shared_strings_subset <- function(xlsx_path, ss_bytes, needed_1based) {
  out_env <- new.env(hash = TRUE, parent = emptyenv())
  needed_set <- sort(unique(as.integer(needed_1based[!is.na(needed_1based) & needed_1based > 0L])))
  if (length(needed_set) == 0L || ss_bytes <= 0) return(out_env)
  max_target <- max(needed_set)
  need_env <- new.env(hash = TRUE, parent = emptyenv())
  for (k in needed_set) need_env[[as.character(k)]] <- TRUE

  con <- unz(xlsx_path, "xl/sharedStrings.xml", open = "rb")
  on.exit(close(con), add = TRUE)

  carry <- ""
  si_counter <- 0L
  chunk_sz <- 1048576L # 1 MB streaming buffer

  repeat {
    raw_b <- readBin(con, "raw", n = chunk_sz)
    if (length(raw_b) == 0L) break
    raw_b <- raw_b[raw_b != as.raw(0L)]
    chunk_txt <- iconv(tryCatch(rawToChar(raw_b), error = function(e) ""), from = "UTF-8", to = "UTF-8", sub = "") %||% ""
    buf <- paste0(carry, chunk_txt)

    # O(N) fixed-string scan for last closing </si> (avoids O(N^2) regex backtracking!)
    end_locs <- gregexpr("</si>", buf, fixed = TRUE)[[1L]]
    if (identical(end_locs[[1L]], -1L)) {
      carry <- buf
      next
    }
    last_end <- end_locs[[length(end_locs)]]
    split_pos <- as.integer(last_end) + 4L # nchar("</si>") - 1L
    complete_part <- substr(buf, 1L, split_pos)
    carry <- if (split_pos < nchar(buf)) substr(buf, split_pos + 1L, nchar(buf)) else ""

    si_blocks <- strsplit(complete_part, "</si>", fixed = TRUE)[[1L]]
    if (length(si_blocks) > 0L) {
      for (blk in si_blocks) {
        if (!grepl("<si", blk, fixed = TRUE)) next
        si_counter <- si_counter + 1L
        key_s <- as.character(si_counter)
        if (!is.null(need_env[[key_s]])) {
          t_parts <- regmatches(blk, gregexpr("<t(?:\\s+[^>]*)?>([\\s\\S]*?)</t>", blk, perl = TRUE))[[1L]]
          t_clean <- gsub("^<t(?:\\s+[^>]*)?>|</t>$", "", t_parts, perl = TRUE)
          out_env[[key_s]] <- clean_str(xml_unescape_basic(paste(t_clean, collapse = "")))
        }
        if (si_counter >= max_target) break
      }
    }
    if (si_counter >= max_target) break
  }

  out_env
}

materialize_streamed_rows <- function(raw_head, ss_map) {
  lapply(raw_head$rows, function(r) {
    if (length(r$cells) == 0L) return(list(row_num = r$row_num, values = character(0)))
    max_c <- max(vapply(r$cells, function(cl) cl$col %||% 1L, integer(1L)))
    max_c <- min(max_c, 80L)
    vec <- rep("", max_c)
    for (cl in r$cells) {
      ci <- cl$col
      if (is.na(ci) || ci < 1L || ci > max_c) next
      vec[[ci]] <- if (isTRUE(cl$is_ss)) {
        ss_map[[as.character(cl$ss_idx)]] %||% ""
      } else {
        clean_str(cl$val)
      }
    }
    list(row_num = r$row_num, values = vec)
  })
}

parse_xlsx_supplement_low_rss <- function(xlsx_path, filename_label) {
  info <- inspect_xlsx_zip_manifest(xlsx_path)
  sh_df <- info$sheets

  # Stream first 64 KB of every sheet (milliseconds, zero DOM!)
  raw_heads <- vector("list", nrow(sh_df))
  all_needed_ss <- integer(0)
  for (i in seq_len(nrow(sh_df))) {
    rh <- stream_sheet_head_raw(xlsx_path, sh_df$xml_entry[[i]], sh_df$uncompressed_bytes[[i]], max_bytes = 65536L, max_rows = 35L)
    raw_heads[[i]] <- rh
    all_needed_ss <- c(all_needed_ss, rh$ss_needed)
  }

  # Also include TOC sheet full rows if a Contents/Index sheet exists (usually < 100 KB)
  toc_idx <- which(grepl("^(contents|index|readme|overview|legend|list)", tolower(sh_df$sheet_name)) & sh_df$uncompressed_bytes <= 2e6)
  if (length(toc_idx) > 0L) {
    ti <- toc_idx[[1L]]
    rh_toc <- stream_sheet_head_raw(xlsx_path, sh_df$xml_entry[[ti]], sh_df$uncompressed_bytes[[ti]], max_bytes = 524288L, max_rows = 120L)
    raw_heads[[ti]] <- rh_toc
    all_needed_ss <- c(all_needed_ss, rh_toc$ss_needed)
  }

  ss_map <- stream_shared_strings_subset(xlsx_path, info$shared_strings_bytes, all_needed_ss)
  mat_heads <- lapply(raw_heads, materialize_streamed_rows, ss_map = ss_map)
  names(mat_heads) <- sh_df$sheet_name

  # Build TOC description lookup if TOC sheet present
  toc_descriptions <- list()
  if (length(toc_idx) > 0L) {
    for (r in mat_heads[[toc_idx[[1L]]]]) {
      line <- paste(r$values[nzchar(r$values)], collapse = " | ")
      for (snm in sh_df$sheet_name) {
        if (grepl(paste0("^", snm, "\\b"), line, ignore.case = TRUE)) {
          toc_descriptions[[snm]] <- line
        }
      }
    }
  }

  out <- list(
    sheet_count = nrow(sh_df),
    n_assays = NA_integer_,
    n_unique_proteins = NA_integer_,
    platform_version = "",
    assay_location = "",
    assay_snippet = "",
    has_beta_se = FALSE,
    has_eaf = FALSE,
    pqtl_cols_location = "",
    pqtl_cols_snippet = "",
    multiallelic_bang = FALSE,
    hg38_header_explicit = FALSE,
    mixed_id_build_note = "",
    variant_probes = data.frame(chrom = character(0), pos = integer(0), rsid = character(0), source_loc = character(0), stringsAsFactors = FALSE)
  )

  # 1. Identify Assay Manifest Sheet (e.g., ST3 in Sun 2023; ST01 in Ferkingstad 2021)
  assay_sheet <- ""
  for (i in seq_len(nrow(sh_df))) {
    snm <- sh_df$sheet_name[[i]]
    top_txt <- paste(
      toc_descriptions[[snm]] %||% "",
      paste(unlist(lapply(head(mat_heads[[i]], 6L), function(r) paste(r$values, collapse = " "))), collapse = " ")
    )
    if (grepl("details of proteins on the olink|aptamers included in the somascan|\\bOlink ID\\b.*\\bUniProt\\b|\\bSeqId\\b.*\\bUniProt\\b", top_txt, ignore.case = TRUE, perl = TRUE)) {
      assay_sheet <- snm
      break
    }
  }

  if (nzchar(assay_sheet)) {
    sh_row <- sh_df[sh_df$sheet_name == assay_sheet, , drop = FALSE]
    head_rows <- mat_heads[[assay_sheet]]
    r1_txt <- if (length(head_rows) > 0L) paste(head_rows[[1L]]$values[nzchar(head_rows[[1L]]$values)], collapse = " ") else ""
    if (grepl("Somascan v4", r1_txt, ignore.case = TRUE)) out$platform_version <- "SomaScan v4 (4.9k)"
    if (grepl("Olink Explore 3k|Olink Explore 3072", r1_txt, ignore.case = TRUE)) out$platform_version <- "Olink Explore 3072"

    # Find header row index (first row with UniProt + Gene/SeqId/Olink ID)
    hdr_row_num <- 1L
    for (r in head(head_rows, 10L)) {
      joined <- paste(r$values, collapse = " | ")
      if (grepl("UniProt", joined, ignore.case = TRUE) && grepl("SeqId|Olink ID|Assay Target|Gene", joined, ignore.case = TRUE)) {
        hdr_row_num <- r$row_num
        break
      }
    }

    # Memory gate ([C-4]): ONLY call readxl::read_excel() if worksheet XML <= 15 MB!
    if (sh_row$uncompressed_bytes[[1L]] <= MAX_READXL_SHEET_XML_BYTES) {
      df_a <- suppressMessages(readxl::read_excel(xlsx_path, sheet = assay_sheet, skip = max(0L, hdr_row_num - 1L)))
      nms <- names(df_a)
      inc_col <- grep("^Included", nms, ignore.case = TRUE, value = TRUE)
      if (length(inc_col) > 0L) {
        keep <- !is.na(df_a[[inc_col[[1L]]]]) & tolower(trimws(as.character(df_a[[inc_col[[1L]]]]))) %in% c("yes", "true", "1")
        df_sub <- df_a[keep, , drop = FALSE]
      } else {
        df_sub <- df_a
      }

      if ("Assay Target" %in% nms && "Olink ID" %in% nms) {
        n_oid <- length(unique(stats::na.omit(df_sub$`Olink ID`)))
        n_tgt <- length(unique(stats::na.omit(df_sub$`Assay Target`)))
        out$n_assays <- n_tgt # Catalogue tracks unique protein targets/assays reported in headline (2923; 2941 panel analytes noted)
        out$n_unique_proteins <- n_tgt
        out$assay_location <- sprintf("%s sheet %s (%d rows)", filename_label, assay_sheet, nrow(df_sub))
        out$assay_snippet <- sprintf("%d Olink IDs across 8 panels capturing %d unique Assay Targets (%s)", n_oid, n_tgt, r1_txt)
      } else if ("SeqId" %in% nms) {
        n_seq <- length(unique(stats::na.omit(df_sub$SeqId)))
        up_vals <- as.character(df_sub$UniProt %||% character(0))
        up_vals <- up_vals[!is.na(up_vals) & nzchar(trimws(up_vals)) & up_vals != "NA"]
        n_up <- length(unique(up_vals))
        out$n_assays <- n_seq
        out$n_unique_proteins <- n_up
        out$assay_location <- sprintf("%s sheet %s (Included=='Yes')", filename_label, assay_sheet)
        out$assay_snippet <- sprintf("%d QC-passing human aptamers (from %d total array rows); %d unique non-NA UniProt IDs in table", n_seq, nrow(df_a), n_up)
      }
      rm(df_a, df_sub)
      gc(verbose = FALSE)
    }
  }

  # 2. Inspect pQTL Summary / Replication Sheets purely from our 64 KB streamed rows (ZERO readxl DOM!)
  probe_rows <- list()
  for (i in seq_len(nrow(sh_df))) {
    snm <- sh_df$sheet_name[[i]]
    hrows <- mat_heads[[i]]
    if (length(hrows) < 4L) next

    # Locate a header row mentioning pQTL / variant / Chrom / Pos / BETA / rsid
    hdr_idx <- NA_integer_
    for (k in seq_len(min(10L, length(hrows)))) {
      v <- hrows[[k]]$values
      jtxt <- paste(v[nzchar(v)], collapse = " | ")
      if (grepl("\\b(CHROM|Chrom|chr)\\b", jtxt) && grepl("\\b(GENPOS|Pos|pos|rsID|rsid|variant)\\b", jtxt)) {
        hdr_idx <- k
        break
      }
    }
    if (is.na(hdr_idx)) next

    hdr_vals <- hrows[[hdr_idx]]$values
    hdr_joined <- paste(hdr_vals[nzchar(hdr_vals)], collapse = " | ")

    if (grepl("\\b(BETA|Beta|Effect|DECODE\\s*Beta)\\b", hdr_joined, ignore.case = TRUE) &&
        grepl("\\b(SE|StdErr)\\b", hdr_joined, ignore.case = FALSE)) {
      out$has_beta_se <- TRUE
      if (!nzchar(out$pqtl_cols_location)) {
        out$pqtl_cols_location <- sprintf("%s sheet %s", filename_label, snm)
        out$pqtl_cols_snippet <- substr(hdr_joined, 1, 210)
      }
    }
    if (grepl("\\b(A1FREQ|EAF|effectAlleleFreq)\\b", hdr_joined, ignore.case = TRUE)) {
      out$has_eaf <- TRUE
    }
    if (grepl("GENPOS\\s*\\(hg38\\)", hdr_joined, ignore.case = TRUE)) {
      out$hg38_header_explicit <- TRUE
      if (grepl("hg37", hdr_joined, ignore.case = TRUE)) {
        out$mixed_id_build_note <- "Supplement Variant ID embeds hg37 position while GENPOS is hg38."
      }
    }

    # Check data rows in our 64 KB head for '!' multiallelic marker and (chrom, pos, rsid) sentinel probes
    chr_cols <- which(grepl("^(CHROM|Chrom|chr\\s*\\(var\\.\\)|chr)$", trimws(hdr_vals), ignore.case = TRUE))
    pos_cols <- which(grepl("^(GENPOS.*|Pos|pos\\s*\\(var\\.\\)|position)$", trimws(hdr_vals), ignore.case = TRUE))
    rs_cols <- which(grepl("^(rsID|rsid|SNP|variant)$", trimws(hdr_vals), ignore.case = TRUE))

    for (dr in hrows[(hdr_idx + 1L):length(hrows)]) {
      dv <- dr$values
      if (any(dv == "!")) out$multiallelic_bang <- TRUE
      if (length(chr_cols) > 0L && length(pos_cols) > 0L && length(rs_cols) > 0L && length(probe_rows) < 25L) {
        cc <- chr_cols[[length(chr_cols)]] # prefer variant chr over protein chr if both exist
        pc <- pos_cols[[1L]]
        rc <- rs_cols[[1L]]
        if (max(cc, pc, rc) <= length(dv)) {
          rs_v <- trimws(dv[[rc]])
          if (grepl("^rs[1-9][0-9]+$", rs_v)) {
            probe_rows[[length(probe_rows) + 1L]] <- data.frame(
              chrom = dv[[cc]],
              pos = suppressWarnings(as.integer(as.numeric(dv[[pc]]))),
              rsid = rs_v,
              source_loc = sprintf("%s sheet %s", filename_label, snm),
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }

  if (length(probe_rows) > 0L) {
    out$variant_probes <- do.call(rbind, probe_rows)
  }
  out
}

run_tier3_supplements <- function(t1, t2, session, repo_root = ".") {
  leg_pref <- legacy_prefix_for_study(t1$pmid, t1$doi)
  slug_base <- t1$study_id

  out <- list(
    sample_matrix = list(value = "", ambiguous = TRUE, snippet = "", location = ""),
    xlsx_parsed = NULL,
    build_resolved = list(resolved = FALSE, genome_build = "Unknown", evidence_snippet = "", location = "", mixed_build_warning = ""),
    bioc_supp_notes = character(0)
  )

  # 3a. NCBI supplmat.cgi BioC JSON for PMC OA PDF supplements (ONLY valid PMC[0-9]+ IDs — [C-6])
  if (grepl("^PMC[0-9]+$", t1$pmcid) && identical(t1$has_suppl, "Y")) {
    list_url <- sprintf("https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/supplmat.cgi/BioC_JSON/%s/LIST", t1$pmcid)
    leg_sm_l <- if (nzchar(leg_pref)) paste0(leg_pref, "_ncbi_supplmat_list.txt") else character(0)
    r_sml <- perform_cached_request(
      session, "Tier 3", "NCBI_supplmat_LIST", list_url,
      cache_key = paste0(slug_base, "_ncbi_supplmat_list.json"),
      legacy_candidates = leg_sm_l
    )
    if (r_sml$ok && file.exists(r_sml$path)) {
      sm_items <- tryCatch(jsonlite::fromJSON(r_sml$path, simplifyVector = FALSE), error = function(e) list())
      # Fetch first 2 PDF BioC JSON supplements (typically Supplementary Information & Reporting Summary)
      for (it in head(sm_items, 2L)) {
        idx <- clean_str(it$num %||% "")
        if (!nzchar(idx)) next
        sm_u <- sprintf("https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/supplmat.cgi/BioC_json/%s/%s", t1$pmcid, idx)
        leg_smi <- if (nzchar(leg_pref)) paste0(leg_pref, "_supplmat_bioc_", idx, ".json") else character(0)
        r_smi <- perform_cached_request(
          session, "Tier 3", paste0("NCBI_supplmat_BioC_", idx), sm_u,
          cache_key = paste0(slug_base, "_supplmat_bioc_", idx, ".json"),
          accept = "application/json",
          legacy_candidates = leg_smi
        )
        if (r_smi$ok && file.exists(r_smi$path)) {
          sm_txt <- paste(readLines(r_smi$path, warn = FALSE), collapse = " ")
          orig_fn <- basename(clean_str(it$original_file %||% paste0("supplmat_", idx)))
          if (!nzchar(out$sample_matrix$value) || isTRUE(t2$sample_matrix$ambiguous)) {
            m_sm <- detect_sample_matrix_12(sm_txt, sprintf("NCBI supplmat BioC (%s)", orig_fn))
            if (nzchar(m_sm$value) && !isTRUE(m_sm$ambiguous)) {
              out$sample_matrix <- m_sm
            }
          }
        }
      }
    }
  }

  # 3b. Discover & download .xlsx supplement via BioStudies (S-EPMC* ONLY — never S-EPPR* [C-6]) or Springer MOESM
  xlsx_urls <- character(0)
  xlsx_names <- character(0)

  if (grepl("^S-EPMC[0-9]+$", t1$biostudies_acc)) {
    bs_url <- paste0("https://www.ebi.ac.uk/biostudies/api/v1/studies/", t1$biostudies_acc)
    leg_bs <- if (nzchar(leg_pref)) paste0(leg_pref, "_biostudies.json") else character(0)
    r_bs <- perform_cached_request(
      session, "Tier 3", "BioStudies_API", bs_url,
      cache_key = paste0(slug_base, "_biostudies.json"),
      accept = "application/json",
      legacy_candidates = leg_bs
    )
    if (r_bs$ok && file.exists(r_bs$path)) {
      jbs <- tryCatch(jsonlite::fromJSON(r_bs$path, simplifyVector = FALSE), error = function(e) NULL)
      fls <- jbs$section$files %||% list()
      for (fl in fls) {
        fp <- clean_str(fl$path %||% "")
        fsz <- as.numeric(fl$size %||% 0)
        if (grepl("\\.xlsx$", fp, ignore.case = TRUE) && fsz > 0 && fsz <= 85e6) {
          xlsx_urls <- c(xlsx_urls, paste0("https://www.ebi.ac.uk/biostudies/files/", t1$biostudies_acc, "/", fp))
          xlsx_names <- c(xlsx_names, fp)
        }
      }
    }
  }

  if (length(xlsx_urls) == 0L && length(t2$public_moesm_urls) > 0L) {
    m_xl <- grep("\\.xlsx$", t2$public_moesm_urls, ignore.case = TRUE, value = TRUE)
    for (u in m_xl) {
      xlsx_urls <- c(xlsx_urls, u)
      xlsx_names <- c(xlsx_names, basename(u))
    }
  }

  # Parse first discovered .xlsx workbook using Zero-DOM Low-RSS parser ([C-4])
  if (length(xlsx_urls) > 0L) {
    xu <- xlsx_urls[[1L]]
    xnm <- xlsx_names[[1L]]
    leg_xl <- if (nzchar(leg_pref)) paste0(leg_pref, "_", xnm) else character(0)
    r_xl <- perform_cached_request(
      session, "Tier 3", "Supplement_XLSX", xu,
      cache_key = paste0(slug_base, "_", xnm),
      timeout_sec = 180,
      expected_data_file = TRUE,
      legacy_candidates = leg_xl
    )
    if (r_xl$ok && file.exists(r_xl$path)) {
      out$xlsx_parsed <- tryCatch(
        parse_xlsx_supplement_low_rss(r_xl$path, xnm),
        error = function(e) {
          warning("Low-RSS XLSX parser error on ", xnm, ": ", conditionMessage(e))
          NULL
        }
      )
    }
  }

  # 3c. Run Deterministic Coordinate Genome-Build Resolver if variant probes found
  if (!is.null(out$xlsx_parsed) && nrow(out$xlsx_parsed$variant_probes) > 0L) {
    out$build_resolved <- resolve_genome_build_coords(out$xlsx_parsed$variant_probes, session, repo_root = repo_root)
  }

  out
}
