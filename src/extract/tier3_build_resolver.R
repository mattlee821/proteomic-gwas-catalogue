# Tier 3e: Deterministic Coordinate Genome-Build Resolver
# Resolves GRCh36 / GRCh37 / GRCh38 / T2T-CHM13 from supplement (chrom, pos, rsid) rows
# using 1) offline bundled catalogue/schema/sentinel_build_anchors.tsv (< 1 ms) and
# 2) NCBI Variation Services SPDI API (/variation/v0/refsnp/{id}) fallback.

normalize_chrom <- function(x) {
  c_str <- toupper(clean_str(x))
  c_str <- sub("^CHR", "", c_str)
  if (identical(c_str, "23")) "X" else c_str
}

load_sentinel_anchors <- function(repo_root = ".") {
  p <- file.path(repo_root, "catalogue", "schema", "sentinel_build_anchors.tsv")
  if (!file.exists(p)) return(NULL)
  utils::read.delim(p, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
}

resolve_genome_build_coords <- function(variant_df, session, repo_root = ".") {
  res_empty <- list(
    resolved = FALSE,
    genome_build = "Unknown",
    matched_count = 0L,
    tested_count = 0L,
    evidence_snippet = "",
    location = "",
    mixed_build_warning = ""
  )
  if (is.null(variant_df) || nrow(variant_df) == 0L) return(res_empty)

  # Clean candidate rows: require valid rsID and positive integer position
  variant_df$rsid_clean <- tolower(trimws(as.character(variant_df$rsid)))
  variant_df$pos_int <- suppressWarnings(as.integer(as.numeric(variant_df$pos)))
  variant_df$chr_clean <- vapply(variant_df$chrom, normalize_chrom, character(1L))
  valid <- !is.na(variant_df$pos_int) & variant_df$pos_int > 0L & grepl("^rs[1-9][0-9]+$", variant_df$rsid_clean)
  vdf <- variant_df[valid, , drop = FALSE]
  vdf <- vdf[!duplicated(vdf$rsid_clean), , drop = FALSE]
  if (nrow(vdf) == 0L) return(res_empty)

  counts <- c(GRCh36 = 0L, GRCh37 = 0L, GRCh38 = 0L, `T2T-CHM13` = 0L)
  ev_parts <- character(0)
  tested <- 0L

  # Pass 1: Offline lookup in sentinel_build_anchors.tsv (< 1 ms, zero network I/O)
  anchors <- load_sentinel_anchors(repo_root)
  if (!is.null(anchors) && nrow(anchors) > 0L) {
    m_idx <- match(vdf$rsid_clean, tolower(anchors$rsid))
    hit_rows <- which(!is.na(m_idx))
    for (i in head(hit_rows, 8L)) {
      arow <- anchors[m_idx[[i]], , drop = FALSE]
      if (nzchar(vdf$chr_clean[[i]]) && nzchar(arow$chrom) && !identical(vdf$chr_clean[[i]], normalize_chrom(arow$chrom))) next
      p_str <- as.character(vdf$pos_int[[i]])
      tested <- tested + 1L
      matched_b <- character(0)
      if (nzchar(arow$pos_grch36) && identical(p_str, arow$pos_grch36)) { counts["GRCh36"] <- counts["GRCh36"] + 1L; matched_b <- c(matched_b, "GRCh36") }
      if (nzchar(arow$pos_grch37) && identical(p_str, arow$pos_grch37)) { counts["GRCh37"] <- counts["GRCh37"] + 1L; matched_b <- c(matched_b, "GRCh37") }
      if (nzchar(arow$pos_grch38) && identical(p_str, arow$pos_grch38)) { counts["GRCh38"] <- counts["GRCh38"] + 1L; matched_b <- c(matched_b, "GRCh38") }
      if (nzchar(arow$pos_t2t_chm13) && identical(p_str, arow$pos_t2t_chm13)) { counts["T2T-CHM13"] <- counts["T2T-CHM13"] + 1L; matched_b <- c(matched_b, "T2T-CHM13") }
      if (length(matched_b) > 0L) {
        ev_parts <- c(ev_parts, sprintf("%s=chr%s:%s (%s)", vdf$rsid_clean[[i]], vdf$chr_clean[[i]], p_str, paste(matched_b, collapse = "/")))
      }
    }
  }

  # Pass 2: Online fallback via NCBI Variation Services SPDI API if < 2 matches from offline table
  if (max(counts) < 2L) {
    rem_df <- if (!is.null(anchors)) vdf[!(vdf$rsid_clean %in% tolower(anchors$rsid)), , drop = FALSE] else vdf
    probe_df <- head(rem_df, 5L)
    for (i in seq_len(nrow(probe_df))) {
      rs <- probe_df$rsid_clean[[i]]
      rs_num <- sub("^rs", "", rs)
      u <- paste0("https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/", rs_num)
      r_v <- perform_cached_request(
        session, "Tier 3", "NCBI_Variation_RefSNP", u,
        cache_key = paste0("refsnp_", rs_num, ".json"),
        accept = "application/json",
        timeout_sec = 15
      )
      if (!r_v$ok) next
      jv <- tryCatch(jsonlite::fromJSON(r_v$path, simplifyVector = FALSE), error = function(e) NULL)
      if (is.null(jv)) next
      tested <- tested + 1L
      p_obs <- probe_df$pos_int[[i]]
      pls <- jv$primary_snapshot_data$placements_with_allele %||% list()
      matched_b <- character(0)
      for (pl in pls) {
        seq_id <- pl$seq_id %||% ""
        if (!grepl("^NC_(0000|0609)", seq_id)) next
        asms <- pl$placement_annot$seq_id_traits_by_assembly %||% list()
        if (length(asms) == 0L || length(pl$alleles) == 0L) next
        asm_nm <- asms[[1L]]$assembly_name %||% ""
        spdi_pos <- pl$alleles[[1L]]$allele$spdi$position
        if (is.null(spdi_pos)) next
        p_ref <- as.integer(spdi_pos) + 1L
        if (identical(p_obs, p_ref)) {
          if (grepl("^NCBI36", asm_nm)) { counts["GRCh36"] <- counts["GRCh36"] + 1L; matched_b <- c(matched_b, "GRCh36") }
          if (grepl("^GRCh37", asm_nm)) { counts["GRCh37"] <- counts["GRCh37"] + 1L; matched_b <- c(matched_b, "GRCh37") }
          if (grepl("^GRCh38", asm_nm)) { counts["GRCh38"] <- counts["GRCh38"] + 1L; matched_b <- c(matched_b, "GRCh38") }
          if (grepl("T2T-CHM13", asm_nm)) { counts["T2T-CHM13"] <- counts["T2T-CHM13"] + 1L; matched_b <- c(matched_b, "T2T-CHM13") }
        }
      }
      matched_b <- unique(matched_b)
      if (length(matched_b) > 0L) {
        ev_parts <- c(ev_parts, sprintf("%s=chr%s:%d (%s)", rs, probe_df$chr_clean[[i]], p_obs, paste(matched_b, collapse = "/")))
      }
      if (max(counts) >= 3L) break
    }
  }

  best_build <- names(which.max(counts))
  best_n <- as.integer(max(counts))
  nonzero_builds <- names(counts[counts > 0L])

  if (best_n >= 2L && length(nonzero_builds) == 1L) {
    list(
      resolved = TRUE,
      genome_build = best_build,
      matched_count = best_n,
      tested_count = tested,
      evidence_snippet = paste(head(ev_parts, 4L), collapse = "; "),
      location = clean_str(vdf$source_loc[[1L]] %||% "Supplement pQTL variant coordinates"),
      mixed_build_warning = ""
    )
  } else if (best_n >= 1L && length(nonzero_builds) > 1L) {
    list(
      resolved = FALSE,
      genome_build = "Unknown",
      matched_count = best_n,
      tested_count = tested,
      evidence_snippet = paste(head(ev_parts, 4L), collapse = "; "),
      location = clean_str(vdf$source_loc[[1L]] %||% "Supplement pQTL variant coordinates"),
      mixed_build_warning = paste0("Conflicting coordinate builds detected across variants: ", paste(nonzero_builds, collapse = " vs "))
    )
  } else {
    res_empty
  }
}
