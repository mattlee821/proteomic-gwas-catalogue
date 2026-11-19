# Tier 2: Legal Full-Text Retrieval Cascade, Section-Aware JATS/BioC Parser,
# 12-Value Sample Matrix Classifier ([C-2]), and 4-Rule Sample-Size Disambiguator ([C-5]).

ALLOWED_SAMPLE_MATRICES <- c(
  "EDTA plasma", "Citrate plasma", "Heparin plasma", "Plasma (unspecified)",
  "Serum", "Whole blood / DBS", "CSF", "Brain tissue",
  "Urine", "Saliva", "Cell line / LCL", "Other tissue"
)

classify_section_bucket <- function(title_txt, sec_type = "") {
  ttl <- tolower(clean_str(title_txt))
  st <- tolower(clean_str(sec_type))
  if (grepl("data-availability|associated-data", st) || grepl("data availability|code availability", ttl)) {
    return("DATA_AVAILABILITY")
  }
  if (grepl("subsampl|downsampl|power simulation|projection", ttl)) {
    return("SUBSAMPLING_SIMULATION")
  }
  if (grepl("supplementary|extended data|reporting summary", ttl) || grepl("supplementary", st)) {
    return("SUPP_CAPTIONS")
  }
  if (grepl("genomic data|genotyping|imputation|whole-genome sequenc|variant calling", ttl)) {
    return("GENOMIC_IMPUTATION")
  }
  if (grepl("proteomic measurement|protein measurement|plasma samples and protein|olink|somascan|aptamer|data processing and quality control|assay", ttl)) {
    return("PROTEOMIC_QC")
  }
  if (grepl("(?<!non-)genetic association|genome-wide association|pqtl mapping|discovery of pqtl|replication of pqtl|non-european pqtl|linear mixed|heritability|fine-mapping|fine mapping|definition and refinement", ttl, perl = TRUE)) {
    return("GWAS_ASSOCIATION")
  }
  if (grepl("participant|sample selection|study population|cohort|overview of", ttl)) {
    return("COHORT_SAMPLES")
  }
  if (grepl("^abstract$|^main$|^summary$", ttl)) {
    return("ABSTRACT_OVERVIEW")
  }
  "OTHER"
}

parse_jats_sections <- function(xml_path) {
  doc <- xml2::read_xml(xml_path)
  buckets <- list()
  add_b <- function(b, heading, txt) {
    txt_c <- clean_str(txt)
    if (!nzchar(txt_c)) return()
    buckets[[b]] <<- c(buckets[[b]] %||% character(0), paste0("[", clean_str(heading), "] ", txt_c))
  }

  abs_nodes <- xml2::xml_find_all(doc, "//front//abstract")
  for (an in abs_nodes) add_b("ABSTRACT_OVERVIEW", "Abstract", xml2::xml_text(an))

  # Extract leaf or meaningful sections without double-counting parent + child text
  sec_nodes <- xml2::xml_find_all(doc, "//body//sec | //back//sec | //back//notes")
  for (sn in sec_nodes) {
    ttl_nd <- xml2::xml_find_first(sn, "title")
    ttl <- if (!inherits(ttl_nd, "xml_missing")) xml2::xml_text(ttl_nd) else ""
    st <- xml2::xml_attr(sn, "sec-type") %||% ""
    # Take direct <p> children of this <sec> so nested <sec> paragraphs are classified under their own more specific title!
    p_nodes <- xml2::xml_find_all(sn, "p")
    if (length(p_nodes) > 0L) {
      p_txt <- paste(xml2::xml_text(p_nodes), collapse = " ")
      b <- classify_section_bucket(ttl, st)
      add_b(b, ttl, p_txt)
    }
  }

  # Extract <supplementary-material> and <media> inventory
  supp_items <- list()
  sm_nodes <- xml2::xml_find_all(doc, "//supplementary-material")
  for (sm in sm_nodes) {
    sm_id <- clean_str(xml2::xml_attr(sm, "id") %||% "")
    lbl <- clean_str(xml2::xml_text(xml2::xml_find_first(sm, ".//label")))
    cap <- clean_str(xml2::xml_text(xml2::xml_find_first(sm, ".//caption")))
    med <- xml2::xml_find_first(sm, ".//media")
    href <- ""
    mime <- ""
    if (!inherits(med, "xml_missing")) {
      href <- clean_str(xml2::xml_attr(med, "href") %||% xml2::xml_attr(med, "xlink:href") %||% "")
      if (!nzchar(lbl)) lbl <- clean_str(xml2::xml_text(xml2::xml_find_first(med, "label")))
      if (!nzchar(cap)) cap <- clean_str(xml2::xml_text(xml2::xml_find_first(med, "caption")))
      mime <- paste0(clean_str(xml2::xml_attr(med, "mimetype") %||% ""), "/", clean_str(xml2::xml_attr(med, "mime-subtype") %||% ""))
    } else {
      href <- clean_str(xml2::xml_attr(sm, "href") %||% xml2::xml_attr(sm, "xlink:href") %||% "")
    }
    if (nzchar(href)) {
      supp_items[[href]] <- list(id = sm_id, filename = basename(href), href = href, mime = mime, label = lbl, caption = cap)
    }
  }

  list(buckets = buckets, jats_supplements = unname(supp_items), n_tables = length(xml2::xml_find_all(doc, "//table-wrap")))
}

parse_bioc_sections <- function(json_path) {
  j <- jsonlite::fromJSON(json_path, simplifyVector = FALSE)
  docs <- if (!is.null(j$documents)) j$documents else if (is.list(j) && length(j) > 0L && !is.null(j[[1L]]$documents)) j[[1L]]$documents else list()
  buckets <- list()
  if (length(docs) > 0L) {
    for (psg in (docs[[1L]]$passages %||% list())) {
      st <- clean_str(psg$infons$section_type %||% psg$infons$type %||% "")
      txt <- clean_str(psg$text %||% "")
      if (!nzchar(txt)) next
      b <- classify_section_bucket(st, st)
      buckets[[b]] <- c(buckets[[b]] %||% character(0), paste0("[", st, "] ", txt))
    }
  }
  list(buckets = buckets, jats_supplements = list(), n_tables = 0L)
}

detect_cohorts_from_text <- function(txt, title_txt = "") {
  comb <- paste(title_txt, txt)
  cohort_map <- list(
    "UK Biobank" = "\\bUK Biobank\\b|\\bUKB-PPP\\b|\\bUKB participants\\b",
    "deCODE Genetics" = "\\bdeCODE\\b|\\bIcelanders\\b|\\bIcelandic Cancer Project\\b",
    "INTERVAL" = "\\bINTERVAL study\\b|\\bINTERVAL cohort\\b",
    "Fenland" = "\\bFenland\\b",
    "KORA" = "\\bKORA\\b",
    "AGES-Reykjavik" = "\\bAGES-Reykjavik\\b",
    "ARIC" = "\\bAtherosclerosis Risk in Communities\\b|\\bARIC\\b",
    "Framingham Heart Study" = "\\bFramingham Heart Study\\b",
    "SCALLOP" = "\\bSCALLOP consortium\\b",
    "Generation Scotland" = "\\bGeneration Scotland\\b",
    "China Kadoorie Biobank" = "\\bChina Kadoorie Biobank\\b",
    "BioBank Japan" = "\\bBioBank Japan\\b",
    "ROSMAP" = "\\bROSMAP\\b",
    "Knight ADRC" = "\\bKnight ADRC\\b",
    "DIAN" = "\\bDIAN\\b",
    "GCKD" = "\\bGCKD\\b",
    "QMDiab" = "\\bQMDiab\\b",
    "CHRIS" = "\\bCHRIS study\\b",
    "IMPROVE" = "\\bIMPROVE\\b",
    "ORCADES" = "\\bORCADES\\b"
  )
  # Only match primary study cohorts (avoid matching replication comparator mentions like "reported by the INTERVAL study")
  primary_txt <- gsub("reported (in|by) the [A-Za-z0-9 -]+ (study|consortium)|compared (with|to) [A-Za-z0-9 -]+", " ", comb, ignore.case = TRUE)
  hits <- character(0)
  for (nm in names(cohort_map)) {
    if (grepl(cohort_map[[nm]], primary_txt, ignore.case = FALSE, perl = TRUE)) {
      hits <- c(hits, nm)
    }
  }
  unique(hits)
}

detect_sample_matrix_12 <- function(txt, loc_label) {
  if (!nzchar(txt)) return(list(value = "", ambiguous = TRUE, snippet = "", location = ""))
  grab_snip <- function(pat) {
    m <- regmatches(txt, regexpr(paste0(".{0,55}(", pat, ").{0,65}"), txt, ignore.case = TRUE, perl = TRUE))
    if (length(m) == 1L && nzchar(m)) clean_str(m) else ""
  }

  rules <- list(
    list(val = "EDTA plasma", amb = FALSE, pat = "\\bEDTA[\\s-]+(blood\\s+|tubes?\\s+)?plasma\\b|\\bblood\\s+in\\s+EDTA\\s+tubes\\b|\\bK2?-?EDTA\\b"),
    list(val = "Citrate plasma", amb = FALSE, pat = "\\bcitrate[d]?[\\s-]+plasma\\b|\\bsodium\\s+citrate\\s+tubes\\b"),
    list(val = "Heparin plasma", amb = FALSE, pat = "\\b(lithium\\s+)?heparin[\\s-]+plasma\\b"),
    list(val = "Whole blood / DBS", amb = FALSE, pat = "\\bdried\\s+blood\\s+spots?\\b|\\bDBS\\s+proteom|\\bwhole\\s+blood\\s+proteom"),
    list(val = "CSF", amb = FALSE, pat = "\\bcerebrospinal\\s+fluid\\b|\\bCSF\\s+(proteom|protein|samples|levels)"),
    list(val = "Brain tissue", amb = FALSE, pat = "\\bbrain\\s+(tissue|homogenate|proteom)|\\bdorsolateral\\s+prefrontal\\s+cortex\\b|\\bDLPFC\\b"),
    list(val = "Urine", amb = FALSE, pat = "\\burin(e|ary)\\s+(proteom|protein|samples|levels)\\b"),
    list(val = "Saliva", amb = FALSE, pat = "\\bsaliv(a|ary)\\s+(proteom|protein|samples|levels)\\b"),
    list(val = "Cell line / LCL", amb = FALSE, pat = "\\blymphoblastoid\\s+cell\\s+lines?\\b|\\bLCL\\s+proteom"),
    list(val = "Serum", amb = FALSE, pat = "\\b(blood\\s+)?serum\\s+(proteom|protein|samples|levels|profiles)\\b"),
    list(val = "Plasma (unspecified)", amb = TRUE, pat = "\\b(blood\\s+)?plasma\\s+(proteom\\w*|protein\\s+levels|proteins?|samples|profiles)\\b")
  )
  for (rl in rules) {
    if (grepl(rl$pat, txt, ignore.case = TRUE, perl = TRUE)) {
      return(list(value = rl$val, ambiguous = rl$amb, snippet = grab_snip(rl$pat), location = loc_label))
    }
  }
  list(value = "", ambiguous = TRUE, snippet = "", location = "")
}

detect_platform_details <- function(txt, loc_label) {
  if (!nzchar(txt)) return(list(category = "", version = "", snippet = "", location = ""))
  grab_snip <- function(pat) {
    m <- regmatches(txt, regexpr(paste0(".{0,45}(", pat, ").{0,60}"), txt, ignore.case = TRUE, perl = TRUE))
    if (length(m) == 1L && nzchar(m)) clean_str(m) else ""
  }

  ver_rules <- list(
    list(cat = "Olink", ver = "Olink Explore HT", pat = "Olink\\s+Explore\\s+HT"),
    list(cat = "Olink", ver = "Olink Explore 3072", pat = "Olink\\s+Explore\\s+3072"),
    list(cat = "Olink", ver = "Olink Explore 1536", pat = "Olink\\s+Explore\\s+1536"),
    list(cat = "Olink", ver = "Olink Target 96", pat = "Olink\\s+(Target\\s+96|CVD-I|Inflammation)"),
    list(cat = "SomaScan", ver = "SomaScan 11k", pat = "SomaScan\\s*(11k|v5)"),
    list(cat = "SomaScan", ver = "SomaScan v4.1 (7k)", pat = "SomaScan\\s*(7k|v4\\.1|version\\s*4\\.1)"),
    list(cat = "SomaScan", ver = "SomaScan v4 (4.9k)", pat = "SomaScan\\s*(v4|version\\s*4)\\b|4,907\\s+aptamers|5,284\\s+aptamers"),
    list(cat = "SomaScan", ver = "SomaScan v3", pat = "SomaScan\\s*(v3|1\\.3k|1\\.1k)"),
    list(cat = "NULISA", ver = "Alamar NULISAseq", pat = "NULISAseq|NULISA"),
    list(cat = "Mass Spectrometry", ver = "Seer Proteograph DIA-MS", pat = "Proteograph"),
    list(cat = "Mass Spectrometry", ver = "SWATH / DIA-MS", pat = "SWATH|DIA-MS|data-independent\\s+acquisition"),
    list(cat = "Mass Spectrometry", ver = "TMT-MS", pat = "TMT-MS|tandem\\s+mass\\s+tag"),
    list(cat = "Luminex", ver = "Luminex xMAP", pat = "Luminex"),
    list(cat = "Olink", ver = "Olink (version unspecified)", pat = "\\bOlink\\b|proximity\\s+extension\\s+assay"),
    list(cat = "SomaScan", ver = "SomaScan (version unspecified)", pat = "\\bSomaScan\\b|\\bSomaLogic\\b|\\baptamers\\b"),
    list(cat = "Mass Spectrometry", ver = "Mass Spectrometry", pat = "mass\\s+spectrometry")
  )
  for (rl in ver_rules) {
    if (grepl(rl$pat, txt, ignore.case = TRUE, perl = TRUE)) {
      return(list(category = rl$cat, version = rl$ver, snippet = grab_snip(rl$pat), location = loc_label))
    }
  }
  list(category = "", version = "", snippet = "", location = "")
}

parse_int_clean <- function(s) {
  suppressWarnings(as.integer(gsub("[^0-9]", "", s)))
}

# [C-5] 4-Rule Sample-Size & Multi-Ancestry Strata Disambiguator
extract_strata_and_n <- function(abstract_title_txt, buckets, study_id, has_fulltext) {
  # Filter out parent-biobank recruitment sentences ("cohort of approximately 500,000 participants"),
  # subsampling simulations, and repeat-visit counts.
  filter_valid_sentences <- function(txt_vec) {
    if (length(txt_vec) == 0L) return(character(0))
    sents <- unlist(strsplit(paste(txt_vec, collapse = " "), "(?<=\\.)\\s+", perl = TRUE))
    keep <- !grepl(
      paste0(
        "is a population-based cohort of (approximately |~|over )?[0-9,]+|",
        "500,000 participants|",
        "subsampl|downsampl|visit 2\\s*\\(|visit 3\\s*\\("
      ),
      sents, ignore.case = TRUE, perl = TRUE
    )
    sents[keep]
  }

  ov_sents <- filter_valid_sentences(c(abstract_title_txt, buckets$ABSTRACT_OVERVIEW, buckets$COHORT_SAMPLES))
  gw_sents <- filter_valid_sentences(buckets$GWAS_ASSOCIATION)
  ov_txt <- paste(ov_sents, collapse = " ")
  gw_txt <- paste(gw_sents, collapse = " ")

  # 1. Track total proteomic-profiled N (e.g., 54,219 UKB participants)
  n_profiled <- NA_integer_
  prof_snip <- ""
  m_prof <- regmatches(ov_txt, regexpr(
    "(proteomic profiling[^.]{0,60}?from\\s+|proteomic profiles of\\s+|in\\s+)([1-9][0-9,]{1,6})\\s+(UK Biobank|UKB|Icelanders|participants|individuals|subjects)",
    ov_txt, ignore.case = TRUE, perl = TRUE
  ))
  if (length(m_prof) == 1L && nzchar(m_prof)) {
    num_m <- regmatches(m_prof, regexpr("[1-9][0-9,]{1,6}", m_prof))
    n_profiled <- parse_int_clean(num_m)
    prof_snip <- clean_str(m_prof)
  }

  strata_list <- list()
  add_st <- function(suffix, anc_grp, anc_det, n_val, loc, snip, n_note = "") {
    if (is.na(n_val) || n_val < 1L) return()
    strata_list[[suffix]] <<- list(
      stratum_suffix = suffix,
      ancestry_group = anc_grp,
      ancestry_detail = anc_det,
      sex = "Combined",
      n_samples = n_val,
      location = loc,
      snippet = clean_str(snip),
      n_note = n_note
    )
  }

  if (has_fulltext && nzchar(gw_txt)) {
    # Combined / Full-cohort GWAS N (e.g. "pQTL mapping for the full UKB-PPP cohort (n = 52,363)")
    m_comb <- regmatches(gw_txt, regexpr(
      "pQTL mapping for the full[^.]{0,50}?\\(n\\s*=\\s*([1-9][0-9,]+)\\)",
      gw_txt, ignore.case = TRUE, perl = TRUE
    ))
    if (length(m_comb) == 1L && nzchar(m_comb)) {
      nc <- parse_int_clean(regmatches(m_comb, regexpr("(?<=n\\s=\\s|n=)[1-9][0-9,]+", m_comb, perl = TRUE)))
      if (is.na(nc)) nc <- parse_int_clean(regmatches(m_comb, regexpr("[1-9][0-9,]{3,}", m_comb)))
      note_qc <- ""
      if (!is.na(n_profiled) && nc < n_profiled && (nc / n_profiled) >= 0.70) {
        note_qc <- sprintf("Post-QC GWAS N=%d (genotype/baseline QC attrition from %d proteomic-profiled participants)", nc, n_profiled)
      }
      add_st("ALL", "Multi-ancestry", "Combined multi-ancestry cohort", nc, "JATS GWAS_ASSOCIATION", m_comb, note_qc)
    }

    # Discovery European GWAS N (e.g. "Discovery pQTL analyses were performed in participants of European ancestry... (n = 34,557)")
    m_eur <- regmatches(gw_txt, regexpr(
      "Discovery pQTL analyses[^.]{0,110}?European ancestry[^.]{0,80}?\\(n\\s*=\\s*([1-9][0-9,]+)\\)",
      gw_txt, ignore.case = TRUE, perl = TRUE
    ))
    if (length(m_eur) == 1L && nzchar(m_eur)) {
      ne <- parse_int_clean(regmatches(m_eur, regexpr("[1-9][0-9,]{3,}", m_eur)))
      add_st("EUR", "EUR", "European (discovery)", ne, "JATS GWAS_ASSOCIATION", m_eur)
    }

    # Explicit non-European ancestry-specific clauses: African (AFR, n = 931), CSA, MID, EAS, AMR
    anc_specs <- list(
      list(suf = "AFR", grp = "AFR", det = "African", pat = "African\\s*\\(\\s*AFR\\s*,\\s*n\\s*=\\s*([1-9][0-9,]*)\\s*\\)"),
      list(suf = "CSA", grp = "CSA", det = "Central/South Asian", pat = "Central/South Asian\\s*\\(\\s*CSA\\s*,\\s*n\\s*=\\s*([1-9][0-9,]*)\\s*\\)"),
      list(suf = "MID", grp = "MID", det = "Middle Eastern", pat = "Middle Eastern\\s*\\(\\s*MID\\s*,\\s*n\\s*=\\s*([1-9][0-9,]*)\\s*\\)"),
      list(suf = "EAS", grp = "EAS", det = "East Asian", pat = "East Asian\\s*\\(\\s*EAS\\s*,\\s*n\\s*=\\s*([1-9][0-9,]*)\\s*\\)"),
      list(suf = "AMR", grp = "AMR", det = "Admixed American", pat = "admixed American\\s*\\(\\s*AMR\\s*,\\s*n\\s*=\\s*([1-9][0-9,]*)\\s*\\)")
    )
    for (aspec in anc_specs) {
      mm <- regmatches(gw_txt, regexpr(aspec$pat, gw_txt, ignore.case = TRUE, perl = TRUE))
      if (length(mm) == 1L && nzchar(mm)) {
        nv <- parse_int_clean(mm)
        add_st(aspec$suf, aspec$grp, aspec$det, nv, "JATS GWAS_ASSOCIATION [Non-European pQTL mapping]", mm)
      }
    }
  }

  # Fallback to Abstract / Title headline GWAS N when full text is unavailable or single-stratum
  if (length(strata_list) == 0L) {
    # Check "in 35,559 Icelanders" / "in 54,306 UK Biobank participants" / "of 54,219 UK Biobank participants"
    m_abs <- regmatches(abstract_title_txt, regexpr(
      "(in|of)\\s+([1-9][0-9,]{1,6})\\s+(Icelanders|Europeans|Chinese\\s+adults|Scottish\\s+individuals|UK\\s+Biobank\\s+participants|participants|individuals)",
      abstract_title_txt, ignore.case = TRUE, perl = TRUE
    ))
    if (length(m_abs) == 1L && nzchar(m_abs)) {
      na_val <- parse_int_clean(regmatches(m_abs, regexpr("[1-9][0-9,]{1,6}", m_abs)))
      if (grepl("Icelanders", m_abs, ignore.case = TRUE)) {
        add_st("EUR", "EUR", "European (Icelandic)", na_val, "Tier 1 Abstract/Title", m_abs)
      } else if (grepl("Chinese", m_abs, ignore.case = TRUE)) {
        add_st("EAS", "EAS", "East Asian (Chinese)", na_val, "Tier 1 Abstract/Title", m_abs)
      } else if (grepl("Europeans|Scottish|British", m_abs, ignore.case = TRUE)) {
        add_st("EUR", "EUR", "European", na_val, "Tier 1 Abstract/Title", m_abs)
      } else {
        # Check if abstract mentions non-European ancestry-specific mapping
        is_multi <- grepl("non-European|multi-ancestry|ancestry-specific", abstract_title_txt, ignore.case = TRUE)
        grp <- if (is_multi) "Multi-ancestry" else "EUR"
        det <- if (is_multi) "Multi-ancestry (exact per-stratum N requires full text)" else "Ancestry inferred from cohort (verify in full text)"
        add_st(if (is_multi) "ALL" else "EUR", grp, det, na_val, "Tier 1 Abstract/Title", m_abs,
               if (!has_fulltext) "Headline N from abstract/title; verify post-QC per-stratum GWAS N" else "")
      }
    }
  }

  list(strata = strata_list, n_profiled = n_profiled, profiled_snippet = prof_snip)
}

extract_methods_fields <- function(abstract_title_txt, buckets, has_fulltext) {
  prot_txt <- paste(c(abstract_title_txt, buckets$ABSTRACT_OVERVIEW, buckets$PROTEOMIC_QC, buckets$GWAS_ASSOCIATION), collapse = " ")
  geno_txt <- paste(c(buckets$GENOMIC_IMPUTATION, buckets$GWAS_ASSOCIATION), collapse = " ")
  gwas_txt <- paste(c(buckets$GWAS_ASSOCIATION, buckets$PROTEOMIC_QC), collapse = " ")

  # Protein & analyte counts
  n_assays <- NA_integer_
  n_unique_proteins <- NA_integer_
  prot_snip <- ""

  # Pattern 1: "measuring 2,941 protein analytes and capturing 2,923 unique proteins"
  m_both <- regmatches(prot_txt, regexpr(
    "measuring\\s+([1-9][0-9,]+)\\s+protein\\s+analytes\\s+and\\s+capturing\\s+([1-9][0-9,]+)\\s+unique\\s+proteins",
    prot_txt, ignore.case = TRUE, perl = TRUE
  ))
  if (length(m_both) == 1L && nzchar(m_both)) {
    nums <- regmatches(m_both, gregexpr("[1-9][0-9,]+", m_both))[[1L]]
    n_assays <- parse_int_clean(nums[[1L]])
    n_unique_proteins <- parse_int_clean(nums[[2L]])
    prot_snip <- clean_str(m_both)
  } else {
    m_apt <- regmatches(prot_txt, regexpr("([1-9][0-9,]+)\\s+aptamers", prot_txt, ignore.case = TRUE, perl = TRUE))
    if (length(m_apt) == 1L && nzchar(m_apt)) {
      n_assays <- parse_int_clean(m_apt)
      prot_snip <- clean_str(m_apt)
    }
    m_up <- regmatches(prot_txt, regexpr("([1-9][0-9,]+)\\s+(unique\\s+)?proteins", prot_txt, ignore.case = TRUE, perl = TRUE))
    if (length(m_up) == 1L && nzchar(m_up)) {
      up_v <- parse_int_clean(m_up)
      if (is.na(n_assays)) n_assays <- up_v
      n_unique_proteins <- up_v
      if (!nzchar(prot_snip)) prot_snip <- clean_str(m_up)
    }
  }

  # Genome build from prose
  genome_build <- ""
  build_snip <- ""
  if (grepl("\\bGRCh38\\b|\\bhg38\\b|\\bbuild\\s*38\\b", geno_txt, ignore.case = TRUE)) {
    genome_build <- "GRCh38"
    build_snip <- clean_str(regmatches(geno_txt, regexpr(".{0,45}(GRCh38|hg38|build\\s*38).{0,45}", geno_txt, ignore.case = TRUE)))
  } else if (grepl("\\bGRCh37\\b|\\bhg19\\b|\\bbuild\\s*37\\b", geno_txt, ignore.case = TRUE)) {
    genome_build <- "GRCh37"
    build_snip <- clean_str(regmatches(geno_txt, regexpr(".{0,45}(GRCh37|hg19|build\\s*37).{0,45}", geno_txt, ignore.case = TRUE)))
  }

  # Imputation panel
  imp_panel <- ""
  if (grepl("TOPMed", geno_txt, ignore.case = TRUE) && grepl("HRC|Haplotype Reference Consortium", geno_txt, ignore.case = TRUE)) {
    imp_panel <- "HRC + UK10K + TOPMed"
  } else if (grepl("TOPMed", geno_txt, ignore.case = TRUE)) {
    imp_panel <- "TOPMed"
  } else if (grepl("Haplotype Reference Consortium|\\bHRC\\b", geno_txt, ignore.case = TRUE)) {
    imp_panel <- "HRC"
  } else if (grepl("1000 Genomes", geno_txt, ignore.case = TRUE)) {
    imp_panel <- "1000 Genomes"
  }

  # GWAS software
  gwas_sw <- ""
  sw_snip <- ""
  for (sw in c("REGENIE", "BOLT-LMM", "fastGWA", "PLINK2", "PLINK", "SNPTEST", "SAIGE", "METAL", "MatrixEQTL")) {
    if (grepl(paste0("\\b", sw, "\\b"), gwas_txt, ignore.case = TRUE, perl = TRUE)) {
      m_sw <- regmatches(gwas_txt, regexpr(paste0(".{0,30}\\b", sw, "[^.]{0,45}"), gwas_txt, ignore.case = TRUE, perl = TRUE))
      gwas_sw <- sw
      sw_snip <- clean_str(m_sw)
      break
    }
  }

  # Transformation & Effect Scale
  transformation <- ""
  effect_scale <- ""
  tr_snip <- ""
  if (grepl("inverse[- ]rank normal|rank[- ]based inverse normal|rank[- ]inverse normal|\\bIRNT\\b", gwas_txt, ignore.case = TRUE, perl = TRUE)) {
    transformation <- if (grepl("\\bNPX\\b", gwas_txt)) "IRNT of NPX" else "Rank-based inverse normal transformation (IRNT)"
    effect_scale <- "SD_IRNT"
    tr_snip <- clean_str(regmatches(gwas_txt, regexpr(".{0,45}(inverse[- ]rank normal|rank[- ]inverse normal|IRNT).{0,55}", gwas_txt, ignore.case = TRUE, perl = TRUE)))
  } else if (grepl("Z-standardized|standardized residuals", gwas_txt, ignore.case = TRUE)) {
    transformation <- "Z-standardized"
    effect_scale <- "SD_Z"
  }

  # Covariates
  covariates <- ""
  cov_snip <- ""
  m_cov <- regmatches(gwas_txt, regexpr(
    "(covariates:\\s*|adjusting for\\s+|adjusted[^.]{0,35}?for\\s+)([^.]{10,220})",
    gwas_txt, ignore.case = TRUE, perl = TRUE
  ))
  if (length(m_cov) == 1L && nzchar(m_cov)) {
    cov_snip <- clean_str(m_cov)
    covariates <- sub("^(covariates:\\s*|adjusting for\\s+|adjusted[^.]{0,35}?for\\s+)", "", cov_snip, ignore.case = TRUE, perl = TRUE)
  }

  # Coverage & cis_window_kb
  coverage <- ""
  cis_kb <- NA_integer_
  if (has_fulltext || grepl("genome-wide association|\\bpQTL\\b", prot_txt, ignore.case = TRUE)) {
    if (grepl("cis-only|restricted to cis|within\\s+\\+/-\\s*(500\\s*kb|1\\s*Mb)\\s+only", gwas_txt, ignore.case = TRUE)) {
      coverage <- "cis_only"
      cis_kb <- if (grepl("500\\s*kb", gwas_txt, ignore.case = TRUE)) 500L else 1000L
    } else if (grepl("genome-wide|trans\\s+associations|cis\\s+and\\s+trans|trans-pQTL", prot_txt, ignore.case = TRUE)) {
      coverage <- "genome_wide"
    }
  }

  list(
    n_assays = n_assays,
    n_unique_proteins = n_unique_proteins,
    prot_snippet = prot_snip,
    genome_build = genome_build,
    build_snippet = build_snip,
    imputation_panel = imp_panel,
    gwas_software = gwas_sw,
    sw_snippet = sw_snip,
    transformation = transformation,
    effect_scale = effect_scale,
    tr_snippet = tr_snip,
    covariates = covariates,
    cov_snippet = cov_snip,
    coverage = coverage,
    cis_window_kb = cis_kb
  )
}

run_tier2_fulltext <- function(t1, session, local_pdf = "", local_xml = "") {
  leg_pref <- legacy_prefix_for_study(t1$pmid, t1$doi)
  slug_base <- t1$study_id

  ft_source <- "NONE"
  verified_full_text <- FALSE
  parsed <- list(buckets = list(), jats_supplements = list(), n_tables = 0L)
  public_da_text <- ""
  public_moesm_urls <- character(0)
  ft_gate_note <- ""

  # 0. Local curator XML/PDF override (--local-xml / --local-pdf)
  if (nzchar(local_xml) && file.exists(local_xml)) {
    parsed <- parse_jats_sections(local_xml)
    ft_source <- "LOCAL_CURATOR_XML"
    verified_full_text <- TRUE
  } else if (nzchar(local_pdf) && file.exists(local_pdf) && requireNamespace("pdftools", quietly = TRUE)) {
    pdf_pages <- pdftools::pdf_text(local_pdf)
    parsed$buckets$GWAS_ASSOCIATION <- paste(pdf_pages, collapse = "\n")
    ft_source <- "LOCAL_CURATOR_PDF"
    verified_full_text <- TRUE
  }

  # 2a. Europe PMC JATS /fullTextXML — ONLY for valid PMC[0-9]+ IDs (NEVER PPR* or EMS* — [C-6])
  if (!verified_full_text && grepl("^PMC[0-9]+$", t1$pmcid) && (identical(t1$in_epmc, "Y") || identical(t1$is_open_access_epmc, "Y"))) {
    ft_url <- sprintf("https://www.ebi.ac.uk/europepmc/webservices/rest/%s/fullTextXML", t1$pmcid)
    leg_ft <- if (nzchar(leg_pref)) paste0(leg_pref, "_epmc_fullText.xml") else character(0)
    r_ft <- perform_cached_request(
      session, "Tier 2", "EuropePMC_fullTextXML", ft_url,
      cache_key = paste0(slug_base, "_epmc_fullText.xml"),
      accept = "application/xml",
      legacy_candidates = leg_ft
    )
    if (r_ft$ok && grepl("<article|<!DOCTYPE\\s+article|<pmc-articleset", r_ft$gates$head_text, ignore.case = TRUE)) {
      parsed <- parse_jats_sections(r_ft$path)
      ft_source <- paste0("EuropePMC_fullTextXML (", t1$pmcid, ")")
      verified_full_text <- TRUE
    }
  }

  # 2b. NCBI PMC BioC JSON fallback (with strict [Error] rejection)
  if (!verified_full_text && (grepl("^PMC[0-9]+$", t1$pmcid) || identical(t1$is_open_access_epmc, "Y")) && !t1$is_preprint) {
    bioc_id <- if (nzchar(t1$pmcid)) t1$pmcid else t1$pmid
    if (nzchar(bioc_id)) {
      bioc_url <- sprintf("https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/%s/unicode", bioc_id)
      leg_bioc <- if (nzchar(leg_pref)) c(paste0(leg_pref, "_bioc_pmcid.json"), paste0(leg_pref, "_bioc_pmid.json")) else character(0)
      r_bioc <- perform_cached_request(
        session, "Tier 2", "NCBI_PMC_BioC", bioc_url,
        cache_key = paste0(slug_base, "_bioc.json"),
        accept = "application/json",
        legacy_candidates = leg_bioc
      )
      if (r_bioc$ok) {
        parsed <- parse_bioc_sections(r_bioc$path)
        if (length(parsed$buckets) > 0L) {
          ft_source <- paste0("NCBI_PMC_BioC (", bioc_id, ")")
          verified_full_text <- TRUE
        }
      }
    }
  }

  # 2c. Preprint JATS XML single polite attempt ([C-6]: NEVER retry or bypass 403/429 Cloudflare challenge!)
  if (!verified_full_text && t1$is_preprint && nzchar(t1$biorxiv_jats_url)) {
    r_bx_xml <- perform_cached_request(
      session, "Tier 2", "bioRxiv_jatsxml", t1$biorxiv_jats_url,
      cache_key = paste0(slug_base, "_biorxiv_source.xml"),
      accept = "application/xml",
      timeout_sec = 12,
      allow_transient_retry = FALSE
    )
    if (r_bx_xml$ok && grepl("<article", r_bx_xml$gates$head_text, ignore.case = TRUE)) {
      parsed <- parse_jats_sections(r_bx_xml$path)
      ft_source <- "bioRxiv_JATS_XML"
      verified_full_text <- TRUE
    } else {
      ft_gate_note <- sprintf(
        "Preprint JATS XML unavailable via automated HTTP (status=%s%s); degraded cleanly to Tier 1 + datalinks without bypass",
        as.character(r_bx_xml$status %||% "ERR"),
        if (isTRUE(r_bx_xml$gates$antibot)) ", Cloudflare bot gate" else ""
      )
    }
  }

  # 2d. Closed-license opportunistic public metadata check ([M-8]: ONLY #data-availability & public MOESM links, NEVER closed prose!)
  if (!verified_full_text && grepl("^10\\.1038/", t1$doi)) {
    ft_gate_note <- sprintf(
      "Closed license (EuropePMC isOpenAccess=%s, Unpaywall oa_status=%s, license=%s); full-text Methods prose skipped in strict_oa_only mode",
      t1$is_open_access_epmc, t1$oa_status, ifelse(nzchar(t1$license_url), t1$license_url, "unspecified")
    )
    nat_url <- paste0("https://www.nature.com/articles/", sub("^10\\.1038/", "", t1$doi))
    leg_html <- if (nzchar(leg_pref)) paste0(leg_pref, "_nature_landing.html") else character(0)
    r_nat <- perform_cached_request(
      session, "Tier 2", "SpringerNature_public_meta", nat_url,
      cache_key = paste0(slug_base, "_nature_landing.html"),
      accept = "text/html",
      timeout_sec = 12,
      allow_transient_retry = FALSE,
      legacy_candidates = leg_html
    )
    if (r_nat$ok && file.exists(r_nat$path)) {
      html_txt <- paste(readLines(r_nat$path, warn = FALSE), collapse = "\n")
      da_m <- regmatches(html_txt, regexpr("id=\"data-availability\"[\\s\\S]{0,3000}?</section>", html_txt, perl = TRUE))
      if (length(da_m) == 1L && nzchar(da_m)) {
        public_da_text <- clean_str(gsub("<[^>]+>", " ", da_m))
      }
      moesm_m <- unique(regmatches(html_txt, gregexpr(
        "https://static-content\\.springer\\.com/esm/[^\"'<>[:space:]]+",
        html_txt, perl = TRUE
      ))[[1L]])
      public_moesm_urls <- moesm_m
    }
  } else if (!verified_full_text && !nzchar(ft_gate_note)) {
    ft_gate_note <- sprintf("No open-access full-text XML available (isOpenAccess=%s, oa_status=%s)", t1$is_open_access_epmc, t1$oa_status)
  }

  # Run section-aware extractors
  abs_title <- paste(t1$title, t1$abstract_text)
  all_ft_txt <- if (verified_full_text) paste(unlist(parsed$buckets), collapse = " ") else ""

  cohorts <- detect_cohorts_from_text(paste(abs_title, paste(c(parsed$buckets$ABSTRACT_OVERVIEW, parsed$buckets$COHORT_SAMPLES), collapse = " ")), t1$title)

  # Sample matrix: check Tier 2 COHORT_SAMPLES + PROTEOMIC_QC first if verified_full_text, else Tier 1 Abstract
  mat_res <- if (verified_full_text) {
    m2 <- detect_sample_matrix_12(paste(c(parsed$buckets$COHORT_SAMPLES, parsed$buckets$PROTEOMIC_QC, parsed$buckets$ABSTRACT_OVERVIEW), collapse = " "), "Tier 2 Full Text")
    if (nzchar(m2$value)) m2 else detect_sample_matrix_12(abs_title, "Tier 1 Abstract")
  } else {
    detect_sample_matrix_12(abs_title, "Tier 1 Abstract")
  }

  plat_res <- if (verified_full_text) {
    p2 <- detect_platform_details(paste(c(parsed$buckets$ABSTRACT_OVERVIEW, parsed$buckets$PROTEOMIC_QC), collapse = " "), "Tier 2 Full Text")
    if (nzchar(p2$category)) p2 else detect_platform_details(abs_title, "Tier 1 Abstract/Title")
  } else {
    detect_platform_details(abs_title, "Tier 1 Abstract/Title")
  }

  strata_res <- extract_strata_and_n(abs_title, parsed$buckets, t1$study_id, verified_full_text)
  meth_res <- extract_methods_fields(abs_title, parsed$buckets, verified_full_text)

  da_combined <- paste(c(parsed$buckets$DATA_AVAILABILITY %||% character(0), public_da_text), collapse = " ")

  list(
    verified_full_text = verified_full_text,
    ft_source = ft_source,
    ft_gate_note = ft_gate_note,
    buckets = parsed$buckets,
    jats_supplements = parsed$jats_supplements,
    public_moesm_urls = public_moesm_urls,
    data_availability_text = da_combined,
    cohorts = cohorts,
    sample_matrix = mat_res,
    platform = plat_res,
    strata_info = strata_res,
    methods = meth_res
  )
}
