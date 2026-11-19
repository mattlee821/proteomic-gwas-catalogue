# Canonical lexicons, controlled vocabularies, query compilers, ID normalizers,
# and Two-Stage Dual-Axis Triage Scorer (v2) for proteomic-gwas-catalogue.
#
# Single source of truth shared across discovery sources so PubMed E-utilities,
# Europe PMC REST API, NHGRI-EBI GWAS Catalog v2, Citation Snowballing, and
# Triage Scoring remain 100% synchronized.

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x

# ---------------------------------------------------------------------------
# 1. Canonical Controlled Vocabularies & Queue Schema Contract (C-1, C-2, M-1)
# ---------------------------------------------------------------------------

CANONICAL_SAMPLE_MATRICES <- c(
  "EDTA plasma",
  "Citrate plasma",
  "Heparin plasma",
  "Plasma (unspecified)",
  "Serum",
  "Whole blood / DBS",
  "CSF",
  "Brain tissue",
  "Urine",
  "Saliva",
  "Cell line / LCL",
  "Other tissue"
)

CANONICAL_PLATFORM_CATEGORIES <- c(
  "SomaScan",
  "Olink",
  "Mass Spectrometry",
  "NULISA",
  "Luminex",
  "Other Immunoassay"
)

CANONICAL_CANDIDATE_QUEUE_COLS <- c(
  # Original 9 UI-facing columns (order preserved for site/updates.qmd)
  "candidate_id",
  "run_id",
  "source",
  "identifier",
  "detected_platform",
  "detected_n",
  "title",
  "auto_classification",
  "curation_status",
  # 8 v2 Extraction Tier-1 machine-handoff columns (C-1)
  "pmid",
  "doi",
  "pmcid",
  "epmc_id",
  "matched_sources",
  "p_score",
  "s_score",
  "preprint_lineage_note"
)

# Fixed, transparent 6-paper landmark seed set for Citation Snowballing (Source 4)
LANDMARK_PGWAS_SEEDS <- data.frame(
  seed_key = c(
    "suhre_2017",
    "sun_2018",
    "folkersen_2020",
    "ferkingstad_2021",
    "sun_2023",
    "eldjarn_2023"
  ),
  pmid = c(
    "28240269",
    "29875488",
    "33067605",
    "34857953",
    "37794186",
    "37794188"
  ),
  doi = c(
    "10.1038/ncomms14357",
    "10.1038/s41586-018-0175-2",
    "10.1038/s42255-020-00287-2",
    "10.1038/s41588-021-00978-w",
    "10.1038/s41586-023-06592-6",
    "10.1038/s41586-023-06563-x"
  ),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# 2. Canonical Search Phrase Vectors (C-3, M-4)
# ---------------------------------------------------------------------------

# Block A1: Specific proteomic platforms, reagents, consortia, MS & multiplex panels
LEX_BLOCK_A1_PLATFORMS <- c(
  # Aptamer affinity
  "SomaScan", "SomaLogic", "SOMAmer", "SOMAmers", "aptamer", "aptamers",
  # Proximity extension & nucleic-acid immunoassays
  "Olink", "proximity extension", "SCALLOP", "UKB-PPP", "Pharma Proteomics Project",
  "NULISA", "NULISAseq", "Alamar",
  # Untargeted, targeted & nanoparticle mass spectrometry
  "Proteograph", "DIA-MS", "SWATH", "SWATH-MS", "data-independent acquisition",
  "tandem mass tag", "tandem mass tags", "TMT", "iTRAQ", "SILAC",
  "LC-MS/MS", "tandem mass spectrometry",
  "multiple reaction monitoring", "MRM", "parallel reaction monitoring", "PRM",
  "selected reaction monitoring", "SRM", "label-free quantification",
  "shotgun proteomics", "quantitative proteomics", "mass spectrometry proteomics",
  # Multiplex bead, electrochemiluminescence & microarray immunoassays
  "Luminex", "xMAP", "Bio-Plex", "Milliplex",
  "Rules-Based Medicine", "Myriad RBM", "Meso Scale Discovery",
  "electrochemiluminescence",
  "multiplex immunoassay", "multiplex immunoassays",
  "multiplexed immunoassay", "multiplexed immunoassays",
  "multiplex panel", "multiplexed panel",
  "protein microarray", "protein microarrays"
)

# Block A2: Proteome, biofluid/tissue proteins, cytokine/biomarker panels & MS peptides
LEX_BLOCK_A2_PROTEOME_ANALYTES <- c(
  # Proteome stems
  "proteome", "proteomes", "proteomic", "proteomics",
  "proteogenomic", "proteogenomics", "proteo-genomic", "proteo-genomics", "secretome",
  # Biofluid & tissue protein terms
  "plasma protein", "plasma proteins",
  "serum protein", "serum proteins",
  "circulating protein", "circulating proteins",
  "blood protein", "blood proteins",
  "CSF protein", "CSF proteins",
  "cerebrospinal fluid protein", "cerebrospinal fluid proteins",
  "urinary protein", "urinary proteins",
  "urine protein", "urine proteins",
  "brain protein", "brain proteins",
  "salivary protein", "salivary proteins",
  "protein abundance", "protein abundances",
  # Cytokine, chemokine & inflammatory/cardiovascular multiplex panels (C-3)
  "circulating cytokine", "circulating cytokines",
  "serum cytokine", "serum cytokines",
  "plasma cytokine", "plasma cytokines",
  "cytokine levels", "cytokine concentrations",
  "cytokines and chemokines", "cytokines and growth factors",
  "inflammatory proteins", "inflammatory biomarkers",
  "cardiovascular proteins", "cardiovascular biomarkers",
  "protein biomarker", "protein biomarkers",
  # Peptide-level MS terms (C-3)
  "plasma peptides", "serum peptides", "tryptic peptides"
)

LEX_MESH_PROTEOMICS <- c("Proteomics", "Proteome")

# Block B: Germline genetic association, GWAS, exome/WGS, linkage/scans & QTL terms (C-3, M-4)
LEX_BLOCK_B_GENETICS <- c(
  "genome-wide association", "genome wide association",
  "GWAS", "GWASs", "GWA study", "pGWAS",
  "genome-wide scan", "genome wide scan", "GWA scan", "genome-wide linkage",
  "exome-wide association", "exome sequencing", "whole-exome", "whole-genome sequencing",
  "rare variant association", "rare variant associations",
  "pQTL", "pQTLs", "cis-pQTL", "cis-pQTLs", "trans-pQTL", "trans-pQTLs",
  "protein quantitative trait", "quantitative trait loci", "quantitative trait locus",
  "genetic architecture", "genetic regulation",
  "genetic determinant", "genetic determinants", "genetic control",
  "genomic atlas", "genetic associations",
  "genetic variants influencing", "genetic variants associated",
  "genetic loci", "genetic loci associated", "heritability and genetic",
  "SNP array", "genotyping array", "single nucleotide polymorphisms"
)

LEX_MESH_GENETICS <- c("Genome-Wide Association Study", "Quantitative Trait Loci")

# Block C: Self-sufficient explicit pQTL / pGWAS anchor
LEX_BLOCK_C_PQTL_ANCHORS <- c(
  "pQTL", "pQTLs", "cis-pQTL", "cis-pQTLs", "trans-pQTL", "trans-pQTLs",
  "pGWAS", "protein quantitative trait loci", "protein quantitative trait locus"
)

# Ultra-broad snowballs & Stage-1 relevance guards (C-3: prevents pure metabolomics/MRI/eQTL leak)
LEX_SNOWBALL_BROAD_PROTEIN_TERMS <- c(
  "protein", "proteins", "proteome", "proteomes", "proteomic", "proteomics", "secretome",
  "peptide", "peptides", "cytokine", "cytokines", "chemokine", "chemokines",
  "biomarker", "biomarkers", "aptamer", "aptamers", "SOMAmer", "SOMAmers",
  "SomaScan", "SomaLogic", "Olink", "NULISA", "NULISAseq", "Luminex", "Bio-Plex", "Milliplex",
  "immunoassay", "immunoassays", "mass spectrometry", "DIA-MS", "SWATH", "SWATH-MS", "LC-MS/MS",
  "Proteograph", "pQTL", "pQTLs", "cis-pQTL", "cis-pQTLs", "trans-pQTL", "trans-pQTLs", "pGWAS"
)

LEX_SNOWBALL_BROAD_GENETICS_TERMS <- c(
  "GWAS", "GWASs", "pGWAS", "genome-wide", "genome wide", "exome-wide",
  "whole-genome", "whole-exome", "exome sequencing",
  "pQTL", "pQTLs", "cis-pQTL", "cis-pQTLs", "trans-pQTL", "trans-pQTLs",
  "quantitative trait loci", "quantitative trait locus",
  "SNP", "SNPs", "single nucleotide polymorphism", "single nucleotide polymorphisms",
  "genetic architecture", "genetic regulation", "genetic determinant", "genetic determinants",
  "genetic variants", "genetic associations", "genetic loci", "heritability",
  "genotyped", "imputed", "imputation"
)

BROAD_PROTEIN_REGEX <- paste0(
  "\\b(proteins?|proteom[a-z]*|secretom[a-z]*|peptid[a-z]*|cytokin[a-z]*|chemokin[a-z]*|",
  "interleukin[a-z0-9-]*|growth factors?|biomarkers?|aptamers?|somamers?|somascan|somalogic|",
  "olink|nulisa[a-z]*|alamar|luminex|xmap|bio-plex|milliplex|meso scale|immunoassay[a-z]*|",
  "mass spectrometry|dia-ms|swath[a-z-]*|lc-ms|tmt|itraq|silac|proteograph|",
  "pqtls?|cis-pqtls?|trans-pqtls?|pgwas)\\b"
)

BROAD_GENETICS_REGEX <- paste0(
  "\\b(gwas|gwass|pgwas|genome-wide|genome wide|exome-wide|whole-genome|whole-exome|",
  "exome sequenc[a-z]*|pqtls?|cis-pqtls?|trans-pqtls?|quantitative trait loc[a-z]*|",
  "snps?|single nucleotide polymorphism[a-z]*|genetic (architecture|regulation|determinants?|",
  "control|variants?|associations?|loci|mapping)|heritability|genotyped|imput(ed|ation))\\b"
)

# ---------------------------------------------------------------------------
# 3. Programmatic Query Compilers (PubMed & Europe PMC)
# ---------------------------------------------------------------------------

format_pubmed_or <- function(terms, tag = "tiab") {
  paste0("(", paste(sprintf('"%s"[%s]', terms, tag), collapse = " OR "), ")")
}

format_epmc_or <- function(terms) {
  paste0("(", paste(sprintf('"%s"', terms), collapse = " OR "), ")")
}

build_pubmed_query <- function(since_date = NULL, until_date = NULL,
                               date_field = c("edat", "dp"), verify_pmids = NULL) {
  date_field <- match.arg(date_field)
  b_a1 <- format_pubmed_or(LEX_BLOCK_A1_PLATFORMS, "tiab")
  b_a2 <- paste0(
    "(", paste(c(
      sprintf('"%s"[mh]', LEX_MESH_PROTEOMICS),
      sprintf('"%s"[tiab]', LEX_BLOCK_A2_PROTEOME_ANALYTES)
    ), collapse = " OR "), ")"
  )
  b_b <- paste0(
    "(", paste(c(
      sprintf('"%s"[mh]', LEX_MESH_GENETICS),
      sprintf('"%s"[tiab]', LEX_BLOCK_B_GENETICS)
    ), collapse = " OR "), ")"
  )
  # Targeted SNP-MeSH clause paired strictly with explicit proteomics/platform (C-3)
  b_snp_proteome <- paste0(
    '("Polymorphism, Single Nucleotide"[mh] AND (',
    b_a1, ' OR "Proteomics"[mh] OR "Proteome"[mh] OR "proteome"[tiab] OR "proteomic"[tiab] OR "proteomics"[tiab]))'
  )
  b_c <- format_pubmed_or(LEX_BLOCK_C_PQTL_ANCHORS, "tiab")

  # Note (M-2): "Published Erratum"[pt] is intentionally NOT excluded here so G1 can catch
  # errata/corrections linked to already-catalogued studies.
  core <- paste0(
    "(((((", b_a1, " OR ", b_a2, ") AND ", b_b, ") OR ", b_snp_proteome, " OR ", b_c, ") ",
    'NOT ("Editorial"[pt] OR "Comment"[pt] OR "News"[pt] OR "Retracted Publication"[pt])))'
  )

  if (!is.null(verify_pmids) && length(verify_pmids) > 0L) {
    uid_clause <- paste0("(", paste(sprintf("%s[uid]", verify_pmids), collapse = " OR "), ")")
    return(paste0(core, " AND ", uid_clause))
  }

  if (!is.null(since_date) && !is.null(until_date) && nzchar(since_date) && nzchar(until_date)) {
    d_from <- gsub("-", "/", since_date)
    d_to   <- gsub("-", "/", until_date)
    core <- paste0(core, sprintf(' AND ("%s"[%s] : "%s"[%s])', d_from, date_field, d_to, date_field))
  }
  core
}

build_epmc_date_clause <- function(since_date = NULL, until_date = NULL, include_idate = TRUE) {
  if (is.null(since_date) || is.null(until_date) || !nzchar(since_date) || !nzchar(until_date)) {
    return("")
  }
  if (isTRUE(include_idate)) {
    sprintf("(FIRST_PDATE:[%s TO %s] OR FIRST_IDATE:[%s TO %s])",
            since_date, until_date, since_date, until_date)
  } else {
    sprintf("(FIRST_PDATE:[%s TO %s])", since_date, until_date)
  }
}

build_epmc_prong2a_query <- function(since_date = NULL, until_date = NULL,
                                     include_idate = TRUE, verify_pmids = NULL) {
  prot_or <- format_epmc_or(unique(c(LEX_BLOCK_A1_PLATFORMS, LEX_BLOCK_A2_PROTEOME_ANALYTES)))
  gen_or  <- format_epmc_or(LEX_BLOCK_B_GENETICS)
  pqtl_or <- format_epmc_or(LEX_BLOCK_C_PQTL_ANCHORS)

  prot_fields <- paste0(
    "(TITLE:", prot_or, " OR ABSTRACT:", prot_or, " OR KW:", prot_or,
    ' OR MESH:"Proteomics" OR MESH:"Proteome")'
  )
  gen_fields <- paste0(
    "(TITLE:", gen_or, " OR ABSTRACT:", gen_or, " OR KW:", gen_or,
    ' OR MESH:"Genome-Wide Association Study" OR MESH:"Quantitative Trait Loci")'
  )
  pqtl_fields <- paste0("(TITLE:", pqtl_or, " OR ABSTRACT:", pqtl_or, " OR KW:", pqtl_or, ")")

  # Note (M-2): PUB_TYPE:"Erratum" is intentionally NOT excluded here so G1 can route
  # errata linked to catalogued studies; OPEN_ACCESS:y is NEVER added.
  core <- paste0(
    "((SRC:MED OR SRC:PMC OR SRC:PPR) AND ((",
    prot_fields, " AND ", gen_fields, ") OR ", pqtl_fields,
    ') NOT (PUB_TYPE:"Editorial" OR PUB_TYPE:"Comment" OR PUB_TYPE:"News"))'
  )

  if (!is.null(verify_pmids) && length(verify_pmids) > 0L) {
    ext_clause <- paste0("(", paste(sprintf("EXT_ID:%s", verify_pmids), collapse = " OR "), ")")
    return(paste0(core, " AND ", ext_clause))
  }

  dt <- build_epmc_date_clause(since_date, until_date, include_idate = include_idate)
  if (nzchar(dt)) paste0(core, " AND ", dt) else core
}

build_epmc_prong2b_query <- function(since_date = NULL, until_date = NULL,
                                     include_idate = TRUE, verify_pmids = NULL) {
  # Full-text METHODS rescue for multi-omics OA/preprint studies omitting platform/proteome in abstract
  meth_plat <- '("SomaScan" OR "SomaLogic" OR "Olink Explore" OR "Olink Target" OR "NULISAseq" OR "Proteograph" OR "DIA-MS" OR "SWATH-MS")'
  meth_gen  <- '("pQTL" OR "pQTLs" OR "cis-pQTL" OR "cis-pQTLs" OR "REGENIE" OR "BOLT-LMM" OR "fastGWA" OR "MatrixEQTL" OR "tensorQTL")'
  res_gen   <- '("pQTL" OR "pQTLs" OR "cis-pQTL" OR "cis-pQTLs")'

  core <- paste0(
    "((SRC:PMC OR SRC:PPR) AND METHODS:", meth_plat,
    " AND (METHODS:", meth_gen, " OR RESULTS:", res_gen, "))"
  )

  if (!is.null(verify_pmids) && length(verify_pmids) > 0L) {
    ext_clause <- paste0("(", paste(sprintf("EXT_ID:%s", verify_pmids), collapse = " OR "), ")")
    return(paste0(core, " AND ", ext_clause))
  }

  dt <- build_epmc_date_clause(since_date, until_date, include_idate = include_idate)
  if (nzchar(dt)) paste0(core, " AND ", dt) else core
}

build_epmc_snowball_query <- function(seed_pmids = LANDMARK_PGWAS_SEEDS$pmid,
                                      since_date = NULL, until_date = NULL,
                                      include_idate = TRUE, verify_pmids = NULL) {
  cites_clause <- paste0("(", paste(sprintf("CITES:%s_MED", seed_pmids), collapse = " OR "), ")")
  b_prot <- format_epmc_or(LEX_SNOWBALL_BROAD_PROTEIN_TERMS)
  b_gen  <- format_epmc_or(LEX_SNOWBALL_BROAD_GENETICS_TERMS)

  # C-3 Fix: Require BOTH broad protein/assay context AND broad genetics context
  prot_clause <- paste0(
    "(TITLE:", b_prot, " OR ABSTRACT:", b_prot, " OR KW:", b_prot,
    ' OR MESH:"Proteomics" OR MESH:"Proteome" OR MESH:"Blood Proteins")'
  )
  gen_clause <- paste0(
    "(TITLE:", b_gen, " OR ABSTRACT:", b_gen, " OR KW:", b_gen,
    ' OR MESH:"Genome-Wide Association Study" OR MESH:"Quantitative Trait Loci")'
  )

  core <- paste0("(", cites_clause, " AND ", prot_clause, " AND ", gen_clause, ")")

  if (!is.null(verify_pmids) && length(verify_pmids) > 0L) {
    ext_clause <- paste0("(", paste(sprintf("EXT_ID:%s", verify_pmids), collapse = " OR "), ")")
    return(paste0(core, " AND ", ext_clause))
  }

  dt <- build_epmc_date_clause(since_date, until_date, include_idate = include_idate)
  if (nzchar(dt)) paste0(core, " AND ", dt) else core
}

# ---------------------------------------------------------------------------
# 4. Shared Identifier Normalizers & ASCII Transliteration (M-5)
# ---------------------------------------------------------------------------

ascii_fold <- function(x) {
  if (is.null(x) || length(x) == 0L) return(character(0))
  x_chr <- as.character(x)
  x_chr[is.na(x_chr)] <- ""
  # Explicit Nordic/Slavic/Germanic letters that BSD/GNU iconv handle inconsistently
  special_from <- c("ø", "Ø", "ð", "Ð", "þ", "Þ", "ß", "æ", "Æ", "ł", "Ł", "đ", "Đ")
  special_to   <- c("o", "O", "d", "D", "th", "TH", "ss", "ae", "AE", "l", "L", "d", "D")
  for (k in seq_along(special_from)) {
    x_chr <- gsub(special_from[k], special_to[k], x_chr, fixed = TRUE)
  }
  # Transliterate only non-ASCII spans so BSD libiconv accent prefixes ('o, "u) are stripped
  # without disturbing genuine ASCII apostrophes (e.g. O'Connor) or punctuation.
  m <- gregexpr("[^[:ascii:]]+", x_chr, perl = TRUE)
  regmatches(x_chr, m) <- lapply(regmatches(x_chr, m), function(spans) {
    if (length(spans) == 0L) return(spans)
    tr <- iconv(spans, from = "UTF-8", to = "ASCII//TRANSLIT")
    tr[is.na(tr)] <- ""
    gsub("[^A-Za-z0-9]", "", tr)
  })
  x_chr
}

normalize_doi_one <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[1])) return("")
  d <- trimws(as.character(x[1]))
  if (!nzchar(d)) return("")
  d <- tryCatch(utils::URLdecode(d), error = function(e) d)
  d <- tolower(trimws(d))
  d <- sub("^https?://(dx\\.)?doi\\.org/", "", d)
  d <- sub("^doi\\s*:\\s*", "", d)
  # Iteratively strip trailing punctuation, .full/.pdf/.abstract/.short, and preprint version suffixes
  repeat {
    prev <- d
    d <- sub("[.,;:)>\\]\"']+$", "", d, perl = TRUE)
    d <- sub("(\\.full(\\.pdf)?|\\.pdf|\\.abstract|\\.short|/full|/abstract)$", "", d, perl = TRUE)
    if (grepl("^10\\.(1101|21203|26434|48550)/", d)) {
      d <- sub("(/v[0-9]+|v[0-9]+)$", "", d, perl = TRUE)
    }
    if (identical(d, prev)) break
  }
  if (!grepl("^10\\.[0-9]{4,9}/\\S+$", d, perl = TRUE)) return("")
  d
}

extract_dois_from_cell <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[1])) return(character(0))
  raw <- tryCatch(utils::URLdecode(as.character(x[1])), error = function(e) as.character(x[1]))
  if (!nzchar(trimws(raw))) return(character(0))
  m <- regmatches(raw, gregexpr("10\\.[0-9]{4,9}/[^\\s;,\"')\\]>]+", raw, perl = TRUE))[[1]]
  if (length(m) == 0L) return(character(0))
  normed <- vapply(m, normalize_doi_one, character(1), USE.NAMES = FALSE)
  unique(normed[nzchar(normed)])
}

extract_pmids_from_cell <- function(x, strip_dois_first = TRUE) {
  if (is.null(x) || length(x) == 0L || is.na(x[1])) return(character(0))
  raw <- as.character(x[1])
  if (!nzchar(trimws(raw))) return(character(0))
  if (isTRUE(strip_dois_first)) {
    raw <- gsub("10\\.[0-9]{4,9}/[^\\s;,\"')\\]>]+", " ", raw, perl = TRUE)
    raw <- gsub("PMC[0-9]+", " ", raw, ignore.case = TRUE, perl = TRUE)
  }
  m <- regmatches(raw, gregexpr("\\b[1-9][0-9]{5,8}\\b", raw, perl = TRUE))[[1]]
  unique(m[nzchar(m)])
}

normalize_pmid_one <- function(x) {
  ids <- extract_pmids_from_cell(x, strip_dois_first = TRUE)
  if (length(ids) == 0L) "" else ids[1L]
}

normalize_pmcid_one <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[1])) return("")
  m <- regmatches(toupper(as.character(x[1])), regexpr("PMC[0-9]+", toupper(as.character(x[1])), perl = TRUE))
  if (length(m) == 1L && nzchar(m)) m else ""
}

normalize_author_surname <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[1])) return("")
  a <- tolower(trimws(ascii_fold(x[1])))
  if (!nzchar(a)) return("")
  a <- sub("\\bet\\s+al\\.?.*$", "", a)
  # Split on semicolon or comma if "Surname, Initials" or "FirstAuthor; SecondAuthor"
  first_part <- trimws(strsplit(a, "[;,]")[[1]][1])
  # Remove standalone 1-2 letter initials with dots (e.g. "b. b. sun" -> "sun", "suhre k" -> "suhre")
  tokens <- unlist(strsplit(gsub("[^a-z0-9 -]", " ", first_part), "\\s+"))
  tokens <- tokens[nzchar(tokens)]
  if (length(tokens) == 0L) return("")
  non_init <- tokens[nchar(gsub("-", "", tokens)) >= 2L]
  if (length(non_init) == 0L) return(tokens[1L])
  # In PubMed/EuropePMC authorString ("Sun BB" or "Ferkingstad E"), surname is token 1 unless token 1 is initials
  if (grepl(",", a)) {
    non_init[1L]
  } else if (nchar(tokens[1L]) == 1L) {
    non_init[length(non_init)]
  } else {
    non_init[1L]
  }
}

STOPWORDS_EN <- c(
  "the", "and", "for", "with", "from", "that", "this", "into", "through",
  "across", "using", "between", "among", "study", "studies", "analysis",
  "analyses", "human", "humans", "identifies", "identify", "reveals", "reveal",
  "novel", "large", "scale", "preprint"
)

title_content_tokens <- function(title) {
  t_clean <- tolower(ascii_fold(title %||% ""))
  t_clean <- gsub("\\(preprint\\)", " ", t_clean)
  toks <- unlist(strsplit(gsub("[^a-z0-9]+", " ", t_clean), "\\s+"))
  toks <- unique(toks[nchar(toks) >= 3L & !(toks %in% STOPWORDS_EN)])
  toks
}

title_jaccard <- function(t1, t2) {
  w1 <- title_content_tokens(t1)
  w2 <- title_content_tokens(t2)
  if (length(w1) == 0L || length(w2) == 0L) return(0)
  length(intersect(w1, w2)) / length(union(w1, w2))
}

build_title_fingerprint <- function(first_author, year, title) {
  sur <- normalize_author_surname(first_author)
  yr  <- gsub("[^0-9]", "", as.character(year %||% ""))
  t_alnum <- substr(gsub("[^a-z0-9]", "", tolower(ascii_fold(title %||% ""))), 1L, 65L)
  if (!nzchar(t_alnum)) return("")
  paste(sur, yr, t_alnum, sep = "|")
}

# ---------------------------------------------------------------------------
# 5. Platform Detector & Sample-Size Extractor (v2 — Canonical Enum Aligned)
# ---------------------------------------------------------------------------

detect_platform_v2 <- function(text) {
  txt <- tolower(ascii_fold(text %||% ""))
  plats <- character(0)
  if (grepl("somascan|somalogic|\\bsomamers?\\b|\\baptamers?\\b", txt, perl = TRUE)) {
    plats <- c(plats, "SomaScan")
  }
  if (grepl("\\bolink\\b|proximity extension|pharma proteomics project|\\bukb-ppp\\b|\\bscallop\\b", txt, perl = TRUE)) {
    plats <- c(plats, "Olink")
  }
  if (grepl("mass spectrometry|\\bdia-ms\\b|\\bswath(-ms)?\\b|data-independent acquisition|\\btmt\\b|tandem mass tag|\\bitraq\\b|\\bsilac\\b|lc-ms/ms|\\bmrm\\b|\\bprm\\b|\\bsrm\\b|multiple reaction monitoring|parallel reaction monitoring|proteograph|shotgun proteom", txt, perl = TRUE)) {
    plats <- c(plats, "Mass Spectrometry")
  }
  if (grepl("\\bnulisa(seq)?\\b|\\balamar\\b", txt, perl = TRUE)) {
    plats <- c(plats, "NULISA")
  }
  if (grepl("luminex|\\bxmap\\b|bio-plex|milliplex", txt, perl = TRUE)) {
    plats <- c(plats, "Luminex")
  }
  if (grepl("rules-based medicine|myriad rbm|meso scale|electrochemiluminescence|multiplex(ed)? immunoassay|protein microarray", txt, perl = TRUE)) {
    plats <- c(plats, "Other Immunoassay")
  }
  if (length(plats) == 0L) "Unspecified Proteomics" else paste(unique(plats), collapse = " / ")
}

split_sentences_safe <- function(text) {
  txt <- gsub("\\s+", " ", trimws(text %||% ""))
  if (!nzchar(txt)) return(character(0))
  unlist(strsplit(txt, "(?<=[.!?])\\s+|;\\s+", perl = TRUE))
}

is_public_reuse_sentence <- function(sent_lower) {
  grepl(paste0(
    "\\b(obtained|extracted|derived|downloaded|retrieved|sourced|taken|selected|leveraged)\\s+from\\b|",
    "\\b(summary(-level)?\\s+(statistics|data)|gwas\\s+data|pqtl\\s+data)\\s+from\\s+(the\\s+)?",
    "(ukb-ppp|uk\\s+biobank\\s+pharma|decode|fenland|interval|scallop|aric|finngen|gwas\\s+catalog)\\b"
  ), sent_lower, perl = TRUE)
}

extract_sample_size_v2 <- function(text) {
  txt <- tolower(ascii_fold(text %||% ""))
  sents <- split_sentences_safe(txt)
  if (length(sents) == 0L) return("Check full text")

  # Keep sentences not describing parent biobank 500,000 background or public summary-stat downloads
  keep_sents <- sents[!vapply(sents, is_public_reuse_sentence, logical(1)) &
                      !grepl("\\b(biobank|cohort|study)\\s+(is\\s+a\\s+|of\\s+)?(approximately|~|over|more\\s+than)\\s+[0-9,]+", sents, perl = TRUE)]
  if (length(keep_sents) == 0L) keep_sents <- sents
  clean_txt <- paste(keep_sents, collapse = " . ")

  m1 <- regmatches(clean_txt, gregexpr("\\bn\\s*=\\s*[1-9][0-9,]*", clean_txt, perl = TRUE))[[1]]
  nums1 <- suppressWarnings(as.numeric(gsub("[^0-9]", "", m1)))

  m2 <- regmatches(clean_txt, gregexpr(paste0(
    "\\b([1-9][0-9,]{1,6})\\s+(?:[a-z-]+\\s+){0,2}",
    "(?:participants|individuals|subjects|icelanders|europeans|british|adults|children|twins|patients|donors|volunteers|brains|samples)\\b"
  ), clean_txt, perl = TRUE))[[1]]
  nums2 <- suppressWarnings(as.numeric(gsub("^([0-9,]+).*$", "\\1", gsub(",", "", m2))))

  all_n <- c(nums1, nums2)
  all_n <- all_n[!is.na(all_n) & all_n >= 10]
  if (length(all_n) == 0L) return("Check full text")
  format(max(all_n), big.mark = ",", scientific = FALSE, trim = TRUE)
}

# ---------------------------------------------------------------------------
# 6. Two-Stage Dual-Axis Triage Scorer v2 (C-3, M-2, M-3)
# ---------------------------------------------------------------------------

score_triage_v2 <- function(title,
                            abstract_text = "",
                            pub_types = character(0),
                            mesh_headings = character(0),
                            gwas_catalog_study_count = 0L,
                            erratum_targets_catalogued = FALSE) {
  ttl <- tolower(trimws(ascii_fold(title %||% "")))
  abs <- tolower(trimws(ascii_fold(abstract_text %||% "")))
  txt <- trimws(paste(ttl, abs))
  pts <- tolower(trimws(pub_types %||% character(0)))
  msh <- tolower(trimws(mesh_headings %||% character(0)))

  # -------------------------------------------------------------------------
  # STAGE 1: Hard Structural & Domain-Relevance Gates
  # -------------------------------------------------------------------------

  # G1a: Erratum / Correction routing (M-2)
  is_erratum <- any(pts %in% c("published erratum", "erratum", "correction")) ||
    grepl("^(author correction|publisher correction|correction|erratum|corrigendum)\\s*:", ttl, perl = TRUE)
  if (is_erratum) {
    if (isTRUE(erratum_targets_catalogued)) {
      return(list(
        stage = "Stage1_ErratumCatalogued", P_score = 7L, S_score = 0L,
        p_rules = "G1_CataloguedErratum(+7)", s_rules = character(0),
        classification = "Needs Human Check: Erratum on Catalogued Study"
      ))
    }
    return(list(
      stage = "Stage1_HardGate", P_score = 0L, S_score = 0L,
      p_rules = character(0), s_rules = "G1_UncataloguedErratum",
      classification = "Auto-Deprioritised: Non-Original (Erratum)"
    ))
  }

  # G1b: Editorial / Comment / News / Retraction
  if (any(pts %in% c("editorial", "comment", "news", "retracted publication"))) {
    return(list(
      stage = "Stage1_HardGate", P_score = 0L, S_score = 0L,
      p_rules = character(0), s_rules = "G1_NonOriginalPubType",
      classification = "Auto-Deprioritised: Non-Original (Editorial/Comment)"
    ))
  }

  # G0: Broad Proteomic/Protein/Peptide/Cytokine Signal Gate (C-3: blocks pure metabolomics/MRI/eQTL)
  has_protein_signal <- grepl(BROAD_PROTEIN_REGEX, txt, perl = TRUE) ||
    any(msh %in% c("proteomics", "proteome", "blood proteins")) ||
    (!is.na(gwas_catalog_study_count) && gwas_catalog_study_count >= 5L)
  if (nzchar(txt) && !has_protein_signal) {
    return(list(
      stage = "Stage1_HardGate", P_score = 0L, S_score = 0L,
      p_rules = character(0), s_rules = "G0_NoProteinSignal",
      classification = "Auto-Deprioritised: Non-Proteomic Study"
    ))
  }

  # Pre-compute core primary signals needed as safety guards in G2
  sents <- split_sentences_safe(txt)
  non_dl_sents <- sents[!vapply(sents, is_public_reuse_sentence, logical(1))]

  p2_hit <- grepl(paste0(
    "\\bwe (conducted|performed|carried out|undertook|report|present|describe) ([a-z-]+ ){0,3}",
    "(genome-wide association|genome-wide scan|gwa scan|exome-wide association|gwas|gwass|pgwas|pqtl mapping|",
    "protein quantitative trait loc(us|i) \\(?pqtl\\)? mapping|genetic mapping)\\b|",
    "\\b(genome-wide association stud(y|ies)|genome-wide scan|gwas|gwass|pqtl mapping|exome-wide association) ",
    "(of|for) ([0-9,]+ )?(human |blood |plasma |serum |csf |cerebrospinal |urinary |urine |brain |salivary |circulating |inflammatory |cardiovascular )?",
    "(protein|proteom|aptamer|somamer|peptid|cytokin|chemokin|biomarker)|",
    "\\b(genome-wide association stud(y|ies)|genome-wide scan|gwas|pqtl mapping) ([a-z0-9, -]{0,40})?(was|were) (performed|conducted|carried out)\\b|",
    "\\bgenetic variants influencing ([a-z-]+ ){0,3}(proteom|proteins?|peptides?|cytokines?)"
  ), txt, perl = TRUE)

  prot_count_pat <- paste0(
    "\\b([5-9]|[1-9][0-9]{1,4})\\s+([a-z-]+\\s+){0,2}",
    "(proteins?|aptamers?|somamers?|olink\\s+assays?|peptides?|cytokines?|chemokines?|biomarkers?|",
    "protein\\s+biomarkers?|protein\\s+targets?)\\b"
  )
  cohort_n_pat <- paste0(
    "\\b(in|across|of|among|from)\\s+[1-9][0-9,]{1,6}\\s+([a-z-]+\\s+){0,2}",
    "(participants|individuals|subjects|donors|volunteers|icelanders|europeans|british|scottish|finnish|",
    "adults|children|pediatric|neonates|infants|twins|brains|patients|samples|men|women|lcls?|cell\\s+lines)\\b"
  )
  p4_hit <- length(non_dl_sents) > 0L &&
    any(grepl(prot_count_pat, non_dl_sents, perl = TRUE)) &&
    any(grepl(cohort_n_pat, non_dl_sents, perl = TRUE))

  # G2: Non-Human Model Organism Gate (M-2 expanded human cohort/donor/LCL vocabulary + P2/P4 guard)
  has_nonhuman <- grepl("\\b(mice|mouse|murine|rat|rats|drosophila|zebrafish|arabidopsis|yeast|saccharomyces|caenorhabditis|bovine|porcine|canine|macaque)\\b", txt, perl = TRUE)
  has_human <- grepl(paste0(
    "\\b(human|humans|participant|participants|patient|patients|individuals|subjects|donors|volunteers|",
    "cases|controls|icelanders|europeans|british|scottish|finnish|swedish|chinese|japanese|african|hispanic|",
    "men|women|adults|children|pediatric|neonates|infants|newborns|mothers|offspring|families|twins|",
    "biobank|cohort|ancestry|lymphoblastoid|lcl|lcls|hapmap|1000 genomes|kora|framingham|aric|ages|interval|fenland|ukb|decode)\\b"
  ), txt, perl = TRUE) || ("humans" %in% msh)

  if (has_nonhuman && !has_human && !p2_hit && !p4_hit) {
    return(list(
      stage = "Stage1_HardGate", P_score = 0L, S_score = 0L,
      p_rules = character(0), s_rules = "G2_NonHumanOnly",
      classification = "Auto-Deprioritised: Non-Human Study"
    ))
  }

  # G3: Pure Epigenome-Wide / DNA Methylation Study without Germline pGWAS (M-2 fixed regex)
  has_ewas <- grepl("epigenome-wide|dna methylation|\\bewas\\b|cpg site|methylation quantitative trait", txt, perl = TRUE)
  has_germline_pgwas <- grepl(paste0(
    "\\b(pqtls?|cis-pqtls?|trans-pqtls?|protein quantitative trait|pgwas|",
    "(genome-wide association stud(y|ies)|gwas|gwass|qtl mapping|cis-qtl) (of|for) ([0-9,a-z -]+ )?(protein|proteom|aptamer|somamer|cytokin|peptid)|",
    "genotyped|whole-genome sequenc|whole-exome sequenc|imput(ed|ation))\\b"
  ), txt, perl = TRUE)
  if (has_ewas && !has_germline_pgwas) {
    return(list(
      stage = "Stage1_HardGate", P_score = 0L, S_score = 0L,
      p_rules = character(0), s_rules = "G3_EWAS_NoGermlinepGWAS",
      classification = "Auto-Deprioritised: Non-GWAS (EWAS)"
    ))
  }

  # -------------------------------------------------------------------------
  # STAGE 2: Dual-Axis Point Rubric (P_score vs S_score)
  # -------------------------------------------------------------------------
  P <- 0L
  S <- 0L
  p_rules <- character(0)
  s_rules <- character(0)

  # P1 (+5): Active or passive primary proteomic measurement / cohort profiling (M-3 expanded)
  if (grepl(paste0(
    "\\bwe (measured|quantified|profiled|assayed|characterized|characterised|generated proteomic)\\b|",
    "\\bprotein(s| levels| concentrations| profiles| abundances)? (were )?(measured|quantified|profiled|assayed) (with|using|by|in|across)\\b|",
    "\\bcharacterizing the [a-z -]*proteomic profiles of [0-9,]+\\b|",
    "\\b(proteomic|olink|somascan|somalogic|nulisa[a-z]*|mass spectrometry|dia-ms|swath[a-z-]*|lc-ms/ms|multiplex)\\s+",
    "([a-z-]+\\s+){0,2}(profiling|measurement|quantification|analysis|characteri[sz]ation|assay[a-z]*)\\s+(was|were)\\s+(performed|conducted|carried out|undertaken|applied)\\b|",
    "\\b(profiled|assayed|quantified|measured)\\s+([5-9]|[1-9][0-9,]+)\\s+([a-z-]+\\s+){0,2}(proteins?|aptamers?|somamers?|peptides?|cytokines?|chemokines?|biomarkers?)\\b"
  ), txt, perl = TRUE)) {
    P <- P + 5L; p_rules <- c(p_rules, "P1_DirectAssay(+5)")
  }

  # P2 (+5): Primary pGWAS / pQTL mapping execution
  if (p2_hit) {
    P <- P + 5L; p_rules <- c(p_rules, "P2_PrimaryMappingExec(+5)")
  }

  # P3 (+5): Quantitative pQTL discovery yield (C-3: requires explicit pQTL or protein-sentence context)
  p3_explicit <- grepl(paste0(
    "\\b(identif(ied|ies|y|ying)|found|detected|discovered|mapped|report(ed)?)\\s+[0-9,]+\\s+([a-z-]+\\s+){0,3}(pqtls?|cis-?pqtls?|trans-?pqtls?)\\b|",
    "\\b[0-9,]+\\s+(cis-?pqtls?|trans-?pqtls?|pqtls?|sentinel variants|genetic loci|loci|primary genetic associations|associations)\\s+",
    "(for|across|in|associated with|between sequence variants and levels of)\\s+([0-9,]+\\s+)?([a-z-]+\\s+){0,2}",
    "(proteins?|proteome|aptamers?|somamers?|peptides?|cytokines?|biomarkers?|targets?|assays?)\\b"
  ), txt, perl = TRUE)
  p3_sent_ctx <- length(non_dl_sents) > 0L && any(
    grepl("\\b(pqtls?|proteins?|proteom[a-z]*|aptamers?|somamers?|peptides?|cytokines?)\\b", non_dl_sents, perl = TRUE) &
    !grepl("\\b(metabolites?|lipids?|glycans?|methylation sites|cpgs?|mri)\\b", non_dl_sents, perl = TRUE) &
    grepl("\\b(identif(ied|ies|y|ying)|found|detected|discovered|mapped)\\s+[0-9,]+\\s+([a-z-]+\\s+){0,2}(primary\\s+)?(genetic\\s+)?(associations|loci|sentinel\\s+variants)\\b", non_dl_sents, perl = TRUE)
  )
  if (p3_explicit || p3_sent_ctx) {
    P <- P + 5L; p_rules <- c(p_rules, "P3_DiscoveryYield(+5)")
  }

  # P4 (+4): Multi-protein assay scale (>= 5) + human cohort N in non-public-download sentences (C-3, M-3)
  if (p4_hit) {
    P <- P + 4L; p_rules <- c(p_rules, "P4_AssayScaleAndN(+4)")
  }

  # P5 (+3): Genotyping / SNP array / WGS / WES / Imputation / Variant QC / Heritability / Cis-Trans (C-3)
  if (grepl(paste0(
    "\\b(genotyped|genotyping array|snp array|\\bsnps?\\b|single nucleotide polymorphism|genome-wide scan|",
    "whole-genome sequenc|wgs|whole-exome sequenc|wes|imputed|imputation|topmed|haplotype reference consortium|1000 genomes|",
    "sequence variants|rare variants|minor allele frequency|regenie|bolt-lmm|fastgwa|snptest|plink|saige|matrixeqtl|tensorqtl|",
    "heritability|epitope|protein-altering variant|trans\\s*-?pqtls?|cis\\s*-?pqtls?)\\b"
  ), txt, perl = TRUE)) {
    P <- P + 3L; p_rules <- c(p_rules, "P5_GenomicArchitecture(+3)")
  }

  # P6 (+2): High-priority biofluid/tissue, MS/NULISA platform, small-N/pediatric, or non-EUR ancestry
  if (grepl(paste0(
    "\\b(cerebrospinal fluid|\\bcsf\\b|dorsolateral prefrontal|brain tissue|post-mortem|urinary|urine|salivary|saliva|dried blood spot|",
    "lymphoblastoid|\\blcls?\\b|nanoparticle|proteograph|nulisa|nulisaseq|dia-ms|swath|tmt|lc-ms/ms|",
    "non-european|african ancestry|african american|east asian|chinese|japanese|south asian|hispanic|admixed|pediatric|children|neonates|infants|twins)\\b"
  ), txt, perl = TRUE)) {
    P <- P + 2L; p_rules <- c(p_rules, "P6_PriorityMatrixAncestryTech(+2)")
  }

  # P7 (+5): EBI GWAS Catalog multi-protein accession count (>= 5 GCST studies)
  if (!is.na(gwas_catalog_study_count) && gwas_catalog_study_count >= 5L) {
    P <- P + 5L; p_rules <- c(p_rules, "P7_GWASCatalogMultiTrait(+5)")
  }

  # P8 (+2): Open-access resource / atlas deposition statement
  if (grepl("\\b(open-access proteomics resource|valuable resource|summary statistics (are|have been) ([a-z]+ )?(available|deposited)|genomic atlas of the [a-z -]*proteome)\\b", txt, perl = TRUE)) {
    P <- P + 2L; p_rules <- c(p_rules, "P8_ResourceDeposition(+2)")
  }

  # P9 (+5): Primary consortium / multi-cohort pGWAS meta-analysis protection
  if (grepl("\\b(meta-analys(is|ed|ing)|combined analysis) of ([5-9]|[1-9][0-9,]+) ([a-z-]+ )?(proteins?|biomarkers?|cytokines?|peptides?)\\b|\\bscallop consortium\\b|\\b[0-9,]+ cardiovascular proteins in [0-9,]+ individuals\\b", txt, perl = TRUE) &&
      !grepl("\\btwo-sample mendelian\\b", txt, perl = TRUE)) {
    P <- P + 5L; p_rules <- c(p_rules, "P9_PrimaryConsortiumMeta(+5)")
  }

  # --- Secondary Summary-Stat Consumer Signals (S_score) ---
  if (grepl(paste0(
    "\\b(summary(-level)? (statistics|data)|genetic instruments|instrumental variables|pqtl data|pqtls?|proteomic data) ",
    "(were )?(obtained|extracted|derived|downloaded|retrieved|sourced|taken|leveraged|selected) from\\b|",
    "\\busing (publicly available |published |summary-level )?(gwas |pqtl |proteomic )?(summary )?(data|statistics) from (the )?",
    "(uk biobank pharma|ukb-ppp|decode|fenland|interval|scallop|aric|finngen)\\b"
  ), txt, perl = TRUE)) {
    S <- S + 7L; s_rules <- c(s_rules, "S1_PublicSumStatsReuse(+7)")
  }

  if (grepl("\\b(inverse-variance weighted|\\bivw\\b|mr-egger|weighted median|mr-presso|wald ratio|cochran'?s q|egger intercept|leave-one-out|summary-data-based mendelian|smr and heidi|\\bheidi\\b|coloc\\.abf|pp\\.h4\\s*[>=])\\b", txt, perl = TRUE)) {
    S <- S + 6L; s_rules <- c(s_rules, "S2_2SMREstimatorBoilerplate(+6)")
  }

  if (grepl("\\b(mendelian randomi[sz]ation|two-sample|proteome-wide association stud(y|ies)|\\bpwas\\b)\\b", ttl, perl = TRUE) ||
      grepl("\\b(therapeutic|druggable|drug|causal|potential|novel) targets? for [a-z0-9 ,/-]+\\b", ttl, perl = TRUE) ||
      grepl("\\bcausal (role|effect|relationship|association|link)s? (of|between) ([a-z]+ )?proteins? (on|in|and|with)\\b", ttl, perl = TRUE)) {
    S <- S + 4L; s_rules <- c(s_rules, "S3_SecondaryTitleTemplate(+4)")
  }

  pub_cohorts <- sum(vapply(
    c("\\bukb-ppp\\b|pharma proteomics project", "\\bdecode\\b", "\\bfenland\\b", "\\baric\\b", "\\binterval\\b", "\\bscallop\\b"),
    function(pat) grepl(pat, txt, perl = TRUE),
    logical(1)
  ))
  has_outcome_db <- grepl("\\b(finngen|ieu open gwas|opengwas|megastroke|diagram consortium|cardiogram|pgc|igap)\\b", txt, perl = TRUE)
  if (pub_cohorts >= 2L && has_outcome_db) {
    S <- S + 4L; s_rules <- c(s_rules, "S4_PublicCohortOutcomePair(+4)")
  }

  if (grepl("\\b(fusion|predixcan|s-predixcan|utmost|jti)\\b", txt, perl = TRUE) &&
      grepl("\\b(pwas|proteome-wide association|genetically predicted protein)\\b", txt, perl = TRUE) &&
      !grepl("\\b(train(ed|ing)|build(ing)?|models for proteome-wide|cis-heritability)\\b", txt, perl = TRUE)) {
    S <- S + 5L; s_rules <- c(s_rules, "S5_PWASWeightConsumer(+5)")
  }

  # G4: Narrative Review Gate (only if P_score < 5)
  if ((any(pts == "review") || grepl("^(a |systematic )?review of\\b|\\bnarrative review\\b", ttl, perl = TRUE)) && P < 5L) {
    return(list(
      stage = "Stage1_ReviewGate", P_score = P, S_score = S,
      p_rules = p_rules, s_rules = c(s_rules, "G4_ReviewLowPrimary"),
      classification = "Auto-Deprioritised: Review Article"
    ))
  }

  # -------------------------------------------------------------------------
  # STAGE 2 DECISION MATRIX (Asymmetric Primary-Protective Thresholds)
  # -------------------------------------------------------------------------
  cls <- if (P >= 7L && (S <= 6L || P >= S + 2L)) {
    "Likely Primary pGWAS"
  } else if ((P >= 5L && S > 6L && P < S + 2L) ||
             (P >= 3L && P < 7L && S < 5L) ||
             (!nzchar(abs) && S < 4L)) {
    "Needs Human Check: Possible Primary / Hybrid"
  } else if (S >= 4L && P < 5L) {
    "Auto-Deprioritised: Secondary 2SMR / PWAS"
  } else {
    "Auto-Deprioritised: Low-Signal / Non-pGWAS"
  }

  list(
    stage = "Stage2_DualAxis",
    P_score = as.integer(P),
    S_score = as.integer(S),
    p_rules = p_rules,
    s_rules = s_rules,
    classification = cls
  )
}
