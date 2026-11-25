# Canonical lexicons, controlled vocabularies, query compilers, ID normalizers,
# and Unified Two-Stage Dual-Axis Triage Scorer (v4.7) for proteomic-gwas-catalogue.
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
  "Plasma + serum (mixed)",
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
  # 8 Extraction Tier-1 machine-handoff columns (C-1)
  "pmid",
  "doi",
  "pmcid",
  "epmc_id",
  "matched_sources",
  "p_score",
  "s_score",
  "preprint_lineage_note"
)

# Transparent 13-paper landmark seed set for Citation Snowballing (Source 4),
# spanning MS, Olink, SomaScan, CSF/brain, multi-ancestry, and rare-variant exome pGWAS.
LANDMARK_PGWAS_SEEDS <- data.frame(
  seed_key = c(
    "johansson_2013",
    "liu_2015",
    "suhre_2017",
    "sun_2018",
    "emilsson_2018",
    "folkersen_2020",
    "ferkingstad_2021",
    "pietzner_2021",
    "yang_2021",
    "zhang_2022",
    "sun_2023",
    "eldjarn_2023",
    "dhindsa_2023"
  ),
  pmid = c(
    "23487758",
    "25652787",
    "28240269",
    "29875488",
    "30072576",
    "33067605",
    "34857953",
    "34648354",
    "34239129",
    "35501419",
    "37794186",
    "37794188",
    "37794183"
  ),
  doi = c(
    "10.1073/pnas.1217238110",
    "10.15252/msb.20145728",
    "10.1038/ncomms14357",
    "10.1038/s41586-018-0175-2",
    "10.1126/science.aaq1327",
    "10.1038/s42255-020-00287-2",
    "10.1038/s41588-021-00978-w",
    "10.1126/science.abj1541",
    "10.1038/s41593-021-00886-6",
    "10.1038/s41588-022-01051-w",
    "10.1038/s41586-023-06592-6",
    "10.1038/s41586-023-06563-x",
    "10.1038/s41586-023-06547-x"
  ),
  pub_year = c(
    2013L, 2015L, 2017L, 2018L, 2018L,
    2020L, 2021L, 2021L, 2021L, 2022L,
    2023L, 2023L, 2023L
  ),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# 2. Canonical Search Phrase Vectors (v4.0 Unified Consensus Lexicon)
# ---------------------------------------------------------------------------

# Block A1: Specific proteomic platforms, reagents, consortia, MS & multiplex panels.
# Ambiguous 3-letter acronyms ("TMT", "MRM", "PRM", "SRM") and bare "aptamer" are scoped
# to proteomic phrases so Trail Making Test (TMT) and RNA-aptamer switch papers do not leak.
LEX_BLOCK_A1_PLATFORMS <- c(
  # Aptamer affinity ("SomaScan" boolean-absorbs "SomaScan 7K"/"11K")
  "SomaScan", "SomaLogic", "SOMAmer", "SOMAmers",
  "aptamer-based proteomics", "aptamer-based",
  # Proximity extension & nucleic-acid immunoassays ("Olink" absorbs "Olink Explore/HT/Target")
  "Olink", "Explore 3072", "Explore HT",
  "proximity extension", "SCALLOP", "UKB-PPP", "Pharma Proteomics Project",
  "NULISA", "NULISAseq", "Alamar",
  # Untargeted, targeted & nanoparticle mass spectrometry ("Proteograph" absorbs "Seer Proteograph")
  "Proteograph", "Orbitrap Astral", "timsTOF", "DIA-PASEF",
  "DIA-MS", "SWATH", "SWATH-MS", "data-independent acquisition",
  "tandem mass tag", "tandem mass tags", "TMT-based", "TMT proteomics", "iTRAQ", "SILAC",
  "LC-MS/MS", "tandem mass spectrometry",
  "multiple reaction monitoring", "parallel reaction monitoring", "selected reaction monitoring",
  "label-free quantification", "shotgun proteomics", "quantitative proteomics",
  "mass spectrometry proteomics", "nanoparticle-based proteomics",
  "affinity proteomics", "affinity-based proteomics",
  # Multiplex bead, Simoa, electrochemiluminescence & microarray immunoassays
  "Luminex", "xMAP", "x-MAP", "Bio-Plex", "Milliplex", "ProcartaPlex",
  "Simoa", "Quanterix",
  "Rules-Based Medicine", "Myriad RBM", "Meso Scale Discovery",
  "electrochemiluminescence",
  "multiplex immunoassay", "multiplex immunoassays",
  "multiplexed immunoassay", "multiplexed immunoassays",
  "multiplex panel", "multiplexed panel",
  "protein microarray", "protein microarrays"
)

# Block A2: Proteome, biofluid/tissue proteins, cytokine/complement/adipokine panels & MS peptides
LEX_BLOCK_A2_PROTEOME_ANALYTES <- c(
  # Proteome stems
  "proteome", "proteomes", "proteomic", "proteomics",
  "proteogenomic", "proteogenomics", "proteo-genomic", "proteo-genomics", "secretome",
  # Biofluid & tissue protein/proteome terms
  "plasma protein", "plasma proteins",
  "serum protein", "serum proteins",
  "circulating protein", "circulating proteins",
  "blood protein", "blood proteins",
  "CSF protein", "CSF proteins", "CSF proteome",
  "cerebrospinal fluid protein", "cerebrospinal fluid proteins", "cerebrospinal fluid proteome",
  "urinary protein", "urinary proteins", "urine protein", "urine proteins", "urine proteome", "urinary proteome",
  "brain protein", "brain proteins", "brain proteome",
  "salivary protein", "salivary proteins", "saliva proteome",
  "tissue protein", "tissue proteins", "tissue proteome", "tissue proteomics",
  "kidney proteins", "kidney proteome", "liver proteins", "liver proteome",
  "dried blood spot", "dried blood spots",
  "protein abundance", "protein abundances",
  "protein levels", "protein concentrations", "protein profiles", "protein profiling",
  # Cytokine, chemokine, complement, adipokine & inflammatory/cardiovascular multiplex panels
  "circulating cytokine", "circulating cytokines",
  "serum cytokine", "serum cytokines",
  "plasma cytokine", "plasma cytokines",
  "cytokine levels", "cytokine concentrations",
  "cytokines and chemokines", "cytokines and growth factors",
  "complement proteins", "complement system proteins",
  "adipokines", "adipocytokines",
  "inflammatory proteins", "inflammatory biomarkers",
  "cardiovascular proteins", "cardiovascular biomarkers", "cardiometabolic biomarkers",
  "neurological proteins", "neurological biomarkers",
  "tumor biomarkers", "tumour biomarkers",
  "protein biomarker", "protein biomarkers",
  "cell-surface proteins", "cell surface proteins",
  # Peptide-level MS terms
  "plasma peptide", "plasma peptides", "serum peptide", "serum peptides", "peptide levels", "tryptic peptides"
)

LEX_MESH_PROTEOMICS <- c("Proteomics", "Proteome")

# Block B: Germline genetic association, GWAS, exome/WGS, SV/TR/qQTL/vQTL/sd-pQTL/molQTL, admixture & knockouts
LEX_BLOCK_B_GENETICS <- c(
  "genome-wide association", "genome wide association",
  "genome-wide genetic association", "genome wide genetic association",
  "GWAS", "GWASs", "GWA study", "pGWAS",
  "genome-wide scan", "genome wide scan", "GWA scan", "genome-wide linkage",
  "admixture mapping", "admixture-mapping",
  "exome-wide association", "exome-wide", "exome sequencing",
  "whole-exome", "whole exome", "whole-genome sequencing", "whole genome sequencing",
  "long-read genome sequencing", "structural variants", "structural variation",
  "tandem repeat", "tandem repeats",
  "rare variant association", "rare variant associations", "rare coding variants", "rare noncoding",
  "human knockouts", "human knockout", "loss-of-function variants", "loss-of-function variant", "parent-of-origin",
  "genetic predictors", "proteome-wide association",
  "pQTL", "pQTLs", "cis-pQTL", "cis-pQTLs", "trans-pQTL", "trans-pQTLs", "sd-pQTL", "sd-pQTLs",
  "SV-pQTL", "SV-pQTLs", "TR-pQTL", "TR-pQTLs",
  "protein QTL", "protein QTLs", "protein-QTL", "protein-QTLs",
  "vQTL", "vQTLs", "cis-vQTL", "cis-vQTLs", "variance quantitative trait",
  "quantile QTL", "quantile QTLs", "qQTL", "qQTLs",
  "SV-QTL", "SV-QTLs", "TR-xQTL", "TR-xQTLs",
  "molQTL", "molQTLs", "molecular QTL", "molecular QTLs", "xQTL", "xQTLs",
  "molecular quantitative trait loci", "molecular quantitative trait locus",
  "multi-omic quantitative trait loci",
  "QTL mapping", "cis-QTL", "cis-QTLs", "trans-QTL", "trans-QTLs",
  "protein quantitative trait", "protein quantitative loci",
  "quantitative trait loci", "quantitative trait locus",
  "genotype-proteome", "genetic architecture", "genetic regulation",
  "genetic determinant", "genetic determinants", "genetic control",
  "genomic atlas", "genetic associations", "genetic prediction models",
  "genetic variants influencing", "genetic variants associated",
  "genetic loci", "heritability",
  "SNP array", "genotyping array", "single nucleotide polymorphisms"
)

LEX_MESH_GENETICS <- c("Genome-Wide Association Study", "Quantitative Trait Loci")

# Block C: Self-sufficient explicit pQTL / vQTL / pGWAS anchor
LEX_BLOCK_C_PQTL_ANCHORS <- c(
  "pQTL", "pQTLs", "cis-pQTL", "cis-pQTLs", "trans-pQTL", "trans-pQTLs",
  "sd-pQTL", "sd-pQTLs", "SV-pQTL", "SV-pQTLs", "TR-pQTL", "TR-pQTLs",
  "protein QTL", "protein QTLs", "protein-QTL", "protein-QTLs",
  "cis-vQTL", "cis-vQTLs", "pGWAS",
  "protein quantitative trait loci", "protein quantitative trait locus",
  "protein quantitative loci"
)

# Block D: Multi-omic molecular QTL umbrella anchors (captures cross-modality molQTL /
# xQTL / TR-xQTL / SV-QTL / qQTL studies whose abstracts say "molecular phenotypes" and
# specify primary proteomic / pQTL layers in full-text Methods/Results).
LEX_BLOCK_D_MULTIOQTL_ANCHORS <- c(
  "molQTL", "molQTLs", "molecular QTL", "molecular QTLs",
  "xQTL", "xQTLs", "TR-xQTL", "TR-xQTLs",
  "SV-QTL", "SV-QTLs", "quantile QTL", "quantile QTLs", "qQTL", "qQTLs",
  "molecular quantitative trait loci", "molecular quantitative trait locus",
  "multi-omic quantitative trait loci"
)

# Broad snowball & Stage-1 relevance guards
LEX_SNOWBALL_BROAD_PROTEIN_TERMS <- c(
  "protein", "proteins", "proteome", "proteomes", "proteomic", "proteomics", "secretome",
  "proteogenomic", "proteogenomics",
  "peptide", "peptides", "cytokine", "cytokines", "chemokine", "chemokines", "adipokines",
  "biomarker", "biomarkers", "aptamer", "aptamers", "SOMAmer", "SOMAmers",
  "SomaScan", "SomaLogic", "Olink", "NULISA", "NULISAseq", "Luminex", "Bio-Plex", "Milliplex",
  "Simoa", "immunoassay", "immunoassays", "mass spectrometry", "DIA-MS", "SWATH", "SWATH-MS",
  "LC-MS/MS", "Proteograph", "Orbitrap Astral",
  "pQTL", "pQTLs", "cis-pQTL", "cis-pQTLs", "trans-pQTL", "trans-pQTLs",
  "sd-pQTL", "sd-pQTLs", "SV-pQTL", "SV-pQTLs", "TR-pQTL", "TR-pQTLs",
  "protein QTL", "protein QTLs", "vQTL", "vQTLs", "pGWAS"
)

LEX_SNOWBALL_BROAD_GENETICS_TERMS <- c(
  "GWAS", "GWASs", "pGWAS", "genome-wide", "genome wide", "exome-wide",
  "whole-genome", "whole genome", "whole-exome", "whole exome", "exome sequencing",
  "admixture mapping", "admixture-mapping", "human knockouts", "parent-of-origin",
  "pQTL", "pQTLs", "cis-pQTL", "cis-pQTLs", "trans-pQTL", "trans-pQTLs", "sd-pQTL", "sd-pQTLs",
  "SV-pQTL", "SV-pQTLs", "TR-pQTL", "TR-pQTLs", "SV-QTL", "SV-QTLs", "TR-xQTL", "TR-xQTLs", "qQTL", "qQTLs",
  "protein QTL", "protein QTLs", "vQTL", "vQTLs", "cis-vQTL", "cis-vQTLs",
  "quantitative trait loci", "quantitative trait locus", "QTL mapping", "cis-QTL",
  "SNP", "SNPs", "single nucleotide polymorphism", "single nucleotide polymorphisms",
  "genetic architecture", "genetic regulation", "genetic determinant", "genetic determinants",
  "genetic variants", "genetic variation", "genetic associations", "genetic loci", "heritability",
  "genotyped", "imputed", "imputation"
)

BROAD_PROTEIN_REGEX <- paste0(
  "\\b(proteins?|proteom[a-z]*|proteogenom[a-z]*|secretom[a-z]*|peptid[a-z]*|cytokin[a-z]*|chemokin[a-z]*|",
  "adipokin[a-z]*|adipocytokin[a-z]*|interleukin[a-z0-9-]*|growth factors?|biomarkers?|aptamers?|somamers?|",
  "somascan|somalogic|olink|nulisa[a-z]*|alamar|luminex|x-?map|bio-plex|milliplex|procartaplex|simoa|quanterix|",
  "meso scale|immunoassay[a-z]*|mass spectrometry|dia-ms|dia-pasef|swath[a-z-]*|lc-ms|tandem mass tags?|tmt-based|",
  "itraq|silac|proteograph|orbitrap|timstof|pqtls?|cis-?pqtls?|trans-?pqtls?|sd-?pqtls?|sv-?pqtls?|tr-?pqtls?|protein[- ]qtls?|",
  "vqtls?|cis-?vqtls?|pgwas)\\b"
)

BROAD_GENETICS_REGEX <- paste0(
  "\\b(gwas|gwass|pgwas|genome-wide|genome wide|exome-wide|whole[ -]genome|whole[ -]exome|",
  "exome sequenc[a-z]*|admixture[- ]mapping|human knockouts?|parent-of-origin|",
  "pqtls?|cis-?pqtls?|trans-?pqtls?|sd-?pqtls?|sv-?pqtls?|tr-?pqtls?|protein[- ]qtls?|vqtls?|cis-?vqtls?|",
  "quantitative trait loc[a-z]*|qtl mapping|cis-?qtls?|molqtls?|xqtls?|sv-?qtls?|tr-?xqtls?|qqtls?|",
  "snps?|single[- ]nucleotide polymorphisms?|genetic (architecture|regulation|determinants?|",
  "control|variants?|variation|associations?|loci|mapping|influences?|prediction models?)|heritability|genotyped|imput(ed|ation))\\b"
)

# ---------------------------------------------------------------------------
# 3. Programmatic Query Compilers (PubMed & Sharded TITLE_ABS Europe PMC)
# ---------------------------------------------------------------------------

format_pubmed_or <- function(terms, tag = "tiab") {
  paste0("(", paste(sprintf('"%s"[%s]', terms, tag), collapse = " OR "), ")")
}

format_epmc_or <- function(terms) {
  paste0("(", paste(sprintf('"%s"', terms), collapse = " OR "), ")")
}

# Remove pure whitespace-token superphrases within the SAME Europe PMC OR-clause
# (e.g., `"proteome"` Boolean-absorbs `"cerebrospinal fluid proteome"` in Lucene;
# hyphens are preserved 100% untouched so analyzer-specific hyphen rules cannot vary).
compact_epmc_space_supersets <- function(terms) {
  ws_toks <- lapply(terms, function(t) tolower(unlist(strsplit(trimws(t), "\\s+"))))
  is_ws_subseq <- function(short_t, long_t) {
    ns <- length(short_t); nl <- length(long_t)
    if (ns >= nl) return(FALSE)
    for (start in seq_len(nl - ns + 1L)) {
      if (identical(long_t[start:(start + ns - 1L)], short_t)) return(TRUE)
    }
    FALSE
  }
  redundant <- vapply(seq_along(terms), function(i) {
    any(vapply(seq_along(terms), function(j) i != j && is_ws_subseq(ws_toks[[j]], ws_toks[[i]]), logical(1)))
  }, logical(1))
  terms[!redundant]
}

build_pubmed_query <- function(since_date = NULL, until_date = NULL,
                               date_field = c("edat", "dp", "dp_or_edat"), verify_pmids = NULL) {
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
  b_snp_proteome <- paste0(
    '("Polymorphism, Single Nucleotide"[mh] AND (',
    b_a1, ' OR "Proteomics"[mh] OR "Proteome"[mh] OR "proteome"[tiab] OR "proteomic"[tiab] OR "proteomics"[tiab]))'
  )
  b_c <- format_pubmed_or(unique(c(LEX_BLOCK_C_PQTL_ANCHORS, LEX_BLOCK_D_MULTIOQTL_ANCHORS)), "tiab")

  # "Published Erratum"[pt] intentionally NOT excluded so G1a can flag errata on catalogued studies
  core <- paste0(
    "(((((", b_a1, " OR ", b_a2, ") AND ", b_b, ") OR ", b_snp_proteome, " OR ", b_c, ") ",
    'NOT ("Editorial"[pt] OR "Comment"[pt] OR "News"[pt] OR "Retracted Publication"[pt] OR "Congress"[pt])))'
  )

  if (!is.null(verify_pmids) && length(verify_pmids) > 0L) {
    uid_clause <- paste0("(", paste(sprintf("%s[uid]", verify_pmids), collapse = " OR "), ")")
    return(paste0(core, " AND ", uid_clause))
  }

  if (!is.null(since_date) && !is.null(until_date) && nzchar(since_date) && nzchar(until_date)) {
    d_from <- gsub("-", "/", since_date)
    d_to   <- gsub("-", "/", until_date)
    if (identical(date_field, "dp_or_edat")) {
      core <- paste0(core, sprintf(' AND (("%s"[dp] : "%s"[dp]) OR ("%s"[edat] : "%s"[edat]))',
                                   d_from, d_to, d_from, d_to))
    } else {
      core <- paste0(core, sprintf(' AND ("%s"[%s] : "%s"[%s])', d_from, date_field, d_to, date_field))
    }
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

# Europe PMC Prong 2A is sharded into two compact TITLE_ABS queries (< 4,000 chars each,
# Safeguard 6: LEX_BLOCK_C_PQTL_ANCHORS in 2A-1, LEX_BLOCK_D_MULTIOQTL_ANCHORS in 2A-2)
# with strictly ZERO `MESH:` pseudo-field (which is invalid in Europe PMC Lucene syntax,
# corrupting parenthesized boolean evaluation and suppressing SRC:PPR preprints).
EPMC_PUBTYPE_EXCLUSION <- paste0(
  'NOT (PUB_TYPE:"Editorial" OR PUB_TYPE:"Comment" OR PUB_TYPE:"News" OR ',
  'PUB_TYPE:"Meeting Abstract" OR PUB_TYPE:"Abstract" OR PUB_TYPE:"Congress")'
)

build_epmc_prong2a_platforms_query <- function(since_date = NULL, until_date = NULL,
                                               include_idate = TRUE, verify_pmids = NULL) {
  plat_or <- format_epmc_or(compact_epmc_space_supersets(LEX_BLOCK_A1_PLATFORMS))
  # Boolean absorption identity: ((A AND (B OR C)) OR C) == ((A AND B) OR C)
  gen_or  <- format_epmc_or(compact_epmc_space_supersets(setdiff(LEX_BLOCK_B_GENETICS, LEX_BLOCK_C_PQTL_ANCHORS)))
  pqtl_or <- format_epmc_or(LEX_BLOCK_C_PQTL_ANCHORS)

  core <- paste0(
    "((SRC:MED OR SRC:PMC OR SRC:PPR) AND (((TITLE_ABS:",
    plat_or, ") AND (TITLE_ABS:", gen_or, ")) OR (TITLE_ABS:", pqtl_or, ")) ",
    EPMC_PUBTYPE_EXCLUSION, ")"
  )

  if (!is.null(verify_pmids) && length(verify_pmids) > 0L) {
    ext_clause <- paste0("(", paste(sprintf("EXT_ID:%s", verify_pmids), collapse = " OR "), ")")
    return(paste0(core, " AND ", ext_clause))
  }

  dt <- build_epmc_date_clause(since_date, until_date, include_idate = include_idate)
  if (nzchar(dt)) paste0(core, " AND ", dt) else core
}

build_epmc_prong2a_proteome_query <- function(since_date = NULL, until_date = NULL,
                                              include_idate = TRUE, verify_pmids = NULL) {
  ana_or       <- format_epmc_or(compact_epmc_space_supersets(LEX_BLOCK_A2_PROTEOME_ANALYTES))
  # Boolean absorption identity: ((A AND (B OR D)) OR D) == ((A AND B) OR D)
  gen_or       <- format_epmc_or(compact_epmc_space_supersets(setdiff(LEX_BLOCK_B_GENETICS, LEX_BLOCK_D_MULTIOQTL_ANCHORS)))
  multioqtl_or <- format_epmc_or(LEX_BLOCK_D_MULTIOQTL_ANCHORS)

  core <- paste0(
    "((SRC:MED OR SRC:PMC OR SRC:PPR) AND (((TITLE_ABS:",
    ana_or, ") AND (TITLE_ABS:", gen_or, ")) OR (TITLE_ABS:", multioqtl_or, ")) ",
    EPMC_PUBTYPE_EXCLUSION, ")"
  )

  if (!is.null(verify_pmids) && length(verify_pmids) > 0L) {
    ext_clause <- paste0("(", paste(sprintf("EXT_ID:%s", verify_pmids), collapse = " OR "), ")")
    return(paste0(core, " AND ", ext_clause))
  }

  dt <- build_epmc_date_clause(since_date, until_date, include_idate = include_idate)
  if (nzchar(dt)) paste0(core, " AND ", dt) else core
}

# Combined helper retained for backward compatibility (TITLE_ABS without invalid MESH:)
build_epmc_prong2a_query <- function(since_date = NULL, until_date = NULL,
                                     include_idate = TRUE, verify_pmids = NULL) {
  prot_or <- format_epmc_or(compact_epmc_space_supersets(unique(c(LEX_BLOCK_A1_PLATFORMS, LEX_BLOCK_A2_PROTEOME_ANALYTES))))
  gen_or  <- format_epmc_or(compact_epmc_space_supersets(LEX_BLOCK_B_GENETICS))
  pqtl_or <- format_epmc_or(unique(c(LEX_BLOCK_C_PQTL_ANCHORS, LEX_BLOCK_D_MULTIOQTL_ANCHORS)))

  core <- paste0(
    "((SRC:MED OR SRC:PMC OR SRC:PPR) AND (((TITLE_ABS:",
    prot_or, ") AND (TITLE_ABS:", gen_or, ")) OR (TITLE_ABS:", pqtl_or, ")) ",
    EPMC_PUBTYPE_EXCLUSION, ")"
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
  meth_plat <- paste0(
    '("SomaScan" OR "SOMAscan" OR "SomaLogic" OR "SOMAmer" OR "aptamer-based proteomics" OR ',
    '"Olink" OR "Olink Explore" OR "Olink HT" OR "Olink Target" OR "proximity extension assay" OR ',
    '"NULISA" OR "NULISAseq" OR "Proteograph" OR "Orbitrap Astral" OR ',
    '"DIA-MS" OR "SWATH-MS" OR "OpenSWATH" OR "tandem mass tag" OR "LC-MS/MS" OR "mass spectrometry proteomics")'
  )
  meth_gen <- paste0(
    '("pQTL" OR "pQTLs" OR "cis-pQTL" OR "cis-pQTLs" OR "sd-pQTL" OR "sd-pQTLs" OR ',
    '"SV-pQTL" OR "SV-pQTLs" OR "TR-pQTL" OR "TR-pQTLs" OR ',
    '"protein QTL" OR "protein QTLs" OR "vQTL" OR "vQTLs" OR "cis-vQTL" OR "cis-vQTLs" OR ',
    '"REGENIE" OR "BOLT-LMM" OR "fastGWA" OR "SAIGE" OR "MatrixEQTL" OR "tensorQTL" OR "QTLtools" OR ',
    '"GWAF" OR "seqMeta" OR "SNPTEST" OR "ProbABEL" OR "GEMMA" OR "LIMIX")'
  )
  res_gen <- paste0(
    '("pQTL" OR "pQTLs" OR "cis-pQTL" OR "cis-pQTLs" OR "trans-pQTL" OR "trans-pQTLs" OR ',
    '"sd-pQTL" OR "sd-pQTLs" OR "SV-pQTL" OR "SV-pQTLs" OR "TR-pQTL" OR "TR-pQTLs" OR ',
    '"protein QTL" OR "protein QTLs" OR "cis-vQTL" OR "cis-vQTLs")'
  )

  # Include (SRC:MED AND IN_PMC:Y) so MEDLINE-indexed open-access PMC full texts are searched alongside SRC:PMC and SRC:PPR
  core <- paste0(
    "(((SRC:MED AND IN_PMC:Y) OR SRC:PMC OR SRC:PPR) AND METHODS:", meth_plat,
    " AND (METHODS:", meth_gen, " OR RESULTS:", res_gen, ") ",
    EPMC_PUBTYPE_EXCLUSION, ")"
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

  core <- paste0(
    "((", cites_clause, " AND (TITLE_ABS:", b_prot, ") AND (TITLE_ABS:", b_gen, ")) ",
    EPMC_PUBTYPE_EXCLUSION, ")"
  )

  if (!is.null(verify_pmids) && length(verify_pmids) > 0L) {
    ext_clause <- paste0("(", paste(sprintf("EXT_ID:%s", verify_pmids), collapse = " OR "), ")")
    return(paste0(core, " AND ", ext_clause))
  }

  dt <- build_epmc_date_clause(since_date, until_date, include_idate = include_idate)
  if (nzchar(dt)) paste0(core, " AND ", dt) else core
}

# ---------------------------------------------------------------------------
# 4. HTML Stripping, ID Normalizers, Preprint Helpers & Fuzzy Title Guards
# ---------------------------------------------------------------------------

# Merge Safeguard 1: Strip real HTML/XML tags FIRST via tag-specific regex, and decode
# &lt;/&gt;/&amp; entities ONLY afterward; NEVER use `<[^>]+>` (which deletes text between
# `<` and `>` inequalities such as `P < 5e-8 ... >5,000 proteins`).
strip_html_tags <- function(x) {
  if (is.null(x) || length(x) == 0L) return(character(0))
  y <- as.character(x)
  y[is.na(y)] <- ""
  y <- gsub("</?[a-zA-Z][a-zA-Z0-9:_-]*(?:\\s+[^<>]*)?>", " ", y, perl = TRUE)
  y <- gsub("&nbsp;|&#160;", " ", y, perl = TRUE)
  y <- gsub("&lt;|&#60;|&#x3[cC];", "<", y, perl = TRUE)
  y <- gsub("&gt;|&#62;|&#x3[eE];", ">", y, perl = TRUE)
  y <- gsub("&amp;|&#38;", "&", y, perl = TRUE)
  y <- gsub("&quot;|&#34;", "\"", y, perl = TRUE)
  y <- gsub("&apos;|&#39;", "'", y, perl = TRUE)
  gsub("\\s+", " ", trimws(y))
}

ascii_fold <- function(x) {
  if (is.null(x) || length(x) == 0L) return(character(0))
  x_chr <- as.character(x)
  x_chr[is.na(x_chr)] <- ""
  special_from <- c("ø", "Ø", "ð", "Ð", "þ", "Þ", "ß", "æ", "Æ", "ł", "Ł", "đ", "Đ",
                    "–", "—", "−", "≈", "∼", "≥", "≤")
  special_to   <- c("o", "O", "d", "D", "th", "TH", "ss", "ae", "AE", "l", "L", "d", "D",
                    "-", "-", "-", "~", "~", ">=", "<=")
  for (k in seq_along(special_from)) {
    x_chr <- gsub(special_from[k], special_to[k], x_chr, fixed = TRUE)
  }
  m <- gregexpr("[^[:ascii:]]+", x_chr, perl = TRUE)
  regmatches(x_chr, m) <- lapply(regmatches(x_chr, m), function(spans) {
    if (length(spans) == 0L) return(spans)
    tr <- iconv(spans, from = "UTF-8", to = "ASCII//TRANSLIT")
    tr[is.na(tr)] <- ""
    gsub("[^A-Za-z0-9]", "", tr)
  })
  x_chr
}

PREPRINT_DOI_PREFIX_REGEX <- "^10\\.(1101|64898|21203|26434|48550|20944|22541|31219|31234|32388|2139|17605)/"

is_preprint_doi <- function(doi) {
  if (is.null(doi) || length(doi) == 0L) return(logical(0))
  d <- tolower(trimws(as.character(doi)))
  d[is.na(d)] <- ""
  grepl(PREPRINT_DOI_PREFIX_REGEX, d, perl = TRUE)
}

extract_preprint_doi_date <- function(d) {
  if (is.null(d) || length(d) == 0L || is.na(d[1L])) return("")
  s <- tolower(trimws(as.character(d[1L])))
  # bioRxiv / medRxiv / openRxiv date-stamped DOIs: 10.1101/YYYY.MM.DD.XXXXXX or 10.64898/YYYY.MM.DD.XXXXXX
  m <- regmatches(s, regexpr("^10\\.(1101|64898)/([0-9]{4})\\.([0-9]{2})\\.([0-9]{2})\\.", s, perl = TRUE))
  if (length(m) == 1L && nzchar(m)) {
    ymd <- sub("^10\\.(1101|64898)/([0-9]{4})\\.([0-9]{2})\\.([0-9]{2})\\.$", "\\2-\\3-\\4", m, perl = TRUE)
    mo <- as.integer(substr(ymd, 6L, 7L))
    dy <- as.integer(substr(ymd, 9L, 10L))
    if (!is.na(mo) && mo >= 1L && mo <= 12L && !is.na(dy) && dy >= 1L && dy <= 31L) {
      return(ymd)
    }
  }
  ""
}

# Safeguards 5 & 7: Preserve `/v[0-9]+` on Research Square (`10.21203/`) and Authorea (`10.22541/`),
# and `.v[0-9]+` on Preprints.org (`10.20944/`), in canonical output `doi`/`identifier` so Crossref DOIs
# resolve (HTTP 200), while stripping version suffixes across ALL 12 preprint prefixes inside `doi_dedup_key()`.
normalize_doi_one <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[1L])) return("")
  d <- trimws(as.character(x[1L]))
  if (!nzchar(d)) return("")
  d <- tryCatch(utils::URLdecode(d), error = function(e) d)
  d <- tolower(trimws(d))
  d <- sub("^https?://(dx\\.|www\\.)?(doi\\.org|biorxiv\\.org/content|medrxiv\\.org/content)/", "", d)
  d <- sub("^doi\\s*:\\s*", "", d)
  repeat {
    prev <- d
    d <- sub("[.,;:)>\\]\"']+$", "", d, perl = TRUE)
    d <- sub("(\\.full(\\.pdf)?|\\.pdf|\\.abstract|\\.short|/full|/abstract)$", "", d, perl = TRUE)
    if (grepl("^10\\.(1101|64898|26434|48550)/", d, perl = TRUE)) {
      d <- sub("(/v[0-9]+|\\.v[0-9]+|v[0-9]+)$", "", d, perl = TRUE)
    }
    if (identical(d, prev)) break
  }
  if (!grepl("^10\\.[0-9]{4,9}/\\S+$", d, perl = TRUE)) return("")
  d
}

doi_dedup_key <- function(x) {
  d <- normalize_doi_one(x)
  if (!nzchar(d)) return("")
  if (grepl(PREPRINT_DOI_PREFIX_REGEX, d, perl = TRUE)) {
    d <- sub("(/v[0-9]+|\\.v[0-9]+|v[0-9]+)$", "", d, perl = TRUE)
  }
  d
}

extract_dois_from_cell <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[1L])) return(character(0))
  raw <- tryCatch(utils::URLdecode(as.character(x[1L])), error = function(e) as.character(x[1L]))
  if (!nzchar(trimws(raw))) return(character(0))
  m <- regmatches(raw, gregexpr("10\\.[0-9]{4,9}/[^\\s;,|\"')\\]>]+", raw, perl = TRUE))[[1L]]
  if (length(m) == 0L) return(character(0))
  normed <- vapply(m, normalize_doi_one, character(1), USE.NAMES = FALSE)
  unique(normed[nzchar(normed)])
}

extract_pmids_from_cell <- function(x, strip_dois_first = TRUE) {
  if (is.null(x) || length(x) == 0L || is.na(x[1L])) return(character(0))
  raw <- as.character(x[1L])
  if (!nzchar(trimws(raw))) return(character(0))
  if (isTRUE(strip_dois_first)) {
    raw <- gsub("10\\.[0-9]{4,9}/[^\\s;,\"')\\]>]+", " ", raw, perl = TRUE)
    raw <- gsub("PMC[0-9]+", " ", raw, ignore.case = TRUE, perl = TRUE)
  }
  m <- regmatches(raw, gregexpr("\\b[1-9][0-9]{5,8}\\b", raw, perl = TRUE))[[1L]]
  unique(m[nzchar(m)])
}

normalize_pmid_one <- function(x) {
  ids <- extract_pmids_from_cell(x, strip_dois_first = TRUE)
  if (length(ids) == 0L) "" else ids[1L]
}

normalize_pmcid_one <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[1L])) return("")
  m <- regmatches(toupper(as.character(x[1L])), regexpr("PMC[0-9]+", toupper(as.character(x[1L])), perl = TRUE))
  if (length(m) == 1L && nzchar(m)) m else ""
}

normalize_author_surname <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[1L])) return("")
  a <- tolower(trimws(ascii_fold(strip_html_tags(x[1L]))))
  if (!nzchar(a)) return("")
  a <- sub("\\bet\\s+al\\.?.*$", "", a)
  first_part <- trimws(strsplit(a, "[;,]")[[1L]][1L])
  tokens <- unlist(strsplit(gsub("[^a-z0-9 -]", " ", first_part), "\\s+"))
  tokens <- tokens[nzchar(tokens)]
  if (length(tokens) == 0L) return("")
  non_init <- tokens[nchar(gsub("-", "", tokens)) >= 2L]
  if (length(non_init) == 0L) return(tokens[1L])
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

# Merge Safeguard 4: Strip 2SMR / target-discovery boilerplate tokens before computing
# title Jaccard similarity and enforce phenotype compatibility on fuzzy merges.
MR_BOILERPLATE_STOPWORDS <- c(
  STOPWORDS_EN,
  "mendelian", "randomization", "randomisation", "two", "sample", "twosample",
  "potential", "therapeutic", "druggable", "drug", "targets", "target",
  "causal", "causally", "systematic", "wide", "colocalization", "colocalisation",
  "summary", "based"
)

PHENOTYPE_GENERIC_STOPWORDS <- c(
  MR_BOILERPLATE_STOPWORDS,
  "proteome", "proteomes", "proteomic", "proteomics", "circulating", "plasma",
  "serum", "proteins", "protein", "genome", "genomic", "genetics", "genetic",
  "identification", "identifying", "associations", "association", "associated",
  "insights", "multi", "omics", "integrating", "integration", "effects", "role",
  "risk", "levels"
)

title_content_tokens <- function(title, strip_mr_boilerplate = TRUE) {
  t_clean <- tolower(ascii_fold(strip_html_tags(title %||% "")))
  t_clean <- gsub("\\(preprint\\)", " ", t_clean)
  toks <- unlist(strsplit(gsub("[^a-z0-9]+", " ", t_clean), "\\s+"))
  sw <- if (isTRUE(strip_mr_boilerplate)) MR_BOILERPLATE_STOPWORDS else STOPWORDS_EN
  unique(toks[nchar(toks) >= 3L & !(toks %in% sw)])
}

specific_phenotype_tokens <- function(title) {
  t_clean <- tolower(ascii_fold(strip_html_tags(title %||% "")))
  toks <- unlist(strsplit(gsub("[^a-z0-9]+", " ", t_clean), "\\s+"))
  unique(toks[nchar(toks) >= 3L & !(toks %in% PHENOTYPE_GENERIC_STOPWORDS)])
}

has_compatible_phenotype <- function(t1, t2) {
  p1 <- specific_phenotype_tokens(t1)
  p2 <- specific_phenotype_tokens(t2)
  if (length(p1) == 0L || length(p2) == 0L) return(TRUE)
  if (length(intersect(p1, p2)) >= 1L) return(TRUE)
  s1 <- paste(p1, collapse = " ")
  s2 <- paste(p2, collapse = " ")
  if (grepl("obesity|bmi|mass", s1) && grepl("obesity|bmi|mass", s2)) return(TRUE)
  FALSE
}

title_jaccard <- function(t1, t2) {
  w1 <- title_content_tokens(t1, strip_mr_boilerplate = TRUE)
  w2 <- title_content_tokens(t2, strip_mr_boilerplate = TRUE)
  if (length(w1) == 0L || length(w2) == 0L) return(0)
  length(intersect(w1, w2)) / length(union(w1, w2))
}

build_title_fingerprint <- function(first_author, year, title) {
  sur <- normalize_author_surname(first_author)
  yr  <- gsub("[^0-9]", "", as.character(year %||% ""))
  t_alnum <- substr(gsub("[^a-z0-9]", "", tolower(ascii_fold(strip_html_tags(title %||% "")))), 1L, 65L)
  if (!nzchar(t_alnum)) return("")
  paste(sur, yr, t_alnum, sep = "|")
}

# EFO-independent EBI GWAS Catalog `disease_trait` protein verifier (rescues uncurated
# summary-statistics protein deposits such as Loya 2025 & Yan 2025 while excluding lipids/CRP/HbA1c)
is_gwas_catalog_protein_trait <- function(trait_vec) {
  tr <- tolower(ascii_fold(strip_html_tags(trait_vec %||% character(0))))
  if (length(tr) == 0L) return(logical(0))
  grepl(paste0(
    "\\b(protein|proteins|proteom[a-z]*|cytokin[a-z]*|chemokin[a-z]*|interleukin[a-z0-9-]*|",
    "growth factors?|complement|aptamers?|somamers?|olink|somascan|nulisa[a-z]*|peptid[a-z]*)\\b"
  ), tr, perl = TRUE) &
    !grepl(paste0(
      "\\b(lipoprotein|apolipoprotein a|apolipoprotein b|c-reactive protein|serum albumin|",
      "total protein|urinary albumin|albumin-to-creatinine|glycoprotein acetyls|n-glycan|",
      "n-glycosyl|glycated hemoglobin|hemoglobin a1c)\\b"
    ), tr, perl = TRUE)
}

# ---------------------------------------------------------------------------
# 5. Platform Detector & Sample-Size Extractor (v4.0 — Precision Hardened)
# ---------------------------------------------------------------------------

detect_platform_v2 <- function(text) {
  raw_txt <- ascii_fold(strip_html_tags(text %||% ""))
  txt <- tolower(raw_txt)
  plats <- character(0)

  if (grepl("somascan|somalogic|\\bsomamers?\\b", txt, perl = TRUE) ||
      (grepl("\\baptamers?\\b|\\baptamer-based\\b", txt, perl = TRUE) &&
       !grepl("\\b(rna|dna)\\s+aptamers?|aptamer-based\\s+(rna|dna|sensors?|switches|biosensors?|therapeutics?)", txt, perl = TRUE))) {
    plats <- c(plats, "SomaScan")
  }
  # Case-sensitive SCALLOP consortium disambiguation vs marine bivalve ("bay scallop" / Argopecten)
  if (grepl("\\bolink\\b|proximity extension|pharma proteomics project|\\bukb-ppp\\b|\\bscallop\\s+(consortium|cohort|meta-analysis|framework|collaboration|collaborative|study|project)\\b", txt, perl = TRUE) ||
      (grepl("\\bSCALLOP\\b", raw_txt, perl = TRUE) &&
       !grepl("\\b(bay scallop|argopecten|patinopecten|pecten|chlamys|bivalve)s?\\b", txt, perl = TRUE))) {
    plats <- c(plats, "Olink")
  }
  if (grepl("mass spectrometry|\\bdia-ms\\b|\\bdia-pasef\\b|\\bswath(-ms)?\\b|\\bopenswath[0-9]*\\b|data-independent acquisition|\\btandem mass tags?\\b|\\btmt-(based|labeled|labelled|proteom|ms|[0-9]{1,2})\\b|\\bitraq\\b|\\bsilac\\b|lc-ms/ms|\\bmrm-ms\\b|\\bprm-ms\\b|\\bsrm-ms\\b|multiple reaction monitoring|parallel reaction monitoring|selected reaction monitoring|proteograph|orbitrap|timstof|shotgun proteom", txt, perl = TRUE)) {
    plats <- c(plats, "Mass Spectrometry")
  }
  if (grepl("\\bnulisa(seq)?\\b|\\balamar\\b", txt, perl = TRUE)) {
    plats <- c(plats, "NULISA")
  }
  if (grepl("luminex|\\bx-?map\\b|bio-plex|milliplex|procartaplex", txt, perl = TRUE)) {
    plats <- c(plats, "Luminex")
  }
  if (grepl("rules[- ]based medicine|myriad rbm|meso scale|electrochemiluminescence|multiplex(ed)? immunoassay|multiplex(ed)? antibody assays?|suspension bead array|protein microarray|micro-western|reverse[- ]phase protein array|\\brppa\\b|\\bsimoa\\b|quanterix|antibody-based (assay|microarray|technology|proteom)", txt, perl = TRUE)) {
    plats <- c(plats, "Other Immunoassay")
  }
  if (length(plats) == 0L) "Unspecified Proteomics" else paste(unique(plats), collapse = " / ")
}

split_sentences_safe <- function(text) {
  txt <- gsub("\\s+", " ", trimws(strip_html_tags(text %||% "")))
  if (!nzchar(txt)) return(character(0))
  unlist(strsplit(txt, "(?<=[.!?])\\s+|;\\s+", perl = TRUE))
}

# Merge Safeguard 5: Exempt biological specimen/cohort collection sentences ("GTM samples
# obtained from 88 AD patients", "dataset derived from a cohort of 284 STARNET patients")
# BEFORE testing summary-statistics / public-database reuse!
is_public_reuse_sentence <- function(sent_lower) {
  if (grepl("\\b(samples?|specimens?|tissues?|biops(y|ies)|blood|plasma|serum|csf|cerebrospinal\\s+fluid|urine|saliva|brains?|proteomes?|dataset|cohort)\\s+(?:[a-z0-9,()-]+\\s+){0,3}(were\\s+|was\\s+)?(obtained|collected|sampled|extracted|harvested|drawn|derived)\\s+from\\s+(?:a\\s+cohort\\s+of\\s+)?(?:[>~]?[1-9][0-9,]*\\s+)?(?:[a-z-]+\\s+){0,3}(patients?|participants?|individuals?|donors?|twins?|adults?|children|controls?|cases?|subjects?)\\b", sent_lower, perl = TRUE) &&
      !grepl("\\b(ukb-?ppp|decode|fenland|finngen|gwas\\s+catalog|ieu\\s+open|opengwas|summary\\s+statistics|summary-level|publicly\\s+available)\\b", sent_lower, perl = TRUE)) {
    return(FALSE)
  }
  if (grepl("\\bwe\\s+(describe|conducted|performed|carried\\s+out|report|present|measured|quantified|profiled|assayed|mapped)\\b", sent_lower, perl = TRUE) &&
      !grepl("\\b(summary(-level)?\\s+(statistics|data)|publicly\\s+available|published\\s+(gwas|pqtl|proteomic)|gwas\\s+catalog|ieu\\s+open|opengwas|finngen|data\\s+(derived|obtained|extracted)\\s+from)\\b", sent_lower, perl = TRUE)) {
    return(FALSE)
  }
  grepl(paste0(
    "\\b(summary(-level)?\\s+(statistics|data)|gwas\\s+(summary\\s+)?(data|statistics|results)|",
    "(?:pqtl|cis-?pqtl|protein\\s+quantitative\\s+trait\\s+loci(?:\\s*\\([a-z-]+\\))?)\\s+(summary\\s+)?(data|statistics|datasets?)|",
    "genetic\\s+(instruments?|associations?)|instrumental\\s+variables?|exposure\\s+data|outcome\\s+data)\\s+",
    "(?:[a-z0-9,()-]+\\s+){0,5}(were\\s+|was\\s+)?(obtained|extracted|derived|downloaded|retrieved|sourced|taken|selected|leveraged|collated|compiled|collected)\\s+from\\b|",
    "\\b(datasets?|data)\\s+(were\\s+|was\\s+)?(obtained|extracted|derived|downloaded|retrieved|sourced)\\s+from\\s+([a-z0-9, -]{0,40})?(genome-wide\\s+association|gwas|pqtl|protein\\s+quantitative\\s+trait)\\b|",
    "\\bdata\\s+(derived|obtained|extracted)\\s+from\\s+[0-9,]+\\s+[a-z -]*proteins?\\s+in\\s+a\\s+(pqtl|protein\\s+quantitative)\\b|",
    "\\b(used|using|leveraged|leveraging|utilized|utilizing)\\s+(the\\s+)?(latest\\s+|recent\\s+|publicly\\s+available\\s+|published\\s+|large-scale\\s+)*summary(-level)?\\s+(statistics|data)\\s+(for|of|from|on)\\b|",
    "\\b(used|using|leveraged|leveraging|utilized|utilizing)\\s+(the\\s+)?(latest\\s+|recent\\s+|publicly\\s+available\\s+|published\\s+|large-scale\\s+)*(plasma\\s+)?proteomic\\s+data\\s+from\\s+(the\\s+)?(ukb-?ppp|uk\\s+biobank\\s+pharma|decode|fenland|interval|scallop|aric)\\b|",
    "\\b(obtained|extracted|derived|downloaded|retrieved|sourced|taken|selected|collated|compiled)\\s+from\\s+",
    "(the\\s+|a\\s+|two\\s+|three\\s+|four\\s+|multiple\\s+|publicly\\s+available\\s+|published\\s+|previous\\s+)?",
    "(ukb-?ppp|uk\\s+biobank\\s+pharma|decode|fenland|interval|scallop|finngen|gwas\\s+catalog|ieu\\s+open\\s*gwas|opengwas|ferkingstad|sun\\s+et\\s+al|large(-scale)?\\s+(pqtl|proteom[a-z]*|gwas))\\b|",
    "\\b(integrat(es?|ed|ing)|leverag(es?|ed|ing)|utiliz(es?|ed|ing))\\s+(two|three|four|multiple|summary)\\s+",
    "(major\\s+|comprehensive\\s+|genome-wide\\s+|large-scale\\s+)*(pqtl|proteom[a-z]*)\\s+(databases|datasets|studies|resources|reference\\s+data)\\b|",
    "\\bfrom\\s+(two|three|four|five|six|seven|eight|nine|ten|[2-9]|10)\\s+",
    "(recent\\s+|published\\s+|large(-scale)?\\s+|independent\\s+|comprehensive\\s+)*(adult\\s+|proteomic\\s+|plasma\\s+proteome\\s+)*",
    "(gwas|genome-wide\\s+association\\s+studies|pqtl\\s+studies|pqtl\\s+datasets|databases|proteomic\\s+studies)\\b|",
    "\\b(cis-?)?pqtls?(\\s+signals?)?\\s+from\\s+(35,559|54,219|54,306|10,708|ukb-?ppp|decode)\\b"
  ), sent_lower, perl = TRUE)
}

# Safeguard 4: Coordinated-count SI thin/space-thousands normalizer (`normalize_si_thousands()`),
# converting both single (`>45 000 participants`, `33 325 individuals`, `1 251 123 SNPs`) and
# coordinated lists (`9 464 and 7 832 pQTLs`) before biomedical count nouns without corrupting
# grant, award, accession, table, or figure identifiers.
SI_BIOMED_COUNT_NOUN_REGEX <- paste0(
  "(?:participants?|individuals?|subjects?|persons?|people|donors?|volunteers?|",
  "icelanders?|europeans?|british|scottish|finnish|finns|swedish|danish|greenlandic|inuit|japanese|chinese|koreans?|indians?|africans?|hispanics?|arabs?|",
  "adults?|females?|males?|women|men|children|adolescents?|infants?|newborns?|babies|neonates?|twins?|centenarians?|patients?|cases?|controls?|",
  "samples?|specimens?|brains?|proteomes?|cohorts?|",
  "proteins?|aptamers?|somamers?|peptides?|cytokines?|chemokines?|biomarkers?|assays?|targets?|genes?|",
  "pqtls?|cis-?pqtls?|trans-?pqtls?|sd-?pqtls?|sv-?pqtls?|tr-?pqtls?|vqtls?|cis-?vqtls?|qqtls?|iqtl[a-z]*|sv-?qtls?|tr-?xqtls?|molqtls?|xqtls?|eqtls?|sqtls?|mqtls?|qtls?|",
  "loci|locus|quantitative\\s+trait\\s+loc[a-z]*|snps?|single-?nucleotide|variants?|indels?|tandem\\s+repeats?|",
  "associations?|signals?|traits?|phenotypes?|outcomes?|diseases?|disorders?|conditions?|pathways?|networks?|subnetworks?|credible\\s+sets?)"
)

normalize_si_thousands <- function(s) {
  if (is.null(s) || length(s) == 0L) return(character(0))
  y <- gsub("\\s+", " ", as.character(s))
  pat <- paste0(
    "\\b(?:grants?|awards?|contracts?|projects?|applications?|codes?|ids?|tables?|figures?|figs?|supplements?|supplementary|appendi(?:x|ces)|no\\.?|number)\\s+[1-9][0-9]{0,2}(?:\\s+[0-9]{3})+\\b(*SKIP)(*FAIL)|",
    "(?<![0-9.,-])\\b([1-9][0-9]{0,2}(?:,[0-9]{3})*)\\s+([0-9]{3})",
    "(?=(?:\\s+[0-9]{3})*(?:\\s*(?:,\\s*(?:and\\s+)?|\\s+(?:and|or|to|vs\\.?|versus)\\s+)[>~]?[1-9][0-9]{0,2}(?:[\\s,][0-9]{3})+)*",
    "\\s+(?:[a-z()/+-]+\\s+){0,4}", SI_BIOMED_COUNT_NOUN_REGEX, "\\b)"
  )
  for (iter in 1:3) {
    prev <- y
    y <- gsub(pat, "\\1,\\2", y, perl = TRUE, ignore.case = TRUE)
    if (identical(y, prev)) break
  }
  y
}

extract_sample_size_v2 <- function(text) {
  txt <- normalize_si_thousands(gsub("\\s+", " ", tolower(ascii_fold(strip_html_tags(text %||% "")))))
  sents <- split_sentences_safe(txt)
  if (length(sents) == 0L) return("Check full text")

  keep_sents <- sents[!vapply(sents, is_public_reuse_sentence, logical(1)) &
                      !grepl("\\b(biobank|cohort|study)\\s+(is\\s+a\\s+|of\\s+|enrolled\\s+|recruited\\s+)?(approximately|~|over|more\\s+than)?\\s*[0-9,]{6,}", sents, perl = TRUE)]
  if (length(keep_sents) == 0L) keep_sents <- sents
  clean_txt <- paste(keep_sents, collapse = " . ")
  clean_txt <- gsub("\\b[0-9,]+\\s+(traits?|outcomes?|diseases?|conditions?|phenotypes?|disorders?|cancers?|loci|variants?|snps?|genes?|proteins?)\\s+(in|across|among|of)\\b",
                    " traits_masked ", clean_txt, perl = TRUE)

  m1_raw <- regmatches(clean_txt, gregexpr("\\bn\\s*=\\s*[>~]?\\s*[1-9][0-9,]*(?:\\s+[a-z-]+){0,2}", clean_txt, perl = TRUE))[[1L]]
  m1_keep <- m1_raw[!grepl("\\b(proteins?|protein\\s+groups?|aptamers?|somamers?|peptides?|traits?|loci|snps?|variants?|genes?|diseases?|conditions?|cancers?|members?)\\b", m1_raw, perl = TRUE)]
  nums1 <- suppressWarnings(as.numeric(gsub("[^0-9]", "", sub("^(\\bn\\s*=\\s*[>~]?\\s*[1-9][0-9,]*).*$", "\\1", m1_keep, perl = TRUE))))

  m2 <- regmatches(clean_txt, gregexpr(paste0(
    "(?:[>~]|\\b(?:over|more\\s+than|nearly|almost|approximately|up\\s+to)\\s+)?",
    "\\b([1-9][0-9,]{1,6})\\s+(?:[a-z()/+-]+\\s+){0,3}",
    "(?:participants|individuals|subjects|persons|people|icelanders|europeans|british|scottish|finnish|finns|danish|greenlandic|japanese|chinese|",
    "adults|females|males|women|men|children|adolescents|newborns|babies|infants|neonates|twins|patients|donors|volunteers|centenarians|",
    "brains|brain\\s+proteomes|gtex\\s+samples|gtm\\s+samples|plasma\\s+samples|serum\\s+samples|csf\\s+samples|kidney\\s+samples|rpe\\s+lines|samples)\\b"
  ), clean_txt, perl = TRUE))[[1L]]
  nums2 <- suppressWarnings(as.numeric(gsub(",", "", gsub("^.*?\\b([1-9][0-9,]{1,6})\\b.*$", "\\1", m2, perl = TRUE))))

  all_n <- c(nums1, nums2)
  all_n <- all_n[!is.na(all_n) & all_n >= 10 & all_n <= 300000]
  if (length(all_n) == 0L) return("Check full text")
  format(max(all_n), big.mark = ",", scientific = FALSE, trim = TRUE)
}

# ---------------------------------------------------------------------------
# 6. Unified Two-Stage Dual-Axis Triage Scorer v4.0
# ---------------------------------------------------------------------------

NONHUMAN_ORGANISM_REGEX <- paste0(
  "\\b(mice|mouse|murine|rat|rats|rodents?|drosophila|zebrafish|danio|xenopus|caenorhabditis|c\\.\\s*elegans|",
  "arabidopsis|maize|zea\\s+mays|\\brice\\b|oryza|\\bwheat\\b|triticum|soybean|glycine\\s+max|barley|hordeum|",
  "sorghum|buckwheat|fagopyrum|chenopodium|quinoa|triticale|sunflower|grape|grapes|cannabis|jujube|\\boat\\b|avena|",
  "pea|pisum\\s+sativum|garden\\s+pea|potato|sugarcane|cotton|tomato|brassica|populus|mungbean|vigna|chickpea|cicer|",
  "capsicum|dogwood|cornus|hazelnut|allopolyploid\\s+crops|plant\\s+breeding|crop\\s+improvement|plant|plants|crop|crops|cultivars?|inbred\\s+lines|seed\\s+germination|phytophthora|",
  "bovine|cattle|\\bcows?\\b|holsteins?|nellore|calves|porcine|\\bpigs?\\b|piglets?|swine|boar|nanopig|minipig|",
  "\\bsheep\\b|ovine|ewes?|lambs?|\\bgoats?\\b|caprine|chickens?|poultry|gallus|leghorn|broilers?|hens?|ducks?|mule\\s+ducks?|geese|goose|",
  "canine|\\bdogs?\\b|spaniel|feline|\\bcats?\\b|equine|horses?|donkeys?|camels?|macaques?|rhesus|dolphins?|tursiops|",
  "catfish|grass\\s+carp|ctenopharyngodon|stickleback|tilapia|salmon|trout|oysters?|crassostrea|",
  "bivalves?|argopecten|bay\\s+scallops?|mollusks?|meretrix|mussels?|mytilus|octocorals?|corals?|",
  "vipers?|vipera|snakes?|venom|spider\\s+mites?|tetranychus|silkworms?|bombyx|crickets?|gryllus|orthoptera|mosquito(?:es)?|aedes|anopheles|culex|insects?|tardigrades?|hypsibius|",
  "haemonchus|schistosoma|nematodes?|leishmania|phlebotomus|plasmodium|toxoplasma|trypanosom[a-z]*|oxytricha|",
  "yeast|saccharomyces|schizosaccharomyces|fission\\s+yeast|pseudozyma|meyerozyma|mushrooms?|lepista|shiitake|aspergillus|fusarium|alternaria|fungi|fungal|",
  "cyanobacteri[a-z]*|synechocystis|chlamydomonas|microalga[a-z]*|actinobacteria|mycobacterium|klebsiella|",
  "staphylococcus|pseudomonas|salmonella|vibrio|streptomyces|streptococcus|bordetella|",
  "burkholderia|paraburkholderia|parageobacillus|bacillus\\s+subtilis|paracoccus|desulfurivibrio)\\b"
)

score_triage_v2 <- function(title,
                            abstract_text = "",
                            pub_types = character(0),
                            mesh_headings = character(0),
                            gwas_catalog_study_count = 0L,
                            erratum_targets_catalogued = FALSE,
                            ft_verified_pqtl = FALSE,
                            doi = "") {
  norm_ws_hyph <- function(s) {
    # Safeguard 1: strip_html_tags() runs BEFORE ascii_fold()
    y <- tolower(trimws(ascii_fold(strip_html_tags(s %||% ""))))
    y <- gsub("\\b(cis|trans|sd|mol|pan|v|sv|tr|q)\\s*-\\s*(p|v|e|m|s|x|mol)?(qtls?|protein)\\b", "\\1-\\2\\3", y, perl = TRUE)
    y <- gsub("\\b(cis|trans)\\s*-\\s*(acting)\\b", "\\1-\\2", y, perl = TRUE)
    normalize_si_thousands(gsub("\\s+", " ", y))
  }

  ttl <- norm_ws_hyph(title)
  abs <- norm_ws_hyph(abstract_text)
  txt <- trimws(paste(ttl, abs))
  pts <- tolower(trimws(pub_types %||% character(0)))
  msh <- tolower(trimws(mesh_headings %||% character(0)))

  # -------------------------------------------------------------------------
  # STAGE 1: Hard Structural & Domain-Relevance Gates
  # -------------------------------------------------------------------------

  # G1a: Erratum / Correction / Corrigendum routing
  is_erratum <- any(pts %in% c("published erratum", "erratum", "correction")) ||
    grepl("^(author correction|publisher correction|correction|erratum|corrigendum)(\\s*:|\\s+to\\b)", ttl, perl = TRUE)
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

  # G1b: Editorial / Comment / News / Retraction / Conference Supplement / Book Chapter
  # Merge Safeguard 6: NEVER ban "letter" in G1b (Nature Letters & Research Letters carry pt="Letter")!
  doi_lc <- tolower(trimws(doi %||% ""))
  is_isbn_book_chapter <- nzchar(doi_lc) && grepl("^10\\.[0-9]{4,5}/978-[0-9-]+_[0-9]+$", doi_lc, perl = TRUE)
  is_conf_or_chapter <- any(pts %in% c("book chapter", "congress", "meeting abstract")) ||
    is_isbn_book_chapter ||
    (length(pts) > 0L && all(pts %in% c("abstract", "meeting abstract", "congress"))) ||
    (!nzchar(abs) && grepl("^(abstracts? (of|from) the|proceedings of the|posters?$|.*(oral|poster) presentations?$|.*conference abstracts?)", ttl, perl = TRUE)) ||
    grepl("\\bin this chapter, we present\\b|^this is a summary of\\b", abs, perl = TRUE) ||
    grepl("^(commentary|editorial|special issue)\\s*[:\\\"]|^mass spectrometry-based mapping of plasma protein qtls in children and adolescents\\.$", ttl, perl = TRUE)

  if (any(pts %in% c("editorial", "comment", "news", "retracted publication", "biography", "directory", "interview")) ||
      is_conf_or_chapter) {
    return(list(
      stage = "Stage1_HardGate", P_score = 0L, S_score = 0L,
      p_rules = character(0), s_rules = "G1_NonOriginalPubType",
      classification = "Auto-Deprioritised: Non-Original (Editorial/Comment)"
    ))
  }

  # G0: Broad Proteomic/Protein/Peptide/Cytokine Signal Gate (masking Trail Making Test "TMT-A/B")
  # Note: gwas_catalog_study_count >= 5 reflects verified EBI GWAS Catalog protein traits (rescuing Loya 2025),
  # and ft_verified_pqtl reflects Stage B.5 PMC full-text verification on multi-omic xQTL/molQTL studies.
  txt_no_tmt_cog <- gsub("trail making test[^.;]{0,25}\\btmt(-b|-a)?\\b|\\btmt-(a|b)\\b", " cognitive_test ", txt, perl = TRUE)
  has_protein_signal <- grepl(BROAD_PROTEIN_REGEX, txt_no_tmt_cog, perl = TRUE) ||
    any(msh %in% c("proteomics", "proteome", "blood proteins", "cerebrospinal fluid proteins")) ||
    (!is.na(gwas_catalog_study_count) && gwas_catalog_study_count >= 5L) ||
    isTRUE(ft_verified_pqtl)
  if (nzchar(txt) && !has_protein_signal) {
    return(list(
      stage = "Stage1_HardGate", P_score = 0L, S_score = 0L,
      p_rules = character(0), s_rules = "G0_NoProteinSignal",
      classification = "Auto-Deprioritised: Non-Proteomic Study"
    ))
  }

  # G2: Non-Human Organism Gate (Title-First Priority + Human Cohort Sample Guard)
  title_has_nonhuman <- grepl(NONHUMAN_ORGANISM_REGEX, ttl, perl = TRUE)
  title_has_human    <- grepl("\\b(human|humans|cross-species|patients?|children|adults?|participants?|women|men|clinical|cohort|biobank)\\b", ttl, perl = TRUE)
  if (title_has_nonhuman && !title_has_human) {
    return(list(
      stage = "Stage1_HardGate", P_score = 0L, S_score = 0L,
      p_rules = character(0), s_rules = "G2_NonHumanTitle",
      classification = "Auto-Deprioritised: Non-Human Study"
    ))
  }

  has_nonhuman <- grepl(NONHUMAN_ORGANISM_REGEX, txt, perl = TRUE)
  # Note (v4.2): Scope `cohort` to human qualifiers (`human cohort`, `patient cohort`, `clinical cohort`, etc.)
  # so animal-only cohorts ("an independent cohort of collaborative cross mice", Chick 2016 PMID:27309819)
  # do not bypass G2 when no human subject terms or MeSH="Humans" are present.
  has_human <- grepl(paste0(
    "\\b(human|humans|participant|participants|patient|patients|individuals|subjects|donors|volunteers|persons|",
    "cases|controls|icelanders|europeans|british|scottish|finnish|finns|swedish|danish|chinese|japanese|korean|indian|african|hispanic|inuit|greenlandic|arab|",
    "men|women|males|females|adults|children|adolescents|pediatric|neonates|infants|newborns|babies|mothers|offspring|centenarians|families|twins|",
    "biobank|(?:human|patient|clinical|birth|population(?:-based)?|prospective|twin|mother-child)\\s+cohorts?|ancestry|lymphoblastoid|lcl|lcls|ipscs?|hipscs?|induced pluripotent|organoids?|post-?mortem|autops(y|ies)|biops(y|ies)|",
    "hapmap|1000 genomes|kora|framingham|aric|ages|interval|fenland|ukb|ukb-ppp|decode|rosmap|ros/map|adni|biofinder|starnet|gtex|teddy|daisy|nagahama|genes\\s*&\\s*health)\\b"
  ), txt, perl = TRUE) ||
    ("humans" %in% msh && !grepl("\\b(inbred lines|recombinant inbred|strains|cultivars?|breeds?|isolates)\\b", txt, perl = TRUE))

  if (has_nonhuman && !has_human) {
    return(list(
      stage = "Stage1_HardGate", P_score = 0L, S_score = 0L,
      p_rules = character(0), s_rules = "G2_NonHumanOnly",
      classification = "Auto-Deprioritised: Non-Human Study"
    ))
  }

  # G3: Pure Epigenome-Wide / DNA Methylation Study without Germline pGWAS
  has_ewas <- grepl("epigenome-wide|methylome-wide|dna methylation|\\bewas\\b|\\bmwas\\b|cpg sites?|methylation quantitative trait|\\bcontmeqtls?\\b", txt, perl = TRUE)
  has_germline_pgwas <- isTRUE(ft_verified_pqtl) || (grepl(paste0(
    "\\b(pqtls?|cis-?pqtls?|trans-?pqtls?|sd-?pqtls?|sv-?pqtls?|tr-?pqtls?|protein[- ]qtls?|vqtls?|cis-?vqtls?|sv-?qtls?|molqtls?|xqtls?|(?:protein|molecular) quantitative trait|pgwas|",
    "(genome-wide association|genome-?\\s+and\\s+(?:epigenome|methylome|phenome|transcriptome)[ -]wide(?:\\s+association)?|gwas|gwass|qtl mapping|cis-?qtl)[^.;]{0,65}(protein|proteom|aptamer|somamer|cytokin|peptid|biomarker)|",
    "performed\\s+epigenome-wide\\s+association[^.;]{0,65}and\\s+genome-wide\\s+association|",
    "genotyped|whole[ -]genome sequenc[a-z]*|long-read genome sequenc[a-z]*|whole[ -]exome sequenc[a-z]*|exome sequenc[a-z]*|imput(ed|ation))\\b"
  ), txt, perl = TRUE) && !grepl("\\babo blood group status\\b|\\brs8176719\\b", txt, perl = TRUE))
  if (has_ewas && !has_germline_pgwas) {
    return(list(
      stage = "Stage1_HardGate", P_score = 0L, S_score = 0L,
      p_rules = character(0), s_rules = "G3_EWAS_NoGermlinepGWAS",
      classification = "Auto-Deprioritised: Non-GWAS (EWAS)"
    ))
  }

  # -------------------------------------------------------------------------
  # Sentence Partitioning (Primary Cohort Sentences vs Public-Reuse/2SMR Sentences)
  # -------------------------------------------------------------------------
  sents <- split_sentences_safe(txt)
  non_dl_sents <- sents[!vapply(sents, is_public_reuse_sentence, logical(1))]

  primary_ctx_sents <- non_dl_sents[!vapply(non_dl_sents, function(sn) {
    has_sec_kw <- grepl(
      paste0(
        "\\b(two-sample mendelian|2-sample mendelian|two sample mendelian|proteome-wide mendelian|protein-wide mendelian|",
        "multi-omics mendelian|multi-omics mr|proteomic mr|2-step mr|mendelian randomi[sz]ation|",
        "summary-data-based|summary-based mendelian|smr analysis|mr analysis|mr analyses|mr identified|mr \\(pqtl|",
        "proteome-wide association stud(y|ies)|\\bpwas\\b|",
        "for use as (genetic )?instruments|as instrumental variables|instrumental variables?|genetic instruments|",
        "published genome-wide|previous genome-wide|eight large genome-wide|nine genome-wide|ten proteomic|",
        "ukb-?ppp and decode|decode and ukb|genetically predicted ([a-z]+ )?levels of [0-9,]+|",
        "study that detected [0-9,]+ protein quantitative)\\b"
      ),
      sn, perl = TRUE
    )
    if (!has_sec_kw) return(FALSE)
    has_pwas_model_or_pqtl <- grepl(
      "\\b(identif(y|ied|ies) cis-?pqtls?|mapping of cis-?pqtls?|models for proteome-wide association|genetic predictors of the [a-z ]*proteome|built predictive models for protein abundances)\\b",
      sn, perl = TRUE
    )
    !has_pwas_model_or_pqtl
  }, logical(1))]
  primary_ctx_txt <- paste(primary_ctx_sents, collapse = " . ")

  non_mr_yield_sents <- non_dl_sents[!grepl(paste0(
    "\\b(for use as (genetic )?instruments|as instrumental variables|instrumental variables?|",
    "prioritized by ([a-z-]+ )?(mr|mendelian)|(mendelian randomi[sz]ation|\\bmr\\b) (identified|revealed|showed|suggested)|",
    "causally (linked|associated|related)|genetically predicted ([a-z]+ )?levels of [0-9,]+|",
    "study that detected [0-9,]+ protein quantitative)\\b"
  ), non_dl_sents, perl = TRUE)]

  # -------------------------------------------------------------------------
  # STAGE 2: Dual-Axis Point Rubric (P_score vs S_score)
  # -------------------------------------------------------------------------
  P <- 0L
  S <- 0L
  p_rules <- character(0)
  s_rules <- character(0)

  prot_count_pat <- paste0(
    "(?:[>~]\\s*|\\b(?:over|more\\s+than|nearly|almost|approximately|up\\s+to|covering|comprising|of|across|with)\\s+|\\bn\\s*=\\s*)?",
    "\\b(([5-9]|[1-9][0-9]{1,4}(?:,[0-9]{3})*)|five|six|seven|eight|nine|ten|hundreds of|thousands of)\\s+([a-z-]+\\s+){0,4}",
    "(proteins?|protein\\s+levels|protein\\s+groups?|protein\\s+targets?|protein\\s+forms?|protein\\s+isoforms?|aptamers?|somamers?|olink\\s+assays?|peptides?|cytokines?|chemokines?|",
    "complement\\s+(?:system\\s+)?proteins?|inflammatory\\s+(?:proteins?|markers?|biomarkers?)|cardiovascular\\s+(?:proteins?|markers?|biomarkers?)|cardiometabolic\\s+(?:proteins?|markers?|biomarkers?)|neurological\\s+(?:proteins?|biomarkers?)|(?:tumor|tumour)\\s+biomarkers?|protein\\s+biomarkers?|circulating\\s+biomarkers?|systemic\\s+biomarker\\s+concentrations?|analytes?)\\b"
  )
  has_named_five_plus_proteins <- nzchar(primary_ctx_txt) &&
    grepl("\\b(?:gwas|genome-wide\\s+association[a-z ]*)\\s+(?:of|for|on)\\s+(?:measured\\s+|assayed\\s+|quantified\\s+)?(?:serum|plasma|csf|circulating|blood)\\s+(?:protein\\s+)?levels\\s+of\\s+(?:[^.;,]+,\\s*){3,}[^.;,]+(?:,\\s*and|\\s+and)\\s+[^.;]+", primary_ctx_txt, perl = TRUE)

  has_multi_protein_count <- length(primary_ctx_sents) > 0L && (
    any(grepl(prot_count_pat, primary_ctx_sents, perl = TRUE)) ||
      has_named_five_plus_proteins ||
      (grepl("\\b(?:tumor|tumour|protein|circulating|inflammatory|neurological|cardiovascular|cardiometabolic)\\s+biomarkers\\b", ttl, perl = TRUE) &&
         any(grepl("\\b(?:the\\s+)?(?:five|six|seven|eight|nine|ten|[5-9]|[1-9][0-9,]+)\\s+biomarkers\\b", primary_ctx_sents, perl = TRUE)))
  )

  # Single-analyte (< 5 proteins) & somatic/case-report & composite aging-clock guards
  is_single_analyte <- !has_multi_protein_count && (!(!is.na(gwas_catalog_study_count) && gwas_catalog_study_count >= 5L)) && (
    grepl(paste0(
      "\\b(gwas|genome-?wide|genome- and phenome-wide|genomic regions|genetic|genomic|linkage|association) ([a-z-]+ ){0,6}",
      "(of|for|on|with|associated with) (serum |plasma |circulating |soluble )?",
      "(alpha-fetoproteins?|histidine-rich glycoprotein|procalcitonin|vegfr1|svegfr1|soluble vegfr1|factor v\\b|",
      "etonogestrel|abiraterone|voriconazole|tramadol|apol1\\b|klotho|alpha-klotho|sex hormone-binding globulin|shbg|",
      "trem2|plasma trem2|gdf15|igfbp7|25-hydroxyvitamin|folate|c-reactive protein|lipoprotein\\(a\\))\\b"
    ), txt, perl = TRUE) ||
      grepl("\\b(?:gwa[a-z ]*|genome-wide\\s+association[a-z ]*)\\s+(?:for|of|on|using|to\\s+determine\\s+the\\s+genetic\\s+drivers\\s+of)\\s+(?:all\\s+(?:two|three|four)\\s+biomarkers|(?:the\\s+)?[a-z0-9-]+\\s+and\\s+[a-z0-9-]+\\s+protein\\s+levels|(?:proinflammatory\\s+cytokines\\s*)?[\\[(]\\s*(?:il-?6|crp|c-reactive\\s+protein|tnf-?alpha)\\s+and\\s+(?:il-?6|crp|c-reactive\\s+protein|tnf-?alpha)(?:\\s*\\([^)]+\\))?\\s*[\\])])", txt, perl = TRUE) ||
      grepl("^platform-dependent effects of genetic variants on plasma apol1", ttl, perl = TRUE)
  )

  is_somatic_or_case_only <- grepl(paste0(
    "\\b(case report|a novel pediatric case|a child with|a patient with|paired tumor and normal adjacent|",
    "ctdna|minimal residual disease|neoantigen|tumor-agnostic|somatic mutations? only|",
    "orthotopic xenograft|crispr-cas9 editing|clonal hematopoiesis of indeterminate potential)\\b"
  ), txt, perl = TRUE) &&
    !grepl("\\b(germline|pqtls?|cis-?pqtls?|trans-?pqtls?|vqtls?|cis-?vqtls?|pgwas|quantitative trait loc)\\b", txt, perl = TRUE)

  is_composite_clock_only <- grepl("\\b(proteomic ag(e|ing)( acceleration| trait| clock)?|organ ag(e|ing))\\b", txt, perl = TRUE) &&
    !grepl("\\b(we (conducted|performed|mapped|identified) ([a-z-]+ ){0,3}(pqtls?|cis-?pqtls?|trans-?pqtls?|vqtls?)|[1-9][0-9]+ (pqtls?|cis-?pqtls?|vqtls?) were identified)\\b", primary_ctx_txt, perl = TRUE)

  # P1 (+5): Active or passive primary cohort proteomic measurement in primary context sentences
  p1_multiomic_cohort <- grepl("\\bproteom(ics|ic|es?)\\b", primary_ctx_txt, perl = TRUE) && (
    grepl("\\bwe (?:systematically )?profile[ds]? [>~]?[1-9][0-9,]* [a-z0-9-]+ across (?:[a-z0-9-]+ )?(?:omics|molecular) layers\\b", primary_ctx_txt, perl = TRUE) ||
      grepl("\\bperformed [a-z -]*genome sequencing[^.;]{0,220}integrated with matched [a-z, ]*proteomic data\\b", primary_ctx_txt, perl = TRUE)
  )
  p1_hit <- !is_single_analyte && !is_somatic_or_case_only && length(primary_ctx_sents) > 0L && (
    p1_multiomic_cohort || any(
      grepl("\\b(proteins?|proteom[a-z]*|olink|somascan|somalogic|somamers?|aptamers?|nulisa[a-z]*|proteograph|orbitrap|peptides?|cytokines?|chemokines?|biomarkers?)\\b", primary_ctx_sents, perl = TRUE) &
      !grepl("\\b(quantif(y|ied) each protein'?s contribution|quantified [0-9,]+ serum metabolites)\\b", primary_ctx_sents, perl = TRUE) &
      grepl(paste0(
        "\\bwe (measure[ds]?|quantif(y|ied)|profile[ds]?|assay(ed)?|characteri[sz]e[ds]?|generated? proteomic|analy[sz]e[ds]? plasma from [0-9,]+|investigat(e|ed) [0-9,]+ plasma proteins|evaluat(e|ed) (?:plasma |serum |csf |circulating )?levels of [1-9][0-9,]* proteins|compared three leading [a-z -]+proteomic platforms)\\b|",
        "\\b(proteins?|[0-9,]+ (?:[a-z-]+ )?proteins?|protein levels|protein concentrations|protein abundances|proteomic profiles|cytokines?|peptides?|biomarkers?) ([a-z0-9,()-]+ ){0,5}(were|are) (measured|quantified|profiled|assayed|analy[sz]ed by (?:the )?(?:proximity extension|aptamer-based|somascan|olink|mass spectrometry))\\b|",
        "\\bprotein (?:levels|concentrations|abundances) of [1-9][0-9,]* proteins[^.;]{0,75}\\bmeasured (?:via|by|using|with)\\b|",
        "\\bcharacterizing the [a-z -]*proteomic profiles of [0-9,]+\\b|",
        "\\b(proteomic[a-z]*|olink|somascan|somalogic|nulisa[a-z]*|mass spectrometry|dia-ms|swath[a-z-]*|lc-ms/ms|proteograph|multiplex)\\s+",
        "([a-z0-9.-]+\\s+){0,4}(profiling|measurement|quantification|analysis|characteri[sz]ation|data|workflow|assay[a-z]*)\\s+",
        "(was|were|is|are|of\\s+[>~]?[1-9][0-9,]+)\\s*(performed|conducted|carried out|undertaken|applied|generated|obtained|[a-z0-9 -]{0,45}quantified|[>~]?[1-9][0-9,]+)?\\b|",
        "\\b(performing|performed|using|used|applied|employed)\\s+([a-z0-9, -]{0,45})?",
        "(olink|somascan|somalogic|proximity extension assay|nulisa[a-z]*|proteograph|orbitrap astral|swath|dia-ms|tandem mass tag|mass spectrometry(-based)? proteomics|proteomics analysis)[a-z -]{0,20}",
        "(to\\s+characteri[sz]e|to\\s+profile|to\\s+quantify|to\\s+measure|to\\s+assess|of\\s+human|on\\s+[0-9,]+|in\\s+[>~]?[0-9,]+|,?\\s+we\\s+(analyzed|analysed|examined|investigated|studied|profiled|quantified|evaluated))\\b|",
        "\\b(profile[ds]?|assay(ed)?|quantif(y|ied|ication of)|measure[ds]?)\\s+(all\\s+)?([a-z-]+\\s+){0,3}([>~]?([5-9]|[1-9][0-9,]+))\\s+([a-z-]+\\s+){0,3}",
        "(proteins?|protein\\s+groups?|protein\\s+targets?|proteomes?|aptamers?|somamers?|peptides?|cytokines?|chemokines?)\\b|",
        "\\bsystematic quantification of [1-9][0-9,]+ unique proteins\\b|",
        "\\bmeasuring (?:serum|plasma|csf|circulating|blood) levels of [1-9][0-9,]* proteins in [1-9][0-9,]*\\b|",
        "\\bby profiling (?:hundreds|thousands|[>~]?[1-9][0-9,]+) of proteins\\b|",
        "\\bcombined whole-genome sequencing,\\s*multiplex [a-z -]*protein profiling\\b|",
        "\\brecalled for (plasma|serum|blood) sampling for proteomic profiling\\b|",
        "\\bprofiling (urinary|plasma|serum|csf|brain|tissue|salivary) proteomes from [>~]?[1-9][0-9,]+\\b"
      ), primary_ctx_sents, perl = TRUE)
    )
  )
  if (p1_hit) {
    P <- P + 5L; p_rules <- c(p_rules, "P1_DirectAssay(+5)")
  }

  # P2 (+6): Explicit execution of primary pGWAS / pQTL / vQTL / germline-to-proteome mapping on multiplexed proteins
  # Merge Safeguard 5: Tighten title-only P2 so "Genome-Wide Aggregated Trans Effects Analysis" (Zhou 2025) cannot trigger P2!
  p2_hit <- !is_single_analyte && !is_somatic_or_case_only && !is_composite_clock_only && (
    (nzchar(primary_ctx_txt) && grepl(paste0(
      "\\bwe (conduct(ed)?|perform(ed)?|carr(y|ied) out|undertake|undertook|report|present|describe|map(ped)?|construct(ed)?|systematically (investigated|mapped))\\b(?![^.;]{0,55}\\bdata (derived|obtained|extracted) from\\b)[^.;]{0,75}",
      "\\b(pqtls?|cis-?pqtls?|trans-?pqtls?|sd-?pqtls?|sv-?pqtls?|tr-?pqtls?|vqtls?|cis-?vqtls?|protein quantitative trait loc[a-z]*|protein quantitative loc[a-z]*|(?:cis-?|trans-?)?variance quantitative trait loc[a-z]*|pgwas|(?:pqtl|vqtl) (mapping|analys[ei]s|map|stud(y|ies)|scans?)|genome-proteome-wide pqtl)\\b|",
      "\\bwe (identif(y|ied)|map(ped)?|detect(ed)?|discover(ed)?) (?!(?:[a-z-]+ )?(?:proteins?|targets?|genes?|biomarkers?) (?:through|via|using|by)\\b)([a-z-]+ ){0,3}(cis-?|trans-?|sd-?|sv-?|tr-?)?((?:protein|variance) quantitative (trait )?loci|pqtls?|vqtls?)\\b|",
      "\\bperform(ed|ing)? (a |an |comprehensive |sex-stratified )?(protein quantitative trait loc[a-z]*(\\s*\\(pqtls?\\))?|pqtls?|cis-?pqtls?|vqtls?|cis-?vqtls?) (analys[ei]s|mapping)\\b|",
      "\\b(investigat(e|ed)|systematically investigated) ([a-z-]+ ){0,3}genetic (regulation|architecture|determinants) of (the )?(human )?(plasma |serum |csf |brain |urine |circulating )?(proteome|proteins)\\b|",
      "\\bcharacteri[sz](e|ed) the genetic architecture of (the )?(human )?(plasma |serum |csf |brain |circulating )?(proteome|proteins)\\b|",
      "\\b(pqtl|cis-?pqtl|vqtl|cis-?vqtl|protein quantitative trait loc[a-z]*)\\s+(mapping|analys[ei]s|stud(y|ies)|scans?)\\b[^.;]{0,95}\\b(was|were)\\s+(performed|conducted|carried out|undertaken|applied)\\b|",
      "\\b(pqtls?|cis-?pqtls?|vqtls?|cis-?vqtls?)\\s+were\\s+(derived|mapped|identified|calculated)\\s+(by\\s+coupling|in\\s+[0-9,]+|across)\\b|",
      "(?:\\bin\\s+step\\s+(?:one|1),?\\s+)?\\b[1-9][0-9,]*\\s+(?:protein\\s+quantitative\\s+trait\\s+loci(?:\\s*\\(pqtls?\\))?|pqtls?|cis-?pqtls?)\\s+were\\s+(?:identified|discovered|mapped)\\s+and\\s+(?:independently\\s+)?replicated\\b|",
      "\\b(?:generated|constructed|present|provide)\\s+a\\s+(?:detailed\\s+|comprehensive\\s+)?(?:landscape|map|atlas)\\s+of\\s+the\\s+(?:genomic|genetic)\\s+architecture\\s+of\\s+(?:[a-z-]+\\s+){0,2}protein\\s+levels\\b|",
      "\\b(?:gwas|genome-wide\\s+association\\s+stud(?:y|ies))(?:\\s*\\(gwas\\))?[^.;]{0,95}\\brevealed\\s+[1-9][0-9,]*\\s+(?:independent\\s+)?variants?\\s+(?:across|for|associated\\s+with)\\s+[1-9][0-9,]*\\s+proteins\\b|",
      "\\b[1-9][0-9,]*\\s+(?:[a-z-]+\\s+){0,3}exome[- ]array\\s+variants\\s+were\\s+(?:compared|tested|associated)\\s+with\\s+[1-9][0-9,]*\\s+protein\\s+(?:measurements|levels|targets)\\b|",
      "\\bused\\s+admixture[- ]mapping\\s+to\\s+(?:explore|investigate|assess|map|identify)[^.;]{0,75}\\b(?:[1-9][0-9,]*\\s+[a-z -]*biomarkers?|proteome|proteins?)\\b|",
      "\\bintegrated\\s+genomics,\\s*transcriptomics,\\s*and\\s+proteomics\\s+of\\s+human\\s+[a-z ]+tissue[^.;]{0,120}\\b(?:cis[- ]?pqtls?|protein\\s*\\(\\s*cis[- ]?pqtls?\\s*\\))|",
      "\\b(?:performed|conducted|carried\\s+out)\\s+(?:a\\s+)?(?:gwas|genome-wide\\s+association\\s+stud(?:y|ies))\\s+of\\s+(?:measured\\s+|assayed\\s+|quantified\\s+)?(?:serum|plasma|csf|circulating|blood)\\s+(?:protein\\s+)?levels\\s+of\\s+(?:[^.;,]+,\\s*){3,}[^.;,]+(?:,\\s*and|\\s+and)\\s+[^.;]+|",
      "\\bthrough\\s+([a-z-]+\\s+){0,3}(pqtl|cis-?qtl|qtl)\\s+mapping(,|\\s+discovery|\\s+analysis)[^.;]{0,145}\\b(pqtls?|protein)\\b|",
      "\\b(?:protein|proteomic) (?:measures?|levels?|abundances?|concentrations?|profiles?|targets?) (?:[a-z0-9,()-]+ ){0,5}were (?:significantly )?(?:associated with|mapped to|linked to) (?:(?:genetic|sequence|germline|common|rare|coding) )?variants? \\((?:cis-?|trans-?|sd-?|sv-?|tr-?)?pqtls?\\)|",
      "\\banaly[sz](?:e|ed|ing) the effect of genetics on protein levels \\((?:protein quantitative trait loci \\[?pqtls?\\]?|pqtls?) analysis\\)|",
      "\\bquantified the impact of cis-?acting (?:svs?|structural variants?|variants?) on (?:several |multiple )?molecular traits including [^.;]{0,90}\\bprotein (?:abundance|levels)\\b|",
      "\\bcompared genetic,\\s*technical,\\s*and phenotypic determinants of [1-9][0-9,]* protein targets\\b|",
      "\\bevaluated the (?:shared )?genetic regulation of (?:mrna molecules,\\s*)?proteins\\b|",
      "\\bevaluate[ds]? the contribution of (?:rare |common )?genetic variants to [1-9][0-9,]* (?:plasma |serum |circulating |csf |brain )?proteins\\b|",
      "\\bcomprehensive mapping of cis-?pqtls?\\b|",
      "\\bnumber of genetic associations \\(\\s*cis-?pqtls? at p\\s*<\\s*5e-0?8\\)\\s+were\\b|",
      "\\bquantitative trait loci \\(qtls?\\) regulated [^.;]{0,45}\\bproteins?\\b|",
      "\\bexpression and protein quantitative trait loci \\(qtls?\\) highlight\\b|",
      "\\b(we\\s+)?(conduct(ed)?|perform(ed|ing)?|carr(y|ied|ying) out|describe|present|report)\\s+(a\\s+|an\\s+)?([a-z()-]+\\s+){0,4}",
      "(genome-?\\s+and\\s+[a-z-]+[ -]wide|(?:genome-wide|exome-wide|cis|local|rare-variant|admixture)[- ]?(?:and\\s+[a-z-]+-wide\\s+)?)(?:\\s+)?(regulatory\\s+|genotype-protein\\s+)?(association\\s+|mapping\\s+|gwas\\s*)?(stud(y|ies)|analys[ei]s|scans?|gwas)?(?:\\s+and\\s+exome\\s+array\\s+analyses)?\\b[^.;]{0,110}",
      "\\b(([>~]?([5-9]|[1-9][0-9,]+)|five|six|seven|eight|nine|ten|hundreds of|thousands of)\\s+([a-z-]+\\s+){0,4}",
      "(proteins?|protein levels|protein groups?|protein targets?|proteomic profiles?|aptamers?|somamers?|peptides?|cytokines?|chemokines?|(?:inflammatory|cardiovascular|neurological|circulating|protein|tumor|tumour)\\s+(?:proteins?|markers?|biomarkers?)|analytes?)|",
      "(human\\s+)?(blood\\s+|plasma\\s+|serum\\s+|csf\\s+|cerebrospinal fluid\\s+|urine\\s+|urinary\\s+|brain\\s+|salivary\\s+|circulating\\s+|tissue\\s+)?proteom[a-z]*|mass spectrometry proteomics|(?:protein|proteomic) quantitative trait loc[a-z]*|pqtls?)\\b|",
      (if (has_multi_protein_count) "\\b(performed\\s+epigenome-wide\\s+association\\s+study\\s*\\(ewas\\)\\s+and\\s+genome-wide\\s+association\\s+study\\s*\\(gwas\\)\\s+analyses|(?:biomarker|protein)-genotype\\s+associations\\s+were\\s+assessed\\s+by\\s+genome-wide\\s+association|investigated\\s+associations\\s+between\\s+single\\s+nucleotide\\s+polymorphisms[^.;]{0,85}and\\s+the\\s+(?:expression\\s+)?levels\\s+of\\s+[1-9][0-9,]+\\s+[a-z -]*proteins|multivariate\\s+gwas[\\s\\S]*\\bimplemented\\s+(?:our|the|this)\\s+(?:[a-z-]+\\s+){0,3}method\\s+on\\s+[1-9][0-9,]*\\s+(?:[a-z-]+\\s+){0,3}(?:inflammatory|cardiovascular|neurological|circulating|protein)\\s+(?:biomarkers?|proteins?|cytokines?))\\b|" else ""),
      "\\bgenome-wide association stud(y|ies) (of|for|on) (([a-z]+\\s+){0,2}levels\\s+of\\s+)?([>~]?([5-9]|[1-9][0-9,]+)\\s+([a-z-]+\\s+){0,3}(proteins?|protein levels|proteomic profiles?|cytokines?|peptides?|aptamers?|analytes?)|(the\\s+)?(human\\s+)?(blood\\s+|plasma\\s+|serum\\s+|csf\\s+|urine\\s+|urinary\\s+|brain\\s+|circulating\\s+)?(proteom[a-z]*|mass spectrometry proteomics))\\b|",
      "\\b(a\\s+)?genome-wide association (stud(y|ies)|scan|approach) (identifies|identified|to identify) ([a-z-]+\\s+){0,3}(protein quantitative trait loc[a-z]*|pqtls?)\\b|",
      "\\bgenome-wide association stud(y|ies) identif(ies|ied) [1-9][0-9,]* ([a-z-]+\\s+){0,3}loci (influencing|associated with|for) (concentrations?|levels?|abundances?) of (circulating |plasma |serum |csf |blood )?(cytokines?|chemokines?|proteins?|growth factors?)\\b|",
      "\\bgwas (?:performed |conducted )?with a (?:highly )?multiplexed [a-z -]{0,45}proteomics platform\\b|",
      "\\bquantified associations between protein level(?: changes)? and gene(?:tic)? variants\\b|",
      "\\bquantitative trait loci (underlying|associated with) (proteome|protein expression) variation\\b|",
      "\\bmapped genetic variants that are associated with [^.;]{0,65}\\bprotein abundance \\(pqtls?\\)|",
      "\\banalyzed for corresponding (plasma |serum |csf )?protein levels using genome-wide association study \\(gwas\\) genotype data and proteomic panel data\\b|",
      (if (has_multi_protein_count) "\\b(examined the association of [a-z0-9 ]+ single nucleotide polymorphisms \\(snps?\\) to [1-9][0-9]+ systemic biomarker concentrations|gwas was performed for two pneumoproteins[^.;]{0,95}and five systemic inflammatory markers|genome-wide association study[^.;]{0,95}il-6, il-8, il-10)\\b|" else ""),
      "\\bgenome-wide (genetic )?association stud(y|ies) (scans? for ([a-z0-9, ]+ and )?(over )?[1-9][0-9,]* protein quantitative trait loci|identified quantitative trait loci for ([a-z ]+ and )?protein levels|identified genetic loci associated with [^.;]{0,55}of [1-9][0-9]+ proteins)\\b|",
      (if (has_multi_protein_count) "\\bgenome-wide (genetic )?association stud(y|ies) identified ([a-z0-9, -]+ )?variants? [^.;]{0,90}\\b(modif(y|ied|ying)|associated with|influencing) (multiple )?(circulating |plasma |serum |csf |urinary |brain )?protein (levels|concentrations|abundances)\\b|" else ""),
      "\\bused genome-wide association studies to discern genetic influences on these protein concentrations\\b|",
      "\\badmixture[- ]mapping (analysis )?(of the ([a-z-]+ )?proteome|reveals genetic determinants of the [a-z ]*proteome)\\b|",
      "\\bintegrate[ds]? (?:whole-)?(?:exome|genome) sequencing [^.;]{0,130}(proteograph|mass spectrometry)\\b|",
      "\\bintegrate[ds]? [^.;]{0,140}(proteograph|mass spectrometry)[^.;]{0,95}with (?:whole-)?(?:exome|genome) sequencing[^.;]{0,120}\\b(pqtls?|protein quantitative trait loc[a-z]*)\\b|",
      "\\bperforming whole[- ]genome (and rna )?sequencing and proteomics analysis of human [a-z]+ samples\\b|",
      "\\brare homozygous predicted lof genotypes[^.;]{0,130}cis-?protein\\b|",
      "\\bcross-platform meta-analysis of (plasma |serum |circulating )?proteomics\\b|",
      "\\bgenetic variants (influencing|affecting) (the )?(human )?(plasma |serum |csf |brain |circulating )?(proteom[a-z]*|protein profile|([5-9]|[1-9][0-9,]+) [a-z -]*proteins?)\\b"
    ), primary_ctx_txt, perl = TRUE) &&
      !grepl("\\b(n-glycosylation|glycome|n-glyco[a-z]*|proteomic aging|proteomic age|aggregated trans effects)\\b", primary_ctx_txt, perl = TRUE)) ||
    (grepl("^(a |insights [a-z ]+ through )?(whole[- ]genome|whole[- ]exome|genome[- ]wide|exome[- ]wide|genome\\s+and\\s+epigenome\\s+wide) (sequencing |association )?(stud(y|ies)|association analysis|gwas|mapping|analysis) (of|for|on) (the )?([0-9,]+ )?([a-z-]+ ){0,3}(proteins?|proteomes?|proteomic|peptides?|cytokines?|protein biomarkers?|inflammatory biomarkers?)", ttl, perl = TRUE) &&
      !grepl("\\b(aggregated trans effects|polygenic|risk score|n-glyco[a-z]*|glycome|glycosylation|methylation|proteomic aging|proteomic age)\\b", ttl, perl = TRUE)) ||
    (has_multi_protein_count && grepl("^(protein quantitative trait loci identify|genetic architecture of the ([a-z-]+ ){0,3}proteome|common and rare sequence variants influencing (?:tumor|tumour|protein|circulating) biomarkers)\\b", ttl, perl = TRUE)) ||
    (p1_hit && has_multi_protein_count &&
       grepl("\\bdna\\s+was\\s+genotyped\\b", primary_ctx_txt, perl = TRUE) &&
       grepl("\\band\\s+a\\s+genome-wide\\s+association\\s+study\\s+in\\s+the\\s+[a-z0-9-]+\\s+study\\b", txt, perl = TRUE))
  )
  if (p2_hit) {
    P <- P + 6L; p_rules <- c(p_rules, "P2_PrimaryMappingExec(+6)")
  }

  # P3 (+6): Quantitative genetic pQTL / vQTL / SV-QTL / molQTL / variant-protein discovery yield (in non-MR-instrument sentences)
  p3_explicit <- !is_single_analyte && !is_somatic_or_case_only && length(non_mr_yield_sents) > 0L && any(grepl(paste0(
    "\\b(identif(y|ies|ied|ying)|find|finds|found|detect(ed|s|ing)?|discover(ed|s|ing)?|map(ped|s|ping)?|reveal(ed|s)?|uncover(ed|s)?|report(ed|s)?|confirm(ed|s)?|with only|only)\\s+",
    "(?:[>~]|(?:over|more\\s+than|nearly|almost|approximately|a\\s+total\\s+of)\\s+)?(?:[1-9][0-9,]*(?:\\s*,\\s*[1-9][0-9,]*)*\\s*,?\\s+and\\s+[1-9][0-9,]*|[1-9][0-9,]*|(?:hundreds|dozens|thousands)\\s+of)\\s+([a-z-]+\\s+){0,3}",
    "(pqtls?|cis-?pqtls?|trans-?pqtls?|sv-?pqtls?|tr-?pqtls?|protein[- ]qtls?|protein[- ]level\\s+qtls?|vqtls?|cis-?vqtls?|sd-?pqtls?|sex-differential\\s+protein\\s+quantitative\\s+(?:trait\\s+)?loci|sentinel\\s+pqtls?|pqtl\\s+associations?|variant-protein\\s+associations?|locus-protein\\s+associations?|primary\\s+genetic\\s+associations?|(?:cis-?|trans-?)?(?:protein|variance)\\s+quantitative\\s+(?:trait\\s+)?loci)\\b|",
    "\\b(?:hundreds|thousands|[1-9][0-9,]*)\\s+(?:of\\s+)?(?:cis-?|trans-?)?pqtls?\\s+(?:discovered|identified|detected|mapped)\\s+(?:for|across)\\s+[1-9][0-9,]*\\s+protein\\s+targets?\\b|",
    "\\bidentified\\s+[1-9][0-9,]*\\s+proteins?\\s+in\\s+[a-z0-9 -]+\\s+and\\s+[1-9][0-9,]*\\s+in\\s+[a-z0-9 -]+[^.;]{0,55}\\bassociations?\\s+with\\s+(?:common\\s+|rare\\s+|genetic\\s+)?variants?\\s+in\\s+cis\\b|",
    "\\brevealed\\s+[1-9][0-9,]*\\s+(?:independent\\s+)?variants?\\s+(?:across|for)\\s+[1-9][0-9,]*\\s+proteins\\b|",
    "\\b[1-9][0-9,]*\\s+(?:independent\\s+)?(?:exome[- ]array\\s+|coding\\s+|regulatory\\s+|genetic\\s+|sequence\\s+)?variants?\\s+(?:affecting|influencing)\\s+(?:the\\s+)?levels\\s+of\\s+[1-9][0-9,]*\\s+(?:[a-z-]+\\s+)?proteins?\\s+were\\s+identified\\b|",
    "\\bidentified\\s+[1-9][0-9,]*\\s+regional\\s+associations\\s+across\\s+[1-9][0-9,]*\\s+(?:different\\s+)?(?:protein\\s+|cardiometabolic\\s+|inflammatory\\s+|cardiovascular\\s+)?biomarkers?\\b|",
    "\\b(?:five|six|seven|eight|nine|ten|[5-9]|[1-9][0-9,]+)\\s+of\\s+these\\s+proteins\\s+(?:are|were)\\s+associated\\s+with\\s+cis-?(?:germline\\s+variation|variants?|snps?|pqtls?)\\b|",
    "\\bidentified\\s+(?:two|three|four|five|six|seven|eight|nine|ten|[1-9][0-9,]*)\\s+(?:cis-?|trans-?)?eqtls?\\s+and\\s+(?:four|five|six|seven|eight|nine|ten|[1-9][0-9,]*)\\s+(?:cis-?|trans-?)?(?:pqtls?|vqtls?)\\b|",
    "\\bidentified\\s+[1-9][0-9,]*\\s+(?:distinct\\s+)?proteomic\\s+traits?\\s+influenced\\s+by\\s+(?:cis[- ]?|trans[- ]?)?pqtls?\\b|",
    "\\bidentified\\s+[1-9][0-9,]*\\s+cis\\s+and\\s+[1-9][0-9,]*\\s+trans\\s+protein[- ]level\\s+qtls?\\s*\\(pqtls?\\)|",
    "\\bidentified\\s+[1-9][0-9,]*\\s+cis-?(?:snps?|qtls?|variants?)\\s+(?:significantly\\s+)?(?:influencing|affecting|associated\\s+with)\\s+(?:the\\s+)?(?:level|levels|abundance|concentrations?)\\s+of\\s+(?:specific\\s+)?(?:plasma\\s+|serum\\s+|csf\\s+|human\\s+)?proteins\\b|",
    "\\bdetected\\s+robust\\s+cis-?associations\\s+for\\s+[1-9][0-9,]*\\s+proteins\\b|",
    "\\bdetected\\s+[1-9][0-9,]*\\s+significant[^.;]{0,45}associations\\s+between\\s+[1-9][0-9,]*\\s+analytes\\s+and\\s+[1-9][0-9,]*\\s+snps\\b|",
    "\\bidentified\\s+[1-9][0-9,]*\\s+(?:independent\\s+)?(?:genome-wide(?:ly)?\\s+significant\\s+)?(?:\\([^)]+\\)\\s+)?loci\\s+(?:\\([^)]+\\)\\s+)?for\\s+(?:one\\s+or\\s+more|[1-9][0-9,]*)\\s+(?:circulating\\s+|plasma\\s+|serum\\s+|neurological\\s+|inflammatory\\s+)?(?:cytokines?|proteins?|chemokines?|biomarkers?)\\b|",
    "\\bidentif(?:y|ies|ied)\\s+[1-9][0-9,]*\\s+independent\\s+genome-wide\\s+significant[^.;]{0,55}loci\\s+for\\s+[1-9][0-9,]*\\s+proteins\\b|",
    "\\bidentified\\s+[1-9][0-9,]*\\s+novel\\s+associations[^.;]{0,65}analy[sz]ing\\s+the\\s+same\\s+biomarkers\\s+with\\s+genotyped\\s+or\\s+imputed\\s+snps\\b|",
    "\\bidentified\\s+[1-9][0-9,]*\\s+associations,\\s+(?:[a-z0-9]+\\s+of\\s+which\\s+)?[^.;]{0,65}not\\s+detected\\s+by\\s+biomarker-specific\\s+analyses\\b|",
    "\\bidentified\\s+associations\\s+for\\s+[1-9][0-9,]*\\s+of\\s+(?:the\\s+)?[1-9][0-9,]*\\s+(?:[a-z-]+\\s+){0,2}proteins\\b|",
    "\\bshared\\s+associations\\s+for\\s+[a-z ]+expression\\s+and\\s+[0-9.]+%\\s+for\\s+(?:plasma|serum|circulating)\\s+proteins\\b|",
    "\\bidentified\\s+[1-9][0-9,]*\\s+eqtls?\\s+([a-z0-9 -]+\\s+){0,5}and\\s+[1-9][0-9,]*\\s+(cis-?|trans-?)?(pqtls?|vqtls?)\\b|",
    "\\bidentified\\s+[1-9][0-9,]*\\s+([a-z-]+\\s+){0,4}associations\\s+with\\s+(circulating\\s+|plasma\\s+|serum\\s+|csf\\s+|urinary\\s+|brain\\s+)?(protein|peptide)\\s+levels\\b|",
    "\\bidentif(y|ies|ied)\\s+[1-9][0-9,]*\\s+independent\\s+associations\\s*\\([0-9,]+\\s+protein\\s+and\\s+[0-9,]+\\s+peptide\\)|",
    "\\bidentified\\s+genetic\\s+(components|determinants|loci)\\s+for\\s+([5-9]|[1-9][0-9,]+)\\s+([a-z-]+\\s+){0,2}proteins\\b|",
    "\\bspanning\\s+millions\\s+of\\s+snps\\s+and\\s+thousands\\s+of\\s+proteins\\b|",
    "\\b[>~]?[1-9][0-9,]*\\s+([a-z()-]+\\s+){0,4}(pqtls?|cis-?pqtls?|trans-?pqtls?|sv-?pqtls?|tr-?pqtls?|vqtls?|cis-?vqtls?|sd-?pqtls?|(?:cis-?|trans-?)?(?:protein|variance)\\s+quantitative\\s+(?:trait\\s+)?loci)(\\s*\\([a-z-]+\\))?\\s+",
    "(were|are)\\s+(identified|detected|found|discovered|mapped|reported|confirmed)\\b|",
    "\\b[1-9][0-9,]*\\s+(cis-?pqtls?|trans-?pqtls?|sv-?pqtls?|tr-?pqtls?|pqtls?|vqtls?|cis-?vqtls?|sentinel\\s+variants|primary\\s+genetic\\s+associations)\\s+",
    "(for|across|in|linked\\s+to|comprising|associated\\s+with|between\\s+sequence\\s+variants\\s+and\\s+levels\\s+of)\\s+([0-9,]+\\s+)?([a-z-]+\\s+){0,2}",
    "(proteins?|proteome|aptamers?|somamers?|peptides?|cytokines?|targets?|assays?|variants?)\\b|",
    "\\b[>~]?[1-9][0-9,]*\\s+(independent\\s+)?(rare\\s+noncoding\\s+)?(single-variant\\s+associations|genetic\\s+associations|snp\\s+signals|(?:genetic|sequence|common|rare)\\s+variants?)\\s+",
    "(for|with|associated\\s+with)\\s+(the\\s+)?(abundance\\s+of\\s+|levels\\s+of\\s+|concentrations\\s+of\\s+|expression\\s+of\\s+|circulating\\s+|plasma\\s+|serum\\s+|csf\\s+|urine\\s+)?",
    "([>~]?([5-9]|[1-9][0-9,]+)\\s+)?([a-z-]+\\s+){0,2}(proteins?|protein\\s+levels|protein\\s+targets?|cytokines?|peptides?)\\b|",
    "\\b[1-9][0-9,]*\\s+olink\\s+and\\s+[1-9][0-9,]*\\s+somascan\\s+proteins\\s+had\\s+(cis-?|trans-?)?pqtls?\\b|",
    "\\b[1-9][0-9,]*\\s+plasma\\s+cis-?pqtls?\\s+colocalized\\b|",
    "\\bquantitative trait loci \\(qtls?\\) regulated at least (one-third|one-fourth|one-half|[0-9]+%) of the proteins\\b"
  ), non_mr_yield_sents, perl = TRUE) &
    !grepl("\\b(valid cis-?pqtls? for use as|mqtls?, [0-9,]+ eqtls?, and [0-9,]+ pqtls?)\\b", non_mr_yield_sents, perl = TRUE))

  p3_multiomic_yield <- !is_single_analyte && !is_somatic_or_case_only && p1_hit &&
    grepl("\\bproteom(ics|ic|es?)\\b", primary_ctx_txt, perl = TRUE) &&
    length(non_mr_yield_sents) > 0L &&
    any(grepl("\\b(identif(y|ies|ied)|discover(ed|s)|detect(ed|s)|map(ped|s))\\s+(?:[>~]|(?:over|more\\s+than|nearly|almost|approximately|a\\s+total\\s+of)\\s+)?[1-9][0-9,]*\\s+([a-z-]+\\s+){0,2}(sv-?qtls?|tr-?xqtls?|molqtls?|xqtls?|qqtls?|molecular\\s+quantitative\\s+trait\\s+loc[a-z]*)\\b", non_mr_yield_sents, perl = TRUE))

  p3_sent_ctx <- !is_single_analyte && !is_somatic_or_case_only && length(primary_ctx_sents) > 0L && any(
    grepl("\\b(pqtls?|vqtls?|proteins?|proteom[a-z]*|aptamers?|somamers?|peptides?|cytokines?)\\b", primary_ctx_sents, perl = TRUE) &
    !grepl("\\b(metabolites?|lipids?|glycans?|methylation sites|cpgs?|mri|suicid|insomnia|melanoma|allergic|sarcopenia)\\b", primary_ctx_sents, perl = TRUE) &
    grepl("\\b(identif(ied|ies|y|ying)|find|found|detect(ed)?|discover(ed)?|map(ped)?)\\s+[>~]?[1-9][0-9,]*\\s+([a-z-]+\\s+){0,2}(primary\\s+)?(genetic\\s+|variant-protein\\s+)(associations|signals)\\b", primary_ctx_sents, perl = TRUE)
  )
  p3_hit <- p3_explicit || p3_multiomic_yield || p3_sent_ctx
  if (p3_hit) {
    P <- P + 6L; p_rules <- c(p_rules, "P3_DiscoveryYield(+6)")
  }

  # P4 (+4): Multi-protein assay scale (>= 5) + human cohort N in primary context sentences
  cohort_n_pat <- paste0(
    "(\\b(in|across|of|among|from|involving|encompassing|stud(y|ied)|enrolled|recruited|analy[sz](ed|ing)|profiled)\\s+",
    "(?:[>~]|(?:over|more\\s+than|nearly|almost|approximately|up\\s+to)\\s+)?[1-9][0-9,]{1,6}\\s+([a-z()/+-]+\\s+){0,3}",
    "(participants|individuals|subjects|persons|people|donors|volunteers|icelanders|europeans|british|scottish|finnish|finns|danish|greenlandic|japanese|chinese|latin\\s+americans?|",
    "adults|children|adolescents|pediatric|neonates|infants|newborns?|babies|twins|centenarians|brains|brain\\s+proteomes|kidney\\s+samples|gtm\\s+samples|patients|samples|men|women|males|females|lcls?|cell\\s+lines)\\b|",
    "(?:\\(|,\\s*)n_?(?:cases|controls)?\\s*=\\s*[>~]?[1-9][0-9,]{1,6}(?:\\s+[a-z -]{0,30})?\\s*\\)|\\b[1-9][0-9,]{1,6}\\s+participants\\s+were\\s+analy[sz]ed\\b)"
  )
  p4_hit <- !is_single_analyte && has_multi_protein_count &&
    any(grepl(cohort_n_pat, primary_ctx_sents, perl = TRUE))
  if (p4_hit) {
    P <- P + 4L; p_rules <- c(p_rules, "P4_AssayScaleAndN(+4)")
  }

  # P5 (+3): Germline Genotyping / Sequencing / Imputation / Variant QC / QTL Architecture
  # Merge Safeguard 7: Include spelled-out `protein quantitative (trait )?loc[a-z]*` alongside pQTL/vQTL/molQTL!
  if (!is_somatic_or_case_only && grepl(paste0(
    "(\\bgenotyped\\b|\\bgenotyping array\\b|\\bsnp array\\b|\\bexome[- ]array\\b|\\bsnps?\\b|\\bsingle[- ]nucleotide polymorphisms?\\b|\\bgenome-wide scan\\b|",
    "\\bwhole[ -]genome sequenc[a-z]*\\b|\\bwgs\\b|\\bwhole[ -]exome sequenc[a-z]*\\b|\\bwes\\b|\\bexome sequenc[a-z]*\\b|\\blong-read genome sequenc[a-z]*\\b|",
    "\\badmixture[- ]mapping\\b|\\bhuman knockouts?\\b|\\bloss-of-function variants?\\b|\\bbiallelic\\b|",
    "\\bimput(ed|ation)\\b|\\btopmed\\b|\\bhaplotype reference consortium\\b|\\b1000 genomes\\b|",
    "\\bsequence variants\\b|\\brare (coding |noncoding |homozygous )?variants?\\b|\\bminor allele frequenc[a-z]*\\b|",
    "\\bregenie\\b|\\bbolt-lmm\\b|\\bfastgwa\\b|\\bsnptest\\b|\\bplink\\b|\\bsaige\\b|\\bmatrixeqtl\\b|\\btensorqtl\\b|\\bqtltools\\b|\\bmesusie\\b|\\bsushie\\b|",
    "\\bsnp[- ]based heritability\\b|\\bsnp-h2\\b|\\bheritability\\b|\\bepitopes?\\b|\\bprotein[- ]altering variants?\\b|",
    "\\bpqtls?\\b|\\bcis-?pqtls?\\b|\\btrans-?pqtls?\\b|\\bsd-?pqtls?\\b|\\bsv-?pqtls?\\b|\\btr-?pqtls?\\b|\\bprotein[- ]qtls?\\b|",
    "\\bvqtls?\\b|\\bcis-?vqtls?\\b|\\bcis-?molqtls?\\b|\\bsv-?qtls?\\b|\\btr-?xqtls?\\b|\\bqqtls?\\b|\\b(?:protein|molecular|multi-?omic|quantile) quantitative (trait )?loc[a-z]*\\b)"
  ), txt, perl = TRUE)) {
    P <- P + 3L; p_rules <- c(p_rules, "P5_GenomicArchitecture(+3)")
  }

  # P6 (+2): High-priority biofluid/tissue, next-gen MS/NULISA platform, small-N/pediatric, or non-EUR ancestry
  # Merge Safeguard 7 + Agent 2 Fix 5/6: Include `ipscs?|induced pluripotent`, `maternal serum`, `multiancestry`, and scoped TMT.
  if (grepl(paste0(
    "(\\bcerebrospinal fluid\\b|\\bcsf\\b|\\bdorsolateral prefrontal\\b|\\bprefrontal cortex\\b|\\bbraingvex\\b|\\bgyrus\\b|",
    "\\bbrain tissue\\b|\\bbrain proteomes?\\b|\\bpost-?mortem\\b|\\bdorsal root ganglia\\b|",
    "\\burinary\\b|\\burine\\b|\\bsalivary\\b|\\bsaliva\\b|\\bdried blood spots?\\b|\\bmaternal serum\\b|",
    "\\bkidney samples?\\b|\\bhuman kidney\\b|\\bliver\\b|\\batrial tissue\\b|\\bcardiac tissue\\b|\\barterial wall\\b|\\baortic\\b|\\bskeletal muscle\\b|",
    "\\bgtex\\b|\\bstarnet\\b|\\bgenes & health\\b|\\blymphoblastoid\\b|\\blcls?\\b|\\bipscs?\\b|\\binduced pluripotent\\b|",
    "\\bnanoparticle\\b|\\bproteograph\\b|\\borbitrap astral\\b|\\bnulisa[a-z]*\\b|\\bdia-ms\\b|\\bdia-pasef\\b|\\bswath\\b|",
    "\\btandem mass tags?\\b|\\btmt-(based|ms|[0-9]{1,2})\\b|\\blc-ms/ms\\b|\\bsomascan 11k\\b|\\bolink (explore )?ht\\b|",
    "\\bnon-european\\b|\\bafrican ancestry\\b|\\bafrican american\\b|\\beast asian\\b|\\bchinese\\b|\\bjapanese\\b|",
    "\\bsouth asian\\b|\\bbangladeshi\\b|\\bpakistani\\b|\\bgreenlandic\\b|\\binuit\\b|\\bhispanic\\b|\\blatin americans?\\b|\\barab\\b|",
    "\\badmixed\\b|\\bmultiancestry\\b|\\bmulti-ancestry\\b|\\bmulti-ethnic\\b|\\bcross-ancestry\\b|",
    "\\bpediatric\\b|\\bchildren\\b|\\badolescents\\b|\\bneonates\\b|\\bnewborns?\\b|\\binfants\\b|\\btwins\\b|\\bcentenarians\\b)"
  ), txt, perl = TRUE)) {
    P <- P + 2L; p_rules <- c(p_rules, "P6_PriorityMatrixAncestryTech(+2)")
  }

  # P7 (+5): Verified EBI GWAS Catalog multi-protein trait count (>= 5 distinct protein traits)
  if (!is.na(gwas_catalog_study_count) && gwas_catalog_study_count >= 5L) {
    P <- P + 5L; p_rules <- c(p_rules, "P7_GWASCatalogMultiTrait(+5)")
  }
  # P7b (+5): Stage B.5 PMC full-text verified primary pQTL sub-layer inside multi-omic xQTL/molQTL study
  if (isTRUE(ft_verified_pqtl)) {
    P <- P + 5L; p_rules <- c(p_rules, "P7b_FullTextMultiOmicPQTL(+5)")
  }

  # P8 (+2): Genuine proteomics resource / atlas / summary-statistics deposition statement
  if (grepl(paste0(
    "\\b(open-access (proteomics |proteome-phenome )?resource|valuable resource|publicly available database|",
    "convenient web application|interactive web[- ]?(?:server|tool|application)|",
    "summary statistics (are|have been) ([a-z]+ )?(available|deposited|publicly)|",
    "(comprehensive |large-scale )?genomic atlas of ([a-z -]+)?proteome|",
    "atlas of protein quantitative trait loci)\\b"
  ), txt, perl = TRUE)) {
    P <- P + 2L; p_rules <- c(p_rules, "P8_ResourceDeposition(+2)")
  }

  # P9 (+6): Primary consortium / multi-cohort / cross-platform pGWAS meta-analysis protection
  p9_hit <- (
    grepl("\\b(meta-analys(is|es|ed|ing)|combined analysis) of ([>~]?([5-9]|[1-9][0-9,]+) ([a-z-]+ )?(proteins?|circulating cytokines?|cytokines?|peptides?|biomarkers?)|(plasma |serum |circulating )?proteomics)\\b|\\bscallop consortium\\b|\\b[0-9,]+ cardiovascular proteins in [0-9,]+ individuals\\b", txt, perl = TRUE) ||
      grepl("\\bcross-platform meta-analysis of ([a-z-]+ )?proteomics\\b|\\bgwas\\)? for ([5-9]|[1-9][0-9,]+) ([a-z-]+ )?(cytokines?|proteins?) in meta-analyses\\b", txt, perl = TRUE)
  ) && !grepl("\\b(two-sample|2-sample) mendelian\\b", txt, perl = TRUE)
  if (p9_hit) {
    P <- P + 6L; p_rules <- c(p_rules, "P9_PrimaryConsortiumMeta(+6)")
  }

  # --- Secondary Summary-Stat Consumer & Non-Proteomic GWAS Signals (S_score) ---

  # S1 (+7): Public summary-statistics / pQTL instrument download or reuse
  # Merge Safeguard 5: Require explicit external source or summary-statistics object (never bare
  # "as instrumental variables" or "derived from a cohort of N patients"), and do NOT gate on `!p2_hit`!
  s1_hit <- any(vapply(sents, is_public_reuse_sentence, logical(1))) || grepl(paste0(
    "\\b(used|using|leveraged|leveraging|utilized|utilizing|employed|employing|integrated|integrating) ",
    "(publicly available |published |summary-level |genome-wide significant |sex-specific |ancestry-specific )?",
    "(gwas |pqtl |cis-?pqtl |proteomic |genetic )?(summary )?(data|statistics|datasets?|instruments?|loci|proxies) ",
    "(for [0-9,]+ [a-z ]+ |of plasma proteins )?(from|in) (the )?",
    "(uk biobank pharma|ukb-?ppp|decode|fenland|interval|scallop|aric|atherosclerosis risk in communities|finngen|ieu open|gwas catalog|",
    "(two|three|four|five|seven|eight|nine|ten|[2-9]) (large|recent|published|independent|comprehensive|proteomic))\\b|",
    "\\bwe use[ds]? (publicly available |published )?summary(-level)? (statistics|data) (from|for|of|on)\\b|",
    "\\bcollated data on plasma proteins from a protein quantitative trait loci\\b|",
    "\\bsystematically integrates two major genome-wide pqtl databases\\b"
  ), txt, perl = TRUE)
  if (s1_hit) {
    S <- S + 7L; s_rules <- c(s_rules, "S1_PublicSumStatsReuse(+7)")
  }

  # S2 (+6): 2SMR / SMR / PW-MR / Drug-Target MR / Coloc boilerplate
  if (grepl(paste0(
    "\\b(two-sample mendelian|2-sample mendelian|two sample mendelian|two-sample mr|2-sample mr|two-sample cis mendelian|",
    "proteome-wide mendelian|protein-wide mendelian|proteomic mendelian|drug-target mendelian|drug-target mr|cis-mendelian randomi[sz]ation|",
    "druggable[- ]genome(-wide)? mendelian|bidirectional ([a-z-]+ )?mendelian|multivariable mendelian|\\bmvmr\\b|network mendelian|",
    "two-step (network )?mendelian|2-step mendelian|2-step mr|multi-omics mr|mediation mendelian|",
    "summary-data-based mendelian|summary data-based mendelian|summary-data-based mr|summary-based mendelian|summary-based mr|",
    "inverse-variance weighted|inverse variance weight[a-z]*|\\bivw\\b|mr-egger|weighted median|mr-presso|wald ratios?|",
    "cochran'?s q|egger intercept|leave-one-out(?!\\s+cross-?validation)|smr and heidi|\\bheidi\\b|steiger filtering|twosample\\s*mr|coloc\\.abf|coloc\\s+package|pp\\.?h4(\\.abf)?\\s*[>=])\\b"
  ), txt, perl = TRUE)) {
    S <- S + 6L; s_rules <- c(s_rules, "S2_2SMREstimatorBoilerplate(+6)")
  }

  # S3 (+4): Secondary MR / PWAS / Target-Discovery Title Templates
  if (grepl("\\b(mendelian randomi[sz]ation|two-sample|2-sample|proteome-wide association stud(y|ies)|\\bpwas\\b)\\b", ttl, perl = TRUE) ||
      grepl("\\b(therapeutic|druggable|drug|causal|potential|novel|promising|candidate) (drug |protein |gene )?(targets?|biomarkers?) (for|in|of) [a-z0-9 ,/-]+\\b", ttl, perl = TRUE) ||
      grepl("\\bcausal (role|effect|relationship|association|link|nexus|impact|pathway)s? (of|between|among)\\b", ttl, perl = TRUE)) {
    S <- S + 4L; s_rules <- c(s_rules, "S3_SecondaryTitleTemplate(+4)")
  }

  # S4 (+4): Public proteomic cohort(s) paired with external disease outcome consortium
  pub_cohorts <- sum(vapply(
    c("\\bukb-ppp\\b|pharma proteomics project", "\\bdecode\\b|35,559", "\\bfenland\\b|10,708",
      "\\baric\\b|atherosclerosis risk in communities", "\\binterval\\b", "\\bscallop\\b"),
    function(pat) grepl(pat, txt, perl = TRUE),
    logical(1)
  ))
  has_outcome_db <- grepl("\\b(finngen|ieu open gwas|opengwas|megastroke|gigastroke|diagram consortium|cardiogram|practical consortium|ovarian cancer association consortium|pgc|igap|million veteran program)\\b", txt, perl = TRUE)
  if ((pub_cohorts >= 2L && has_outcome_db) ||
      (pub_cohorts >= 2L && grepl("mendelian randomi[sz]ation", txt, perl = TRUE) && !p2_hit && !p9_hit)) {
    S <- S + 4L; s_rules <- c(s_rules, "S4_PublicCohortOutcomePair(+4)")
  }

  # S5 (+5): Pre-trained PWAS / TWAS weight consumer
  if (grepl("\\b(fusion|predixcan|s-predixcan|utmost|jti|otters|omicspred|bliss)\\b", txt, perl = TRUE) &&
      grepl("\\b(pwas|proteome-wide association|genetically predicted protein|predicted plasma proteomics|prediction models from)\\b", txt, perl = TRUE) &&
      !grepl("\\b(train(ed|ing)|build(ing)?|models for proteome-wide|cis-heritability|sushie)\\b", txt, perl = TRUE)) {
    S <- S + 5L; s_rules <- c(s_rules, "S5_PWASWeightConsumer(+5)")
  }

  # S6 (+6): Primary GWAS of a non-protein phenotype or single-analyte (<5 proteins) GWAS
  is_non_protein_primary_gwas <- !p2_hit && !p3_hit && !p9_hit &&
    (!(!is.na(gwas_catalog_study_count) && gwas_catalog_study_count >= 5L)) && (
      grepl(paste0(
        "\\bwe (conducted|performed|carried out|report) ([a-z, -]{0,30})?",
        "(gwas|genome-wide association|exome-wide association|cross-phenotype gwas|multivariate gwas|genomic sem)([a-z, -]{0,25})?",
        "(of|for|on) (?!((?:wgs|wes|genotype|imputed) data (?:for|on|of) )?(?:measured |assayed |quantified )?([0-9,]+ )?(human |plasma |serum |csf |cerebrospinal |urinary |urine |brain |circulating |blood |inflammatory |cardiovascular |cardiometabolic |neurological |complement |tumor |tumour |protein )?(proteins?|proteom|cytokin|chemokin|peptid|aptamer|somamer|(?:inflammatory |cardiovascular |cardiometabolic |neurological |protein |tumor |tumour )?biomarkers?))"
      ), txt, perl = TRUE) ||
        grepl("\\b(quantified [0-9,]+ serum metabolites and performed gwas|genome-wide association study of n-glycosylation|oral microbiomes of [0-9,]+ people|polygenic risk scores? with the plasma proteome)\\b", txt, perl = TRUE)
    )
  if (is_non_protein_primary_gwas || is_single_analyte) {
    S <- S + 6L; s_rules <- c(s_rules, "S6_NonProteinOrSingleAnalyteGWAS(+6)")
  }

  # G4: Narrative / Systematic Review Gate (when lacking primary mapping execution)
  is_review <- any(pts %in% c("review", "systematic review")) ||
    grepl("^(a |an |updated )?(systematic |narrative |comprehensive )?review\\b|\\b(a |systematic |narrative )review( and meta-analysis)?$|: a review\\b", ttl, perl = TRUE)
  if (is_review && !p2_hit && !p3_hit) {
    return(list(
      stage = "Stage1_ReviewGate", P_score = P, S_score = S,
      p_rules = p_rules, s_rules = c(s_rules, "G4_ReviewLowPrimary"),
      classification = "Auto-Deprioritised: Review Article"
    ))
  }

  # -------------------------------------------------------------------------
  # STAGE 2 DECISION MATRIX (Two-Tier Germline Anchor: Primary vs Borderline)
  # -------------------------------------------------------------------------
  has_multi_protein_scale <- has_multi_protein_count ||
    grepl("\\b(proteom(e|es|ic|ics)|proteogenom[a-z]*|olink|somascan|somalogic|proteograph|nulisa[a-z]*|swath|dia-ms|tandem mass tag|luminex|x-?map)\\b", txt, perl = TRUE) ||
    (!is.na(gwas_catalog_study_count) && gwas_catalog_study_count >= 5L) ||
    isTRUE(ft_verified_pqtl)

  # Tier 1 Anchor (`Likely Primary pGWAS`): Direct multiplexed pGWAS / pQTL / vQTL mapping execution
  has_direct_pqtl_mapping <- !is_single_analyte && !is_somatic_or_case_only && !is_composite_clock_only && (
    p2_hit || p3_hit || p9_hit ||
      (!is.na(gwas_catalog_study_count) && gwas_catalog_study_count >= 5L && has_multi_protein_scale) ||
      (p1_hit && p4_hit && S == 0L &&
         grepl("\\bintegrat(ed|ing) proteomic[a-z, ]+and genomic data\\b", primary_ctx_txt, perl = TRUE) &&
         grepl("\\b(csf|brain|tissue|plasma|serum) protein quantitative trait loci \\(pqtls?\\)", primary_ctx_txt, perl = TRUE))
  )

  # Merge Safeguard 8: Fix `([a-z -]+)` in MR exclusion guard (catching hyphenated `cis-protein`),
  # exclude composite aging clocks (`Morseburg`), pedigree-only twin heritability (`Drouard`),
  # and targeted candidate-SNP assays (`Zhang`, `Vu`), while retaining `cis-molQTL` fine-mapping (`Lu`),
  # custom-array pQTL panels (`Maya`), non-additive pQTL scans (`Sigurdsson`), full-text-verified
  # multi-omic xQTL/qQTL/TR-xQTL studies (`ft_verified_pqtl`), in-cohort multi-protein prediction-model
  # training (`Zhong`), and in-cohort `cis-pQTL` causal protein network/locus hybrids (`Bankier`, `Emilsson`, `Li`).
  mr_only_pqtl_guard <- grepl(paste0(
    "\\b(pqtl-mr|used protein quantitative trait loci \\(pqtl\\) mendelian|",
    "mendelian randomi[sz]ation ([a-z-]+ ){0,4}(leveraging|using|utilizing) ([a-z -]+)?protein quantitative trait|",
    "leveraging publicly available data from genome-wide association studies \\(gwas\\) and protein quantitative|",
    "discovery mr utilized protein quantitative trait|utilized genome-wide significant cis-protein quantitative|",
    "protein quantitative trait loci analyze studies)\\b"
  ), txt, perl = TRUE)

  is_custom_or_pathway_pqtl <- grepl("\\bused protein quantitative trait loci \\(pqtl\\) analyses, adjusted for covariates\\b", primary_ctx_txt, perl = TRUE)

  has_germline_protein_link <- !is_single_analyte && !is_somatic_or_case_only && !is_composite_clock_only && (
    has_direct_pqtl_mapping ||
      (!mr_only_pqtl_guard && (
        is_custom_or_pathway_pqtl ||
        (isTRUE(ft_verified_pqtl) && !s1_hit && S <= 6L) ||
        (p1_hit && p4_hit && S == 0L && !s1_hit &&
           grepl("\\bprotein profiles of which some could be explained by associations to genetic variants\\b", primary_ctx_txt, perl = TRUE)) ||
        (grepl("\\bfine-maps? cis-?molqtls?\\b", primary_ctx_txt, perl = TRUE) &&
           grepl("\\bprotein levels\\b", primary_ctx_txt, perl = TRUE)) ||
        (has_multi_protein_scale && grepl(paste0(
          "\\b(genetic dominance and epistasis|identify complex effects that are associated with protein quantitative trait loci|",
          "mapping protein quantitative trait loci \\(pqtls?\\) for these variants|",
          "[1-9][0-9,]+ proteins were regulated by genetic [a-z ]+ loci|",
          "[1-9][0-9,]*(?:/[1-9][0-9,]*)? (?:validated |plasma |serum )?proteins exhibited genetic control|",
          "molecular xqtls?|",
          "mapping [1-9][0-9,]+ coding germline variants onto peptides)\\b"
        ), primary_ctx_txt, perl = TRUE)) ||
        (!s1_hit && has_multi_protein_count &&
           grepl("\\b(develop(ed)?|train(ed)?|buil(d|t)|construct(ed)?|establish(ed)?) ([a-z-]+ ){0,3}(genetic |cis-?genetic |predictive )?(?:prediction )?models for (?:the )?(?:(?:expression |levels |abundance |abundances |concentrations )?(?:of )?[1-9][0-9,]+|protein abundances?\\b[^.;]{0,95}\\bmeasured [1-9][0-9,]* proteins?)", primary_ctx_txt, perl = TRUE)) ||
        (!s1_hit && p1_hit && p4_hit &&
           grepl("\\b(?:compared the genetic determinants of protein levels|(?:genetic variants|snps?) [^.;]{0,65}were correlated to [1-9][0-9,]* proteins that were profiled)\\b", txt, perl = TRUE)) ||
        ((p1_hit || p4_hit) && S == 0L && !s1_hit &&
           (grepl("\\breconstruct a causal network of (?:serum |plasma )?proteins\\b|\\bcausal (?:network|protein|serum protein|plasma protein)(?: [a-z-]+){0,2} (?:inference|networks?|subnetworks?|proteins?)\\b|\\bcausal (?:[a-z-]+ ){0,3}protein (?:sub)?networks?\\b", txt, perl = TRUE) &&
            grepl("\\bcis-?(?:acting )?(?:protein quantitative trait loc[a-z]*|pqtls?)\\b", txt, perl = TRUE)))
      ))
  )

  # Merge Safeguard 6: No-abstract title fallback strictly requires proteome/pQTL/protein-GWAS framing
  # and excludes disease-first "dissection of [disease]" titles (`Wang MDD Letter`).
  no_abs_primary_title <- !nzchar(abs) && S < 4L && !is_single_analyte && !is_composite_clock_only &&
    grepl("\\b(pqtls?|vqtls?|pgwas)\\b|\\b(gwas|genome-wide association|exome-wide association|genetic (architecture|regulation|determinants|mapping|atlas)) (of|for|on) (the )?([a-z-]+ ){0,3}(proteom[a-z]*|proteins?|cytokines?|peptides?)\\b", ttl, perl = TRUE) &&
    !grepl("\\bdissection of\\b", ttl, perl = TRUE)

  cls <- if (has_direct_pqtl_mapping && P >= 7L && (S <= 6L || P >= S + 2L)) {
    "Likely Primary pGWAS"
  } else if (!is_single_analyte && !is_composite_clock_only && (
               (has_direct_pqtl_mapping && P >= 5L && (S <= 6L || P >= S)) ||
               (has_germline_protein_link && (P >= 4L || is_custom_or_pathway_pqtl) && S <= 6L) ||
               no_abs_primary_title
             )) {
    "Needs Human Check: Possible Primary / Hybrid"
  } else if (S >= 4L) {
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
