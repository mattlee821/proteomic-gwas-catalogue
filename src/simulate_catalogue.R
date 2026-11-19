#!/usr/bin/env Rscript
# Generate a realistic simulated proteomic GWAS catalogue (v0.2.0-sim) to drive
# website/app design, schema validation, and release workflow testing.
#
# Covers:
# - Platforms: SomaLogic (1.1k to 11k), Olink (Target 96 to Explore HT 5k),
#   Mass Spectrometry (DIA-MS, TMT-MS, Seer Proteograph), Alamar NULISAseq, Luminex
# - Sample sizes: no floor (N = 64 brain tissue up to N = 54,219 biobank)
# - Publication types: peer-reviewed articles, active preprints, superseded preprints
# - Matrices: EDTA plasma, Citrate plasma, Serum, CSF, Brain tissue, Urine
# - Strata: multi-ancestry (EUR, AFR, EAS, CSA, AMR, Multi-ancestry) & sex-stratified

out_dir <- "catalogue"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write_tsv <- function(df, name) {
  path <- file.path(out_dir, name)
  utils::write.table(
    df,
    file = path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
  message("Wrote ", nrow(df), " rows -> ", path)
}

# 1. STUDIES (1 row per publication / preprint) ====
studies <- data.frame(
  study_id = c(
    "suhre_2017_28240269",
    "folk_2017_28369058",
    "sun_2018_29875488",
    "yao_2018_30148068",
    "folkersen_2020_33067605",
    "ferkingstad_2021_34857953",
    "pietzner_2021_34648354",
    "yang_2021_34239129",
    "gudjonsson_2022_35078996",
    "zhang_2022_35501419",
    "sun_2022_preprint",
    "sun_2023_37794186",
    "eldjarn_2023_37794188",
    "wingo_2023_brain_ms",
    "schlosser_2023_urine",
    "chen_2024_ckb_dual",
    "koprulu_2025_sex_pgwas",
    "campbell_2025_gs_diams",
    "suhre_2025_proteograph",
    "lopez_2026_csf_nulisa_pre",
    "patel_2026_somascan7k_pre",
    "nakamura_2026_olinkht_pre"
  ),
  first_author = c(
    "Suhre", "Folkersen", "Sun", "Yao", "Folkersen",
    "Ferkingstad", "Pietzner", "Yang", "Gudjonsson", "Zhang",
    "Sun", "Sun", "Eldjarn", "Wingo", "Schlosser",
    "Chen", "Koprulu", "Campbell", "Suhre", "Lopez-Garcia",
    "Patel", "Nakamura"
  ),
  year = c(
    2017L, 2017L, 2018L, 2018L, 2020L,
    2021L, 2021L, 2021L, 2022L, 2022L,
    2022L, 2023L, 2023L, 2023L, 2023L,
    2024L, 2025L, 2025L, 2025L, 2026L,
    2026L, 2026L
  ),
  status = c(
    "Published", "Published", "Published", "Published", "Published",
    "Published", "Published", "Published", "Published", "Published",
    "Superseded Preprint", "Published", "Published", "Published", "Published",
    "Published", "Published", "Published", "Published", "Preprint",
    "Preprint", "Preprint"
  ),
  superseded_by = c(
    "", "", "", "", "",
    "", "", "", "", "",
    "sun_2023_37794186", "", "", "", "",
    "", "", "", "", "",
    "", ""
  ),
  pmid = c(
    "28240269", "28369058", "29875488", "30148068", "33067605",
    "34857953", "34648354", "34239129", "35078996", "35501419",
    "", "37794186", "37794188", "36959201", "37812445",
    "38412099", "40360480", "40911204", "41310232", "",
    "", ""
  ),
  doi = c(
    "10.1038/ncomms14357",
    "10.1371/journal.pgen.1006706",
    "10.1038/s41586-018-0175-2",
    "10.1038/s41467-018-05512-x",
    "10.1038/s42255-020-00287-2",
    "10.1038/s41588-021-00978-w",
    "10.1038/s41467-021-25589-7",
    "10.1038/s41593-021-00886-6",
    "10.1038/s41467-021-27850-z",
    "10.1038/s41588-022-01051-w",
    "10.1101/2022.06.17.496443",
    "10.1038/s41586-023-06592-6",
    "10.1038/s41586-023-06563-x",
    "10.1038/s41588-023-01350-0",
    "10.1038/s41588-023-01501-3",
    "10.1038/s41467-024-46100-2",
    "10.1038/s41467-025-56112-9",
    "10.1038/s41588-025-02114-5",
    "10.1038/s41467-025-60911-4",
    "10.1101/2026.04.12.26305118",
    "10.1101/2026.02.19.581022",
    "10.1101/2026.09.03.611408"
  ),
  journal = c(
    "Nature Communications", "PLOS Genetics", "Nature", "Nature Communications", "Nature Metabolism",
    "Nature Genetics", "Nature Communications", "Nature Neuroscience", "Nature Communications", "Nature Genetics",
    "bioRxiv", "Nature", "Nature", "Nature Genetics", "Nature Genetics",
    "Nature Communications", "Nature Communications", "Nature Genetics", "Nature Communications", "medRxiv",
    "bioRxiv", "bioRxiv"
  ),
  title = c(
    "Connecting genetic risk to disease end points through the human blood plasma proteome",
    "Mapping of 79 loci for 83 plasma protein biomarkers in cardiovascular disease",
    "Genomic atlas of the human plasma proteome",
    "Clinical cardiovascular disease risk variants are enriched in plasma protein quantitative trait loci",
    "Genomic and drug target evaluation of 90 cardiovascular proteins in 30,931 individuals",
    "Large-scale integration of the plasma proteome with genetics and disease",
    "Mapping the proteo-genomic convergence of human diseases",
    "Genomic atlas of the proteome from brain, CSF and plasma in sporadic and genetic Alzheimer disease",
    "A genome-wide association study of serum proteins reveals shared loci with common diseases",
    "Plasma proteome analyses in individuals of European and African ancestry identify cis-pQTLs and models for proteome-wide association studies",
    "Genetic regulation of the human plasma proteome in 54,306 UK Biobank participants (Preprint)",
    "Plasma proteomic associations with genetics and health in the UK Biobank",
    "Large-scale plasma proteomics comparisons through genetics and disease associations",
    "Deep brain tissue TMT mass-spectrometry pQTL mapping in dorsolateral prefrontal cortex",
    "Genetic determinants of the urinary proteome in chronic kidney disease",
    "Comparative genetic architecture of 2,168 plasma proteins across affinity platforms in Chinese adults",
    "Sex differences in the genetic regulation of the human plasma proteome",
    "Genome-wide analysis of 439 mass spectrometry-based proteomic profiles in 15,035 Scottish individuals",
    "A genome-wide association study of mass spectrometry proteomics using a nanoparticle enrichment platform",
    "Pediatric cerebrospinal fluid pQTL mapping using attomolar NULISAseq immunoassays",
    "Beyond the classical plasma secretome: expanded SomaScan 7k genetic architecture in 13,445 Europeans",
    "Multi-ancestry pGWAS of 5,400 plasma proteins using Olink Explore HT and SomaScan 11k Pilot"
  ),
  repository = c(
    "Metabolomix / Helmholz", "Zenodo", "EBI GWAS Catalog / Box", "dbGaP", "SCALLOP / Zenodo",
    "deCODE Genetics Portal", "Omicscience Portal", "NIAGADS / WashU Portal", "Zenodo", "ARIC / Zenodo",
    "Synapse (syn51364943)", "Synapse (syn51364943)", "deCODE Genetics Portal", "Synapse (syn2580853)", "EBI GWAS Catalog",
    "CNCB OMIX", "Omicscience / Synapse", "EBI GWAS Catalog", "Figshare / PRIDE", "Zenodo (Preprint)",
    "EBI GWAS Catalog (Staged)", "Synapse (syn64019920)"
  ),
  accession_or_url = c(
    "http://proteomics.gwas.eu/",
    "https://doi.org/10.5281/zenodo.264128",
    "https://www.ebi.ac.uk/gwas/publications/29875488",
    "https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000007",
    "https://doi.org/10.5281/zenodo.2615265",
    "https://www.decode.com/summarydata/",
    "https://omicscience.org/apps/pgwas/",
    "https://www.niagads.org/datasets/ng00102",
    "https://doi.org/10.5281/zenodo.5711426",
    "http://nilanjanchatterjeelab.org/pwas/",
    "https://www.synapse.org/Synapse:syn51364943",
    "https://www.synapse.org/Synapse:syn51364943",
    "https://www.decode.com/summarydata/",
    "https://www.synapse.org/Synapse:syn2580853",
    "https://www.ebi.ac.uk/gwas/publications/37812445",
    "https://ngdc.cncb.ac.cn/omix/",
    "https://omicscience.org/apps/sex_pgwas/",
    "https://www.ebi.ac.uk/gwas/publications/40911204",
    "https://figshare.com/projects/Nanoparticle_MS_pGWAS/198201",
    "https://doi.org/10.5281/zenodo.14991021",
    "https://www.ebi.ac.uk/gwas/downloads/summary-statistics",
    "https://www.synapse.org/Synapse:syn64019920"
  ),
  access_level = c(
    "Open download", "Open download", "Open download", "Controlled access", "Open download",
    "Portal registration", "Open download", "Portal registration", "Upon request", "Open download",
    "Portal registration", "Portal registration", "Portal registration", "Portal registration", "Open download",
    "Portal registration", "Open download", "Open download", "Open download", "Open download",
    "Open download", "Portal registration"
  ),
  added_release = c(
    "v0.1.0-sim", "v0.1.0-sim", "v0.1.0-sim", "v0.1.0-sim", "v0.1.0-sim",
    "v0.1.0-sim", "v0.1.0-sim", "v0.1.0-sim", "v0.1.0-sim", "v0.1.0-sim",
    "v0.1.0-sim", "v0.1.0-sim", "v0.1.0-sim", "v0.1.0-sim", "v0.1.0-sim",
    "v0.1.0-sim", "v0.1.0-sim", "v0.1.0-sim", "v0.1.0-sim", "v0.2.0-sim",
    "v0.2.0-sim", "v0.2.0-sim"
  ),
  verified_full_text = c(
    TRUE, TRUE, TRUE, TRUE, TRUE,
    TRUE, TRUE, TRUE, TRUE, TRUE,
    TRUE, TRUE, TRUE, TRUE, TRUE,
    TRUE, TRUE, TRUE, TRUE, FALSE,
    FALSE, FALSE
  ),
  stringsAsFactors = FALSE
)

# 2. STRATA (1 row per distinct GWAS summary-stat dataset within a study) ====
strata <- data.frame(
  stratum_id = c(
    "suhre_2017_EUR",
    "folk_2017_EUR",
    "sun_2018_EUR",
    "yao_2018_EUR",
    "folkersen_2020_EUR",
    "ferkingstad_2021_EUR",
    "pietzner_2021_EUR",
    "yang_2021_CSF_EUR",
    "yang_2021_BRAIN_EUR",
    "gudjonsson_2022_EUR",
    "zhang_2022_EUR",
    "zhang_2022_AFR",
    "sun_2022_pre_EUR",
    "sun_2023_ALL",
    "sun_2023_EUR",
    "sun_2023_AFR",
    "sun_2023_CSA",
    "sun_2023_EAS",
    "eldjarn_2023_OLINK_EUR",
    "wingo_2023_DLPFC_EUR",
    "schlosser_2023_URINE_EUR",
    "chen_2024_OLINK_EAS",
    "chen_2024_SOMA_EAS",
    "koprulu_2025_FEMALE_EUR",
    "koprulu_2025_MALE_EUR",
    "campbell_2025_MS_EUR",
    "suhre_2025_NPMS_MULTI",
    "lopez_2026_CSF_AMR",
    "patel_2026_SOMA7K_EUR",
    "nakamura_2026_OLINKHT_EAS",
    "nakamura_2026_SOMA11K_EAS"
  ),
  study_id = c(
    "suhre_2017_28240269",
    "folk_2017_28369058",
    "sun_2018_29875488",
    "yao_2018_30148068",
    "folkersen_2020_33067605",
    "ferkingstad_2021_34857953",
    "pietzner_2021_34648354",
    "yang_2021_34239129",
    "yang_2021_34239129",
    "gudjonsson_2022_35078996",
    "zhang_2022_35501419",
    "zhang_2022_35501419",
    "sun_2022_preprint",
    "sun_2023_37794186",
    "sun_2023_37794186",
    "sun_2023_37794186",
    "sun_2023_37794186",
    "sun_2023_37794186",
    "eldjarn_2023_37794188",
    "wingo_2023_brain_ms",
    "schlosser_2023_urine",
    "chen_2024_ckb_dual",
    "chen_2024_ckb_dual",
    "koprulu_2025_sex_pgwas",
    "koprulu_2025_sex_pgwas",
    "campbell_2025_gs_diams",
    "suhre_2025_proteograph",
    "lopez_2026_csf_nulisa_pre",
    "patel_2026_somascan7k_pre",
    "nakamura_2026_olinkht_pre",
    "nakamura_2026_olinkht_pre"
  ),
  cohorts_display = c(
    "KORA F4",
    "IMPROVE",
    "INTERVAL",
    "Framingham Heart Study",
    "SCALLOP Consortium (15 cohorts incl. INTERVAL, KORA, IMPROVE)",
    "deCODE Genetics",
    "Fenland",
    "Knight ADRC; DIAN",
    "Knight ADRC",
    "AGES-Reykjavik",
    "ARIC",
    "ARIC",
    "UK Biobank",
    "UK Biobank",
    "UK Biobank",
    "UK Biobank",
    "UK Biobank",
    "UK Biobank",
    "deCODE Genetics",
    "ROSMAP (DLPFC subset)",
    "GCKD (Urine subset)",
    "China Kadoorie Biobank",
    "China Kadoorie Biobank",
    "UK Biobank; Fenland",
    "UK Biobank; Fenland",
    "Generation Scotland",
    "QMDiab",
    "CHOP Pediatric Neuro Cohort",
    "INTERVAL; CHRIS",
    "BioBank Japan; Nagahama",
    "BioBank Japan (Pilot)"
  ),
  ancestry_group = c(
    "EUR", "EUR", "EUR", "EUR", "EUR",
    "EUR", "EUR", "EUR", "EUR", "EUR",
    "EUR", "AFR", "EUR", "Multi-ancestry", "EUR",
    "AFR", "CSA", "EAS", "EUR", "EUR",
    "EUR", "EAS", "EAS", "EUR", "EUR",
    "EUR", "Multi-ancestry", "AMR", "EUR", "EAS",
    "EAS"
  ),
  ancestry_detail = c(
    "European (German)",
    "European (Multi-centre)",
    "European (British)",
    "European (US Framingham)",
    "European (15 SCALLOP cohorts)",
    "European (Icelandic)",
    "European (British)",
    "European (US ADRC)",
    "European (US ADRC)",
    "European (Icelandic)",
    "European American",
    "African American",
    "European (White British)",
    "Combined (EUR, AFR, CSA, EAS, AMR)",
    "European (British)",
    "African / Afro-Caribbean",
    "Central / South Asian",
    "East Asian",
    "European (Icelandic)",
    "European (US post-mortem)",
    "European (German CKD)",
    "East Asian (Chinese)",
    "East Asian (Chinese)",
    "European (Female only)",
    "European (Male only)",
    "European (Scottish)",
    "Arab, South Asian & Filipino",
    "Admixed American / Hispanic",
    "European (British & Alpine)",
    "East Asian (Japanese)",
    "East Asian (Japanese)"
  ),
  sex = c(
    "Combined", "Combined", "Combined", "Combined", "Combined",
    "Combined", "Combined", "Combined", "Combined", "Combined",
    "Combined", "Combined", "Combined", "Combined", "Combined",
    "Combined", "Combined", "Combined", "Combined", "Combined",
    "Combined", "Combined", "Combined", "Female", "Male",
    "Combined", "Combined", "Combined", "Combined", "Combined",
    "Combined"
  ),
  n_samples = c(
    997L, 3394L, 3301L, 6861L, 30931L,
    35559L, 10708L, 835L, 380L, 5368L,
    7213L, 1871L, 34557L, 54219L, 34557L,
    931L, 920L, 262L, 35892L, 64L,
    280L, 3978L, 3978L, 29110L, 25109L,
    15035L, 1420L, 118L, 13445L, 8420L,
    410L
  ),
  sample_matrix = c(
    "EDTA plasma", "EDTA plasma", "EDTA plasma", "Citrate plasma", "EDTA plasma",
    "EDTA plasma", "EDTA plasma", "CSF", "Brain tissue", "Serum",
    "EDTA plasma", "EDTA plasma", "EDTA plasma", "EDTA plasma", "EDTA plasma",
    "EDTA plasma", "EDTA plasma", "EDTA plasma", "EDTA plasma", "Brain tissue",
    "Urine", "EDTA plasma", "EDTA plasma", "EDTA plasma", "EDTA plasma",
    "EDTA plasma", "EDTA plasma", "CSF", "EDTA plasma", "EDTA plasma",
    "EDTA plasma"
  ),
  platform_category = c(
    "SomaScan", "Olink", "SomaScan", "Luminex", "Olink",
    "SomaScan", "SomaScan", "SomaScan", "SomaScan", "SomaScan",
    "SomaScan", "SomaScan", "Olink", "Olink", "Olink",
    "Olink", "Olink", "Olink", "Olink", "Mass Spectrometry",
    "Olink", "Olink", "SomaScan", "Olink", "Olink",
    "Mass Spectrometry", "Mass Spectrometry", "NULISA", "SomaScan", "Olink",
    "SomaScan"
  ),
  platform_version = c(
    "SomaScan v3 (1.1k)",
    "Olink Target 96 CVD-I",
    "SomaScan v4 Custom (3.6k)",
    "Luminex xMAP CVD (71-plex)",
    "Olink Target 96 CVD-I",
    "SomaScan v4 (4.9k)",
    "SomaScan v4 (5k)",
    "SomaScan v3 (1.3k)",
    "SomaScan v3 (1.3k)",
    "SomaScan v4 Custom (4.8k)",
    "SomaScan v4 (4.9k)",
    "SomaScan v4 (4.9k)",
    "Olink Explore 1536",
    "Olink Explore 3072",
    "Olink Explore 3072",
    "Olink Explore 3072",
    "Olink Explore 3072",
    "Olink Explore 3072",
    "Olink Explore 3072",
    "TMT-MS 16-plex",
    "Olink Target 96 Inflammation/Organ",
    "Olink Explore 3072",
    "SomaScan v4 (4.9k)",
    "Olink Explore 3072",
    "Olink Explore 3072",
    "SWATH / DIA-MS",
    "Seer Proteograph DIA-MS",
    "Alamar NULISAseq CNS 120",
    "SomaScan v4.1 (7k)",
    "Olink Explore HT (5.4k)",
    "SomaScan 11k"
  ),
  n_assays = c(
    1124L, 83L, 3622L, 71L, 90L,
    4907L, 4979L, 1305L, 1305L, 4782L,
    4657L, 4657L, 1463L, 2923L, 2923L,
    2923L, 2923L, 2923L, 2941L, 8412L,
    184L, 2168L, 2168L, 2923L, 2923L,
    439L, 2980L, 124L, 7144L, 5416L,
    10776L
  ),
  n_unique_proteins = c(
    1095L, 83L, 3283L, 71L, 90L,
    4719L, 4775L, 1279L, 1279L, 4137L,
    4435L, 4435L, 1463L, 2923L, 2923L,
    2923L, 2923L, 2923L, 2925L, 6102L,
    181L, 2168L, 2168L, 2923L, 2923L,
    312L, 2410L, 120L, 6267L, 5390L,
    9650L
  ),
  genome_build = c(
    "GRCh37", "GRCh37", "GRCh37", "GRCh37", "GRCh37",
    "GRCh38", "GRCh37", "GRCh38", "GRCh38", "GRCh37",
    "GRCh38", "GRCh38", "GRCh38", "GRCh38", "GRCh38",
    "GRCh38", "GRCh38", "GRCh38", "GRCh38", "GRCh38",
    "GRCh38", "GRCh38", "GRCh38", "GRCh38", "GRCh38",
    "GRCh38", "GRCh38", "GRCh38", "GRCh38", "GRCh38",
    "GRCh38"
  ),
  imputation_panel = c(
    "1000G Phase 1", "1000G Phase 3", "1000G Phase 3 + UK10K", "HRC r1.1", "1000G Phase 3 / HRC",
    "deCODE WGS", "HRC r1.1 + UK10K", "TOPMed r2", "TOPMed r2", "1000G Phase 3",
    "TOPMed r2", "TOPMed r2", "UKB HRC + UK10K", "UKB HRC + UK10K + TOPMed", "UKB HRC + UK10K + TOPMed",
    "UKB HRC + UK10K + TOPMed", "UKB HRC + UK10K + TOPMed", "UKB HRC + UK10K + TOPMed", "deCODE WGS", "WGS",
    "TOPMed r2", "CKB WGS Reference", "CKB WGS Reference", "UKB HRC + UK10K", "UKB HRC + UK10K",
    "HRC r1.1", "TOPMed r3", "TOPMed r3", "TOPMed r2", "BBJ WGS + 1000G",
    "BBJ WGS + 1000G"
  ),
  coverage = c(
    "genome_wide", "genome_wide", "genome_wide", "top_hits_only", "genome_wide",
    "genome_wide", "genome_wide", "genome_wide", "cis_only", "top_hits_only",
    "cis_only", "cis_only", "genome_wide", "genome_wide", "genome_wide",
    "genome_wide", "genome_wide", "genome_wide", "genome_wide", "cis_only",
    "genome_wide", "cis_only", "cis_only", "genome_wide", "genome_wide",
    "genome_wide", "genome_wide", "cis_only", "genome_wide", "genome_wide",
    "cis_only"
  ),
  cis_window_kb = c(
    NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_,
    NA_integer_, NA_integer_, NA_integer_, 1000L, NA_integer_,
    500L, 500L, NA_integer_, NA_integer_, NA_integer_,
    NA_integer_, NA_integer_, NA_integer_, NA_integer_, 1000L,
    NA_integer_, 500L, 500L, NA_integer_, NA_integer_,
    NA_integer_, NA_integer_, 1000L, NA_integer_, NA_integer_,
    500L
  ),
  transformation = c(
    "IRNT of log10(RFU) residuals",
    "IRNT of NPX",
    "IRNT of log(RFU) residuals",
    "IRNT of log(conc)",
    "IRNT of NPX per cohort",
    "IRNT of age/sex-adjusted RFU",
    "IRNT of residuals",
    "log10(RFU) Z-standardized",
    "log10(RFU) Z-standardized",
    "Box-Cox / IRNT of RFU",
    "IRNT of PEER-adjusted log2(RFU)",
    "IRNT of PEER-adjusted log2(RFU)",
    "IRNT of NPX",
    "IRNT of NPX",
    "IRNT of NPX",
    "IRNT of NPX",
    "IRNT of NPX",
    "IRNT of NPX",
    "IRNT of NPX",
    "IRNT of log2(TMT ratio)",
    "IRNT of creatinine-adjusted NPX",
    "IRNT of NPX",
    "IRNT of log2(RFU)",
    "Sex-specific IRNT of NPX",
    "Sex-specific IRNT of NPX",
    "IRNT of log2(MS intensity)",
    "IRNT of peptide LFQ intensity",
    "NPQ (log2 attomolar scale)",
    "IRNT of ANML-normalized RFU",
    "IRNT of Explore HT NPX",
    "IRNT of log2(RFU)"
  ),
  effect_scale = c(
    "SD_IRNT", "SD_IRNT", "SD_IRNT", "SD_IRNT", "SD_IRNT",
    "SD_IRNT", "SD_IRNT", "SD_Z", "SD_Z", "SD_IRNT",
    "SD_IRNT", "SD_IRNT", "SD_IRNT", "SD_IRNT", "SD_IRNT",
    "SD_IRNT", "SD_IRNT", "SD_IRNT", "SD_IRNT", "SD_IRNT",
    "SD_IRNT", "SD_IRNT", "SD_IRNT", "SD_IRNT", "SD_IRNT",
    "SD_IRNT", "SD_IRNT", "NPX_log2", "SD_IRNT", "SD_IRNT",
    "SD_IRNT"
  ),
  covariates = c(
    "Age, sex, BMI, 3 PCs",
    "Age, sex, site, 4 PCs",
    "Age, sex, duration blood draw, 3 PCs, subcohort",
    "Age, sex, cohort generation, 10 PCs",
    "Age, sex, cohort-specific PCs",
    "Age, sex, sample age, kinship LMM",
    "Age, sex, test site, 10 PCs",
    "Age, sex, plate, 2 PCs",
    "Age, sex, PMI, plate, 2 PCs",
    "Age, sex, eGFR, 5 PCs",
    "Age, sex, center, 10 PCs, 90–100 PEER factors",
    "Age, sex, center, 10 PCs, 80 PEER factors",
    "Age, age2, sex, batch, center, array, 20 PCs",
    "Age, age2, sex, age*sex, batch, center, array, time-to-process, 20 PCs",
    "Age, age2, sex, age*sex, batch, center, array, time-to-process, 20 PCs",
    "Age, age2, sex, batch, center, 20 PCs",
    "Age, age2, sex, batch, center, 20 PCs",
    "Age, age2, sex, batch, center, 20 PCs",
    "Age, sex, sample age, kinship LMM",
    "Age at death, sex, PMI, RIN, batch, 3 PCs",
    "Age, sex, eGFR, UACR, 10 PCs",
    "Age, sex, region, fasting time, 10 PCs",
    "Age, sex, region, fasting time, 10 PCs",
    "Age, age2, batch, center, 20 PCs",
    "Age, age2, batch, center, 20 PCs",
    "Age, sex, storage time, GRM (fastGWA)",
    "Age, sex, BMI, study group, 4 PCs",
    "Age, sex, lumbar puncture indication, 5 PCs",
    "Age, sex, season, batch, 10 PCs",
    "Age, sex, biobank center, 10 PCs",
    "Age, sex, biobank center, 10 PCs"
  ),
  gwas_software = c(
    "SNPTEST", "PLINK", "SNPTEST", "FASTGWA", "METAL",
    "BOLT-LMM / deCODE", "BOLT-LMM", "PLINK2", "PLINK2", "PLINK2",
    "PLINK2", "PLINK2", "REGENIE", "REGENIE", "REGENIE",
    "REGENIE", "REGENIE", "REGENIE", "BOLT-LMM / deCODE", "MatrixEQTL / PLINK2",
    "REGENIE", "REGENIE", "REGENIE", "REGENIE / METAL", "REGENIE / METAL",
    "GCTA-fastGWA", "REGENIE", "REGENIE", "REGENIE", "REGENIE",
    "REGENIE"
  ),
  file_format = c(
    "TSV.gz", "TSV.gz", "TSV.gz (per chr)", "CSV (top hits)", "TSV.gz",
    "TXT.gz + Variant Anno", "TXT.gz", "TXT.gz", "TXT.gz", "XLSX (supplement)",
    "PLINK2 .glm.linear", "PLINK2 .glm.linear", "TAR of space-delimited", "TAR of space-delimited", "TAR of space-delimited",
    "TAR of space-delimited", "TAR of space-delimited", "TAR of space-delimited", "TXT.gz + Variant Anno", "TSV.gz",
    "GWAS-SSF TSV.gz", "TSV.gz", "TSV.gz", "GWAS-SSF TSV.gz", "GWAS-SSF TSV.gz",
    "GWAS-SSF TSV.gz", "GWAS-SSF TSV.gz", "TSV.gz", "GWAS-SSF TSV.gz", "TAR of REGENIE .gz",
    "TSV.gz"
  ),
  has_beta_se = c(
    TRUE, TRUE, TRUE, TRUE, TRUE,
    TRUE, TRUE, TRUE, TRUE, TRUE,
    TRUE, TRUE, TRUE, TRUE, TRUE,
    TRUE, TRUE, TRUE, TRUE, TRUE,
    TRUE, TRUE, TRUE, TRUE, TRUE,
    TRUE, TRUE, TRUE, TRUE, TRUE,
    TRUE
  ),
  has_eaf = c(
    TRUE, TRUE, TRUE, TRUE, TRUE,
    TRUE, TRUE, FALSE, FALSE, TRUE,
    TRUE, TRUE, TRUE, TRUE, TRUE,
    TRUE, TRUE, TRUE, TRUE, TRUE,
    TRUE, TRUE, TRUE, TRUE, TRUE,
    TRUE, TRUE, TRUE, TRUE, TRUE,
    TRUE
  ),
  indel_encoding = c(
    "snps_only", "snps_only", "sequence", "snps_only", "D_I_codes",
    "sequence_plus_multiallelic", "D_I_codes", "snps_only", "snps_only", "snps_only",
    "sequence", "sequence", "sequence", "sequence", "sequence",
    "sequence", "sequence", "sequence", "sequence_plus_multiallelic", "sequence",
    "sequence", "sequence", "sequence", "sequence", "sequence",
    "sequence", "sequence", "sequence", "sequence", "sequence",
    "sequence"
  ),
  gwasprocessing_yaml = c(
    "", "", "", "", "",
    "ferkingstad_2021_34857953.yaml", "pietzner_2021_34648354.yaml", "", "", "",
    "zhang_2022_35501419.yaml", "zhang_2022_35501419.yaml", "", "sun_2023_37794186.yaml", "sun_2023_37794186.yaml",
    "sun_2023_37794186.yaml", "sun_2023_37794186.yaml", "sun_2023_37794186.yaml", "", "",
    "", "", "", "", "",
    "", "", "", "", "",
    ""
  ),
  harmonisation_notes = c(
    "Early 1.1k SomaScan panel; hg19 SNPs only.",
    "CVD-I panel across IMPROVE.",
    "Per-chromosome bgzipped files; full genome-wide summary statistics.",
    "Individual genotypes in dbGaP; public supplement contains significant pQTLs only.",
    "D/I indel alleles across 15 SCALLOP meta-analysed cohorts.",
    "Requires joining deCODE assocvariants.annotated.txt.gz for EAF and '!' multiallelic flags.",
    "Indels coded as D/I without sequence; resolve via dbSNP rsID lookup.",
    "MAF provided without signed EAF in public files; check palindromic SNPs carefully.",
    "Cis-window (+/- 1 Mb) brain homogenate pQTLs only.",
    "Full genome-wide serum files available only via collaboration; supplement has 10k+ top hits.",
    "Cis-windows (+/- 500 kb) only; OA derived from REF/ALT + A1.",
    "Cis-windows (+/- 500 kb) in African American ARIC participants; valuable fine-mapping LD.",
    "Superseded by Sun 2023 Explore 3072 release.",
    "Per-protein .tar archives split by chromosome; LOG10P reported instead of raw P.",
    "European discovery + replication subset of UKB-PPP.",
    "African ancestry subset of UKB-PPP; check MAF >= 0.01 filter for rare variants.",
    "Central/South Asian subset of UKB-PPP.",
    "East Asian subset of UKB-PPP (N = 262; no sample-size floor applied).",
    "Same deCODE cohort as Ferkingstad 2021, assayed on Olink Explore 3072.",
    "Small-N post-mortem DLPFC brain TMT-MS pGWAS (N = 64); cis +/- 1 Mb.",
    "Urinary Olink proteomics adjusted for eGFR and creatinine.",
    "Matched CKB cohort assayed on Olink Explore 3072 (cis +/- 500 kb).",
    "Matched CKB cohort assayed on SomaScan v4 (cis +/- 500 kb); ideal for cross-platform PAV check.",
    "Female-only meta-analysis (UKB + Fenland); watch sample overlap with Sun 2023 & Pietzner 2021.",
    "Male-only meta-analysis (UKB + Fenland); watch sample overlap with Sun 2023 & Pietzner 2021.",
    "Untargeted SWATH/DIA mass spectrometry in Generation Scotland; peptide-to-protein rollup.",
    "Nanoparticle corona enrichment MS; measures deep intracellular & secreted plasma proteins.",
    "Pediatric CSF NULISAseq CNS panel (N = 118 preprint); unstandardized log2 NPQ units.",
    "Expanded SomaScan 7k preprint combining INTERVAL + CHRIS; overlaps Sun 2018 INTERVAL.",
    "Olink Explore HT 5.4k preprint in BioBank Japan + Nagahama.",
    "Pilot SomaScan 11k cis-pQTL arm (N = 410) in BioBank Japan."
  ),
  stringsAsFactors = FALSE
)

# 3. COHORT MEMBERSHIP (Normalized table for exact cohort-overlap checking) ====
cohort_membership <- rbind(
  data.frame(study_id = "suhre_2017_28240269",      cohort = "KORA",                    n_cohort = 997L),
  data.frame(study_id = "folk_2017_28369058",       cohort = "IMPROVE",                 n_cohort = 3394L),
  data.frame(study_id = "sun_2018_29875488",        cohort = "INTERVAL",                n_cohort = 3301L),
  data.frame(study_id = "yao_2018_30148068",        cohort = "Framingham Heart Study",  n_cohort = 6861L),
  data.frame(study_id = "folkersen_2020_33067605",  cohort = "INTERVAL",                n_cohort = 4890L),
  data.frame(study_id = "folkersen_2020_33067605",  cohort = "KORA",                    n_cohort = 1064L),
  data.frame(study_id = "folkersen_2020_33067605",  cohort = "IMPROVE",                 n_cohort = 3321L),
  data.frame(study_id = "folkersen_2020_33067605",  cohort = "ORCADES",                 n_cohort = 1980L),
  data.frame(study_id = "folkersen_2020_33067605",  cohort = "Other SCALLOP (11)",      n_cohort = 19676L),
  data.frame(study_id = "ferkingstad_2021_34857953", cohort = "deCODE Genetics",        n_cohort = 35559L),
  data.frame(study_id = "pietzner_2021_34648354",   cohort = "Fenland",                 n_cohort = 10708L),
  data.frame(study_id = "yang_2021_34239129",       cohort = "Knight ADRC",             n_cohort = 980L),
  data.frame(study_id = "yang_2021_34239129",       cohort = "DIAN",                    n_cohort = 235L),
  data.frame(study_id = "gudjonsson_2022_35078996", cohort = "AGES-Reykjavik",          n_cohort = 5368L),
  data.frame(study_id = "zhang_2022_35501419",      cohort = "ARIC",                    n_cohort = 9084L),
  data.frame(study_id = "sun_2022_preprint",        cohort = "UK Biobank",              n_cohort = 34557L),
  data.frame(study_id = "sun_2023_37794186",        cohort = "UK Biobank",              n_cohort = 54219L),
  data.frame(study_id = "eldjarn_2023_37794188",    cohort = "deCODE Genetics",         n_cohort = 35892L),
  data.frame(study_id = "wingo_2023_brain_ms",      cohort = "ROSMAP",                  n_cohort = 64L),
  data.frame(study_id = "schlosser_2023_urine",     cohort = "GCKD",                    n_cohort = 280L),
  data.frame(study_id = "chen_2024_ckb_dual",       cohort = "China Kadoorie Biobank",  n_cohort = 3978L),
  data.frame(study_id = "koprulu_2025_sex_pgwas",   cohort = "UK Biobank",              n_cohort = 43511L),
  data.frame(study_id = "koprulu_2025_sex_pgwas",   cohort = "Fenland",                 n_cohort = 10708L),
  data.frame(study_id = "campbell_2025_gs_diams",   cohort = "Generation Scotland",     n_cohort = 15035L),
  data.frame(study_id = "suhre_2025_proteograph",   cohort = "QMDiab",                  n_cohort = 1420L),
  data.frame(study_id = "lopez_2026_csf_nulisa_pre", cohort = "CHOP Pediatric Neuro",   n_cohort = 118L),
  data.frame(study_id = "patel_2026_somascan7k_pre", cohort = "INTERVAL",               n_cohort = 8890L),
  data.frame(study_id = "patel_2026_somascan7k_pre", cohort = "CHRIS",                  n_cohort = 4555L),
  data.frame(study_id = "nakamura_2026_olinkht_pre", cohort = "BioBank Japan",          n_cohort = 6200L),
  data.frame(study_id = "nakamura_2026_olinkht_pre", cohort = "Nagahama",               n_cohort = 2220L)
)
cohort_membership$study_id <- as.character(cohort_membership$study_id)
cohort_membership$cohort <- as.character(cohort_membership$cohort)

# 4. PROTEIN & ASSAY INDEX (Cross-platform protein lookup simulation) ====
protein_catalog <- data.frame(
  uniprot = c(
    "P05231", "Q8NBP7", "P01375", "P10145", "P02649", "Q9NZC2", "P41159", "P02741",
    "Q99988", "Q99523", "P08887", "P39900", "P09237", "P01579", "P13500", "P01019",
    "P03951", "P00734", "P01009", "Q03405", "P01033", "P19320", "P05362", "P10636",
    "P07196", "P14136", "P37840", "Q15389", "P08519", "P16860", "P01034", "Q96D42"
  ),
  gene_symbol = c(
    "IL6", "PCSK9", "TNF", "CXCL8", "APOE", "TREM2", "LEP", "CRP",
    "GDF15", "SORT1", "IL6R", "MMP12", "MMP7", "IFNG", "CCL2", "AGT",
    "F11", "F2", "SERPINA1", "PLAUR", "TIMP1", "VCAM1", "ICAM1", "MAPT",
    "NEFL", "GFAP", "SNCA", "ANGPTL3", "LPA", "NPPB", "CST3", "HAVCR1"
  ),
  protein_name = c(
    "Interleukin-6",
    "Proprotein convertase subtilisin/kexin type 9",
    "Tumor necrosis factor",
    "Interleukin-8",
    "Apolipoprotein E",
    "Triggering receptor expressed on myeloid cells 2",
    "Leptin",
    "C-reactive protein",
    "Growth/differentiation factor 15",
    "Sortilin",
    "Interleukin-6 receptor subunit alpha",
    "Macrophage metalloelastase",
    "Matrilysin",
    "Interferon gamma",
    "C-C motif chemokine 2",
    "Angiotensinogen",
    "Coagulation factor XI",
    "Prothrombin (Coagulation factor II)",
    "Alpha-1-antitrypsin",
    "Urokinase plasminogen activator surface receptor",
    "Metalloproteinase inhibitor 1",
    "Vascular cell adhesion protein 1",
    "Intercellular adhesion molecule 1",
    "Microtubule-associated protein tau",
    "Neurofilament light polypeptide",
    "Glial fibrillary acidic protein",
    "Alpha-synuclein",
    "Angiopoietin-related protein 3",
    "Apolipoprotein(a)",
    "Natriuretic peptides B (BNP)",
    "Cystatin-C",
    "Hepatitis A virus cellular receptor 1 (KIM-1)"
  ),
  panel_class = c(
    "Inflammation", "Cardiometabolic", "Inflammation", "Inflammation", "Neuro / Lipid", "Neurology", "Cardiometabolic", "Abundant Plasma",
    "Cardiometabolic", "Cardiometabolic", "Inflammation", "Oncology / CVD", "Cardiometabolic", "Inflammation", "Inflammation", "Abundant Plasma",
    "Coagulation", "Abundant Plasma", "Abundant Plasma", "Inflammation", "Cardiometabolic", "Cardiometabolic", "Cardiometabolic", "CNS",
    "CNS", "CNS", "CNS", "Cardiometabolic", "Abundant Plasma", "Cardiometabolic", "Renal / Abundant", "Renal"
  ),
  soma_seqid = c(
    "4673-13", "5231-79;13104-55", "5936-53", "3447-64", "2418-55;2937-10", "11851-21;16300-4", "2575-5", "4337-49",
    "4374-45", "4929-55", "4139-71", "4496-60", "4924-32", "15346-31", "2578-67", "4880-21",
    "2797-56", "4157-2", "3580-25", "2652-15", "2211-9", "2967-8", "4342-10", "4878-3",
    "10082-7", "15509-2", "13666-190", "10391-1", "4990-87", "3723-1", "2609-59", "9283-8"
  ),
  olink_id = c(
    "OID20101;OID20563;OID20911;OID21276", "OID20235", "OID20074;OID20555", "OID20153;OID20631", "OID20918", "OID20844", "OID20190", "",
    "OID20251", "OID20312", "OID20611", "OID21439", "OID20087", "OID20495", "OID21004", "",
    "OID20389", "", "", "OID20690", "OID20219", "OID20301", "OID20305", "OID30912",
    "OID20871", "OID20899", "OID31044", "OID20288", "OID20341", "OID20112", "OID20355", "OID20188"
  ),
  pav_sensitive = c(
    FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE,
    TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE,
    FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE,
    FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE
  ),
  stringsAsFactors = FALSE
)

# Deterministically map proteins across strata according to platform & matrix realism
assay_rows <- list()
idx <- 1L

for (i in seq_len(nrow(strata))) {
  st <- strata[i, ]
  plat <- st$platform_category
  ver <- st$platform_version
  mat <- st$sample_matrix

  for (j in seq_len(nrow(protein_catalog))) {
    pr <- protein_catalog[j, ]
    cls <- pr$panel_class

    included <- FALSE
    assay_ids <- character(0)

    if (plat == "SomaScan") {
      sids <- strsplit(pr$soma_seqid, ";", fixed = TRUE)[[1]]
      if (grepl("1\\.1k|1\\.3k", ver)) {
        included <- cls != "CNS" || pr$gene_symbol %in% c("MAPT", "APOE", "TREM2")
        assay_ids <- sids[1]
      } else if (grepl("3\\.6k", ver)) {
        included <- pr$gene_symbol != "GFAP"
        assay_ids <- sids[1]
      } else {
        included <- TRUE
        assay_ids <- sids
      }
    } else if (plat == "Olink") {
      oids <- strsplit(pr$olink_id, ";", fixed = TRUE)[[1]]
      oids <- oids[nzchar(oids)]
      if (length(oids) > 0) {
        if (grepl("CVD-I", ver)) {
          included <- pr$gene_symbol %in% c("IL6", "CXCL8", "GDF15", "IL6R", "MMP12", "MMP7", "CCL2", "PLAUR", "TIMP1", "NPPB", "HAVCR1", "LEP")
          assay_ids <- oids[1]
        } else if (grepl("Target 96", ver)) {
          included <- pr$gene_symbol %in% c("IL6", "TNF", "CXCL8", "GDF15", "MMP7", "CCL2", "PLAUR", "CST3", "HAVCR1", "VCAM1")
          assay_ids <- oids[1]
        } else if (grepl("1536", ver)) {
          included <- !pr$gene_symbol %in% c("MAPT", "SNCA")
          assay_ids <- oids[1]
        } else if (grepl("3072", ver)) {
          included <- !pr$gene_symbol %in% c("MAPT", "SNCA")
          assay_ids <- oids
        } else if (grepl("HT", ver)) {
          included <- TRUE
          assay_ids <- oids
        }
      } else if (grepl("HT", ver) && pr$gene_symbol %in% c("CRP", "AGT", "F2", "SERPINA1")) {
        included <- TRUE
        assay_ids <- paste0("OID45", sprintf("%03d", j))
      }
    } else if (plat == "Mass Spectrometry") {
      if (grepl("SWATH|DIA-MS", ver) && !grepl("Proteograph", ver)) {
        # Classical neat-plasma MS captures high/medium abundance plasma proteins
        included <- cls %in% c("Abundant Plasma", "Coagulation", "Renal / Abundant") || pr$gene_symbol %in% c("APOE", "PCSK9", "VCAM1", "ICAM1", "TIMP1")
        assay_ids <- paste0("MS_", pr$uniprot)
      } else if (grepl("Proteograph", ver)) {
        # Nanoparticle enrichment reaches lower-abundance + intracellular proteins
        included <- cls != "CNS" || pr$gene_symbol %in% c("APOE", "MAPT", "SNCA")
        assay_ids <- paste0("NPMS_", pr$uniprot)
      } else if (grepl("TMT", ver) && mat == "Brain tissue") {
        included <- cls %in% c("CNS", "Neuro / Lipid", "Neurology") || pr$gene_symbol %in% c("SORT1", "TIMP1", "ICAM1", "GDF15")
        assay_ids <- paste0("TMT_", pr$uniprot)
      }
    } else if (plat == "NULISA") {
      included <- cls %in% c("CNS", "Neuro / Lipid", "Neurology", "Inflammation")
      assay_ids <- paste0("NULISA_", pr$gene_symbol)
    } else if (plat == "Luminex") {
      included <- pr$gene_symbol %in% c("IL6", "CRP", "LEP", "GDF15", "MMP12", "TIMP1", "VCAM1", "ICAM1", "LPA", "NPPB", "CST3")
      assay_ids <- paste0("LUM_", pr$gene_symbol)
    }

    if (included && length(assay_ids) > 0) {
      for (aid in assay_ids) {
        # Simulate realistic cis-pQTL detection probability scaling with sample size
        has_cis <- (st$n_samples >= 1500L && !pr$gene_symbol %in% c("MAPT", "SNCA")) ||
          (st$n_samples >= 250L && pr$gene_symbol %in% c("PCSK9", "IL6R", "SORT1", "APOE", "LPA", "SERPINA1", "F11", "PLAUR", "MMP12", "TREM2", "HAVCR1")) ||
          (mat %in% c("CSF", "Brain tissue") && pr$gene_symbol %in% c("APOE", "TREM2", "MAPT", "GFAP"))

        assay_rows[[idx]] <- data.frame(
          uniprot = pr$uniprot,
          gene_symbol = pr$gene_symbol,
          protein_name = pr$protein_name,
          assay_id = aid,
          stratum_id = st$stratum_id,
          study_id = st$study_id,
          cis_pqtl_reported = has_cis,
          pav_epitope_flag = pr$pav_sensitive && plat %in% c("SomaScan", "Olink", "NULISA", "Luminex"),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
}
proteins_index <- do.call(rbind, assay_rows)

# 5. AUTOMATED SEARCH RUNS LOG (PRISMA & Scheduled Action History) ====
search_runs <- data.frame(
  run_id = c(
    "run_2026_07_01",
    "run_2026_08_01",
    "run_2026_09_01",
    "run_2026_10_01",
    "run_2026_11_15"
  ),
  run_date = c(
    "2026-07-01",
    "2026-08-01",
    "2026-09-01",
    "2026-10-01",
    "2026-11-15"
  ),
  trigger = c(
    "Initial seed + retrospective sweep (2015–2026)",
    "Scheduled GitHub Action (cron)",
    "Scheduled GitHub Action (cron)",
    "Scheduled GitHub Action (cron)",
    "Scheduled GitHub Action (cron)"
  ),
  sources_queried = c(
    "PubMed; Europe PMC; bioRxiv/medRxiv; EBI GWAS Catalog; OpenAlex Snowball",
    "PubMed; Europe PMC; bioRxiv/medRxiv; EBI GWAS Catalog; OpenAlex Snowball",
    "PubMed; Europe PMC; bioRxiv/medRxiv; EBI GWAS Catalog; OpenAlex Snowball",
    "PubMed; Europe PMC; bioRxiv/medRxiv; EBI GWAS Catalog; OpenAlex Snowball",
    "PubMed; Europe PMC; bioRxiv/medRxiv; EBI GWAS Catalog; OpenAlex Snowball"
  ),
  raw_records = c(2840L, 164L, 192L, 178L, 211L),
  deduplicated_new = c(1910L, 68L, 81L, 74L, 89L),
  auto_flagged_secondary_mr = c(1685L, 59L, 70L, 66L, 78L),
  candidates_for_review = c(225L, 9L, 11L, 8L, 11L),
  studies_added = c(19L, 0L, 0L, 0L, 3L),
  preprints_upgraded = c(1L, 0L, 0L, 0L, 0L),
  pr_number = c("PR #1", "PR #4", "PR #7", "PR #11", "PR #14"),
  release_tag = c("v0.1.0-sim", "", "", "", "v0.2.0-sim"),
  status = c("Merged & Released", "Merged (0 primary)", "Merged (0 primary)", "Merged (0 primary)", "Open Triage / v0.2.0-sim"),
  stringsAsFactors = FALSE
)

# 6. CANDIDATE TRIAGE QUEUE (Simulated automated PR candidates) ====
candidate_queue <- data.frame(
  candidate_id = c(
    "cand_2026_11_01",
    "cand_2026_11_02",
    "cand_2026_11_03",
    "cand_2026_11_04",
    "cand_2026_11_05",
    "cand_2026_11_06",
    "cand_2026_11_07",
    "cand_2026_11_08"
  ),
  run_id = rep("run_2026_11_15", 8),
  source = c(
    "bioRxiv", "medRxiv", "bioRxiv", "Nature Genetics", "PubMed", "PubMed", "Europe PMC", "medRxiv"
  ),
  identifier = c(
    "10.1101/2026.02.19.581022",
    "10.1101/2026.04.12.26305118",
    "10.1101/2026.09.03.611408",
    "10.1038/s41588-026-09912-1",
    "PMID:41488603",
    "PMID:41398839",
    "PMID:41782372",
    "10.1101/2026.11.02.26311904"
  ),
  detected_platform = c(
    "SomaScan v4.1 (7k)",
    "Alamar NULISAseq CNS 120",
    "Olink Explore HT / SomaScan 11k",
    "SomaScan v4.1 (7k)",
    "SomaScan / Olink (Secondary)",
    "Olink Explore (Secondary)",
    "SomaScan v4 (EWAS)",
    "DIA-MS (Saliva)"
  ),
  detected_n = c("13,445", "118", "8,420", "13,445", "N/A (2SMR)", "N/A (2SMR)", "2,410", "92"),
  title = c(
    "Beyond the classical plasma secretome: expanded SomaScan 7k genetic architecture in 13,445 Europeans",
    "Pediatric cerebrospinal fluid pQTL mapping using attomolar NULISAseq immunoassays",
    "Multi-ancestry pGWAS of 5,400 plasma proteins using Olink Explore HT and SomaScan 11k Pilot",
    "Genetic architecture of the expanded human plasma secretome in 13,445 individuals",
    "Identification of potential therapeutic targets for stroke and its subtypes by integrating proteomes and genetics",
    "Identification of therapeutic targets for insomnia through genetic and plasma proteomic approaches",
    "Epigenome-wide association study of blood proteome in the Atherosclerosis Risk in Communities study",
    "Pilot pQTL mapping of the human salivary proteome by DIA mass spectrometry in 92 twins"
  ),
  auto_classification = c(
    "Likely Primary pGWAS",
    "Likely Primary pGWAS (Small N)",
    "Likely Primary pGWAS",
    "Preprint -> Published Match (patel_2026_somascan7k_pre)",
    "Auto-Deprioritised: Secondary 2SMR",
    "Auto-Deprioritised: Secondary 2SMR",
    "Auto-Deprioritised: Non-GWAS (EWAS)",
    "Needs Human Check: Novel Primary pGWAS"
  ),
  curation_status = c(
    "Accepted -> v0.2.0-sim",
    "Accepted -> v0.2.0-sim",
    "Accepted -> v0.2.0-sim",
    "Pending Curator Review (PR #14)",
    "Excluded (Secondary MR)",
    "Excluded (Secondary MR)",
    "Excluded (DNA methylation EWAS)",
    "Pending Curator Review (PR #14)"
  ),
  pmid = c("", "", "", "41990112", "41488603", "41398839", "41782372", ""),
  doi = c(
    "10.1101/2026.02.19.581022",
    "10.1101/2026.04.12.26305118",
    "10.1101/2026.09.03.611408",
    "10.1038/s41588-026-09912-1",
    "10.1161/strokeaha.126.041111",
    "10.1038/s41398-025-03110-2",
    "10.1038/s41467-026-44100-1",
    "10.1101/2026.11.02.26311904"
  ),
  pmcid = c("", "", "", "PMC11882011", "", "", "PMC11900122", ""),
  epmc_id = c("PPR811022", "PPR830511", "PPR891408", "41990112", "41488603", "41398839", "41782372", "PPR911904"),
  matched_sources = c(
    "EuropePMC_PPR;bioRxiv",
    "EuropePMC_PPR;medRxiv",
    "EuropePMC_PPR;Snowball",
    "PubMed;EuropePMC_MED;Snowball",
    "PubMed;EuropePMC_MED",
    "PubMed;EuropePMC_MED",
    "PubMed;EuropePMC_MED",
    "EuropePMC_PPR;medRxiv"
  ),
  p_score = c(21L, 16L, 22L, 23L, 0L, 0L, 0L, 14L),
  s_score = c(0L, 0L, 0L, 0L, 17L, 17L, 0L, 0L),
  preprint_lineage_note = c(
    "", "", "", "Matches active preprint patel_2026_somascan7k_pre (10.1101/2026.02.19.581022)",
    "", "", "", ""
  ),
  stringsAsFactors = FALSE
)

# 7. TAGGED RELEASES CHANGELOG (Starts empty until the first real release is cut) ====
releases <- data.frame(
  version = character(0),
  release_date = character(0),
  doi_badge = character(0),
  n_studies = integer(0),
  n_strata = integer(0),
  n_preprints_active = integer(0),
  summary = character(0),
  stringsAsFactors = FALSE
)

# 8. EXCLUDED RECORDS LEDGER (Prevents re-screening known secondary 2SMR / EWAS) ====
excluded_records <- data.frame(
  identifier = c("PMID:41488603", "PMID:41398839", "PMID:41782372"),
  reason = c(
    "Secondary two-sample MR study (stroke)",
    "Secondary two-sample MR study (insomnia)",
    "DNA methylation EWAS, not genetic GWAS"
  ),
  screened_date = c("2026-11-15", "2026-11-15", "2026-11-15"),
  stringsAsFactors = FALSE
)

write_tsv(studies, "studies.tsv")
write_tsv(strata, "strata.tsv")
write_tsv(cohort_membership, "cohort_membership.tsv")
write_tsv(proteins_index, "proteins_index.tsv")
write_tsv(search_runs, "search_runs.tsv")
write_tsv(candidate_queue, "candidate_queue.tsv")
write_tsv(releases, "releases.tsv")
write_tsv(excluded_records, "excluded_records.tsv")
