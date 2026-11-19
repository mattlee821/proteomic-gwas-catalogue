
**Project:** `proteomic-gwas-catalogue`  
**Document Type:** Living Methodological Log (updated at every search/extraction version milestone)  
**Last Updated:** 2026-11-19 (`v0` -> `v1` -> `v2` -> `v3.0` documented; Round-2 `v3` replication in progress)

---

## 1. Purpose, Eligibility Scope, and Methodological Principles

### 1.1 Scientific Objective
Published genome-wide association studies (GWAS) of protein levels — protein quantitative trait locus (`pQTL`) studies — span multiple assay technologies, cohorts, biofluids, journals, preprint servers, and summary-statistics portals. Meanwhile, since 2022, biomedical literature databases have experienced an exponential surge in **secondary two-sample Mendelian randomisation (`2SMR`), summary-data-based MR (`SMR`), and proteome-wide association (`PWAS`) papers** that reuse public summary statistics (such as UKB-PPP, deCODE, Fenland, INTERVAL, ARIC, and SCALLOP) without generating new primary genetic associations. In naive keyword searches, secondary summary-statistics consumers now outnumber primary pGWAS papers by approximately **10:1 to 30:1**.

The goal of our search and extraction methodology is to achieve **near-100% recall of primary human protein GWAS** (including preprints and small-sample tissue/biofluid studies) while compressing curator screening noise by **~85–90%** and preserving strict auditability.

### 1.2 Study Eligibility Criteria
A publication or preprint is eligible for inclusion if it satisfies all four criteria:
1. **Human population (`N >= 1`, no sample-size threshold):** Any human cohort, biobank, clinical series, twin sample, post-mortem tissue bank, or donor-derived primary/iPSC cell panel across any ancestry or sex. Small-sample studies (e.g. `N = 50–150` brain homogenate, CSF, urine, saliva, or iPSC studies) are explicitly retained.
2. **Multiplexed / proteomic measurement (`>= 5` protein targets or explicit proteomic framing):** Covers aptamer affinity (`SomaScan` / `SomaLogic` 1.1k–11k), proximity extension and nucleic-acid immunoassays (`Olink Target`, `SCALLOP`, `Olink Explore 1536/3072/HT`, `Alamar NULISAseq`), untargeted and nanoparticle-enriched mass spectrometry (`DIA-MS`, `SWATH-MS`, `DDA`, `TMT`, `iTRAQ`, `Orbitrap Astral`, `timsTOF`, `Seer Proteograph`, `MRM/PRM`), and multiplex bead/planar immunoassays (`Luminex xMAP`, `Bio-Plex`, `Milliplex`, `Rules-Based Medicine`, `Meso Scale Discovery`, `Simoa`). Single-analyte clinical chemistry GWAS (`< 5` proteins, e.g. standalone clinical CRP or Lp(a) GWAS) are excluded.
3. **All human biological matrices:** EDTA/citrate/heparin/unspecified plasma, serum, whole blood / dried blood spots, cerebrospinal fluid (`CSF`), brain and solid tissues (kidney, liver, heart, vascular), urine, saliva, and human cell lines.
4. **Primary germline genetic association mapping:** Generates new genotype-to-protein association results (array GWAS, whole-exome or whole-genome rare-variant single-variant/burden pGWAS, *cis*-only pQTL mapping, variance-pQTL [`vQTL`] mapping, parent-of-origin or sex-stratified pGWAS, or primary multi-cohort consortium meta-analyses), regardless of whether full summary statistics are open, portal-gated, top-hits-only, or unavailable. Active preprints (`bioRxiv`, `medRxiv`, `Research Square`, `openRxiv`, `Preprints.org`) are included and linked via `superseded_by` once peer-reviewed.

**Exclusions:** Secondary `2SMR`/`SMR`/`PWAS`/colocalisation-only studies that exclusively reuse existing summary statistics or pre-trained weights; pure DNA-methylation `EWAS` without germline pQTLs; pedigree-only twin heritability without genotyping; targeted candidate-SNP (`< 100` SNPs) lookups; non-human-only studies; and narrative reviews/editorials.

### 1.3 Core Validation Principles
1. **Strict Blinding to External Positive Controls:** An external benchmark list (Karsten Suhre's curated proteomics GWAS table on *metabolomix.com*) is held back as an air-gapped positive control (`data/positive_controls.tsv`, git-ignored and programmatically blocked from being read by search scripts via `src/validate_catalogue.R`). Neither human prompts nor subagents are permitted to inspect *metabolomix.com* during strategy design.
2. **Adversarial Multi-Agent Design & Peer Review:** Specialist subagents independently design search and extraction modules, while separate methodological reviewer subagents audit code, edge cases, memory footprints, and live API behavior.
3. **Blinded Dual-Agent Replication & Capture-Recapture Estimation:** Independent search-scientist subagents operating in isolated empty-catalogue sandboxes audit and improve the search strategy over fixed calendar windows (`2025-01-01` to `2025-12-31`), adjudicated by an independent biostatistical referee subagent using ground-truth full-text checks, Cohen's $\kappa$, and Chapman capture-recapture sensitivity estimation.

---

## 2. Prototype Baseline (`v0`) — Initial Single-Source Keyword Search

### 2.1 Architecture of `v0`
During initial repository scaffolding, a minimal prototype script (`src/run_literature_search.R`, `v0`) queried the Europe PMC REST API `GET /search` endpoint using a single unstructured Boolean string:
```text
("SomaScan" OR "SomaLogic" OR "Olink" OR "NULISA" OR "NULISAseq" OR "SCALLOP consortium"
 OR "mass spectrometry proteomics" OR "DIA-MS" OR "plasma proteome" OR "CSF proteome" OR "urinary proteome")
AND ("genome-wide association" OR "GWAS" OR "pQTL" OR "protein quantitative trait")
```
Triage relied on simple title-keyword subtraction: any article with `"Mendelian randomization"`, `"proteome-wide association study"`, or `"drug target"` in its title was auto-deprioritised unless its abstract matched five literal phrases (`"we measured"`, `"genotyped and measured"`, etc.).

### 2.2 Identified Deficiencies in `v0`
1. **Unqualified Full-Text & Reference-List Noise in Europe PMC:** Unqualified keywords in Europe PMC search full text *and* bibliography titles across >9M open-access PMC articles, inflating false positives ~10-fold (as seen in the earlier `analysis/somalogic-search/` pilot where >85% of 64 candidates were secondary MR/PWAS papers, EWASs, or unrelated MS papers).
2. **False Demotion of Translational Primary pGWAS:** Title-keyword subtraction falsely demoted landmark primary pGWAS whose titles mention translational applications, such as **Folkersen et al. 2020** (*"...drug target evaluation of 90 cardiovascular proteins in 30,931 individuals"*) and **Zhang et al. 2022 ARIC** (*"...identify cis-pQTLs and models for proteome-wide association studies"*).
3. **Stealth 2SMR Leakage:** ~40% of recent 2SMR papers omit `"Mendelian randomization"` from their title (*"Proteogenomic integration identifies causal plasma proteins for..."*), bypassing title-only subtraction.
4. **Identifier Loss at Extraction Handoff:** `candidate_queue.tsv` collapsed records into a single `identifier` string, discarding `PMID`, `PMCID`, and Europe PMC `PPR*` preprint links needed for downstream full-text and supplement extraction.

---

## 3. Version 1 (`v1.0`) — Specialist Subagent Architecture & 2-Paper Live Benchmark

### 3.1 Subagent Design Setup
Two specialist subagents were defined and executed in parallel without broad literature sweeps, restricted strictly to targeted empirical verification against two flagship benchmark studies:
- **Benchmark Study A:** **Sun et al. 2023 (UKB-PPP Olink)** — *Nature* (`PMID: 37794186`, `DOI: 10.1038/s41586-023-06592-6`, `PMC10567551`).
- **Benchmark Study B:** **Ferkingstad et al. 2021 (deCODE SomaScan)** — *Nature Genetics* (`PMID: 34857953`, `DOI: 10.1038/s41588-021-00978-w`).

Subagents:
1. **`search-strategy-designer`:** Designed the v1 4-source search strategy, connected-components deduplication, 4-tier preprint-to-published upgrade detector, Two-Stage Dual-Axis (`P_score` vs `S_score`) triage rubric, and blinded positive-control evaluator (`design_search_v1.md`).
2. **`paper-extraction-designer`:** Designed and live-tested the 4-tier cascading retrieval and extraction pipeline (`Structured APIs -> Legal OA Full-Text XML/BioC -> Supplementary Spreadsheets/Docs & Portal Preflight -> Provenance & NEEDS_MANUAL_REVIEW Escalation`, `design_extraction_v1.md`).

### 3.2 Critical Live Empirical Discoveries from the 2-Paper Benchmark
Testing live API responses on `PMID: 37794186` and `PMID: 34857953` uncovered **seven critical real-world behaviors** that refuted standard assumptions:

| # | Empirical Discovery | Evidence on `Sun 2023 (37794186)` | Evidence on `Ferkingstad 2021 (34857953)` | Architectural Consequence |
|---:|---|---|---|---|
| 1 | **Platform brand names & `"GWAS"` absent from flagship abstracts** | Abstract has **`0`** occurrences of `"Olink"` and **`0`** occurrences of `"GWAS"` / `"genome-wide association"` (uses `"Pharma Proteomics Project"` & `"pQTL"`). | Abstract has **`0`** occurrences of `"SomaScan"` / `"SomaLogic"` (uses `"4,907 aptamers"` and plural `"GWASs"`). | Platform-brand queries alone miss both papers; modular proteome/reagent blocks + plural `"GWASs"`/`"pQTLs"` are mandatory. |
| 2 | **`OPEN_ACCESS:y` drops subscription flagships** | `isOpenAccess = "Y"` (`CC-BY 4.0`). | **`isOpenAccess = "N"`**, `inPMC = "N"`, `pmcid = ""`. | Discovery queries **must never** include `OPEN_ACCESS:y` or `HAS_FT:y`. |
| 3 | **EBI GWAS Catalog retired v1 API (`410`) & lacks both flagships (`404`)** | `/gwas/rest/api/...` -> `HTTP 410 Gone`; active `/v2/publications/37794186` -> **`HTTP 404`** (`0` studies). | `/v2/publications/34857953` -> **`HTTP 404`** (`0` studies). | Portal/Synapse-hosted pGWAS are absent from EBI GWAS Catalog; multi-source redundancy is essential. |
| 4 | **Key metadata missing from Open-Access JATS main text** | `PMC10567551` JATS XML (`205 KB`) has **`0` `<table-wrap>` nodes** and **`0` mentions of `syn51364943`, `synapse`, or `EDTA`**. | Closed license (`springer.com/tdm`). | Tier 3 must parse `MOESM1` BioC JSON (`EDTA plasma`), BioStudies `MOESM3_ESM.xlsx` (`2,941` assays), and 1-hop `ukb-ppp.gwas.eu` -> Synapse `syn51364943` (discovering all **7 ancestry strata** incl. `MID` and `AMR`). |
| 5 | **Genome build omitted from entire paper & supplements** | Explicitly states `hg38` in `Methods` & `ST9`. | **`GRCh38`/`hg38`/`GRCh37`/`hg19` occur `0` times across full text and `MOESM1–4`!** | Created Tier 3 **Sentinel `rsID` Coordinate Build Resolver** (`ST03` `rs1891906=chr1:1014863` -> `2/2` exact match on `GRCh38`). |
| 6 | **Interactive JS portal gates masquerade as `HTTP 200`** | Synapse `/file?redirect=false` -> `HTTP 403` (`"Anonymous users have only READ access"`). | `download.decode.is/form/...` -> `HTTP 200 text/html` (`605 B` React SPA license gate). | Tier 3 inspects `Content-Type` and response bodies without bypassing gates, classifying both as `"Portal registration"`. |
| 7 | **NCBI `pmcoa.cgi` returns `HTTP 200 text/html` on errors** | Valid `141 KB` BioC JSON. | `HTTP 200` (`6 KB` HTML: `[Error] : No result can be found.`). | Parser validates JSON MIME/structure before parsing. |

---

## 4. Adversarial Subagent Peer Review & Version 2 (`v2.0`) Implementation

### 4.1 Independent Peer-Review Audit (`peer_review_round1.md`)
Before writing production code, a third subagent (`peer-reviewer`) performed a cross-audit of `design_search_v1.md`, `design_extraction_v1.md`, and existing repository validation rules, identifying **`6 CRITICAL`**, **`8 MAJOR`**, and **`4 MINOR`** issues:

1. **`[C-1]` Search-to-Extraction Queue Handoff:** Expanded `catalogue/candidate_queue.tsv` from 9 to **17 columns** (appending `pmid`, `doi`, `pmcid`, `epmc_id`, `matched_sources`, `p_score`, `s_score`, `preprint_lineage_note` while keeping columns 1–9 unchanged for UI compatibility).
2. **`[C-2]` Vocabulary Unification & Null-Safety:** Unified `sample_matrix` to 12 controlled values (`EDTA plasma`, `Citrate plasma`, `Heparin plasma`, `Plasma (unspecified)`, `Serum`, `Whole blood / DBS`, `CSF`, `Brain tissue`, `Urine`, `Saliva`, `Cell line / LCL`, `Other tissue`), enforced non-`NA` booleans on `has_eaf`/`has_beta_se`, and fixed a simulated typo (`zhang_2022` `Olink` -> `SomaScan`).
3. **`[C-3]` Early Panel/Peptide Lexicon & Snowball Guard:** Added pre-2017 multiplex cytokine/immunoassay (`Bio-Plex`, `Milliplex`, `serum/plasma cytokines`) and targeted MS (`LC-MS/MS`, `MRM/PRM/SRM`, `plasma peptides`) terms, lowered `P4` protein count floor from `10` to `5`, and enforced `BROAD_PROTEIN AND BROAD_GENETICS` on citation snowballing + Stage 1 `G0` to block pure metabolomics/eQTL citations.
4. **`[C-4]` Low-RSS Zero-DOM `.xlsx` Supplement Streaming (`1.75 GB -> 273 MB`):**
   - Profiling `readxl::read_excel()` on Ferkingstad's `44.4 MB` (`302.9 MB` uncompressed XML; `sheet12.xml = 118.9 MB`) `MOESM4_ESM.xlsx` revealed a **`1.75 GB` C++ `rapidxml` heap spike** and a **`26.7 s`** runtime across Sun's 31 sheets.
   - Re-engineered Tier 3 (`src/extract/tier3_supplements.R`) to inspect ZIP headers (`utils::unzip(..., list = TRUE)` + `xl/workbook.xml` in `< 2 ms`), parse the `Contents` TOC sheet first, and stream only the first `64 KB` via base `unz()` + `1 MB`-chunked `sharedStrings.xml` scanning whenever `sheet*.xml > 15 MB` — cutting peak OS RAM **6.1-fold (`1.75 GB` -> `272.9 MB`)** in `9.4 s`.
5. **`[C-5]` Sample-Size Disambiguation:** Prevented parent-biobank recruitment counts (`500,000 UK Biobank`) from overriding GWAS `n_samples`, reconciling proteomic-profiled `N` (`54,219`) vs post-QC genetic GWAS `N` (`52,363`) automatically.
6. **`[C-6]` Live Preprint API Resilience (`PPR508031`):** Verified live that Europe PMC `/rest/PPR508031/fullTextXML` returns `HTTP 500` and `biorxiv.org/...source.xml` returns `HTTP 429` Cloudflare blocks; updated Tier 2c to check NIH Preprint Pilot `PMC*` IDs first, harvest `/rest/PPR/{id}/datalinks`, never bypass `403/429`, and escalate cleanly to `NEEDS_MANUAL_REVIEW`.
7. **`[M-5]` Strictly Blinded Positive-Control Summary Mode:** Defaulted `src/evaluate_positive_controls.R` to `--mode blinded-summary`, reporting aggregate recall, triage survival, per-source unique lift, and 8-class failure counts **without printing missed study titles/IDs** so search tuning remains strictly blinded.

Both specialist subagents implemented these fixes into modular R production scripts (`src/search/lexicons.R`, `src/run_literature_search.R`, `src/evaluate_positive_controls.R`, `src/extract/*.R`, `src/extract_candidate_study.R`), verified green on both macOS and Ubuntu GitHub Actions (`v2.0`).

---

## 5. Round-1 Blinded Dual-Agent 2025 Replication & Referee Adjudication (`v2.0 -> v3.0`)

### 5.1 Experimental Design (Option A — Pure Independent Replication)
To stress-test `v2.0` on real, noisy literature without touching pre-2025 studies (preserving historical positive-control blinding), two fresh, independent subagents (`Independent Search Scientist 1` and `Independent Search Scientist 2`) were given **word-for-word identical instructions** in isolated sandboxes (`agent_1/` and `agent_2/`) initialized with **empty** `catalogue/*.tsv` files over the fixed calendar window **`--since 2025-01-01 --until 2025-12-31`**:
1. Run unmodified `v2.0` baseline on `2025-01-01..2025-12-31`.
2. Audit false positives, date leaks, and false negatives across all triage buckets.
3. Independently re-engineer `lexicons.R` and `run_literature_search.R`, re-run `2025-01-01..2025-12-31`, and document paper-by-paper diffs.
4. A third independent subagent (`Independent Search Referee & Adjudicator`) compared code diffs, linked entities across dumps, computed inter-rater reliability and capture-recapture sensitivity, and adjudicated every discordant study against full abstracts and PMC BioC `Methods` (`analysis/search-benchmark-2025/adjudication_report_2025.md`).

### 5.2 Quantitative Convergence & Concordance Results

| Metric (`2025-01-01` to `2025-12-31`) | `v2.0` Baseline | Agent 1 (`v3_A1`) | Agent 2 (`v3_A2`) | Inter-Agent Concordance |
| :--- | ---: | ---: | ---: | :--- |
| **Total Raw Multi-Source Hits** | 2,451 | 3,993 | 2,435 | — |
| **Deduplicated Works in Dump** | 1,815 *(226 date leaks)* | 1,747 *(0 actionable leaks)* | 1,698 *(0 leaks)* | `1,309` shared (`Jaccard = 0.613`) |
| **`Likely Primary pGWAS`** | 92 *(~22.8% precision)* | **42** *(95.2% precision)* | **43** *(97.7% precision)* | **`37` shared (`99.39%` agreement, $\kappa = 0.899$)** |
| **`Needs Human Check`** | 308 | **7** | **7** | `2` shared + `1` cross-bucket (`LP`/`NHC`) |
| **Total Actionable (`LP + NHC`)** | **400** | **49** *(91.8% true)* | **50** *(92.0% true)* | **`40` shared (`98.78%` agreement, $\kappa = 0.827$)** |
| **Concordant Intersection Precision (`A_1 ∩ A_2`)** | — | — | — | **`40 / 40 = 100.0%` (`0` false positives)** |

### 5.3 Convergent vs. Complementary Discoveries Across Agent 1 and Agent 2
- **15 Convergent Fixes (found independently by both agents):**
  1. HTTP exponential backoff + `Retry-After` handling on transient EBI `HTTP 500/503` responses.
  2. Date-bounding `Snowball_Backward` (`/references`), eliminating all `86` `2001–2023` citation leaks.
  3. Switching EBI GWAS Catalog v2 from dateless `/v2/studies` (8,000-study cap) to dated `/v2/publications`.
  4. Normalizing late-2025 online publications (`epub`) assigned to 2026 print volumes back to `2025`.
  5. Preserving `SRC:PMC` IDs (`11 -> 0` silently dropped PMC rows).
  6. Supporting bioRxiv/medRxiv's Dec-2025 **`10.64898/` (`openRxiv`)** DOI prefix + `Preprints.org`/`Authorea` prefixes.
  7. Intra-batch preprint-to-published graph consolidation (`Hartley`, `Hu`).
  8. Repairing broken trailing `\b` regex stems (`whole-genome sequenc\b` -> `sequenc[a-z]*`, rescuing `Hawkes 2025` *Nat Genet* WGS pGWAS).
  9. Normalizing `<i>cis</i>-pQTL` tag-stripped hyphen spacing (`cis - pQTL` -> `cis-pQTL`).
  10. Adding 2025 platform/method lexicons (`Orbitrap Astral`, `timsTOF`, `Seer Proteograph`, `human knockouts`, `admixture mapping`, `SAIGE`, `QTLtools`).
  11. Disambiguating bare `"TMT"` away from neuropsychological *Trail Making Test (TMT-A/B)* dementia GWASs.
  12. Expanding `G2` non-human organism filtering to 80+ livestock, crop, aquaculture, and microbial taxa with title-first priority (`36 -> 148–176` gated).
  13. Scoping `P2` so non-proteomic disease/N-glycome GWASs cannot score `+6`.
  14. Blocking single-analyte (`< 5` proteins) biomarker GWASs (`procalcitonin`, `AFP`, `HRG`, `sVEGFR1`).
  15. Requiring an explicit germline proteogenomic anchor in Stage 2 (eliminating `~140` pure observational Olink/SomaScan cohort papers).

- **Critical Orthogonal Fixes Unique to `Agent 1` (`8` discoveries):**
  - **Invalid `MESH:` field in Europe PMC & `<3,100`-char `TITLE_ABS:` query sharding:** `MESH:` is not a valid Europe PMC search field (`KW:` indexes MeSH); `OR MESH:"..."` corrupted server-side Lucene Boolean precedence on `SRC:PPR` preprints while leaking `347` unrelated `EuropePMC_MED` rows into `D_2`. Sharding Prong 2A into two compact `TITLE_ABS:` queries boosted unique 2025 preprints **2.4-fold (`240` vs `98`)**.
  - **Variance-QTL (`vQTL` / `cis-vQTL`) lexicon & triage:** Uniquely discovered **`Ben-Eghan et al. 2025`** (`10.1101/2025.11.20.25340688`, UKB multi-ancestry `cis-vQTL` pGWAS, `N = 52,706`).
  - **10-seed diverse citation snowballing** (`715` vs `119` forward hits), **biological-sample reuse exemption** (`"samples/brains obtained from N patients"`), **`<PubmedBookArticle>` hydration**, and `"over/more than [N] proteins"` quantifier support (rescuing **`Wojtas 2025`** ADNI CSF + brain vascular TMT-MS pGWAS).

- **Critical Orthogonal Fixes Unique to `Agent 2` (`7` discoveries):**
  - **Inequality `<[^>]+>` HTML-truncation bug (`strip_html_tags`):** Decoding `&lt;`/`&gt;` before running `gsub("<[^>]+>", " ", ...)` deleted `>= 25` characters across **`92` abstracts (`5.4%`)** in `Agent 1` whenever `P < 5e-8` preceded `>5,000` (deleting `544` chars of `Haug 2025`!). Tag-specific `</?[a-zA-Z]...>` stripping prior to entity decoding eliminated all truncations.
  - **Date-sorted GWAS Catalog `/v2/publications` + EFO-independent `disease_trait` verification:** Uncurated summary-statistics submissions lack `EFO_0004747` mappings; inspecting `disease_trait` strings uniquely rescued **`Loya et al. 2025`** (*Nat Genet*, `PMID: 39789286`, `2,929` UKB Olink pGWAS with **zero** `"protein"` mention in title/abstract!) and **`Yan et al. 2025`** (*AJOG*, `9` maternal serum proteins, `27` GCSTs).
  - **Canonical `first_pub_date` precedence** (blocking `~85` Epub-2024/print-2025 and revised 2023/2024 preprint leaks), `.v1`-preserving `10.20944/` DOIs, `MR_BOILERPLATE_STOPWORDS` fuzzy-merge protection, and expanded phrasing rescuing **`Farinas 2025`** (*Nat Med* CSF-to-plasma protein ratio GWAS), **`Dimas 2025`** (diet-responsive `cis-pQTLs`), and **`Hall 2025`** (`N = 69` donor iPSC-RPE MS `cis-pQTLs`).

### 5.4 Ground-Truth Adjudication & Capture-Recapture Sensitivity (`2025`)
Across all `59` union actionable studies (`40` concordant + `19` discordant):
- **`45` `TRUE_PRIMARY_2025`** (`38` in `A_1 ∩ A_2`, `2` in `A_1 \ A_2`, `5` in `A_2 \ A_1`).
- **`6` `TRUE_BORDERLINE_HUMAN_CHECK`** (`2` in `A_1 ∩ A_2`, `3` in `A_1 \ A_2`, `1` in `A_2 \ A_1`).
- **`8` `FALSE_POSITIVE_SECONDARY_OR_NON_PGWAS`** (`4` in `A_1 \ A_2`, `4` in `A_2 \ A_1`; `0` in `A_1 ∩ A_2`).

Using Chapman's bias-corrected capture-recapture estimator on verified `TRUE_PRIMARY_2025` (`n_1 = 40`, `n_2 = 43`, `m_{12} = 38`):
$$\hat{N}_{\text{Chapman}} = \frac{(40 + 1)(43 + 1)}{38 + 1} - 1 = 45.26 \quad (\text{SE} = 0.54,\ \text{Burnham-Chao 95\% CI: } 45.02\text{--}48.32)$$
- **Agent 1 Sensitivity:** `88.4%` (`40 / 45.26`)
- **Agent 2 Sensitivity:** `95.0%` (`43 / 45.26`)
- **Combined `A_1 ∪ A_2` Sensitivity:** **`99.4%`** (`45 / 45.26`; `91.8%–97.8%` under conservative positive source-dependence bounds $\hat{N} \in [46, 49]$).

---

## 6. Unified Version 3 (`v3.0`) Consensus Merge & Round-2 Blinded Replication (`v3.0 -> v4.0`)

### 6.1 Unified `v3.0` Consensus Architecture & Regression Verification
All 15 convergent fixes, Agent 1's 8 unique improvements, Agent 2's 7 unique improvements, and **8 Referee conflict-resolution safeguards** (tag-first HTML stripping before entity decoding; post-merge journal `epub > ppub` precedence over superseded preprints; internal-only `.v1` stripping on `10.20944/`; boilerplate-stripped fuzzy merger; biological-sample reuse exemption paired with `S1` external-source requirement without `!p2_hit`; preserving MEDLINE `Letter` articles while gating disease-only no-abstract letters; spelled-out `protein quantitative trait loci` in `P5` + `iPSCs` in `P6`; and candidate-SNP `< 100` / twin-ACE / composite aging-clock guards) were merged into `src/search/lexicons.R` and `src/run_literature_search.R` (`v3.0`, commit `a68d4d2`).

Prior to launching Round 2, `v3.0` was verified against a 22-assertion unit suite and the full Round-1 59-study union + 1,698-work cached 2025 corpus:
- **Europe PMC Query Shards (`0` invalid `MESH:` fields):** `2A-1 (Platforms + Anchors)` = `3,377` chars; `2A-2 (Proteome + Analytes)` = `3,410` chars; `2B (METHODS Rescue)` = `873` chars; `Snowball (10 Landmark Seeds)` = `2,106` chars.
- **59-Study Adjudicated Union Regression:** `45 / 45` (`100%`) `TRUE_PRIMARY_2025` classified as `Likely Primary pGWAS`; `6 / 6` (`100%`) `TRUE_BORDERLINE_HUMAN_CHECK` classified as `Needs Human Check`; `8 / 8` (`100%`) `FALSE_POSITIVE` classified as `Auto-Deprioritised`.
- **Full-Corpus Specificity Check (`N = 1,698` cached 2025 works):** `0` newly introduced false-positive actionable studies.
- **Live Positive Benchmark Wiring (`--verify-pmids 37794186,34857953`):** `Sun 2023` (`P = 28, S = 0`, `Olink`, `N = 54,219`, preprint lineage `PPR508031|10.1101/2022.06.17.496443|sun_2022_preprint`) and `Ferkingstad 2021` (`P = 20, S = 0`, `SomaScan`, `N = 35,559`) both passed end-to-end across all four discovery sources.

### 6.2 Round-2 Blinded Dual-Agent Replication (`Agent 3`, `Agent 4`, and `Reviewer 2`) — *In Progress*
To verify whether `v3.0` has reached asymptotic convergence or whether independent eyes uncover subtle second-order recall/precision improvements, a second blinded dual-agent replication (`Agent 3` and `Agent 4` starting from `v3.0` in clean empty-catalogue sandboxes on `--since 2025-01-01 --until 2025-12-31`, adjudicated by a fresh `Reviewer 2`) is currently executing. Results, concordance metrics, and any `v3.0 -> v4.0` refinements will be appended here upon completion.

---

## 7. Planned Out-of-Sample Blinded Positive-Control Validation

Once the search and triage strategy reaches multi-agent convergence on `2025`, out-of-sample historical recall and triage survival will be evaluated against the user-curated, air-gapped positive-control benchmark (`data/positive_controls.tsv`, derived independently from Karsten Suhre's proteomics GWAS compilation) via `src/evaluate_positive_controls.R --mode blinded-summary`, recording final sensitivity prior to executing the full historical baseline search (`2005–present`).
