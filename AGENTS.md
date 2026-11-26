# Workspace Rules: Proteomic GWAS Catalogue

## bioRxiv & medRxiv Preprint `.full.pdf` Retrieval Rule

When retrieving full-text PDFs (or JATS XML / supplementary files) for `bioRxiv` and `medRxiv` preprints (`10.1101/*` and `10.64898/*`):

1. **Direct `.full.pdf` & Official Early-Path URL Resolution:**
   - Every `bioRxiv` / `medRxiv` preprint exposes its PDF by appending `.full.pdf` to the versioned content URL:
     - `https://www.{biorxiv|medrxiv}.org/content/{doi}v{version}.full.pdf` (e.g. `https://www.medrxiv.org/content/10.1101/2023.12.20.23300308v1.full.pdf`).
   - In automated scripts (`src/extract/tier1_apis.R`, `src/extract/tier2_fulltext.R`, and period open-access harvesters), resolve `https://api.biorxiv.org/details/{biorxiv|medrxiv}/{doi}/na/json` first to obtain `version` and `jatsxml` (`https://www.{server}.org/content/early/YYYY/MM/DD/{id}.source.xml`), and construct both:
     - **Early-path PDF (`medRxiv` & `bioRxiv` HTTP 200):** `https://www.{server}.org/content/early/YYYY/MM/DD/{id}.full.pdf` (substituting `.source.xml` -> `.full.pdf`), alongside `.source.xml` and supplementary media at `https://www.{server}.org/content/early/YYYY/MM/DD/{id}/DC1/embed/{filename}`.
     - **Direct versioned PDF (`bioRxiv` HTTP 200 & 1-click manual link):** `https://www.{server}.org/content/{doi}v{version}.full.pdf`.

2. **Polite Rate-Limiting, Payload Verification & Zero Anti-Bot Bypass:**
   - Space requests to `biorxiv.org` and `medrxiv.org` by **`>= 3.0` seconds** to respect Cold Spring Harbor Laboratory rate limits (avoiding `HTTP 429` error `1015`).
   - Verify binary/XML magic headers (`%PDF` for `.full.pdf`, `<article` for `.source.xml`, `PK\x03\x04` for `.xlsx`/`.docx`) before accepting downloaded files, and also verify content authenticity on generic EBI BioStudies `S-EPMC*` `media-*`/`mmc*` files.
   - **Never** spoof browser headers or attempt to bypass any `HTTP 403` / `HTTP 429` Cloudflare Turnstile or anti-bot challenge (`allow_transient_retry = FALSE`). Whenever a preprint PDF or supplement remains gated, escalate cleanly in `search_report.md` (Sections 2 & 3) using the direct `.full.pdf` / `.supplementary-material` URL for manual download.

## Search Window End-Date Rule (`current_search_date - 1 day`)

- All literature search periods **must end on the calendar day immediately before the current search execution date** (`until_date = Sys.Date() - 1L`, e.g. `2026-11-25` when searching on `2026-11-26`).
- Ending on `current_search_date - 1 day` guarantees every queried date is a complete 24-hour indexing day across PubMed, Europe PMC, bioRxiv/medRxiv, Crossref, and GWAS Catalog without partial-day boundary gaps, and subsequent incremental periods begin cleanly on `previous_date_to + 1 day`.
