# proteomic-gwas-catalogue

An all-encompassing, continuously updated catalogue and interactive web application of published and preprint genome-wide association studies (GWAS) of proteins — built so researchers can quickly locate, evaluate, and harmonise summary statistics for meta-analysis, Mendelian randomisation (MR), and colocalisation.

**Current release:** `v0.1.20` (`2000-01-01` to `2026-11-25` curated: `305` studies across `456` strata) · **Website:** https://mattlee821.github.io/proteomic-gwas-catalogue/

## Purpose & Scope

Protein GWAS are spread across dozens of papers, preprints, cohorts, proteomic platforms, biological fluids, and data repositories. Unlike static blog tables that edit entries in place or omit operational harmonisation metadata, this catalogue tracks every primary protein GWAS in a version-controlled, machine-readable repository deployed to **GitHub Pages** with **tagged GitHub Releases**:

- **All proteomic technologies:** SomaLogic (SomaScan 1.1k to 11k), Olink (Target 96 / SCALLOP to Explore HT 5k+), mass spectrometry (DIA/SWATH, TMT, nanoparticle-enriched MS), Alamar NULISAseq, and multiplex immunoassays (Luminex).
- **Preprints included:** bioRxiv, medRxiv, and Research Square preprints are tracked alongside peer-reviewed papers and linked cleanly via `superseded_by` once published.
- **No sample-size floor:** Captures everything from small-N tissue, pediatric, urine, and CSF pQTL studies (`N < 100`) up to biobank-scale plasma studies (`N > 50,000`).
- **Operational readiness for downstream genetics:** Records multi-ancestry and sex-specific strata, genome build, *cis*-only vs genome-wide coverage, phenotype transformation and effect scale, signed EAF and indel encoding, cohort membership for overlap checking, and links to `GWASprocessing` YAML manifests.

## Pipeline: Search to Version Update

Each catalogue update combines automated multi-database search and open-access text/supplement extraction with **human review checkpoints** for paywalled papers, missing supplementary files, study metadata verification, and versioned releases (coloured to match the `AI-extracted` and `Human-reviewed` `CURATION` badges in the Catalogue):

<div align="center">

```mermaid
%%{init: {'themeVariables': {'fontSize': '15px'}, 'flowchart': {'wrappingWidth': 390, 'nodeSpacing': 28, 'rankSpacing': 44, 'subGraphTitleMargin': {'top': 12, 'bottom': 30}, 'padding': 14, 'curve': 'basis'}}}%%
flowchart TD
    S1["`**1. Literature & Preprint Search**
Query PubMed, Europe PMC, GWAS Catalog & citations`"]:::auto
    S2["`**2. Deduplicate & Screen**
Match preprints to papers; filter out secondary MR studies`"]:::auto
    S3["`**3. Automated Extraction**
Parse open-access text, supplements & summary-stat links`"]:::auto

    subgraph MANUAL["Human Review & Verification"]
        M1["`**4.1 Eligibility**
Confirm borderline study inclusions & exclusions`"]:::manual
        M2["`**4.2 Missing Files**
Download paywalled PDFs & supplementary tables`"]:::manual
        M3["`**4.3 Metadata Review**
Verify extracted study, cohort & assay details`"]:::manual
        M4["`**4.4 Recall Check**
Test completeness against positive-control studies`"]:::manual
    end

    S4["`**5. Validate & Update Website**
Check data schema, merge update & rebuild web tables`"]:::auto
    M5["`**6. Version Release**
Publish numbered catalogue release & data archive`"]:::manual

    S1 --> S2 --> S3
    S3 --> M1 & M2 & M3 & M4
    M1 & M2 & M3 & M4 --> S4
    S4 --> M5

    classDef auto fill:#f3e8ff,stroke:#c084fc,stroke-width:1.6px,color:#6b21a8
    classDef manual fill:#e0f2fe,stroke:#38bdf8,stroke-width:2px,color:#075985
    style MANUAL fill:#f0f9ff,stroke:#38bdf8,stroke-width:2px,stroke-dasharray: 6 4,color:#075985
```

</div>

- **Purple (`Automated / AI-extracted` — Steps 1–3 & 5):** Queries literature and preprint databases, links preprints to published journal articles, filters out secondary Mendelian randomisation studies, extracts metadata from open-access and downloaded files, validates `catalogue/*.tsv` schemas, and rebuilds the website.
- **Blue (`Human-reviewed` — Steps 4.1–4.4 & 6):** Human reviewer downloads gated PDFs/supplements for automated parsing (which remain `AI-extracted` unless a curatorial choice is required), resolves genuine Step 4 curatorial choices (marking only those studies `Human-reviewed` with a short inline summary such as `Strata check`), tests recall against positive-control studies, and publishes numbered releases.

## Repository Structure

- `catalogue/`: Canonical Git-tracked TSV tables (`studies.tsv`, `strata.tsv`, `cohort_membership.tsv`, `proteins_index.tsv`, `search_runs.tsv`, `candidate_queue.tsv`, `excluded_records.tsv`, `releases.tsv`).
- `site/`: Quarto + Observable JS (`ojs`) static web application with five pages:
  1. **Catalogue (`site/index.qmd`):** Full-width 11-column reference table (`AUTHOR`, `YEAR`, `COHORT`, `OVERLAP`, `POPULATION`, `SEX`, `N`, `SAMPLE`, `PLATFORM`, `ASSAY`, `DATA`) with inline per-column filter/search inputs, pagination, click-to-expand inline 3-column detail rows (`PHENOTYPE`, `GENETIC`, `OTHER`), and full-table `.tsv` download.
  2. **Lookup (`site/lookup.qmd`):** Protein-centric typeahead lookup by UniProt accession, gene symbol, or protein name, populating summary KPI boxes, a paginated study table, and an interactive `Year × Sample Size` scatter plot.
  3. **Log (`site/log.qmd`):** Date-stamped search update log (newest at top) with inline-expandable run details and `preprint -> published` replacement tracking.
  4. **About (`site/about.qmd`) & Search Strategy (`site/search-strategy.qmd`):** Pipeline overview, flowchart, eligibility criteria, curation tiers, and full methodological history synced from `docs/methods.md`.
- `src/`: R and shell scripts:
  - `src/setup_quarto.sh`: Installs Quarto CLI into `~/.local/share/quarto` and `~/.local/bin/quarto` (outside OneDrive).
  - `src/validate_catalogue.R`: Validates primary/foreign keys, controlled vocabularies, and positive-control blinding guardrails (enforced in CI).
  - `src/build_site_data.R`: Validates `catalogue/*.tsv` and compiles `site/data/catalogue.json`.
  - `src/run_literature_search.R` & `src/search/lexicons.R`: Multi-resource discovery engine (PubMed POST, Europe PMC journal + preprint `searchPOST` + `METHODS:` rescue, EBI GWAS Catalog REST v2, and blinded 6-seed citation snowballing), graph deduplication, 4-tier preprint-to-published upgrade detector, and Two-Stage Dual-Axis (`P_score` vs `S_score`) triage scorer.
  - `src/extract_candidate_study.R` & `src/extract/*.R`: 4-tier cascading paper, low-RSS supplementary `.xlsx`/`.docx`/BioC, repository access-gate preflight, sentinel `rsID` coordinate genome-build resolver (`catalogue/schema/sentinel_build_anchors.tsv`), atomic field provenance, and `NEEDS_MANUAL_REVIEW` escalation checklist generator.
  - `src/evaluate_positive_controls.R`: Strictly air-gapped post-hoc evaluator (`--mode blinded-summary`) for testing recall and triage survival against `data/positive_controls.tsv` (git-ignored).
- `.github/`:
  - `workflows/scheduled-search.yml`: PR-gated search + 4-tier staging pipeline (cron paused until final search strategy sign-off; manual `workflow_dispatch` supports targeted `--verify-pmids` tests).
  - `workflows/validate-and-deploy.yml`: Validates `catalogue/*.tsv`, renders `site/`, and deploys to GitHub Pages when enabled.
  - `workflows/release-catalogue.yml`: Manual/tag-triggered workflow publishing versioned TSV/JSON snapshots and `SHA256SUMS` as a GitHub Release only when explicitly requested.
  - `ISSUE_TEMPLATE/`: Structured forms for community submissions and study updates.

## Quick Start

```bash
# 1. Launch local interactive website preview anytime (Ctrl+C to stop)
cd ~/research/proteomic-gwas-catalogue && ~/.local/bin/quarto preview site

# 2. Targeted 2-paper verification of search + 4-tier extraction (no broad search)
Rscript src/run_literature_search.R --verify-pmids 37794186,34857953 --out-dir output/verify
Rscript src/extract_candidate_study.R --queue-file output/verify/candidate_queue_preview.tsv --batch-pending --out-dir output/verify/staging

# 3. Evaluate search recall & triage survival against a blinded positive-control table
Rscript src/evaluate_positive_controls.R --controls data/positive_controls.tsv --dump data/search_dump_<run_id>.tsv --mode blinded-summary
```

## License & Citation

- **Curated Data (`catalogue/`):** [Creative Commons Attribution 4.0 International (CC BY 4.0)](LICENSE)
- **Source Code & Website (`src/`, `site/`, `.github/`):** [MIT License](LICENSE)
- **Citation:** Metadata is provided in [`CITATION.cff`](CITATION.cff); please also cite the underlying primary protein GWAS studies when using their summary statistics.
