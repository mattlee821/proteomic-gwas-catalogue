# proteomic-gwas-catalogue

An all-encompassing, continuously updated catalogue and interactive web application of published and preprint genome-wide association studies (GWAS) of proteins — built so researchers can quickly locate, evaluate, and harmonise summary statistics for meta-analysis, Mendelian randomisation (MR), and colocalisation.

Status: prototype web application, schema, and GitHub Actions automation built on a simulated dataset (`v0.2.0-sim`); full systematic curation next.

## Purpose & Scope

Protein GWAS are spread across dozens of papers, preprints, cohorts, proteomic platforms, biological fluids, and data repositories. Unlike static blog tables that edit entries in place or omit operational harmonisation metadata, this catalogue tracks every primary protein GWAS in a version-controlled, machine-readable repository deployed to **GitHub Pages** with **tagged GitHub Releases**:

- **All proteomic technologies:** SomaLogic (SomaScan 1.1k to 11k), Olink (Target 96 / SCALLOP to Explore HT 5k+), mass spectrometry (DIA/SWATH, TMT, nanoparticle-enriched MS), Alamar NULISAseq, and multiplex immunoassays (Luminex).
- **Preprints included:** bioRxiv, medRxiv, and Research Square preprints are tracked alongside peer-reviewed papers and linked cleanly via `superseded_by` once published.
- **No sample-size floor:** Captures everything from small-N tissue, pediatric, urine, and CSF pQTL studies (`N < 100`) up to biobank-scale plasma studies (`N > 50,000`).
- **Operational readiness for downstream genetics:** Records multi-ancestry and sex-specific strata, genome build, *cis*-only vs genome-wide coverage, phenotype transformation and effect scale, signed EAF and indel encoding, cohort membership for overlap checking, and links to `GWASprocessing` YAML manifests.

## Repository Structure

- `catalogue/`: Canonical Git-tracked TSV tables (`studies.tsv`, `strata.tsv`, `cohort_membership.tsv`, `proteins_index.tsv`, `search_runs.tsv`, `candidate_queue.tsv`, `releases.tsv`). Currently populated with simulated benchmark + search data (`v0.2.0-sim`).
- `site/`: Quarto + Observable JS (`ojs`) static web application with four interactive views:
  1. **Catalogue Explorer (`site/index.qmd`):** Faceted filtering across strata/studies, timeline bubble plot, row-level harmonisation inspector, and filtered TSV export.
  2. **Protein Lookup (`site/proteins.qmd`):** Instant cross-study search by gene symbol, UniProt accession, SomaScan `SeqId`, or `OlinkID`, with epitope/PAV sensitivity warnings.
  3. **Cohort Overlap & Lineage (`site/overlap.qmd`):** Interactive biobank/cohort overlap matrix for meta-analysis and two-sample MR checks, plus preprint-to-published lineage.
  4. **Search Log & Releases (`site/updates.qmd`):** PRISMA screening metrics, automated search run logs, candidate PR triage preview, and versioned TSV/JSON downloads.
- `src/`: R and shell scripts:
  - `src/setup_quarto.sh`: Installs Quarto CLI into `~/.local/share/quarto` and `~/.local/bin/quarto` (outside OneDrive).
  - `src/simulate_catalogue.R`: Generates the simulated `v0.2.0-sim` dataset.
  - `src/validate_catalogue.R`: Validates primary/foreign keys and controlled vocabularies (enforced in CI).
  - `src/build_site_data.R`: Validates `catalogue/*.tsv` and compiles `site/data/catalogue.json`.
  - `src/run_literature_search.R`: Queries Europe PMC (PubMed + preprints), deduplicates against screened studies, deprioritises secondary 2SMR/EWAS papers, and writes `catalogue/candidate_queue.tsv`.
- `.github/`:
  - `workflows/scheduled-search.yml`: Runs on the 1st and 15th of each month (or manually), opening a Pull Request with new candidate studies and preprint upgrade checks.
  - `workflows/validate-and-deploy.yml`: Validates `catalogue/*.tsv`, renders `site/`, and deploys to GitHub Pages on merge to `main`.
  - `workflows/release-catalogue.yml`: Publishes versioned TSV/JSON snapshots and `SHA256SUMS` as a GitHub Release whenever a `v*` tag is pushed.
  - `ISSUE_TEMPLATE/`: Structured forms for community submissions and study updates.

## Quick Start

```bash
# 1. Ensure Quarto CLI is available in ~/.local/bin
bash src/setup_quarto.sh

# 2. (Prototype only) Regenerate simulated catalogue tables
Rscript src/simulate_catalogue.R

# 3. Validate tables and build site/data/catalogue.json
Rscript src/build_site_data.R

# 4. Render or preview the website locally
~/.local/bin/quarto render site
~/.local/bin/quarto preview site
```
