# ATTEMPT Ancillary Analysis

Code repository for the ancillary analyses in the **ATTEMPT trial** manuscript, submitted to *Science Translational Medicine*.

The ATTEMPT trial (A Trial To Evaluate Mechanisms Preserving Tubuloglomerular feedback in Type 1 diabetes) is a randomized, double-blind, placebo-controlled clinical trial studying the effects of Dapagliflozin (SGLT2 inhibitor) vs. Placebo in individuals with Type 1 Diabetes, incorporating kidney biopsies for single-cell and proteomic molecular profiling.

## Repository Structure

| File | Description |
|------|-------------|
| `attempt_functions.R` | Shared utility functions used across all analysis scripts (differential expression processing, volcano plots, GSEA, pseudotime, proteomics correlation, etc.) |
| `attempt_scrna_hyak.qmd` | Single-cell RNA-seq analysis pipeline: QC, doublet removal, normalization, Harmony integration, KPMP cell type annotation, and NEBULA NBLMM differential expression modeling |
| `attempt_urine_somascan_analysis.qmd` | Urine SomaScan proteomics: difference-in-differences (DiD) limma models, volcano plots, GSEA pathway enrichment, clinical correlations, and scRNA integration |
| `attempt_plasma_somascan_analysis.qmd` | Plasma SomaScan proteomics: parallel analysis to urine with plasma-specific DiD models, pathway analysis, and OrganAge biological age estimation |
| `ATTEMPT_zi_crocodile_vis_kpmp.qmd` | CROCODILE cross-study comparison: concordance analysis between ATTEMPT treatment effects and T1D vs. healthy control differences from KPMP |
| `ATTEMPT_zi_emmeans_vis_kpmp.qmd` | Estimated marginal means visualization, DEG analysis across cell types, IPA pathway analysis, GSEA, and chord diagram summaries |
| `attempt_pseudotime.qmd` | Slingshot-based pseudotime trajectory inference on scRNA-seq data with clinical covariate associations |

## Setup and Configuration

### Environment Variables

All scripts use environment variables for configurable paths, with sensible defaults for local use. Set these before running if your data is not in the default locations:

```bash
# Path to the root data directory containing raw data, Seurat objects, etc.
export ATTEMPT_DATA_PATH="/path/to/data"

# Path to the results/output directory
export ATTEMPT_RESULTS_PATH="/path/to/results"

# Path to the credentials JSON file (for S3/cloud data access)
export ATTEMPT_KEYS_PATH="/path/to/keys.json"

# Path to the Python executable (needed for scCODA and OrganAge)
export ATTEMPT_PYTHON_PATH="/path/to/python3"
```

If environment variables are not set, scripts default to:
- `ATTEMPT_DATA_PATH` → `here::here("data")`
- `ATTEMPT_RESULTS_PATH` → `here::here("results")`
- `ATTEMPT_KEYS_PATH` → `"keys.json"` (in working directory)
- `ATTEMPT_PYTHON_PATH` → `"python3"`

### R Dependencies

Key R packages required (non-exhaustive):

- **Single-cell**: Seurat, nebula, harmony, slingshot, SingleCellExperiment, scater, scran, DoubletFinder
- **Proteomics**: limma, SomaDataIO
- **Pathway analysis**: fgsea, enrichR, msigdbr
- **Visualization**: ggplot2, cowplot, patchwork, ComplexHeatmap, circlize
- **Data wrangling**: dplyr, tidyr, purrr, data.table, openxlsx
- **Modeling**: lmerTest, emmeans, glmmTMB
- **Utilities**: here, jsonlite, aws.s3

### Python Dependencies (for specific analyses)

- scCODA (compositional analysis in `attempt_scrna_hyak.qmd`)
- OrganAge (biological age estimation in `attempt_plasma_somascan_analysis.qmd`)

## Data Availability

The analyses in this repository use data from the ATTEMPT clinical trial. Due to the sensitive nature of patient-level clinical and genomic data, the following access model applies:

- **Single-cell RNA-seq data**: Raw and processed scRNA-seq data are deposited in GEO/SRA (accession number provided in the manuscript).
- **SomaScan proteomics data**: Proteomics datasets are available as described in the manuscript's Data and Materials Availability statement.
- **Clinical data**: De-identified clinical variables are available upon reasonable request to the corresponding author, subject to IRB approval and data use agreements, as required for human subjects research.
- **KPMP reference data**: Kidney Precision Medicine Project data used for cell type annotation and cross-study comparisons are publicly available through the KPMP data portal (https://atlas.kpmp.org).
- **Pathway databases**: KEGG, Reactome, Gene Ontology, and MSigDB Hallmark gene sets are accessed programmatically via the `msigdbr` and `fgsea` R packages.

### Credentials

A `keys.json` file is required for accessing data stored on cloud infrastructure (AWS S3). This file is not included in the repository for security reasons. Contact the corresponding author for data access instructions.

## Reproducibility Notes

- The scRNA-seq analyses (`attempt_scrna_hyak.qmd`) were originally run on the University of Washington Hyak high-performance computing cluster due to memory requirements. These analyses may require 64+ GB RAM.
- NEBULA mixed-model analyses are parallelized and may take several hours per cell type on a standard workstation.
- Quarto (`.qmd`) files can be rendered with `quarto render <filename>.qmd` or run interactively in RStudio.
- All scripts source `attempt_functions.R` at the top, which must be in the project root directory.

## Contact

For questions about the code or data access, please contact the corresponding author as listed in the manuscript.

## License

This code is provided for reproducibility of the analyses described in the associated manuscript. Please cite the manuscript if you use or adapt this code.
