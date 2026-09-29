# IMPACT: immune dysfunction in myeloma precursors

Code to reproduce every figure panel in **"Myeloma precursors are associated with suboptimal
immune responses to vaccination"**.

The IMPACT study profiles immune function in myeloma precursor conditions (MGUS and SMM) using
SARS-CoV-2 vaccination as a controlled in vivo challenge: longitudinal serology in up to 731
individuals, single-cell RNA sequencing, T cell receptor sequencing
and plasma proteomics.

## Conventions

The significance threshold throughout is q (or p) < 0.1. Cross-sectional disease-vs-HD contrasts
are age- and sex-adjusted by rank-based ANCOVA. Paired pre/post comparisons are inherently
controlled for age and sex, since each individual is their own baseline. Ordered trends across the
disease continuum use the Jonckheere-Terpstra test on age+sex residuals (helper in `config.R`).

Folders and scripts carry numeric prefixes matching their manuscript order, and each script is
named for the panel it produces.

## Main figures

| Panel | Script | Description |
|-------|--------|-------------|
| 1B | `01_Figure1/01_Figure1B.R` | Spike IgG titers 2-8 weeks post-2nd dose |
| 1C | `01_Figure1/02_Figure1C.R` | Spike IgG titers 2-4 months post-2nd dose |
| 1D | `01_Figure1/03_Figure1D.R` | Per-individual serial titer trajectories |
| 1E | `01_Figure1/04_Figure1E.R` | Mixed-effects waning slope (MGUS vs untreated SMM) |
| 1F | `01_Figure1/05_Figure1F.R` | Spike IgG titers after the 3rd dose |
| 2B | `02_Figure2/01_Figure2B.R` | MMR antibody titers, HD vs SMM |
| 2D | `02_Figure2/02_Figure2D.R` | Panel-wide Olink screen, with APRIL as the top hit |
| 2E | `02_Figure2/03_Figure2E.R` | Plasma APRIL (TNFSF13), paired pre/post |
| 2F | `02_Figure2/04_Figure2F.R` | Post-vaccination APRIL vs serum M-spike |
| 2G | `02_Figure2/05_Figure2G.py` | Pre vs post-teclistamab M-spike, sBCMA, APRIL and BAFF (Immuno-PRISM, NCT05469893) |
| 3B | `03_Figure3/01_Figure3B.py` | Full UMAP with Annotation_Level_2 labels |
| 3C | `03_Figure3/02_Figure3C.py` | Myeloid TNFSF13 expression, pre vs post |
| 3D | `03_Figure3/03_Figure3D.py` | APRIL module in non-malignant plasma cells across NBM/MGUS/SMM/NDMM, marrow and blood (SWIFT-seq) |
| 3E | `03_Figure3/04_Figure3E.py` | APRIL module in tumor vs non-malignant plasma cells, paired within sample (SWIFT-seq) |
| 4B | `04_Figure4/01_Figure4B.R` | Spike-specific TCR clonotypes, paired pre/post |
| 4C | `04_Figure4/02_Figure4C.R` | CEF recall clonotypes, the control for 4B |
| 4E | `04_Figure4/03_Figure4E.R` | IFN-gamma ELISPOT, HD vs SMM |
| 5A, 5B | `05_Figure5/01_Figure5AB.R` | Plasma IL-1B and IL-18, paired pre/post |
| 5C | `05_Figure5/02_Figure5C.py` | IL-1B response signature score (scRNA-seq) |
| 5D, 5E, 5F | `05_Figure5/03_Figure5DEF.R` | Plasma DDX58, NUB1 and MMP7, paired pre/post |

## Supplementary figures

| Panel | Script | Description |
|-------|--------|-------------|
| S1 | `06_SupplementalFigures/01_SupFig1.py` | Vaccine response by SMM 20/2/20 risk tier, treatment-naive SMM |
| S2 | `06_SupplementalFigures/02_SupFig2.py` | Per-lineage annotation UMAPs and marker heatmaps |
| S3A | `06_SupplementalFigures/03_SupFig3A.py` | Bone-marrow myeloid TNFSF13 (Zavidij, GSE124310) |
| S3B | `06_SupplementalFigures/03_SupFig3B.py` | Figure 3C on shipped samples only |
| S4A, S4B | `06_SupplementalFigures/04_SupFig4AB.py` | APRIL signature in stimulated plasmablasts (GSE173644) |
| S4C | `06_SupplementalFigures/04_SupFig4C.py` | APRIL module in marrow plasma cells (Boiarsky, GSE193531) |
| S5 | `06_SupplementalFigures/05_SupFig5.R` | Shipment control for Figure 4B |
| S6 | `06_SupplementalFigures/06_SupFig6.py` | Shipment control for Figure 5C |

Each folder has its own README with the panel inputs and statistics.

## Data

### Zenodo

Processed data are deposited on Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18989222.svg)](https://doi.org/10.5281/zenodo.18989222)

| File | Description | Size |
|------|-------------|------|
| `scRNAseq_IMPACT_Zenodo.h5ad` | De-identified single-cell object: 1,433,497 cells, 42,090 genes, normalized log1p `X`, int64 `counts` layer, `obsm['X_umap']`, and obs `Annotation_Level_1`, `Annotation_Level_2`, `Diagnosis`, `Timepoint`, `TreatmentStatus`. QC-failed cells, doublets and platelets are retained under explicit labels and filtered by the scripts; CLL is not included. Excluding QC-failed cells and doublets leaves 1,109,633 cells, which is what Figure 3B plots; platelets are removed on top of that by the analyses that quantify immune cells. | 9.6 GB |
| `subclusters/` | Five per-lineage objects (B, T, monocyte, dendritic, NK) used by Supplementary Figure 2 | 100 MB |
| Tabular files | De-identified serology, Olink, TCR, metadata and external-cohort inputs | tens of MB |

`README_DATA.md` in the deposit is the authoritative file list.

### Raw sequencing data

Raw sequencing data are deposited at GEO under accession GSEXXXXXX.

### External datasets

| Accession | Reference | Used in |
|-----------|-----------|---------|
| [phs003855](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs003855) (dbGaP, controlled access) | Lightbody et al., *Nat Cancer* 2025 (SWIFT-seq) | Figures 3D and 3E |
| [GSE124310](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124310) | Zavidij et al., *Nat Cancer* 2020 | Supplementary Figure 3C |
| [GSE173644](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173644) | Stephenson et al., *J Immunol* 2022 | Supplementary Figure 4A |
| [GSE193531](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193531) | Boiarsky et al., *Nat Commun* 2022 | Supplementary Figure 4B |

The phs003855 per-sample summary, the GSE193531 metadata and count matrix, and the GSE173644 time
course all ship in the Zenodo deposit. The GSE124310 sample matrices are the one exception:
download the GSE124310 supplementary file from GEO and extract it to
`data/external/zavidij_bm/matrices/`, one directory per sample.

### Local layout

Download the deposit and place the tabular files under `data/`:

```
data/
├── elisa/          spike and MMR titers, cohort demographics, ELISPOT counts
├── olink/          paired pre/post panel, cytokines, tumor-burden summary
├── tcr/            ClusTCR clonotype proportions, shipping status
├── metadata/       Supplementary_Table_4_scRNAseq_sample_list.csv, metadata_deidentified.csv
├── external/       SWIFT-seq summary, GSE193531, GSE173644, teclistamab Olink, zavidij_bm/
├── hvg_2678_genes.txt
└── il1b_response_genes_human.csv
```

The single-cell objects go anywhere you like; point `SCRNA_DIR` at their directory.

## Setup

```bash
conda create -n impact python=3.12
conda activate impact
pip install -r requirements.txt      # exact pinned versions, see the note below
```

```r
source("environment.R")              # CRAN packages used by the R scripts, pinned
```

Both dependency files pin **exact versions**, matching the environment the published figures
were generated in. This matters beyond tidiness: `scanpy`'s `sc.tl.score_genes()` chooses its
control gene bins in a version-sensitive way, so an unpinned install changes the module scores
in Figure 3E and Figure 5C.

Then fetch the data:

```bash
python3 tools/fetch_zenodo_data.py                  # everything, about 10.4 GB
python3 tools/fetch_zenodo_data.py --tabular-only   # about 190 MB, skips the single-cell objects
```

This downloads the deposit and writes each file to the path the scripts expect. The record
stores files flat, with the directory encoded in the filename (`elisa_elisa_mmr.csv`), and the
scripts read a nested tree (`data/elisa/elisa_mmr.csv`), so the script unflattens as it goes.
Checksums are verified, files already present are left alone, and an interrupted run resumes
where it stopped. `--tabular-only` is enough for every panel except Figure 3B, 3C and 5C and
Supplementary Figures 2, 3B and 6, which read the single-cell objects.

Then set the data paths:

```python
# config.py
SCRNA_DIR = Path("/path/to/scrnaseq_data")
```

```r
# config.R
SCRNA_DIR <- "/path/to/scrnaseq_data"
```

## Running

Scripts are independent, with one exception: `05_Figure5/03_Figure5DEF.R` sources
`01_Figure5AB.R`, so run 5AB first.

These read the single-cell object and need the memory to match: `03_Figure3/01_Figure3B.py`,
`03_Figure3/02_Figure3C.py`, `05_Figure5/02_Figure5C.py`,
`06_SupplementalFigures/02_SupFig2.py`, `03_SupFig3B.py` and `06_SupFig6.py`. Every other script
runs from the tabular inputs.

## Reproducibility check

`.github/workflows/run-reproducibility-check.yml` runs `tools/check_reproducibility.py` weekly and
appends the score to `reproducibility.md`. It scores repository structure, not scientific
correctness. Locally:

```bash
python3 tools/check_reproducibility.py .
```
