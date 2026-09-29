# tools/

Helper scripts that are not figure panels.

| Script | Output |
|--------|--------|
| `build_supp_table_1_demographics.py` | `tables/Supplementary_Table_1_Demographics.xlsx`, the demographics of the 731-participant serology cohort by disease group, with each group tested against HD |
| `build_supp_table_il1b_genes.py` | `tables/Supplementary_Table_7_IL1B-response-genes.csv`, the IL-1B response gene list behind Figure 5C, flagged for presence in the object and for survival into the analysis universe |
| `check_reproducibility.py` | A 0-100 repository structure score, run weekly by CI |
| `fetch_zenodo_data.py` | A populated `data/` directory, downloaded from the Zenodo deposit and unflattened into the layout the figure scripts expect |

```bash
python3 tools/build_supp_table_1_demographics.py
python3 tools/build_supp_table_il1b_genes.py
```

## fetch_zenodo_data.py

Downloads the deposited data and places every file where the scripts look for it. The record
stores files flat with the directory encoded in the filename (`elisa_elisa_mmr.csv`), while the
scripts read a nested tree (`data/elisa/elisa_mmr.csv`); the script maps one to the other, so
nothing has to be renamed or moved by hand.

The record is resolved from the concept DOI `10.5281/zenodo.18989222`, which always points at
the latest version, so the script does not go stale when a new version is published.

```bash
python3 tools/fetch_zenodo_data.py                  # everything, 40 files, about 10.4 GB
python3 tools/fetch_zenodo_data.py --tabular-only   # about 190 MB, no single-cell objects
python3 tools/fetch_zenodo_data.py --dry-run        # list what is missing, without downloading
python3 tools/fetch_zenodo_data.py --dest /data/x   # somewhere other than ./data
python3 tools/fetch_zenodo_data.py --record 1234567 # a specific version instead of the latest
```

MD5 checksums are verified against the record. A file already present with the right size and
checksum is skipped, so the script is safe to re-run and resumes an interrupted download.
`--tabular-only` covers every panel except Figure 3B, 3C and 5C and Supplementary Figures 2, 3B
and 6, which need the single-cell objects.

While the deposit is under embargo its files are restricted and a share token is required, given
as `--token` or in the `ZENODO_TOKEN` environment variable. No token is needed once the embargo
lifts.

## check_reproducibility.py

Scores the repository on structure and hygiene, not on scientific correctness. Standard library
only, Python 3.10+.

| Category | Weight | Checks |
|----------|--------|--------|
| Step ordering | 20 | Scripts and folders carry numeric prefixes |
| Documentation | 25 | Header blocks in scripts, a README in each folder |
| Path hygiene | 20 | No hardcoded local paths |
| Data handling | 15 | `.gitignore` covers data files, READMEs say where the data live |
| Naming | 10 | No spaces or special characters in file names |
| PHI and credentials | 10 | No identifiers, passwords or API keys in code |

```bash
python3 tools/check_reproducibility.py            # current directory
python3 tools/check_reproducibility.py /path      # a specific directory
python3 tools/check_reproducibility.py --json     # write reproducibility_score.json
python3 tools/check_reproducibility.py --min-score 80   # exit 1 below the threshold
```

A score of 80 or above is good; below 60 means something structural is wrong.
