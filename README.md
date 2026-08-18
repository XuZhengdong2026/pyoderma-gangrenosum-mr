# Analysis code and data for "Genetic liability to Crohn's disease, but not ulcerative colitis, is associated with increased risk of pyoderma gangrenosum"

This repository contains the analysis code and result tables for the Mendelian randomization (MR), mediation, specificity, transcriptomic (including the independent GSE280220 validation), in silico modelling and reproducibility analyses reported in the manuscript.

**Version:** v1.1.0 (2026-08-18). See `CHANGELOG.md` for changes since the
initial release (v1.0.0, DOI: 10.5281/zenodo.21947642).

## Contents

- `code/` — R and Python scripts used for all analyses and figures
- `results/` — key result tables in CSV format
- `data_manifest/` — input data sources and accessions
- `.zenodo.json` — metadata for Zenodo/GitHub archiving

## Reproducibility notes

- The independent validation script `validate_gse280220.py` reads the
  GSE280220 processed matrix (download from NCBI GEO) and writes outputs under
  `PG_ROOT` (default: the original machine's `outputs` directory); set
  `PG_ROOT` to your own output directory to run elsewhere.
- `save_mrpresso_r5.R` uses a fixed random seed (20260814); the archived
  MR-PRESSO R5 output is in `results/CD_UC_R5_mrpresso.csv/.txt/.rds` and
  `Table_S15_mrpresso_R5.csv`.
- `make_reproducibility_tables.py` assembles the ODE summary, full Enrichr
  tables and MR-PRESSO R5 table; set `TX_DATA_OUT` to the directory containing
  the raw Enrichr output if it differs from the default.

## Software requirements

- R >= 4.6.1 with packages: `TwoSampleMR` (0.7.9), `ieugwasr`, `MRPRESSO`, `dplyr`, `coloc`, `hugene10sttranscriptcluster.db`
- Python >= 3.12 with packages: `numpy`, `scipy`, `matplotlib`, `pandas`, `pydeseq2` (0.5.4), `Pillow`, `pypdf`

An OpenGWAS JWT token is required for some API calls (`Sys.setenv(OPENGWAS_JWT = "<token>")` or the `OPENGWAS_JWT` environment variable; see https://api.opengwas.io).

## Input data

| Dataset | Accession / ID | Use |
|---|---|---|
| IBD GWAS | IEU OpenGWAS `ieu-a-294` (Liu et al. 2015) | IBD exposure |
| CD / UC GWAS | `ebi-a-GCST004132` / `ebi-a-GCST004133` (de Lange et al. 2017) | CD/UC exposure |
| PG outcome | FinnGen R5 `finn-b-L12_PYODERMA`; FinnGen R12 `L12_PYODERMA` | MR outcome |
| Blood-cell traits | GWAS Catalog `GCST90002394`-`GCST90002407` | Mediators |
| Immune-cell subsets | `ebi-a-GCST90001391`-`ebi-a-GCST90001560` (Orrù et al. 2020, selected subsets) | Mediators |
| Dermatological outcomes | `ebi-a-GCST90019017` (psoriasis); `ebi-a-GCST90027161` (atopic dermatitis); FinnGen R5 vitiligo, hidradenitis suppurativa, alopecia areata, erythema nodosum | Specificity analyses |
| PG skin transcriptome | NCBI GEO `GSE298908` | Skin DEG / programmes |
| PG skin transcriptome (validation) | NCBI GEO `GSE280220` (NanoString nCounter) | Independent key-gene and programme validation |
| IBD gut transcriptome | NCBI GEO `GSE75214` (GPL6244) | Gut-skin comparison |
| Plasma pQTL | `prot-a-749`, `prot-a-1504`, `prot-a-3029`, `ebi-a-GCST90011994`, `ebi-a-GCST90019460` | Drug-target MR |

## Suggested run order

1. `mr_pipeline.R` / `mr_run_local_clump.R` — instrument selection and main MR
2. `mr_r12_sensitivity.R` — FinnGen R12 replication
3. `twostep_mr.R` / `cd_mediation_skin_specificity.R` — two-step mediation and specificity analyses
4. `immune_subsets_twostep.R` — immune-cell subset mediation
5. `run_coloc.R` / `run_susie_coloc.R` — colocalization
6. `transcriptome_prep.py` → `map_gpl6244_probes.R` → `gut_skin_compare.py` — transcriptomic analyses
7. `pg_network_knockout.py`, `ode_calibrate.py`, `ode_graded_ablation.py`, `ode_uncertainty.py` — in silico model
8. `validate_gse280220.py` — independent PG skin transcriptome validation
9. `save_mrpresso_r5.R`, `make_reproducibility_tables.py` — reproducibility outputs (MR-PRESSO R5, ODE summary, full Enrichr tables)
10. Figure scripts (`fig*.py`, `make_ibd_figures.py`)

## License

MIT (see `LICENSE`). All input datasets remain subject to their original data-use policies (OpenGWAS, FinnGen, GWAS Catalog, NCBI GEO).

## Citation

Please cite the associated manuscript once published, and this repository by its GitHub URL and DOI: https://github.com/XuZhengdong2026/pyoderma-gangrenosum-mr (version v1.1.0; DOI: `10.5281/zenodo.21991327`).

Previous version: 10.5281/zenodo.21947642 (v1.0.0).
