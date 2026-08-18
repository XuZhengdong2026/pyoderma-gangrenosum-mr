# Changelog

## v1.1.1 (2026-08-18, post-archive maintenance; not published on Zenodo)

### Fixed
- `results/MR_results_summary.csv` replaced with the final strict-LD-clumped
  estimates (IBD/CD/UC, FinnGen R5 and R12); added
  `results/MR_instrument_F_stats.csv` (mean/min F statistics).
- Added `results/S01_R12_sensitivity_reverse_MR.csv` (R12 sensitivity MR plus
  relaxed-threshold reverse MR, Steiger-filtered) and the raw
  `results/PG_to_*_steiger_mr.csv` outputs.
- Added `results/Table_S15_mrpresso.csv` (MR-PRESSO for CD and UC against
  FinnGen R5 and R12).
- Manuscript v5.1 errata: Table 1 GSE75214 sample counts, Table 3 white-blood-cell
  step-2 SNP count, GSE298908 GEO reference title, NanoString panel naming,
  Figure 1 working-model legend, author ORCIDs and word-count labels.

## v1.1.0 (2026-08-18)

Archived at Zenodo (DOI: 10.5281/zenodo.21991327).

### Added
- Independent PG skin transcriptome validation in **GSE280220** (NanoString
  nCounter Human Inflammation panel v2; 5 PG lesions, 2 Sweet syndrome lesions,
  8 healthy skin samples):
  - `code/validate_gse280220.py`
  - `results/gse280220_gene_validation.csv`
  - `results/gse280220_key_genes` content in `Table_S12a_gse280220_key_genes.csv`
  - `results/gse280220_programme_validation.csv` (`Table_S12b_gse280220_programmes.csv`)
  - `results/gse280220_validation_summary.txt`
- Reproducibility outputs:
  - `Table_S13_ode_ablation_summary.csv` (deterministic, calibrated,
    parameter-uncertainty and graded-ablation ODE results)
  - `Table_S14_enrichr_full_up.csv` / `Table_S14_enrichr_full_down.csv`
    (full raw Enrichr output)
  - `Table_S15_mrpresso_R5.csv`, `CD_UC_R5_mrpresso.csv`,
    `CD_UC_R5_mrpresso.txt`, `CD_UC_R5_mrpresso.rds` (seed-fixed MR-PRESSO,
    random seed 20260814)
  - `code/save_mrpresso_r5.R`, `code/make_reproducibility_tables.py`
- `data_manifest/input_data_manifest.csv` now lists GSE280220.
- `.zenodo.json` metadata for Zenodo/GitHub archiving.

### Changed
- README title and run order updated for the validation and reproducibility
  scripts; new scripts read input/output paths from environment variables
  (`PG_ROOT`, `TX_DATA_OUT`) so they can be run outside the original machine.
- Results wording in the manuscript v5.1: MR-PRESSO R5 global P uses the
  seed-fixed value (0.010) consistent with `Table_S15_mrpresso_R5.csv`.

## v1.0.0 (2026-08-15)

Initial archived release with MR, mediation, specificity, transcriptomic,
gut-skin comparison and in silico modelling code and results
(DOI: 10.5281/zenodo.21947642).
