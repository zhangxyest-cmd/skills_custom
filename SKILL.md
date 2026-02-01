---
name: lovo-pipeline
description: "Run and summarize gene-based LOVO (leave-one-variant-out) analyses for regenie outputs. Use when reproducing the DR_FU/DR_FU_years LOVO pipeline, matching tests to thresholded gene-based results, generating driven/not-driven calls, and producing final gene-level tables with expression annotations."
---

# LOVO Pipeline (Regenie Gene-Based)

Use this skill to reproduce the full LOVO workflow used in this project. Follow the steps in order and keep outputs under `results_gene/lovo/`.

## Inputs (project-specific)

- Regenie log + results folder: `/mnt/project/WGS_XY/Results_gene/251224_dr_mac20`
- LOVO jobs list (generated): `results_gene/lovo/lovo_jobs.tsv`
- LOVO runner (from log): `scripts/run_lovo_from_log.py`
- Log parser: `scripts/parse_regenie_log.py`
- Significant-mask extractor: `scripts/extract_sig_masks.py`
- LOVO job builder: `scripts/build_lovo_jobs.py`
- Gene-level table builder: `scripts/build_gene_test_table.py`
- Expression sources:
  - Tabula Sapiens Eye: `/mnt/project/Analysis_XY/scRNASeq/TabulaSapiens/Eye_TSP1_30_version2d_10X_smartseq_scvi_Nov122024_updated.h5ad`
  - Retina CD31 (mouse): `/mnt/project/Analysis_XY/scRNASeq/Retina_EC_Mice/retina_CD31_annotated.h5ad`

## Critical rules

1. **Test alignment**: Only judge LOVO for tests that reached the threshold in the original gene-based results.
2. **GENE_P**: Record it, but **do not** do LOVO-driven calls for GENE_P.
3. **Joint tests**: `--joint` is **not** allowed with `--mask-lovo`. Any tests that only come from joint ACAT (e.g., `ADD-ACATV-ACAT`, `ADD-BURDEN-ACAT`) are recorded but **not** LOVO-run.
4. **Test → vc-tests mapping** (minimal):
   - `ADD` → no `--vc-tests`
   - `ADD-SKAT` → `--vc-tests skat`
   - `ADD-SKATO` → `--vc-tests skato`
   - `ADD-SKATO-ACAT` → `--vc-tests skato-acat`
   - `ADD-ACATV` → `--vc-tests acatv`
   - `ADD-ACATO` → `--vc-tests acato`
5. **AAF handling**: `--aaf-bins`, `--vc-maxAAF`, and `--mask-lovo` must match the AAF from the significant mask. If `aaf == singleton`, use `--mask-lovo setid,mask,singleton` and keep `--aaf-bins/--vc-maxAAF` at `0.001` (per original run).
6. **Keep list**: only use `--keep /mnt/project/PRS/Ancestry_Inference/cohort_EUR.txt` (other keep entries in logs are ignored).
7. **Trait flags**: use `--bt` or `--t2e` exactly as in the original log.
8. **Driven definition**: `single_variant_driven = (max LOVO P > 1e-4)` per test.

## Workflow

### A) Extract significant masks from regenie results

Default thresholds for this project:
- primary max-p: 1e-5
- reporting thresholds: 2.5e-6 and 7.9e-7 (computed in output table)

```bash
scripts/extract_sig_masks.py \
  --results-dir /mnt/project/WGS_XY/Results_gene/251224_dr_mac20 \
  --trait DR_FU \
  --max-p 1e-5 \
  --out /opt/notebooks/results_gene/lovo/sig_masks_DR_FU.tsv \
  --out-all-thresholds
```

### B) Build LOVO jobs (based on significant tests)

```bash
scripts/build_lovo_jobs.py \
  --sig /opt/notebooks/results_gene/lovo/sig_masks_DR_FU.tsv \
  --out /opt/notebooks/results_gene/lovo/lovo_jobs.tsv
```

### C) Run LOVO using original log settings

```bash
scripts/run_lovo_from_log.py \
  --log /mnt/project/WGS_XY/Results_gene/251224_dr_mac20/dr_gene_t2d_chr10_eur_mac20_firth.log \
  --jobs /opt/notebooks/results_gene/lovo/lovo_jobs.tsv \
  --lovo-dir /opt/notebooks/results_gene/lovo \
  --regenie /opt/notebooks/regenie_v4.1.gz_x86_64_Linux_mkl
```

### D) Build gene-level + test table (with expression)

```bash
scripts/build_gene_test_table.py \
  --gene-results /opt/notebooks/results_gene/lovo/sig_masks_DR_FU.tsv \
  --lovo-jobs /opt/notebooks/results_gene/lovo/lovo_jobs.tsv \
  --lovo-dir /opt/notebooks/results_gene/lovo \
  --ts /mnt/project/Analysis_XY/scRNASeq/TabulaSapiens/Eye_TSP1_30_version2d_10X_smartseq_scvi_Nov122024_updated.h5ad \
  --cd31 /mnt/project/Analysis_XY/scRNASeq/Retina_EC_Mice/retina_CD31_annotated.h5ad \
  --out /opt/notebooks/results_gene/lovo/lovo_gene_expression_table_by_test.tsv
```

## Outputs

- Final table (gene-level + test):
  - `results_gene/lovo/lovo_gene_expression_table_by_test.tsv`
- Significant masks (from regenie):
  - `results_gene/lovo/sig_masks_DR_FU.tsv`
- LOVO jobs:
  - `results_gene/lovo/lovo_jobs.tsv`
