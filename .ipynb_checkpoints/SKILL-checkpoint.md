---
name: lovo-pipeline
description: "Run and summarize gene-based LOVO (leave-one-variant-out) analyses for regenie outputs. Use when reproducing the DR_FU/DR_FU_years LOVO pipeline, matching tests to thresholded gene-based results, generating driven/not-driven calls, and producing final gene-level tables with expression annotations."
---

# LOVO Pipeline (Regenie Gene-Based)

Use this skill to reproduce the full LOVO workflow used in this project. Follow the steps in order and keep outputs under `results_gene/lovo/`.

## Inputs (project-specific)

- Gene-based results (DR_FU): `results_gene/DR_FU_p1e-5.csv`
- LOVO jobs list: `results_gene/lovo_jobs.tsv`
- LOVO runner script (DR_FU, full tests): `results_gene/run_lovo_dr_fu_fulltests.sh`
- DR_FU_years LOVO already done (burden-only)
- Expression sources:
  - Tabula Sapiens Eye: `/mnt/project/Analysis_XY/scRNASeq/TabulaSapiens/Eye_TSP1_30_version2d_10X_smartseq_scvi_Nov122024_updated.h5ad`
  - Retina CD31 (mouse): `/mnt/project/Analysis_XY/scRNASeq/Retina_EC_Mice/retina_CD31_annotated.h5ad`

## Critical rules

1. **Test alignment**: Only judge LOVO for tests that reached the threshold in the original gene-based results.
2. **GENE_P**: Record it, but **do not** do LOVO-driven calls for GENE_P.
3. **DR_FU test set**: Use `--vc-tests skat,skato,skato-acat,acatv,acato` with `--mask-lovo`.
   - This produces tests: ADD, ADD-ACATV, ADD-ACATO, ADD-SKAT, ADD-SKATO, ADD-SKATO-ACAT.
   - Do **not** use `--joint` with `--mask-lovo` (regenie disallows it).
4. **Driven definition**: `single_variant_driven = (max LOVO P > 1e-4)` per test.

## Workflow summary

### A) Build/confirm LOVO jobs

- Jobs come from significant mask IDs (ID = `setid.mask.aaf`).
- For DR_FU use `results_gene/lovo_jobs.tsv` (already prepared).

### B) Run DR_FU LOVO (full tests)

- Use `results_gene/run_lovo_dr_fu_fulltests.sh`.
- This script runs all DR_FU LOVO jobs with `--vc-tests skat,skato,skato-acat,acatv,acato`.
- Outputs go to `results_gene/lovo/*.regenie.gz`.

### C) Parse results and build gene-level+test table

- Parse original DR_FU results, **split ID with last 2 dots**:
  - `setid = ID.split('.')[0]`
  - `mask = ID.split('.')[1]`
  - `aaf = '.'.join(ID.split('.')[2:])`
- Determine per gene+test **least_p_mask** and **original_p**.
- Map to LOVO output and compute **lovo_max_p** for the same TEST.
- Apply driven rule (skip for GENE_P).
- Add beta.

### D) Add expression annotations

Compute from the two h5ad files:

**Tabula Sapiens (Eye)**
- mean_expression_retina: mean across all cells
- max_expression_cell_TS: cell_ontology_class with highest mean
- retinal_endothelial_level: mean in `retinal blood vessel endothelial cell`

**CD31 Retina (mouse)**
- mean_expression_CD31_retina: mean across all cells
- max_expression_cell_CD31: cell_type with highest mean
- endothelial_level: mean across cell_type containing "EC"
- DM/Ctrl means: by obs['group']
- Wilcoxon p-value: ranksums(DM, Ctrl)

## Outputs

- Final table (gene-level + test):
  - `results_gene/lovo/lovo_gene_expression_table_by_test.tsv`
- Comparison vs old driven calls:
  - `results_gene/lovo/compare_old_vs_new_driven.tsv`
  - `results_gene/lovo/compare_old_vs_new_driven_changed.tsv`

