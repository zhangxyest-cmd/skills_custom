---
name: pyscenic
description: "Run pySCENIC pipelines (GRN→regulons→AUCell) and summarize regulon activity per cluster/condition for scRNA-seq. Use when running pySCENIC on AnnData/loom, generating GRN/regularons/AUC, or producing DM vs Ctrl (or other group) reports with plots and tables."
---

# pySCENIC workflow skill

## Quick start

Use the bundled scripts to run a full cluster pipeline and generate reports.

1) Run a cluster (full genes) pipeline:

```bash
python skills/pyscenic/scripts/run_cluster_full.py \
  --h5ad scRNASeq/retina_ec_annotated_noRgs5.h5ad \
  --cluster-col EC_subtype \
  --cluster-name IFN \
  --out-dir pySCENIC/IFN_full \
  --tf-list scenic/resources/allTFs_mm.txt \
  --db scenic/db/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
  --db scenic/db/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
  --motif-annot scenic/resources/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl \
  --workers 4
```

2) Generate summary tables and HTML report:

```bash
python skills/pyscenic/scripts/make_report.py \
  --auc-loom pySCENIC/IFN_full/auc_mtx.loom \
  --out-dir pySCENIC/IFN_full \
  --group-col group \
  --title "IFN pySCENIC Report (Full Genes)"
```

## Key outputs

- `adj.tsv` (GRN links)
- `regulons.gmt` (regulons)
- `auc_mtx.loom` (AUCell matrix)
- `*_regulon_mean_auc.tsv` (mean AUC by regulon)
- `*_dm_vs_ctrl.tsv` (or group comparison table)
- `report.html`

## Notes

- AUCell binarization is patched in the script for NumPy/Pandas deprecations.
- Use `--workers` to control parallelism.
- For multiple clusters, run the pipeline per cluster and then generate reports.

## Bundled scripts

- `scripts/run_cluster_full.py`: full pySCENIC pipeline for a single cluster
- `scripts/make_report.py`: create TSV summaries and HTML report from `auc_mtx.loom`
