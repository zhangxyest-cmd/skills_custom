#!/usr/bin/env python
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import loompy


def parse_args():
    p = argparse.ArgumentParser(description='Generate TSV summaries and HTML report from auc_mtx.loom')
    p.add_argument('--auc-loom', required=True)
    p.add_argument('--out-dir', required=True)
    p.add_argument('--group-col', default='group')
    p.add_argument('--title', default='pySCENIC Report')
    p.add_argument('--prefix', default=None, help='Prefix for output TSVs')
    return p.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = out_dir / 'figures'
    fig_dir.mkdir(exist_ok=True)

    lf = loompy.connect(str(args.auc_loom), validate=False)
    auc = lf.ca['RegulonsAUC']
    cell_ids = lf.ca['CellID'].astype(str)
    group = lf.ca[args.group_col].astype(str)
    lf.close()

    auc_df = pd.DataFrame.from_records(auc, index=cell_ids)
    regulon_names = auc_df.columns.tolist()

    mean_auc = auc_df.mean().sort_values(ascending=False)

    prefix = args.prefix or out_dir.name.replace('_full', '')
    mean_path = out_dir / f'{prefix}_regulon_mean_auc.tsv'
    top_path = out_dir / f'{prefix}_top_regulons.tsv'
    dm_path = out_dir / f'{prefix}_dm_vs_ctrl.tsv'

    mean_auc.to_csv(mean_path, sep='\t')
    top_df = mean_auc.head(50).reset_index().rename(columns={'index': 'regulon', 0: 'mean_auc'})
    top_df.to_csv(top_path, sep='\t', index=False)

    from scipy.stats import mannwhitneyu

    dm_idx = group == 'DM'
    ctrl_idx = group == 'Ctrl'

    rows = []
    for reg in regulon_names:
        x = auc_df.loc[dm_idx, reg].values
        y = auc_df.loc[ctrl_idx, reg].values
        try:
            stat, p = mannwhitneyu(x, y, alternative='two-sided')
        except Exception:
            p = np.nan
        rows.append({
            'regulon': reg,
            'mean_DM': float(np.mean(x)),
            'mean_Ctrl': float(np.mean(y)),
            'delta_DM_minus_Ctrl': float(np.mean(x) - np.mean(y)),
            'pval': float(p),
            'n_DM': int(len(x)),
            'n_Ctrl': int(len(y)),
        })

    df = pd.DataFrame(rows)
    # BH FDR
    pvals = df['pval'].to_numpy()
    order = np.argsort(pvals)
    rank = np.empty_like(order)
    rank[order] = np.arange(1, len(pvals) + 1)
    qvals = pvals * len(pvals) / rank
    qvals = np.minimum.accumulate(qvals[order][::-1])[::-1]
    q = np.empty_like(qvals)
    q[order] = np.minimum(qvals, 1.0)
    df['qval'] = q

    # Cohen d
    cohen = []
    for reg in regulon_names:
        x = auc_df.loc[dm_idx, reg].values
        y = auc_df.loc[ctrl_idx, reg].values
        nx, ny = len(x), len(y)
        if nx < 2 or ny < 2:
            cohen.append(np.nan)
            continue
        vx = x.var(ddof=1)
        vy = y.var(ddof=1)
        pooled = ((nx - 1) * vx + (ny - 1) * vy) / (nx + ny - 2)
        if pooled == 0:
            cohen.append(0.0)
        else:
            cohen.append((x.mean() - y.mean()) / np.sqrt(pooled))
    df['cohen_d'] = cohen

    df.sort_values(['qval', 'pval']).to_csv(dm_path, sep='\t', index=False)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns

    plt.figure(figsize=(10, 6))
    sns.boxplot(data=auc_df[top_df['regulon'].head(10)], orient='h')
    plt.title(f'{prefix} top 10 regulon AUC distribution')
    plt.tight_layout()
    plt.savefig(fig_dir / f'{prefix}_top10_boxplot.png', dpi=200)
    plt.close()

    html = f"""
<!DOCTYPE html>
<html>
<head>
<meta charset='utf-8'>
<title>{args.title}</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 24px; }}
img {{ max-width: 100%; height: auto; }}
code {{ background: #f3f3f3; padding: 2px 4px; }}
table {{ border-collapse: collapse; margin-bottom: 24px; }}
th, td {{ border: 1px solid #ccc; padding: 4px 8px; }}
</style>
</head>
<body>
<h1>{args.title}</h1>
<ul>
  <li>Cells: {len(cell_ids)}</li>
  <li>Regulons: {len(regulon_names)}</li>
  <li>DM: {int(dm_idx.sum())}, Ctrl: {int(ctrl_idx.sum())}</li>
</ul>

<h2>Top regulons (mean AUC)</h2>
{top_df.head(20).to_html(index=False)}

<h2>DM vs Ctrl (top 20 by q-value)</h2>
{df.sort_values(['qval','pval']).head(20).to_html(index=False)}

<h2>Top 10 regulon AUC distribution</h2>
<img src='figures/{prefix}_top10_boxplot.png' />

<h2>Output files</h2>
<ul>
  <li><code>{mean_path.name}</code></li>
  <li><code>{top_path.name}</code></li>
  <li><code>{dm_path.name}</code></li>
  <li><code>auc_mtx.loom</code></li>
  <li><code>regulons.gmt</code></li>
</ul>
</body>
</html>
"""

    (out_dir / 'report.html').write_text(html)


if __name__ == '__main__':
    main()
