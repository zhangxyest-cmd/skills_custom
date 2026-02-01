#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import ranksums


def parse_id(id_str: str):
    parts = id_str.split('.')
    if len(parts) >= 3:
        setid = parts[0]
        mask = parts[1]
        aaf = '.'.join(parts[2:])
        return setid, mask, aaf
    return None, None, None


def setid_to_gene(setid: str) -> str:
    parts = setid.split('_')
    if len(parts) >= 3 and parts[-1].startswith('ENST'):
        return '_'.join(parts[1:-1])
    return setid


def load_lovo(path: Path, cache: dict):
    if path in cache:
        return cache[path]
    df = pd.read_csv(path, sep=r'\s+', comment='#')
    df['LOG10P'] = pd.to_numeric(df['LOG10P'], errors='coerce')
    df['P'] = 10 ** (-df['LOG10P'])
    cache[path] = df
    return df


def compute_expression(ts_h5ad: Path, cd31_h5ad: Path, genes):
    # Tabula Sapiens (Eye)
    ad_ts = sc.read_h5ad(ts_h5ad, backed='r')
    var_ts = pd.Index(ad_ts.var_names)
    ts_upper = {g.upper(): g for g in var_ts}

    ts_gene_map = {g: ts_upper.get(g.upper()) for g in genes}
    mapped_ts = [g for g in genes if ts_gene_map[g] is not None]

    ts_group_means = pd.DataFrame()
    ts_overall_means = pd.Series(dtype=float)
    if mapped_ts:
        idx = [var_ts.get_loc(ts_gene_map[g]) for g in mapped_ts]
        X = ad_ts[:, idx].X
        try:
            X = X.toarray()
        except Exception:
            X = np.array(X)
        X = X.astype(float)
        ts_df = pd.DataFrame(X, columns=mapped_ts)
        ts_overall_means = ts_df.mean(axis=0)
        ts_group_means = ts_df.groupby(ad_ts.obs['cell_ontology_class'].values).mean()

    # CD31 retina
    ad_cd = sc.read_h5ad(cd31_h5ad, backed='r')
    var_cd = pd.Index(ad_cd.var_names)
    cd_upper = {g.upper(): g for g in var_cd}

    cd_gene_map = {g: cd_upper.get(g.upper()) for g in genes}
    mapped_cd = [g for g in genes if cd_gene_map[g] is not None]

    cd_group_means = pd.DataFrame()
    cd_overall_means = pd.Series(dtype=float)
    cd_dm_means = pd.Series(dtype=float)
    cd_ctrl_means = pd.Series(dtype=float)

    if mapped_cd:
        idx = [var_cd.get_loc(cd_gene_map[g]) for g in mapped_cd]
        X = ad_cd[:, idx].X
        try:
            X = X.toarray()
        except Exception:
            X = np.array(X)
        X = X.astype(float)
        cd_df = pd.DataFrame(X, columns=mapped_cd)
        cd_overall_means = cd_df.mean(axis=0)
        cd_group_means = cd_df.groupby(ad_cd.obs['cell_type'].values).mean()
        grp = ad_cd.obs['group'].values
        cd_dm_means = cd_df[grp == 'DM'].mean(axis=0)
        cd_ctrl_means = cd_df[grp == 'Ctrl'].mean(axis=0)

    retinal_endothelial_label = 'retinal blood vessel endothelial cell'

    records = []
    for g in genes:
        rec = {'gene': g}
        # TS
        if g in mapped_ts:
            rec['expression_in_TS'] = 'T' if ts_overall_means[g] > 0 else 'F'
            rec['mean_expression_retina'] = float(ts_overall_means[g])
            if not ts_group_means.empty:
                rec['max_expression_cell_TS'] = ts_group_means[g].idxmax()
                rec['expression_in_retinal_endothelial_cell_level'] = (
                    float(ts_group_means.loc[retinal_endothelial_label, g])
                    if retinal_endothelial_label in ts_group_means.index else np.nan
                )
            else:
                rec['max_expression_cell_TS'] = np.nan
                rec['expression_in_retinal_endothelial_cell_level'] = np.nan
        else:
            rec['expression_in_TS'] = 'F'
            rec['mean_expression_retina'] = np.nan
            rec['max_expression_cell_TS'] = np.nan
            rec['expression_in_retinal_endothelial_cell_level'] = np.nan

        # CD31
        if g in mapped_cd:
            rec['expression_in_CD31_mice'] = 'T' if cd_overall_means[g] > 0 else 'F'
            rec['mean_expression_CD31_retina'] = float(cd_overall_means[g])
            if not cd_group_means.empty:
                rec['max_expression_cell_CD31'] = cd_group_means[g].idxmax()
                ec_types = [c for c in cd_group_means.index if 'EC' in str(c)]
                rec['expression_in_endothelial_cell_level'] = (
                    float(cd_group_means.loc[ec_types, g].mean()) if ec_types else np.nan
                )
            else:
                rec['max_expression_cell_CD31'] = np.nan
                rec['expression_in_endothelial_cell_level'] = np.nan

            dm = cd_dm_means.get(g, np.nan)
            ctrl = cd_ctrl_means.get(g, np.nan)
            rec['expression_in_DM'] = float(dm) if pd.notna(dm) else np.nan
            rec['expression_in_Ctrl'] = float(ctrl) if pd.notna(ctrl) else np.nan

            eps = 1e-9
            if pd.notna(dm) and pd.notna(ctrl):
                rec['log2fc_DM_vs_Ctrl'] = float(np.log2((dm + eps) / (ctrl + eps)))
            else:
                rec['log2fc_DM_vs_Ctrl'] = np.nan

            # ranksums
            try:
                # reload vectors lazily for p-value per gene
                # This avoids storing full matrix in memory for all genes
                x = ad_cd[:, var_cd.get_loc(cd_gene_map[g])].X
                try:
                    x = x.toarray().ravel()
                except Exception:
                    x = np.array(x).ravel()
                grp = ad_cd.obs['group'].values
                dm_vals = x[grp == 'DM']
                ctrl_vals = x[grp == 'Ctrl']
                rec['wilcoxon_p_value'] = float(ranksums(dm_vals, ctrl_vals).pvalue)
            except Exception:
                rec['wilcoxon_p_value'] = np.nan
        else:
            rec['expression_in_CD31_mice'] = 'F'
            rec['mean_expression_CD31_retina'] = np.nan
            rec['max_expression_cell_CD31'] = np.nan
            rec['expression_in_endothelial_cell_level'] = np.nan
            rec['expression_in_DM'] = np.nan
            rec['expression_in_Ctrl'] = np.nan
            rec['log2fc_DM_vs_Ctrl'] = np.nan
            rec['wilcoxon_p_value'] = np.nan

        records.append(rec)

    return pd.DataFrame(records)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--gene-results', help='sig_masks_DR_FU.tsv (from extract_sig_masks)')
    p.add_argument('--dr-fu', help='legacy DR_FU_p1e-5.csv')
    p.add_argument('--lovo-jobs', required=True, help='lovo_jobs.tsv')
    p.add_argument('--lovo-dir', required=True, help='results_gene/lovo')
    p.add_argument('--ts', required=True, help='Tabula Sapiens h5ad')
    p.add_argument('--cd31', required=True, help='CD31 h5ad')
    p.add_argument('--out', required=True, help='output TSV')
    args = p.parse_args()

    if args.gene_results:
        orig = pd.read_csv(args.gene_results, sep='\t')
        if 'P' not in orig.columns and 'LOG10P' in orig.columns:
            orig['P'] = 10 ** (-pd.to_numeric(orig['LOG10P'], errors='coerce'))
        if 'setid' not in orig.columns or 'mask' not in orig.columns or 'aaf' not in orig.columns:
            orig[['setid','mask','aaf']] = pd.DataFrame(orig['ID'].apply(parse_id).tolist())
        if 'gene' not in orig.columns:
            orig['gene'] = orig['setid'].map(setid_to_gene)
        mask_df = orig[orig['setid'].notna()].copy()
    else:
        if not args.dr_fu:
            raise SystemExit("Provide --gene-results or --dr-fu")
        orig = pd.read_csv(args.dr_fu)
        orig[['setid','mask','aaf']] = pd.DataFrame(orig['ID'].apply(parse_id).tolist())
        mask_df = orig[orig['setid'].notna()].copy()
        mask_df['gene'] = mask_df['setid'].map(setid_to_gene)

    mask_df['pass_2.5e-6'] = mask_df['P'] <= 2.5e-6
    mask_df['pass_7.9e-7'] = mask_df['P'] <= 7.9e-7
    mask_df['pass_1e-5'] = mask_df['P'] <= 1e-5

    mask_df['mask_label'] = mask_df.apply(lambda r: f"{r['setid']}.{r['mask']}.{r['aaf']}", axis=1)
    idx = mask_df.groupby(['gene','TEST'])['P'].idxmin()
    cols = ['gene','TEST','mask_label','P','BETA','setid','mask','aaf','pass_2.5e-6','pass_7.9e-7','pass_1e-5']
    for c in cols:
        if c not in mask_df.columns:
            mask_df[c] = np.nan
    least = mask_df.loc[idx, cols].copy()
    least.rename(columns={'mask_label':'least_p_mask','P':'original_p','BETA':'beta'}, inplace=True)

    jobs = pd.read_csv(args.lovo_jobs, sep='\t')
    jobs = jobs[jobs['trait']=='DR_FU']
    merge_cols = ['setid','mask','aaf','TEST']
    jobs_cols = ['setid','mask','aaf','TEST','outprefix']
    for c in merge_cols:
        if c not in jobs.columns:
            jobs[c] = np.nan
    least = least.merge(jobs[jobs_cols], on=merge_cols, how='left')

    lovo_dir = Path(args.lovo_dir)
    cache = {}

    max_p_list = []
    for _, r in least.iterrows():
        test = r['TEST']
        if test == 'GENE_P':
            max_p_list.append(np.nan)
            continue
        outprefix = r['outprefix']
        if pd.isna(outprefix):
            max_p_list.append(np.nan)
            continue
        path = lovo_dir / f"{outprefix}_DR_FU.regenie.gz"
        if not path.exists():
            max_p_list.append(np.nan)
            continue
        lovo = load_lovo(path, cache)
        sub = lovo[lovo['TEST'] == test]
        max_p_list.append(float(sub['P'].max()) if not sub.empty else np.nan)

    least['lovo_max_p'] = max_p_list
    least['single_variant_driven'] = np.where(
        (least['TEST'] != 'GENE_P') & (least['lovo_max_p'] > 1e-4), 'T',
        np.where(least['TEST']=='GENE_P','NA','F')
    )

    expr = compute_expression(Path(args.ts), Path(args.cd31), sorted(least['gene'].unique()))

    out = least.merge(expr, on='gene', how='left')
    cols = [
        'gene', 'TEST', 'least_p_mask', 'original_p', 'beta', 'lovo_max_p', 'single_variant_driven',
        'pass_2.5e-6','pass_7.9e-7','pass_1e-5',
        'expression_in_TS','max_expression_cell_TS','expression_in_retinal_endothelial_cell_level','mean_expression_retina',
        'expression_in_CD31_mice','max_expression_cell_CD31','expression_in_endothelial_cell_level','mean_expression_CD31_retina',
        'expression_in_DM','expression_in_Ctrl','log2fc_DM_vs_Ctrl','wilcoxon_p_value'
    ]
    for c in cols:
        if c not in out.columns:
            out[c] = np.nan
    out = out[cols]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep='\t', index=False)


if __name__ == '__main__':
    main()
