#!/usr/bin/env python
import argparse
import os
from pathlib import Path

os.environ.setdefault('NUMBA_CACHE_DIR', '/tmp/numba_cache')
os.environ.setdefault('JOBLIB_TEMP_FOLDER', '/tmp')
Path(os.environ['NUMBA_CACHE_DIR']).mkdir(parents=True, exist_ok=True)

import numpy as np
import pandas as pd
import anndata as ad
import loompy
from scipy.sparse import csr_matrix


def run_cmd(cmd, log_path):
    import subprocess
    with open(log_path, 'a') as lf:
        lf.write('\n$ ' + ' '.join(cmd) + '\n')
        lf.flush()
        proc = subprocess.run(cmd, stdout=lf, stderr=lf)
        if proc.returncode != 0:
            raise RuntimeError(f"Command failed: {' '.join(cmd)}")


def parse_args():
    p = argparse.ArgumentParser(description='Run pySCENIC full pipeline for a single cluster')
    p.add_argument('--h5ad', required=True)
    p.add_argument('--cluster-col', required=True)
    p.add_argument('--cluster-name', required=True)
    p.add_argument('--out-dir', required=True)
    p.add_argument('--tf-list', required=True)
    p.add_argument('--db', action='append', required=True)
    p.add_argument('--motif-annot', required=True)
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--seed', type=int, default=777)
    return p.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    loom_path = out_dir / f"{args.cluster_name}_counts.loom"
    tf_list_exp = out_dir / 'tfs_in_genes.txt'
    adj = out_dir / 'adj.tsv'
    reg_gmt = out_dir / 'regulons.gmt'
    auc_loom = out_dir / 'auc_mtx.loom'

    if not loom_path.exists() or loom_path.stat().st_size == 0:
        adata = ad.read_h5ad(args.h5ad)
        adata = adata[adata.obs[args.cluster_col].astype(str) == args.cluster_name].copy()
        if 'counts' in adata.layers:
            counts = adata.layers['counts']
        else:
            counts = adata.X
        if not hasattr(counts, 'tocsr'):
            counts = csr_matrix(counts)
        genes = adata.var_names.astype(str).values
        cells = adata.obs_names.astype(str).values
        row_attrs = {'Gene': genes}
        col_attrs = {
            'CellID': cells,
            args.cluster_col: adata.obs[args.cluster_col].astype(str).values,
            'group': adata.obs['group'].astype(str).values,
        }
        loompy.create(str(loom_path), counts.T, row_attrs, col_attrs)

    if not tf_list_exp.exists() or tf_list_exp.stat().st_size == 0:
        lf = loompy.connect(str(loom_path), validate=False)
        genes = set(map(str, lf.ra['Gene']))
        lf.close()
        tfs = [line.strip() for line in open(args.tf_list) if line.strip() in genes]
        with open(tf_list_exp, 'w') as f:
            f.write('\n'.join(tfs) + ('\n' if tfs else ''))

    if not adj.exists() or adj.stat().st_size == 0:
        cmd = [
            'pyscenic', 'grn', str(loom_path), str(tf_list_exp),
            '-o', str(adj),
            '--method', 'grnboost2',
            '--num_workers', str(args.workers),
            '--seed', str(args.seed),
            '--gene_attribute', 'Gene',
            '--cell_id_attribute', 'CellID'
        ]
        run_cmd(cmd, out_dir / 'pyscenic_grn.log')

    if not reg_gmt.exists() or reg_gmt.stat().st_size == 0:
        cmd = [
            'pyscenic', 'ctx', str(adj),
            *args.db,
            '--annotations_fname', str(args.motif_annot),
            '--expression_mtx_fname', str(loom_path),
            '--mode', 'custom_multiprocessing',
            '--num_workers', str(args.workers),
            '--gene_attribute', 'Gene',
            '--cell_id_attribute', 'CellID',
            '-o', str(reg_gmt)
        ]
        run_cmd(cmd, out_dir / 'pyscenic_ctx.log')

    if not auc_loom.exists() or auc_loom.stat().st_size == 0:
        if not hasattr(np, 'msort'):
            np.msort = np.sort
        if not hasattr(np, 'float'):
            np.float = float
        if not hasattr(pd.Series, 'iteritems'):
            pd.Series.iteritems = pd.Series.items

        from pyscenic.cli.pyscenic import aucell_command

        class FakeArg:
            def __init__(self, name):
                self.name = name

        args2 = type('Args', (), {})()
        args2.expression_mtx_fname = FakeArg(str(loom_path))
        args2.signatures_fname = FakeArg(str(reg_gmt))
        args2.output = FakeArg(str(auc_loom))
        args2.transpose = 'no'
        args2.weights = 'no'
        args2.num_workers = args.workers
        args2.seed = args.seed
        args2.rank_threshold = 5000
        args2.auc_threshold = 0.05
        args2.nes_threshold = 3.0
        args2.cell_id_attribute = 'CellID'
        args2.gene_attribute = 'Gene'
        args2.sparse = True

        aucell_command(args2)


if __name__ == '__main__':
    main()
