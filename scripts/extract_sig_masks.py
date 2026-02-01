#!/usr/bin/env python3
import argparse
import gzip
from pathlib import Path
import pandas as pd
import numpy as np


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


def iter_regenie(path: Path, chunksize: int = 500000):
    return pd.read_csv(path, sep=r'\s+', comment='#', chunksize=chunksize)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--results-dir', required=True)
    p.add_argument('--trait', required=True)
    p.add_argument('--max-p', type=float, default=1e-5)
    p.add_argument('--out', required=True)
    p.add_argument('--out-all-thresholds', action='store_true',
                   help='auto-generate 2.5e-6 and 7.9e-7 outputs based on --out')
    p.add_argument('--out-2p5e-6', help='optional output for P<=2.5e-6')
    p.add_argument('--out-7p9e-7', help='optional output for P<=7.9e-7')
    p.add_argument('--tests', help='comma-separated list of TEST to include (optional)')
    args = p.parse_args()

    res_dir = Path(args.results_dir)
    files = sorted(res_dir.glob(f"*_{args.trait}.regenie.gz"))
    if not files:
        raise SystemExit(f"No regenie.gz files for trait {args.trait} in {res_dir}")

    test_filter = None
    if args.tests:
        test_filter = set([t.strip() for t in args.tests.split(',') if t.strip()])

    out_rows = []
    for path in files:
        for chunk in iter_regenie(path):
            if 'LOG10P' not in chunk.columns:
                continue
            chunk['LOG10P'] = pd.to_numeric(chunk['LOG10P'], errors='coerce')
            chunk['P'] = 10 ** (-chunk['LOG10P'])
            if test_filter is not None:
                chunk = chunk[chunk['TEST'].isin(test_filter)]
            chunk = chunk[chunk['P'] <= args.max_p]
            if chunk.empty:
                continue
            chunk[['setid','mask','aaf']] = pd.DataFrame(chunk['ID'].apply(parse_id).tolist(), index=chunk.index)
            chunk['gene'] = chunk['setid'].map(setid_to_gene)
            chunk['trait'] = args.trait
            chunk['source_file'] = path.name
            out_rows.append(chunk)

    if out_rows:
        out = pd.concat(out_rows, ignore_index=True)
    else:
        out = pd.DataFrame(columns=['trait','CHROM','ID','TEST','BETA','LOG10P','P','setid','mask','aaf','gene','source_file'])

    out['pass_2.5e-6'] = out['P'] <= 2.5e-6
    out['pass_7.9e-7'] = out['P'] <= 7.9e-7
    out['pass_1e-5'] = out['P'] <= 1e-5

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep='\t', index=False)

    out_2p5 = args.out_2p5e_6
    out_7p9 = args.out_7p9e_7
    if args.out_all_thresholds:
        base = Path(args.out)
        if base.suffix:
            stem = base.with_suffix('')
            out_2p5 = str(stem) + '_p2.5e-6' + base.suffix
            out_7p9 = str(stem) + '_p7.9e-7' + base.suffix
        else:
            out_2p5 = str(base) + '_p2.5e-6.tsv'
            out_7p9 = str(base) + '_p7.9e-7.tsv'

    if out_2p5:
        out2 = out[out['P'] <= 2.5e-6].copy()
        Path(out_2p5).parent.mkdir(parents=True, exist_ok=True)
        out2.to_csv(out_2p5, sep='\\t', index=False)

    if out_7p9:
        out7 = out[out['P'] <= 7.9e-7].copy()
        Path(out_7p9).parent.mkdir(parents=True, exist_ok=True)
        out7.to_csv(out_7p9, sep='\\t', index=False)


if __name__ == '__main__':
    main()
