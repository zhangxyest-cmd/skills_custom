#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

SUPPORTED_TESTS = {
    'ADD',
    'ADD-SKAT',
    'ADD-SKATO',
    'ADD-SKATO-ACAT',
    'ADD-ACATV',
    'ADD-ACATO',
}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--sig', required=True, help='sig_masks_DR_FU.tsv')
    p.add_argument('--out', required=True, help='lovo_jobs.tsv')
    p.add_argument('--include-tests', help='comma-separated TEST list to keep (optional)')
    args = p.parse_args()

    df = pd.read_csv(args.sig, sep='\t')
    if df.empty:
        raise SystemExit("No significant masks found")

    if args.include_tests:
        allow = set([t.strip() for t in args.include_tests.split(',') if t.strip()])
        df = df[df['TEST'].isin(allow)]

    # build jobs for supported tests only
    df['lovo_supported'] = df['TEST'].isin(SUPPORTED_TESTS)

    # use best (min P) per unique mask+test
    grp = df.groupby(['trait','TEST','setid','mask','aaf'], as_index=False)
    df_min = grp.apply(lambda x: x.loc[x['P'].idxmin()]).reset_index(drop=True)

    jobs = df_min[df_min['lovo_supported']].copy()
    jobs['CHR'] = jobs['CHROM'] if 'CHROM' in jobs.columns else jobs.get('CHR')

    def outprefix_row(r):
        aaf = str(r['aaf']).replace('/', '_')
        return f"{r['trait']}_lovo_{r['setid']}_{r['mask']}_{aaf}_{r['TEST']}"

    jobs['outprefix'] = jobs.apply(outprefix_row, axis=1)

    cols = ['trait','setid','mask','aaf','TEST','P','CHR','outprefix']
    jobs = jobs[cols]

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    jobs.to_csv(args.out, sep='\t', index=False)


if __name__ == '__main__':
    main()
