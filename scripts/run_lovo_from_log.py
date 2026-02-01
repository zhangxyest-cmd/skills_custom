#!/usr/bin/env python3
import argparse
import json
import os
import re
import subprocess
from pathlib import Path
import pandas as pd

TEST_TO_VC = {
    'ADD': None,
    'ADD-SKAT': 'skat',
    'ADD-SKATO': 'skato',
    'ADD-SKATO-ACAT': 'skato-acat',
    'ADD-ACATV': 'acatv',
    'ADD-ACATO': 'acato',
}


def parse_log(path: Path):
    opts = {}
    flags = set()
    keep_list = []
    with open(path, 'r') as f:
        for line in f:
            s = line.strip()
            if not s.startswith('--'):
                continue
            s = s.rstrip('\\').strip()
            parts = s.split(None, 1)
            key = parts[0][2:]
            if len(parts) == 1:
                flags.add(key)
            else:
                val = parts[1].strip()
                if key == 'keep':
                    keep_list.append(val)
                else:
                    opts[key] = val
    opts['flags'] = flags
    opts['keep'] = keep_list
    return opts


def make_chr_template(path_str: str) -> str:
    return re.sub(r'chr(\d+)', 'chr{chr}', path_str)


def resolve_path(path_str: str, roots):
    if path_str is None:
        return None
    if os.path.isabs(path_str):
        return path_str
    for r in roots:
        cand = Path(r) / path_str
        if cand.exists():
            return str(cand)
    return path_str


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--log', required=True)
    p.add_argument('--jobs', required=True)
    p.add_argument('--lovo-dir', required=True)
    p.add_argument('--regenie', required=True)
    p.add_argument('--search-root', action='append', default=[])
    p.add_argument('--dry-run', action='store_true')
    args = p.parse_args()

    log_path = Path(args.log)
    opts = parse_log(log_path)

    roots = [log_path.parent, Path('/opt/notebooks'), Path('/mnt/project/WGS_XY'), Path('/mnt/project')]
    roots += [Path(r) for r in args.search_root]

    keep = None
    if opts.get('keep'):
        keep = opts['keep'][0]

    bgen_template = make_chr_template(opts.get('bgen',''))
    sample_template = make_chr_template(opts.get('sample',''))
    anno_template = make_chr_template(opts.get('anno-file',''))
    setlist_template = make_chr_template(opts.get('set-list',''))
    cond_template = make_chr_template(opts.get('condition-list','')) if 'condition-list' in opts else None

    jobs = pd.read_csv(args.jobs, sep='\t')
    lovo_dir = Path(args.lovo_dir)
    lovo_dir.mkdir(parents=True, exist_ok=True)

    for _, r in jobs.iterrows():
        trait = r['trait']
        test = r['TEST']
        setid = r['setid']
        mask = r['mask']
        aaf = str(r['aaf'])
        chr_ = str(int(r['CHR'])) if not pd.isna(r['CHR']) else None

        vc = TEST_TO_VC.get(test)
        if vc is None and test != 'ADD':
            print(f"[skip] unsupported test {test} for {setid}")
            continue

        if not chr_:
            print(f"[skip] missing chr for {setid}")
            continue

        outbase = lovo_dir / r['outprefix']
        outcheck = Path(str(outbase) + f"_{trait}.regenie.gz")
        if outcheck.exists():
            continue

        aaf_bins = aaf
        vc_max = aaf
        mask_lovo_aaf = aaf
        if aaf == 'singleton':
            aaf_bins = '0.001'
            vc_max = '0.001'
            mask_lovo_aaf = 'singleton'

        cmd = [
            args.regenie,
            '--step', '2',
            '--pred', resolve_path(opts.get('pred'), roots),
            '--bgen', resolve_path(bgen_template.format(chr=chr_), roots),
            '--sample', resolve_path(sample_template.format(chr=chr_), roots),
            '--ref-first',
        ]
        if keep:
            cmd += ['--keep', resolve_path(keep, roots)]

        cmd += [
            '--phenoFile', resolve_path(opts.get('phenoFile'), roots),
            f"--phenoCol={trait}",
        ]

        if 'bt' in opts['flags']:
            cmd += ['--bt']
        if 't2e' in opts['flags']:
            cmd += ['--t2e']
        if 'firth' in opts['flags']:
            cmd += ['--firth']
        if 'firth-se' in opts['flags']:
            cmd += ['--firth-se']
        if 'approx' in opts['flags']:
            cmd += ['--approx']

        if vc and test != 'ADD':
            cmd += ['--vc-tests', vc]

        if 'minMAC' in opts:
            cmd += ['--minMAC', str(opts['minMAC'])]

        if 'covarFile' in opts:
            cmd += ['--covarFile', resolve_path(opts['covarFile'], roots)]
        if 'covarCol' in opts:
            cmd += [f"--covarCol={opts['covarCol']}"]
        if 'catCovarList' in opts:
            cmd += [f"--catCovarList={opts['catCovarList']}"]
        if 'maxCatLevels' in opts:
            cmd += ['--maxCatLevels', str(opts['maxCatLevels'])]

        cmd += [
            '--anno-file', resolve_path(anno_template.format(chr=chr_), roots),
            '--mask-def', resolve_path(opts.get('mask-def'), roots),
            '--set-list', resolve_path(setlist_template.format(chr=chr_), roots),
            '--aaf-bins', str(aaf_bins),
            '--vc-maxAAF', str(vc_max),
            '--extract-setlist', str(setid),
            '--mask-lovo', f"{setid},{mask},{mask_lovo_aaf}",
        ]

        if cond_template:
            cond_path = resolve_path(cond_template.format(chr=chr_), roots)
            if Path(cond_path).exists():
                cmd += ['--condition-list', cond_path]

        if 'bsize' in opts:
            cmd += ['--bsize', str(opts['bsize'])]

        cmd += ['--out', str(outbase), '--gz']

        if args.dry_run:
            print(' '.join(cmd))
            continue

        print(f"[run] {trait} {test} {setid} {mask} {aaf} chr{chr_}")
        subprocess.run(cmd, check=True)


if __name__ == '__main__':
    main()
