#!/usr/bin/env python3
import argparse
import json
import re
from pathlib import Path


def make_chr_template(path_str: str) -> str:
    # replace chr<number> with chr{chr}
    return re.sub(r'chr(\d+)', 'chr{chr}', path_str)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--log', required=True)
    p.add_argument('--out', required=True)
    args = p.parse_args()

    opts = {}
    flags = set()
    keep_list = []

    with open(args.log, 'r') as f:
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

    opts['flags'] = sorted(flags)
    opts['keep'] = keep_list

    for k in ['bgen','sample','anno-file','set-list','condition-list']:
        if k in opts:
            opts[k + '_template'] = make_chr_template(opts[k])

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(opts, f, indent=2)


if __name__ == '__main__':
    main()
