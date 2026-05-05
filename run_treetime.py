"""
Run TreeTime on Gubbins divergence trees for temporal signal check and outlier detection.

This step sits between Gubbins and BEAST in the pipeline:
  1. Gubbins produces a divergence tree + filtered alignment
  2. TreeTime (this script) checks clock signal and flags outlier tips
  3. Filtered alignment (outliers removed) goes into BEAUti/BEAST

Usage:
    python run_treetime.py --st st131
    python run_treetime.py  # runs all STs found in gubbins/
"""

import argparse
import glob
import os
import sys

GUBBINS_DIR = './gubbins'
DATES_DIR = './tables'
TREETIME_OUTDIR = './treetime'


def find_gubbins_tree(st_name: str) -> str | None:
    candidates = [
        os.path.join(GUBBINS_DIR, f'{st_name}_norm.nwk'),
        os.path.join(GUBBINS_DIR, f'{st_name}_down100.node_labelled.final_tree.tre'),
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    return None


def find_dates_file(st_name: str) -> str | None:
    path = os.path.join(DATES_DIR, f'{st_name}_dates.txt')
    return path if os.path.exists(path) else None


def run_treetime(st_name: str, clock_rate: float | None = None) -> None:
    tree = find_gubbins_tree(st_name)
    if tree is None:
        print(f'[{st_name}] No Gubbins tree found in {GUBBINS_DIR}/ — skipping.')
        return

    dates = find_dates_file(st_name)
    if dates is None:
        print(f'[{st_name}] No dates file found at {DATES_DIR}/{st_name}_dates.txt — skipping.')
        return

    outdir = os.path.join(TREETIME_OUTDIR, st_name)
    os.makedirs(outdir, exist_ok=True)

    clock_flag = f'--clock-rate {clock_rate}' if clock_rate else ''

    # root-to-tip regression + outlier detection
    # --name-column / --date-column are explicit so TreeTime doesn't try to
    # auto-detect the separator (auto-detection can misread TSV as single-column)
    cmd = (
        f'treetime '
        f'--tree "{tree}" '
        f'--dates "{dates}" '
        f'--name-column name '
        f'--date-column date '
        f'--outdir "{outdir}" '
        f'--reroot best '
        f'--aln {GUBBINS_DIR}/{st_name}_down100.filtered_polymorphic_sites.fasta'
        f'{clock_flag}'
    )
    print(f'[{st_name}] Running: {cmd}')
    ret = os.system(cmd)

    if ret != 0:
        print(f'[{st_name}] TreeTime exited with code {ret}. Check {outdir} for partial output.')
        return

    outliers_path = os.path.join(outdir, 'outliers.tsv')
    if os.path.exists(outliers_path):
        with open(outliers_path) as fh:
            lines = [l for l in fh if not l.startswith('#') and l.strip()]
        # header line + data lines
        outlier_samples = [l.split('\t')[0] for l in lines[1:]] if len(lines) > 1 else []
        print(f'[{st_name}] {len(outlier_samples)} outlier(s) detected — see {outliers_path}')
        if outlier_samples:
            print('  Outliers:', ', '.join(outlier_samples[:10]),
                  '...' if len(outlier_samples) > 10 else '')
    else:
        print(f'[{st_name}] No outliers.tsv produced — inspect {outdir} manually.')

    print(f'[{st_name}] Done. Root-to-tip plot: {outdir}/root_to_tip.pdf')


def main():
    parser = argparse.ArgumentParser(description='TreeTime outlier detection for each ST')
    parser.add_argument('--st', help='Specific ST to process (e.g. st131). '
                                     'If omitted, all STs in gubbins/ are processed.')
    parser.add_argument('--clock-rate', type=float, default=None,
                        help='Fixed substitution rate (per site per year). '
                             'Omit to let TreeTime infer it from the data.')
    args = parser.parse_args()

    if args.st:
        targets = [args.st]
    else:
        targets = [
            os.path.basename(d)
            for d in glob.glob(os.path.join(GUBBINS_DIR, '*'))
            if os.path.isdir(d)
        ]
        if not targets:
            sys.exit(f'No subdirectories found in {GUBBINS_DIR}. '
                     'Run align_with_gubbins.py first.')

    for st in sorted(targets):
        run_treetime(st, clock_rate=args.clock_rate)

    print(f'\nAll done. Outputs in {TREETIME_OUTDIR}/')
    print('Next step: review root_to_tip.pdf and outliers.tsv, '
          'then pass cleaned alignment + dates to prepare_beauti_inputs.py')


if __name__ == '__main__':
    main()