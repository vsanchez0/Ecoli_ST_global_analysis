"""
Split the FastTree/Gubbins/TreeTime tree into sub-clades using a
pairwise tip-distance approach with a sparse connectivity matrix.

  1. Compute all pairwise patristic distances between tips.
  2. Build a sparse adjacency matrix: pair (i,j) = 1 if distance <= threshold,
     0 otherwise.
  3. Find connected components of that graph — each component is one clade.

The threshold is:  mutation_rate (subs/site/year) × max_years

Usage:
    python split_st131_tree.py --tree fasttree/st131/st131.nwk
    python split_st131_tree.py --tree fasttree/st131/st131.nwk --max-years 50
    python split_st131_tree.py --tree fasttree/st131/st131.nwk --cutoff 2.5e-5
"""

import argparse
import glob
import os
import re
import sys

try:
    import dendropy
except ImportError:
    sys.exit('dendropy is required: conda install -c conda-forge dendropy')

try:
    import numpy as np
except ImportError:
    sys.exit('numpy is required: pip install numpy')

try:
    from scipy.sparse import lil_matrix
    from scipy.sparse.csgraph import connected_components
except ImportError:
    sys.exit('scipy is required: pip install scipy')


TREETIME_DIR   = './treetime/st131'
GUBBINS_DIR    = './gubbins/st131'
OUT_DIR        = './st131_clades'
ASSEMBLIES_DIR = './assemblies'
MAX_YEARS      = 50

ECOLI_MUTATION_RATE = 7e-7


def load_tree(tree_path: str) -> dendropy.Tree:
    schema = 'nexus' if tree_path.endswith('.nexus') else 'newick'
    return dendropy.Tree.get(path=tree_path, schema=schema,
                             preserve_underscores=True)

def _regex_dates_from_nexus(nexus_path: str, known_taxa: set) -> dict:
    tip_dates = {}
    pattern = re.compile(r'([\w.]+):\s*[\d.]+\s*\[&[^\]]*\bdate=([\d.]+)')
    try:
        with open(nexus_path) as fh:
            content = fh.read()
        for m in pattern.finditer(content):
            label = m.group(1)
            if label in known_taxa:
                tip_dates[label] = float(m.group(2))
    except Exception as e:
        print(f'  Warning: regex date extraction failed on {nexus_path}: {e}')
    return tip_dates


def _annotation_date(leaf):
    for source in [getattr(leaf, 'annotations', []),
                   getattr(getattr(leaf, 'edge', None), 'annotations', [])]:
        for ann in (source or []):
            if getattr(ann, 'name', None) == 'date':
                try:
                    return float(ann.value)
                except (TypeError, ValueError):
                    pass
    return None


def load_dates_file(dates_path: str) -> dict:
    import csv as _csv
    result = {}
    with open(dates_path) as fh:
        sample = fh.read(1024)
        fh.seek(0)
        dialect = _csv.Sniffer().sniff(sample, delimiters=',\t')
        reader = _csv.DictReader(fh, dialect=dialect)
        for row in reader:
            name = row.get('name', '').strip()
            date = row.get('date', '').strip()
            if not name or not date:
                continue
            parts = date.split('-')
            try:
                yr = int(parts[0])
                mo = int(parts[1]) if len(parts) > 1 else 6
                da = int(parts[2]) if len(parts) > 2 else 15
                result[name] = yr + (mo - 1) / 12.0 + (da - 1) / 365.0
            except (ValueError, IndexError):
                pass
    return result


def tip_years(tree: dendropy.Tree, tree_path: str = '',
              dates_file: str = '') -> dict:
    years = {}
    if dates_file and os.path.exists(dates_file):
        years = load_dates_file(dates_file)
        print(f'  Loaded dates for {len(years)} tips from {dates_file}')
        return years
    if tree_path.endswith('.nexus'):
        known_taxa = {lf.taxon.label for lf in tree.leaf_node_iter() if lf.taxon}
        years = _regex_dates_from_nexus(tree_path, known_taxa)
        print(f'  Regex extracted dates for {len(years)} / {len(known_taxa)} tips.')
    for leaf in tree.leaf_node_iter():
        label = leaf.taxon.label if leaf.taxon else ''
        if label not in years:
            val = _annotation_date(leaf)
            if val is not None:
                years[label] = val
    if not years:
        for leaf in tree.leaf_node_iter():
            label = leaf.taxon.label if leaf.taxon else ''
            if '|' in label:
                date_str = label.split('|')[-1]
                parts = date_str.split('-')
                try:
                    yr = int(parts[0])
                    mo = int(parts[1]) if len(parts) > 1 else 6
                    years[label] = yr + (mo - 1) / 12.0
                except (ValueError, IndexError):
                    pass
    return years


def cluster_by_pairwise_distance(
        tree: dendropy.Tree,
        threshold: float,
        min_clade_size: int = 10,
) -> list[set]:
    """
    Cluster tips by pairwise patristic distance using a sparse adjacency graph.

    For every pair of tips (i, j):
        - if patristic_distance(i, j) <= threshold  →  connect them (edge = 1)
        - otherwise                                  →  no edge (0)

    Parameters
    ----------
    tree            : dendropy.Tree (rooted or unrooted — both fine)
    threshold       : max patristic distance for two tips to be in the same clade
    min_clade_size  : clusters smaller than this are reported but NOT merged;
                      they remain as separate (likely outlier) clusters

    Returns
    -------
    list[set[str]]  : each set is a clade's tip labels, sorted largest first
    """
    print(f'\n  Computing pairwise patristic distances for all tips...')
    print(f'  Threshold: {threshold:.3e} subs/site')

    pdm = tree.phylogenetic_distance_matrix()
    taxa = [t for t in pdm.taxon_namespace if t is not None]
    n    = len(taxa)
    print(f'  Tips: {n}  →  up to {n*(n-1)//2:,} pairs to evaluate')

    adj = lil_matrix((n, n), dtype=np.int8)

    connected_pairs = 0
    for i in range(n):
        for j in range(i + 1, n):
            d = pdm(taxa[i], taxa[j])
            if d <= threshold:
                adj[i, j] = 1
                adj[j, i] = 1
                connected_pairs += 1

    print(f'  Connected pairs (distance <= threshold): {connected_pairs:,}')

    n_components, labels = connected_components(
        adj.tocsr(), directed=False, return_labels=True
    )
    print(f'  Connected components found: {n_components}')

    clusters: dict[int, set] = {}
    for i, taxon in enumerate(taxa):
        comp = int(labels[i])
        clusters.setdefault(comp, set()).add(taxon.label)

    clade_list = sorted(clusters.values(), key=len, reverse=True)

    sizes = sorted((len(c) for c in clade_list), reverse=True)
    small = sum(1 for s in sizes if s < min_clade_size)
    print(f'\n  Clade size summary:')
    print(f'    Largest clade:  {sizes[0]} tips')
    print(f'    Median size:    {int(np.median(sizes))} tips')
    print(f'    Smallest clade: {sizes[-1]} tips')
    print(f'    Clades >= {min_clade_size} tips: {len(sizes) - small}')
    print(f'    Clades <  {min_clade_size} tips: {small}  '
          f'(kept as individual clades — do not pass to TreeTime)')

    return clade_list



def _iter_leaves(node):
    if node.is_leaf():
        yield node
    else:
        for child in node.child_node_iter():
            yield from _iter_leaves(child)


def find_assembly(isolate: str, assemblies_dir: str) -> str | None:
    pattern = os.path.join(assemblies_dir, '*', f'{isolate}.fa.gz')
    hits = glob.glob(pattern)
    if hits:
        return os.path.normpath(hits[0])
    return None


def write_clade(tips_set: set, tree: dendropy.Tree, tip_year_map: dict,
                idx: int, out_dir: str, assemblies_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)

    tips = sorted(tips_set)
    n_tips = len(tips)
    years = [tip_year_map[t] for t in tips if t in tip_year_map]
    span  = f'{min(years):.1f}–{max(years):.1f}' if years else 'unknown'

    print(f'\n  Clade {idx}: {n_tips} tips, date range {span}')

    # Isolate TSV
    rows, missing = [], []
    for tip in tips:
        isolate = tip.split('|')[0] if '|' in tip else tip
        fa_path = find_assembly(isolate, assemblies_dir)
        if fa_path is None:
            missing.append(isolate)
            fa_path = 'NOT_FOUND'
        rows.append((isolate, fa_path))

    if missing:
        print(f'    WARNING: .fa.gz not found for {len(missing)} isolate(s): '
              f'{missing[:5]}{"..." if len(missing) > 5 else ""}')

    list_path = os.path.join(out_dir, f'clade_{idx}_isolates.tsv')
    with open(list_path, 'w') as fh:
        for isolate, fa_path in rows:
            fh.write(f'{isolate}\t./{fa_path}\n')

    # Sub-tree newick
    sub_tree = tree.extract_tree_with_taxa_labels(
        labels=set(tips), suppress_unifurcations=True
    )
    nwk_path = os.path.join(out_dir, f'clade_{idx}.nwk')
    sub_tree.write(path=nwk_path, schema='newick')

    print(f'    Wrote {list_path}')
    print(f'    Wrote {nwk_path}')


def main():
    parser = argparse.ArgumentParser(
        description=(
            'Split an ST131 divergence tree into clades using pairwise tip '
            'distances and a temporal threshold.  Works correctly on unrooted '
            'FastTree output.'
        )
    )
    parser.add_argument('--tree', default=None,
                        help='Path to Newick or Nexus tree file. '
                             'Defaults to TreeTime output, then Gubbins tree.')
    parser.add_argument('--dates', default=None,
                        help='Path to dates CSV/TSV (columns: name, date). '
                             'Used only for reporting date ranges per clade — '
                             'NOT for splitting. Optional.')
    parser.add_argument('--max-years', type=float, default=MAX_YEARS,
                        help=f'Maximum pairwise distance in years for two tips '
                             f'to be placed in the same clade. '
                             f'Converted to branch-length threshold as: '
                             f'mutation-rate × max-years  (default: {MAX_YEARS})')
    parser.add_argument('--mutation-rate', type=float, default=ECOLI_MUTATION_RATE,
                        help=f'Substitution rate in subs/site/year used to convert '
                             f'--max-years into a branch-length threshold '
                             f'(default: {ECOLI_MUTATION_RATE:.1e}). '
                             f'IMPORTANT: if your SKA alignment covers only variant '
                             f'sites (e.g. 50k sites from a 5Mb genome), inflate this '
                             f'by genome_size / alignment_length, e.g. 5e-7 × 100 = 5e-5.')
    parser.add_argument('--cutoff', type=float, default=None,
                        help='Override the threshold directly in subs/site, '
                             'bypassing --mutation-rate and --max-years.')
    parser.add_argument('--min-clade-size', type=int, default=10,
                        help='Clades smaller than this are written out but flagged '
                             'as too small for TreeTime (default: 10).')
    parser.add_argument('--out-dir', default=OUT_DIR,
                        help=f'Output directory (default: {OUT_DIR})')
    parser.add_argument('--assemblies-dir', default=ASSEMBLIES_DIR,
                        help=f'Root of assembly directories (default: {ASSEMBLIES_DIR})')
    args = parser.parse_args()

    tree_path = args.tree
    if tree_path is None:
        candidates = [
            os.path.join(TREETIME_DIR, 'timetree.nexus'),
            os.path.join(TREETIME_DIR, 'timetree.nwk'),
            os.path.join(GUBBINS_DIR,  'st131.final_tree.tre'),
            os.path.join(GUBBINS_DIR,  'st131.node_labelled.final_tree.tre'),
        ]
        for c in candidates:
            if os.path.exists(c):
                tree_path = c
                break

    if tree_path is None or not os.path.exists(tree_path):
        sys.exit(
            'No tree found. Pass --tree explicitly, or run run_treetime.py first.'
        )

    if args.cutoff is not None:
        threshold = args.cutoff
        print(f'Threshold (user-supplied): {threshold:.3e} subs/site')
    else:
        threshold = args.mutation_rate * args.max_years
        print(f'Threshold = {args.mutation_rate:.2e} subs/site/yr '
              f'× {args.max_years} yr = {threshold:.3e} subs/site')

    schema = 'nexus' if tree_path.endswith('.nexus') else 'newick'
    print(f'Loading tree: {tree_path}  (schema={schema})')
    tree = load_tree(tree_path)
    n_tips = sum(1 for _ in tree.leaf_node_iter())
    print(f'Tips loaded: {n_tips}')

    tip_year_map = tip_years(tree, tree_path=tree_path,
                             dates_file=args.dates or '')
    if tip_year_map:
        all_y = list(tip_year_map.values())
        print(f'Dates loaded for {len(tip_year_map)} tips  '
              f'({min(all_y):.1f}–{max(all_y):.1f})')
    else:
        print('No dates loaded — clade date ranges will show as "unknown". '
              'Pass --dates if you have a dates file.')

    clades = cluster_by_pairwise_distance(
        tree,
        threshold=threshold,
        min_clade_size=args.min_clade_size,
    )

    print(f'\nWriting {len(clades)} clade(s) to {args.out_dir}/')
    for idx, tips_set in enumerate(clades, start=1):
        write_clade(tips_set, tree, tip_year_map, idx,
                    args.out_dir, args.assemblies_dir)

    summary_path = os.path.join(args.out_dir, 'clade_assignments.tsv')
    os.makedirs(args.out_dir, exist_ok=True)
    with open(summary_path, 'w') as fh:
        fh.write('# Pairwise-distance clade assignments\n')
        fh.write(f'# Threshold: {threshold:.3e} subs/site  '
                 f'({args.mutation_rate:.2e} subs/site/yr × {args.max_years} yr)\n')
        fh.write('taxon\tclade\n')
        for idx, tips_set in enumerate(clades, start=1):
            for tip in sorted(tips_set):
                fh.write(f'{tip}\tclade_{idx}\n')
    print(f'\nAssignment summary: {summary_path}')

    print('\nDone.')
    print('For each clade >= min-clade-size, re-run align_with_gubbins.py,')
    print('then run_treetime.py, before proceeding to BEAST.')


if __name__ == '__main__':
    main()