"""
Split the ST131 Gubbins/TreeTime tree into sub-clades whose most recent common
ancestor (MRCA) falls within a given number of years of the most recent tip.

ST131 is known to split very early (two deeply-diverged clades), making a single
BEAST run covering the full tree impractical. This script:
  1. Loads the time-scaled tree produced by TreeTime (or the Gubbins divergence
     tree if a TreeTime tree is not yet available).
  2. Identifies the earliest split(s) that divide the tree into sub-clades each
     spanning <= MAX_CLADE_YEARS years from MRCA to present.
  3. Writes each sub-clade as its own Newick file + a matching isolate TSV
     (isolate name + path to .fa.gz), so align_with_gubbins.py can be
     re-run per clade.

Usage:
    python split_st131_tree.py
    python split_st131_tree.py --max-years 50 --tree-dir treetime/st131
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


TREETIME_DIR   = './treetime/st131'
GUBBINS_DIR    = './gubbins/st131'
OUT_DIR        = './st131_clades'
ASSEMBLIES_DIR = './assemblies'
MAX_YEARS      = 50   # split so each clade MRCA is within this many years of tips


def load_tree(tree_path: str) -> dendropy.Tree:
    schema = 'nexus' if tree_path.endswith('.nexus') else 'newick'
    return dendropy.Tree.get(path=tree_path, schema=schema,
                             preserve_underscores=True)


def _regex_dates_from_nexus(nexus_path: str,
                            known_taxa: set[str]) -> dict[str, float]:
    """
    Parse tip dates directly from a TreeTime nexus file with regex,
    restricted to labels that are known leaf taxa.

    TreeTime annotates every node (tips AND internals) with [&...,date=...].
    Internal node dates are relative to the root (not calendar years), so we
    must filter by the taxon set from the tree to avoid matching them.
    """
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


def _annotation_date(leaf) -> float | None:
    """Try to read the 'date' annotation from a dendropy leaf node."""
    for source in [getattr(leaf, 'annotations', []),
                   getattr(getattr(leaf, 'edge', None), 'annotations', [])]:
        for ann in (source or []):
            if getattr(ann, 'name', None) == 'date':
                try:
                    return float(ann.value)
                except (TypeError, ValueError):
                    pass
    return None


def _iter_leaves(node):
    """Yield all leaf nodes in a subtree (works on any dendropy Node)."""
    if node.is_leaf():
        yield node
    else:
        for child in node.child_node_iter():
            yield from _iter_leaves(child)


def tip_years(tree: dendropy.Tree, tree_path: str = '') -> dict[str, float]:
    """
    Extract decimal years for every tip, trying three strategies:

    1. Regex on the raw nexus file  (primary for TreeTime output — most complete)
    2. dendropy node/edge annotations (supplement for any gaps)
    3. '|YYYY-MM-DD' in the tip label (prepare_beauti_inputs.py format)
    """
    years = {}

    # strategy 1: regex on raw nexus — dendropy may only surface a subset of
    # [&...] annotations, so always parse the file directly for nexus trees.
    # Pass known_taxa so internal node annotations are excluded.
    if tree_path.endswith('.nexus'):
        known_taxa = {lf.taxon.label
                      for lf in tree.leaf_node_iter() if lf.taxon}
        years = _regex_dates_from_nexus(tree_path, known_taxa)
        print(f'  Regex extracted dates for {len(years)} / {len(known_taxa)} '
              f'tips from nexus file.')

    # strategy 2: dendropy annotations — fill gaps not caught by regex
    for leaf in tree.leaf_node_iter():
        label = leaf.taxon.label if leaf.taxon else ''
        if label not in years:
            val = _annotation_date(leaf)
            if val is not None:
                years[label] = val

    # strategy 3: '|date' embedded in tip labels
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


def node_mrca_year(node, tip_year_map: dict) -> float | None:
    """
    For a subtree rooted at *node*, return the estimated year of the MRCA
    using the minimum tip year in that clade as a rough lower bound.
    Works on divergence trees (uses minimum tip date as MRCA proxy) and on
    time-scaled trees (reads node age from distance-to-root if available).
    """
    leaf_years = [
        tip_year_map[lf.taxon.label]
        for lf in _iter_leaves(node)
        if lf.taxon and lf.taxon.label in tip_year_map
    ]
    return min(leaf_years) if leaf_years else None


def clade_span_years(node, tip_year_map: dict) -> float:
    """Max tip year minus estimated MRCA year for this clade."""
    leaf_years = [
        tip_year_map[lf.taxon.label]
        for lf in _iter_leaves(node)
        if lf.taxon and lf.taxon.label in tip_year_map
    ]
    if not leaf_years:
        return 0.0
    return max(leaf_years) - min(leaf_years)


def find_split_nodes(tree: dendropy.Tree, tip_year_map: dict,
                     max_years: float) -> list:
    """
    Traverse from the root. Return the set of nodes where we 'cut' the tree
    so that every resulting sub-clade spans <= max_years.

    If a node's clade spans > max_years we recurse into its children.
    If a node's clade spans <= max_years we keep it as a single clade.
    """
    def _recurse(node):
        span = clade_span_years(node, tip_year_map)
        if span <= max_years or node.is_leaf():
            return [node]
        result = []
        for child in node.child_nodes():
            result.extend(_recurse(child))
        return result

    return _recurse(tree.seed_node)


def find_assembly(isolate: str, assemblies_dir: str) -> str | None:
    """
    Search for <isolate>.fa.gz under any immediate subdirectory of
    assemblies_dir.  Returns the (relative) path if found, else None.
    """
    pattern = os.path.join(assemblies_dir, '*', f'{isolate}.fa.gz')
    hits = glob.glob(pattern)
    if hits:
        return os.path.normpath(hits[0])
    return None


def write_clade(clade_node, tree: dendropy.Tree, tip_year_map: dict,
                idx: int, out_dir: str, assemblies_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)

    # collect tip names
    tips = [
        lf.taxon.label
        for lf in _iter_leaves(clade_node)
        if lf.taxon
    ]
    n_tips = len(tips)
    years = [tip_year_map.get(t) for t in tips if tip_year_map.get(t)]
    span  = f'{min(years):.1f}–{max(years):.1f}' if years else 'unknown'

    print(f'  Clade {idx}: {n_tips} tips, date range {span}')

    # Build TSV rows: isolate_name <TAB> path/to/isolate.fa.gz
    # The isolate name is the part before the first '|' (if dates are encoded).
    rows = []
    missing = []
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

    # write isolate TSV
    list_path = os.path.join(out_dir, f'clade_{idx}_isolates.tsv')
    with open(list_path, 'w') as fh:
        for isolate, fa_path in rows:
            fh.write(f'{isolate}\t./{fa_path}\n')

    # extract and write sub-tree
    sub_tree = tree.extract_tree_with_taxa_labels(
        labels=set(tips), suppress_unifurcations=True
    )
    nwk_path = os.path.join(out_dir, f'clade_{idx}.nwk')
    sub_tree.write(path=nwk_path, schema='newick')

    print(f'    Wrote {list_path}')
    print(f'    Wrote {nwk_path}')
    print(f'    Next: re-run align_with_gubbins.py using {list_path} as isolate input')


def main():
    parser = argparse.ArgumentParser(
        description='Split ST131 tree into sub-clades with MRCA <= MAX_YEARS apart'
    )
    parser.add_argument('--tree', default=None,
                        help='Path to Newick tree. Defaults to TreeTime output, '
                             'then Gubbins tree.')
    parser.add_argument('--max-years', type=float, default=MAX_YEARS,
                        help=f'Maximum clade span in years (default: {MAX_YEARS})')
    parser.add_argument('--out-dir', default=OUT_DIR,
                        help=f'Output directory (default: {OUT_DIR})')
    parser.add_argument('--assemblies-dir', default=ASSEMBLIES_DIR,
                        help=f'Root directory containing per-batch assembly '
                             f'subdirectories (default: {ASSEMBLIES_DIR})')
    args = parser.parse_args()

    # find tree
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
            'No tree found. Run run_treetime.py first, or pass --tree explicitly.\n'
            f'Searched: {[os.path.join(TREETIME_DIR, "timetree.nexus"), os.path.join(GUBBINS_DIR, "st131.final_tree.tre")]}'
        )

    schema = 'nexus' if tree_path.endswith('.nexus') else 'newick'
    print(f'Loading tree: {tree_path} (schema={schema})')
    tree = load_tree(tree_path)

    tip_year_map = tip_years(tree, tree_path=tree_path)
    if not tip_year_map:
        sys.exit(
            'Could not extract dates from this tree.\n'
            'Tried: dendropy annotations, nexus regex, and |date tip labels.\n'
            'Make sure you are passing a TreeTime timetree.nexus or a tree\n'
            'whose tip labels contain |YYYY-MM-DD dates.'
        )

    all_years = list(tip_year_map.values())
    print(f'Tips with dates: {len(tip_year_map)}')
    print(f'Date range: {min(all_years):.1f} – {max(all_years):.1f}')
    print(f'Splitting into clades with span <= {args.max_years} years...')
    print(f'Searching for assemblies under: {args.assemblies_dir}')

    clades = find_split_nodes(tree, tip_year_map, max_years=args.max_years)
    print(f'Found {len(clades)} clade(s)')

    for idx, node in enumerate(clades, start=1):
        write_clade(node, tree, tip_year_map, idx, args.out_dir,
                    assemblies_dir=args.assemblies_dir)

    print(f'\nDone. Sub-clade files written to {args.out_dir}/')
    print('For each clade, re-run align_with_gubbins.py with the clade isolate TSV,')
    print('then re-run run_treetime.py on the per-clade tree before proceeding to BEAST.')


if __name__ == '__main__':
    main()