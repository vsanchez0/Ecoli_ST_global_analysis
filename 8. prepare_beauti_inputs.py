"""
Prepare input files for BEAUti / BEAST DTA analysis.

Produces:

  1. <st>_alignment.fasta
       Gubbins filtered alignment with tip names renamed to isolate|YYYY-MM-DD
       (format BEAUti reads tip dates from).

  2. <st>_amr_traits.fasta
       Presence/absence matrix encoded as a binary FASTA alignment.
       Each sequence is a string of 0/1 characters (one per AMR gene).
       Genes that are fixed (100% prevalence) or absent (0% prevalence) are dropped.
       Import into BEAUti as "Standard" data type to get Lewis Mk model.

  3. <st>_locations.txt
       Tab-delimited HHS region per isolate.
       Columns: name, hhs_region
"""

import argparse
import glob
import os
import re
import sys

import pandas as pd

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
except ImportError:
    sys.exit('biopython is required: conda install -c conda-forge biopython')


GUBBINS_DIR  = './gubbins'
TABLES_DIR   = './tables'
RGI_DIR      = './rgi_results'
BEAUTI_DIR   = './beauti_inputs'

RGI_SAMPLE_COL  = 'Best_Hit_ARO'
RGI_CUT_OFF     = 'Cut_Off'      # "Strict" or "Perfect" hits only
RGI_VALID_CUTS  = {'Strict', 'Perfect'}


def decimal_year(date_str: str) -> float | None:
    """Convert YYYY-MM-DD (or YYYY-MM or YYYY) to decimal year."""
    if not isinstance(date_str, str):
        return None
    parts = date_str.strip().split('-')
    try:
        yr = int(parts[0])
        mo = int(parts[1]) if len(parts) > 1 else 6
        return yr + (mo - 1) / 12.0
    except (ValueError, IndexError):
        return None


def beauti_date_label(sample_id: str, date_str: str) -> str:
    """Return tip name in BEAUti format: sampleID|YYYY-MM-DD."""
    if not isinstance(date_str, str):
        return sample_id
    clean = date_str.strip()[:10]  # keep YYYY-MM-DD at most
    return f'{sample_id}|{clean}'


def load_metadata(st_name: str) -> pd.DataFrame:
    path = os.path.join(TABLES_DIR, f'{st_name}_metadata_with_hhs.csv')
    if not os.path.exists(path):
        path = os.path.join(TABLES_DIR, f'{st_name}_metadata.csv')
    if not os.path.exists(path):
        return pd.DataFrame()
    return pd.read_csv(path)


def load_treetime_dates(path: str) -> dict[str, str]:
    """
    Load dates from a TreeTime dates TSV.
    Expected header: #node  date  numeric_date
    Returns a dict mapping node name -> date string (YYYY-MM-DD or YYYY-MM or YYYY).
    Rows with '--' as the date are included as-is; callers are responsible for
    filtering them out via the excluded_sids mechanism.
    """
    date_map = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            node, date_str = parts[0], parts[1]
            if node and date_str:
                date_map[node] = date_str
    return date_map


def load_rgi_for_st(st_name: str, sample_ids: list[str]) -> pd.DataFrame:
    """
    Load RGI output for the given samples and return a binary presence/absence
    DataFrame (samples x genes).

    Expected RGI output layout: one txt per sample at
        {RGI_DIR}/{sample_id}/{sample_id}.txt   or   {RGI_DIR}/{sample_id}.txt
    """
    records = {}
    for sid in sample_ids:
        candidates = [
            os.path.join(RGI_DIR, sid, f'{sid}.txt'),
            os.path.join(RGI_DIR, f'{sid}.txt'),
            os.path.join(RGI_DIR, st_name, sid, f'{sid}.txt'),
            os.path.join(RGI_DIR, st_name, f'{sid}.txt'),
        ]
        rgi_file = next((c for c in candidates if os.path.exists(c)), None)
        if rgi_file is None:
            records[sid] = {}
            continue
        try:
            df = pd.read_csv(rgi_file, sep='\t', low_memory=False)
            if RGI_CUT_OFF in df.columns:
                df = df[df[RGI_CUT_OFF].isin(RGI_VALID_CUTS)]
            if RGI_SAMPLE_COL in df.columns:
                genes = df[RGI_SAMPLE_COL].dropna().unique().tolist()
                records[sid] = {g: 1 for g in genes}
        except Exception as e:
            print(f'  Warning: could not parse RGI file {rgi_file}: {e}')
            records[sid] = {}

    if not records:
        return pd.DataFrame()

    pa = pd.DataFrame(records).T.fillna(0).astype(int)
    return pa


def filter_amr_matrix(pa: pd.DataFrame) -> pd.DataFrame:
    """Drop genes that are present in all or none of the samples."""
    n = len(pa)
    if n == 0:
        return pa
    prevalence = pa.sum(axis=0) / n
    keep = pa.columns[(prevalence > 0) & (prevalence < 1)]
    dropped = len(pa.columns) - len(keep)
    if dropped:
        print(f'  Dropped {dropped} gene(s) at 0% or 100% prevalence; '
              f'{len(keep)} gene(s) remain.')
    return pa[keep]


def write_amr_fasta(pa: pd.DataFrame, out_path: str, gene_map_path: str) -> None:
    """
    Write AMR presence/absence matrix as a binary FASTA alignment.

    Each sequence is a string of '0'/'1' characters, one per AMR gene column.
    A companion gene map TSV is written alongside so the column order is recoverable.

    Import the FASTA into BEAUti as data type "Standard" (Lewis Mk).
    """
    records = []
    for tip_label, row in pa.iterrows():
        seq_str = ''.join(str(v) for v in row.values)
        records.append(SeqRecord(Seq(seq_str), id=str(tip_label), description=''))

    SeqIO.write(records, out_path, 'fasta')

    # Write a companion file mapping column index -> gene name
    gene_map = pd.DataFrame({
        'col_index': range(len(pa.columns)),
        'gene': pa.columns.tolist()
    })
    gene_map.to_csv(gene_map_path, sep='\t', index=False)


def process_st(st_name: str, clade_isolates_path: str | None = None,
               alignment_path: str | None = None,
               treetime_dates_path: str | None = None) -> None:
    print(f'\n=== {st_name} ===')
    out_dir = os.path.join(BEAUTI_DIR, st_name)
    os.makedirs(out_dir, exist_ok=True)

    if alignment_path:
        aln_path = alignment_path
        if not os.path.exists(aln_path):
            print(f'  Alignment file not found: {aln_path} — skipping.')
            return
    else:
        gubbins_subdir = os.path.join(GUBBINS_DIR, st_name)
        aln_candidates = [
            os.path.join(gubbins_subdir, f'{st_name}.filtered_polymorphic_sites.fasta'),
            os.path.join(gubbins_subdir, f'{st_name}.snp_sites.aligned.fa'),
        ]
        aln_path = next((c for c in aln_candidates if os.path.exists(c)), None)
        if aln_path is None:
            print(f'  No Gubbins alignment found in {gubbins_subdir}/ — skipping.')
            print(f'  Pass --alignment to specify the file explicitly.')
            return

    meta = load_metadata(st_name)
    if meta.empty:
        print(f'  No metadata found for {st_name} — skipping.')
        return

    date_map = dict(zip(meta['sample_accession'], meta['collection_date']))
    hhs_map  = dict(zip(meta['sample_accession'],
                        meta.get('hhs_region', pd.Series(dtype=str))))

    if treetime_dates_path:
        tt_dates = load_treetime_dates(treetime_dates_path)
        if tt_dates:
            date_map.update(tt_dates)
            print(f'  Using TreeTime dates for {len(tt_dates)} node(s) '
                  f'from {treetime_dates_path}')
        else:
            print(f'  Warning: --treetime-dates file {treetime_dates_path} '
                  f'was empty or unreadable — falling back to metadata dates.')

    if treetime_dates_path:
        excluded_sids = {sid for sid, d in date_map.items() if d == '--'}
        if excluded_sids:
            print(f'  Excluding {len(excluded_sids)} sample(s) with "--" date from TreeTime.')
            date_map = {sid: d for sid, d in date_map.items() if d != '--'}
    else:
        excluded_sids = set()

    if clade_isolates_path:
        with open(clade_isolates_path) as fh:
            keep_ids = {l.strip().split('\t')[0] for l in fh if l.strip()}
        print(f'  Restricting to {len(keep_ids)} isolates from {clade_isolates_path}')
    else:
        keep_ids = None

    records_out = []
    sample_ids  = []
    for rec in SeqIO.parse(aln_path, 'fasta'):
        sid = rec.id
        if keep_ids is not None and sid not in keep_ids:
            continue
        if sid in excluded_sids:
            continue
        date_str  = date_map.get(sid, '')
        new_label = beauti_date_label(sid, date_str)
        records_out.append(SeqRecord(Seq(str(rec.seq)), id=new_label,
                                     description=''))
        sample_ids.append(sid)

    if not records_out:
        print('  No sequences remaining after filtering — skipping.')
        return

    aln_out = os.path.join(out_dir, f'{st_name}_alignment.fasta')
    SeqIO.write(records_out, aln_out, 'fasta')
    print(f'  Alignment ({len(records_out)} sequences): {aln_out}')

    # --- AMR traits as binary FASTA ---
    pa = load_rgi_for_st(st_name, sample_ids)
    if pa.empty:
        print('  No RGI data found — AMR trait file not written.')
        print(f'  Expected RGI files under {RGI_DIR}/')
    else:
        pa = filter_amr_matrix(pa)
        if not pa.empty:
            pa.index = [beauti_date_label(sid, date_map.get(sid, ''))
                        for sid in pa.index]
            pa.index.name = 'name'

            amr_fasta_out   = os.path.join(out_dir, f'{st_name}_amr_traits.fasta')
            amr_genemap_out = os.path.join(out_dir, f'{st_name}_amr_gene_map.tsv')
            write_amr_fasta(pa, amr_fasta_out, amr_genemap_out)

            print(f'  AMR traits FASTA ({pa.shape[1]} genes × {pa.shape[0]} isolates): {amr_fasta_out}')
            print(f'  AMR gene map (col index → gene name): {amr_genemap_out}')
            print(f'  → Import {amr_fasta_out} into BEAUti as data type "Standard"')

    loc_rows = []
    for sid in sample_ids:
        hhs = hhs_map.get(sid)
        if hhs and not pd.isna(hhs):
            tip_label = beauti_date_label(sid, date_map.get(sid, ''))
            loc_rows.append({'name': tip_label,
                             'hhs_region': f'HHS{int(hhs)}'})

    if loc_rows:
        loc_df  = pd.DataFrame(loc_rows)
        loc_out = os.path.join(out_dir, f'{st_name}_locations.txt')
        loc_df.to_csv(loc_out, sep='\t', index=False)
        print(f'  HHS locations ({len(loc_rows)} isolates): {loc_out}')
    else:
        print('  No HHS region data available — location file not written.')

    print(f'  BEAUti inputs ready in {out_dir}/')
    print('  Next: open BEAUti, load template XML, import alignment and trait files.')


def main():
    parser = argparse.ArgumentParser(description='Prepare BEAUti/BEAST input files')
    parser.add_argument('--st', default=None,
                        help='ST to process (e.g. st131). Omit to process all.')
    parser.add_argument('--clade-isolates', default=None,
                        help='Path to clade isolate list (from split_st131_tree.py). '
                             'Only those isolates will be included.')
    parser.add_argument('--alignment', default=None,
                        help='Path to alignment FASTA. Overrides the default '
                             'gubbins/{st}/{st}.filtered_polymorphic_sites.fasta lookup.')
    parser.add_argument('--treetime-dates', default=None,
                        help='Path to TreeTime dates TSV with header '
                             '"#node\\tdate\\tnumeric_date". Overrides collection_date '
                             'from metadata. Samples with "--" as date are excluded.')
    args = parser.parse_args()

    if args.st:
        process_st(args.st, clade_isolates_path=args.clade_isolates,
                   alignment_path=args.alignment,
                   treetime_dates_path=args.treetime_dates)
    else:
        gubbins_sts = [
            os.path.basename(d)
            for d in glob.glob(os.path.join(GUBBINS_DIR, '*'))
            if os.path.isdir(d)
        ]
        if not gubbins_sts:
            sys.exit(f'No ST directories found in {GUBBINS_DIR}/. '
                     'Run align_with_gubbins.py first.')
        for st in sorted(gubbins_sts):
            process_st(st, treetime_dates_path=args.treetime_dates)

    print(f'\nAll BEAUti inputs written to {BEAUTI_DIR}/')


if __name__ == '__main__':
    main()