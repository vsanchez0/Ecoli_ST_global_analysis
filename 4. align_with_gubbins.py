import pandas as pd
import glob
import os
import re
import numpy as np
import datetime as _dt

# hhs region mappings for future dta migration analysis
STATE_TO_HHS = {
    # Region 1 – New England
    'CT': 1, 'CONNECTICUT': 1,
    'ME': 1, 'MAINE': 1,
    'MA': 1, 'MASSACHUSETTS': 1,
    'NH': 1, 'NEW HAMPSHIRE': 1,
    'RI': 1, 'RHODE ISLAND': 1,
    'VT': 1, 'VERMONT': 1,
    # Region 2 – NY / NJ
    'NJ': 2, 'NEW JERSEY': 2,
    'NY': 2, 'NEW YORK': 2,
    # Region 3 – Mid-Atlantic
    'DE': 3, 'DELAWARE': 3,
    'MD': 3, 'MARYLAND': 3,
    'PA': 3, 'PENNSYLVANIA': 3,
    'VA': 3, 'VIRGINIA': 3,
    'WV': 3, 'WEST VIRGINIA': 3,
    'DC': 3, 'DISTRICT OF COLUMBIA': 3,
    # Region 4 – Southeast
    'AL': 4, 'ALABAMA': 4,
    'FL': 4, 'FLORIDA': 4,
    'GA': 4, 'GEORGIA': 4,
    'KY': 4, 'KENTUCKY': 4,
    'MS': 4, 'MISSISSIPPI': 4,
    'NC': 4, 'NORTH CAROLINA': 4,
    'SC': 4, 'SOUTH CAROLINA': 4,
    'TN': 4, 'TENNESSEE': 4,
    # Region 5 – Midwest
    'IL': 5, 'ILLINOIS': 5,
    'IN': 5, 'INDIANA': 5,
    'MI': 5, 'MICHIGAN': 5,
    'MN': 5, 'MINNESOTA': 5,
    'OH': 5, 'OHIO': 5,
    'WI': 5, 'WISCONSIN': 5,
    # Region 6 – South-Central
    'AR': 6, 'ARKANSAS': 6,
    'LA': 6, 'LOUISIANA': 6,
    'NM': 6, 'NEW MEXICO': 6,
    'OK': 6, 'OKLAHOMA': 6,
    'TX': 6, 'TEXAS': 6,
    # Region 7 – Heartland
    'IA': 7, 'IOWA': 7,
    'KS': 7, 'KANSAS': 7,
    'MO': 7, 'MISSOURI': 7,
    'NE': 7, 'NEBRASKA': 7,
    # Region 8 – Mountain
    'CO': 8, 'COLORADO': 8,
    'MT': 8, 'MONTANA': 8,
    'ND': 8, 'NORTH DAKOTA': 8,
    'SD': 8, 'SOUTH DAKOTA': 8,
    'UT': 8, 'UTAH': 8,
    'WY': 8, 'WYOMING': 8,
    # Region 9 – Pacific
    'AZ': 9, 'ARIZONA': 9,
    'CA': 9, 'CALIFORNIA': 9,
    'HI': 9, 'HAWAII': 9,
    'NV': 9, 'NEVADA': 9,
    # Region 10 – Northwest
    'AK': 10, 'ALASKA': 10,
    'ID': 10, 'IDAHO': 10,
    'OR': 10, 'OREGON': 10,
    'WA': 10, 'WASHINGTON': 10,
}

CITY_TO_STATE = {
    'NEW YORK CITY': 'NY', 'NYC': 'NY',
    'BOSTON': 'MA',
    'CHICAGO': 'IL',
    'LOS ANGELES': 'CA', 'LA': 'CA',
    'HOUSTON': 'TX',
    'PHOENIX': 'AZ',
    'PHILADELPHIA': 'PA', 'PHILLY': 'PA',
    'SAN ANTONIO': 'TX',
    'SAN DIEGO': 'CA',
    'DALLAS': 'TX',
    'SEATTLE': 'WA',
    'DENVER': 'CO',
    'NASHVILLE': 'TN',
    'BALTIMORE': 'MD',
    'ATLANTA': 'GA',
    'MIAMI': 'FL',
    'MINNEAPOLIS': 'MN',
    'PORTLAND': 'OR',
    'LAS VEGAS': 'NV',
    'MEMPHIS': 'TN',
    'LOUISVILLE': 'KY',
    'MILWAUKEE': 'WI',
    'ALBUQUERQUE': 'NM',
    'TUCSON': 'AZ',
    'FRESNO': 'CA',
    'SACRAMENTO': 'CA',
    'MESA': 'AZ',
    'KANSAS CITY': 'MO',
    'OMAHA': 'NE',
    'RALEIGH': 'NC',
    'COLORADO SPRINGS': 'CO',
    'VIRGINIA BEACH': 'VA',
    'LONG BEACH': 'CA',
    'OLYMPIA': 'WA',
    'RICHMOND': 'VA',
    'HARTFORD': 'CT',
    'CONCORD': 'NH',
    'MONTPELIER': 'VT',
    'PROVIDENCE': 'RI',
    'AUGUSTA': 'ME',
    'ALBANY': 'NY',
    'TRENTON': 'NJ',
    'DOVER': 'DE',
    'ANNAPOLIS': 'MD',
    'CHARLESTON': 'WV',
    'COLUMBIA': 'SC',
    'JACKSON': 'MS',
    'MONTGOMERY': 'AL',
    'TALLAHASSEE': 'FL',
    'FRANKFORT': 'KY',
    'INDIANAPOLIS': 'IN',
    'LANSING': 'MI',
    'COLUMBUS': 'OH',
    'SPRINGFIELD': 'IL',
    'MADISON': 'WI',
    'ST. PAUL': 'MN', 'SAINT PAUL': 'MN',
    'DES MOINES': 'IA',
    'TOPEKA': 'KS',
    'JEFFERSON CITY': 'MO',
    'LINCOLN': 'NE',
    'BISMARCK': 'ND',
    'PIERRE': 'SD',
    'CHEYENNE': 'WY',
    'HELENA': 'MT',
    'BOISE': 'ID',
    'SANTA FE': 'NM',
    'OKLAHOMA CITY': 'OK',
    'LITTLE ROCK': 'AR',
    'BATON ROUGE': 'LA',
    'AUSTIN': 'TX',
    'JUNEAU': 'AK',
    'HONOLULU': 'HI',
    'SALT LAKE CITY': 'UT',
    'SACRAMENTO': 'CA',
    'CARSON CITY': 'NV',
    'CINCINNATI': 'OH',
    'LA JOLLA': 'CA',   'LAJOLLA': 'CA',
    'SAN FRANCISCO': 'CA',
    'ST. LOUIS': 'MO',  'SAINT LOUIS': 'MO',
    'STATE COLLEGE': 'PA',
    'STANFORD': 'CA',
    'GAMBIER': 'OH',
    'ANCHORAGE': 'AK',
    'BETHEL': 'AK',
    'UTQIAGVIK': 'AK',
    'SOLDOTNA': 'AK',
    'LOWER KASILOF RIVER': 'AK',
    'LOWER KENAI RIVER': 'AK',
}

# usa 2020 census state populations for weighted subsampling
STATE_POPULATION = {
    'AL': 5024279, 'AK': 733391, 'AZ': 7151502, 'AR': 3011524,
    'CA': 39538223, 'CO': 5773714, 'CT': 3605944, 'DE': 989948,
    'FL': 21538187, 'GA': 10711908, 'HI': 1455271, 'ID': 1839106,
    'IL': 12812508, 'IN': 6785528, 'IA': 3190369, 'KS': 2937880,
    'KY': 4505836, 'LA': 4657757, 'ME': 1362359, 'MD': 6177224,
    'MA': 7029917, 'MI': 10077331, 'MN': 5706494, 'MS': 2961279,
    'MO': 6154913, 'MT': 1084225, 'NE': 1961504, 'NV': 3104614,
    'NH': 1377529, 'NJ': 9288994, 'NM': 2117522, 'NY': 20201249,
    'NC': 10439388, 'ND': 779094, 'OH': 11799448, 'OK': 3959353,
    'OR': 4237256, 'PA': 13002700, 'RI': 1097379, 'SC': 5118425,
    'SD': 886667, 'TN': 6910840, 'TX': 29145505, 'UT': 3271616,
    'VT': 643077, 'VA': 8631393, 'WA': 7705281, 'WV': 1793716,
    'WI': 5893718, 'WY': 576851, 'DC': 689545,
}


# collects unresolved USA: values during a run, can print after if results look sus
UNMATCHED_USA: set[str] = set()


def parse_usa_state(country_str, _unmatched_log=None):
    """
    Parse a USA country field and return the two-letter state abbreviation,
    or None if it cannot be determined.

    Extended behaviour:
    - Substring matching: if a known state name/abbreviation or city name
      appears *anywhere* in the remainder string, it is accepted.
    - If the string starts with 'USA:' but no state can be resolved, the
      remainder is appended to _unmatched_log (if provided) so callers can
      inspect/print the unknown values.
    """
    if not isinstance(country_str, str):
        return None

    s = country_str.strip().upper()

    if not re.match(r'^USA\b', s):
        return None

    remainder = re.sub(r'^USA\s*:?\s*', '', s).strip()

    if not remainder:
        return None

    def _resolve_token(token):
        token = token.strip()
        if re.fullmatch(r'[A-Z]{2}', token) and token in STATE_TO_HHS:
            return token
        if token in STATE_TO_HHS:
            return _abbrev_from_name(token)
        if token in CITY_TO_STATE:
            return CITY_TO_STATE[token]
        return None

    # 1. Comma-split tokens (e.g. "TX, HOUSTON" or "HOUSTON, TX")
    if ',' in remainder:
        parts = [p.strip() for p in remainder.split(',')]
        for part in reversed(parts):
            result = _resolve_token(part)
            if result:
                return result
        # fall through to substring scan below

    # 2. Exact match on the whole remainder
    result = _resolve_token(remainder)
    if result:
        return result

    # 3. Substring scan – state full names (longest first to avoid partials)
    for name in sorted(
        (k for k in STATE_TO_HHS if len(k) > 2), key=len, reverse=True
    ):
        if re.search(r'\b' + re.escape(name) + r'\b', remainder):
            abbrev = _abbrev_from_name(name)
            if abbrev:
                return abbrev

    # 4. Substring scan – two-letter abbreviations (word-boundary guarded)
    for abbrev in (k for k in STATE_TO_HHS if len(k) == 2):
        if re.search(r'(?<![A-Z])' + abbrev + r'(?![A-Z])', remainder):
            return abbrev

    # 5. Substring scan – city names (longest first to avoid partial matches)
    for city in sorted(CITY_TO_STATE.keys(), key=len, reverse=True):
        if re.search(r'\b' + re.escape(city) + r'\b', remainder):
            return CITY_TO_STATE[city]

    # 6. Unmatched – log the remainder for inspection
    if _unmatched_log is not None:
        _unmatched_log.add(remainder)

    return None

_NAME_TO_ABBREV: dict[str, str] = {}
for _k in STATE_TO_HHS:
    if len(_k) == 2:
        _NAME_TO_ABBREV[_k] = _k  # abbreviation maps to itself

_FULL_NAMES = {k: k for k in STATE_TO_HHS if len(k) == 2}
for _name in list(STATE_TO_HHS.keys()):
    if len(_name) > 2:
        # find the 2-letter key with the same HHS region value
        _region = STATE_TO_HHS[_name]
        for _abbrev, _r in STATE_TO_HHS.items():
            if len(_abbrev) == 2 and _r == _region:
                # only assign if we haven't already (first match wins)
                if _name not in _NAME_TO_ABBREV:
                    _NAME_TO_ABBREV[_name] = _abbrev
 
 
def _abbrev_from_name(full_name: str) -> str | None:
    """Return the 2-letter abbreviation for a full US state name, or None."""
    return _NAME_TO_ABBREV.get(full_name)


LONG_READ_PLATFORMS = {'OXFORD_NANOPORE', 'PACBIO_SMRT'}


def has_long_read(platform_str) -> bool:
    """Return True if the instrument_platform indicates long-read sequencing."""
    if not isinstance(platform_str, str):
        return False
    return platform_str.strip().upper() in LONG_READ_PLATFORMS


def _to_decimal_year(year: int, month: int = 6, day: int = 15) -> float:
    """Convert a calendar date to a decimal year (e.g. 2015.37)."""
    date = _dt.date(year, month, day)
    year_start = _dt.date(year, 1, 1)
    year_end   = _dt.date(year + 1, 1, 1)
    return year + (date - year_start).days / (year_end - year_start).days


def _parse_decimal_year(s: str) -> float | None:
    """Return a decimal year for a YYYY, YYYY-MM, or YYYY-MM-DD string."""
    s = s.strip()[:10]
    try:
        if re.fullmatch(r'\d{4}-\d{2}-\d{2}', s):
            yr, mo, dy = int(s[:4]), int(s[5:7]), int(s[8:10])
            return _to_decimal_year(yr, mo, dy)
        if re.fullmatch(r'\d{4}-\d{2}', s):
            yr, mo = int(s[:4]), int(s[5:7])
            if 1 <= mo <= 12:
                return _to_decimal_year(yr, mo, 15)
        if re.fullmatch(r'\d{4}', s):
            return float(s) + 0.5
    except (ValueError, OverflowError):
        pass
    return None

def normalize_collection_date(date_str):
    """
    Return collection_date as YYYY-MM-DD (or as much as is known).
    Handles DD-Mon-YYYY (e.g. '15-Jul-2009') in addition to standard formats.
    """
    if pd.isna(date_str):
        return date_str
    d = str(date_str).strip()
    m = re.fullmatch(r'(\d{1,2})-([A-Za-z]{3})-(\d{4})', d)
    if m:
        try:
            dt = _dt.datetime.strptime(d, '%d-%b-%Y')
            return dt.strftime('%Y-%m-%d')
        except ValueError:
            pass
    return d


# treetime date format helper
def format_date_treetime(d):
    """
    Convert a collection_date string to a TreeTime-compatible format.

    TreeTime accepts: YYYY-MM-DD, YYYY-MM, YYYY, or a decimal year float.

    Edge cases handled:
    - YYYY-WW  (ISO week notation, e.g. '2015-17'): converted to decimal year
      of the Wednesday of that week.
    - Slash-delimited ranges (e.g. '2015-01/2016-07'): midpoint returned as a
      decimal year. The bracket range format [A:B] triggers a bug in TreeTime
      when dates lack a day component, so we avoid it entirely.
    - YYYY only: returned as decimal year (mid-year).
    """
    if pd.isna(d):
        return 'XX'

    d = str(d).strip()

    # slash-delimited range e.g. "2020-01/2020-03"
    if '/' in d:
        parts = [p.strip() for p in d.split('/')]
        dec_years = [_parse_decimal_year(p) for p in parts]
        dec_years = [y for y in dec_years if y is not None]
        if dec_years:
            return str(round(sum(dec_years) / len(dec_years), 4))
        return 'XX'

    # DD-Mon-YYYY (e.g. "15-Jul-2009")
    m = re.fullmatch(r'(\d{1,2})-([A-Za-z]{3})-(\d{4})', d)
    if m:
        try:
            dt = _dt.datetime.strptime(d.strip(), '%d-%b-%Y')
            return dt.strftime('%Y-%m-%d')
        except ValueError:
            pass

    # full date YYYY-MM-DD
    if re.fullmatch(r'\d{4}-\d{2}-\d{2}', d):
        return d

    # YYYY-MM or YYYY-WW
    if re.fullmatch(r'\d{4}-\d{2}', d):
        num = int(d[5:7])
        year = int(d[:4])
        if 1 <= num <= 12:
            return d  # valid month → pass through as YYYY-MM
        # ISO week number (e.g. 2015-17 means week 17 of 2015)
        try:
            # %W treats week 1 as the week containing the first Monday;
            # use Wednesday (%w=3) as the representative day
            dt = _dt.datetime.strptime(f'{year}-{num:02d}-3', '%Y-%W-%w')
            return str(round(_to_decimal_year(dt.year, dt.month, dt.day), 4))
        except ValueError:
            return str(year)  # fall back to year only

    # YYYY only → mid-year decimal
    if re.fullmatch(r'\d{4}', d):
        return str(float(d) + 0.5)

    # try pandas as a last resort
    try:
        dt = pd.to_datetime(d)
        return dt.strftime('%Y-%m-%d')
    except Exception:
        return 'XX'


# weighted subsampling
def weighted_subsample(df, n, random_state=42):
    """
    Two-stage downsampling to at most *n* rows:

    Stage 1 - temporal allocation (by month):
        Distribute the n slots evenly across all months that have at least one
        eligible sample, giving each month an equal quota. Any remainder slots
        from months that have fewer samples than their quota are redistributed
        to months that still have capacity, repeating until all n slots are
        allocated or no capacity remains.

    Stage 2 - within-month state sampling (by population):
        Within each month, draw the allocated quota from available samples
        using weights proportional to US Census state population. This means
        larger states are more likely to be selected within a given month, but
        every month gets an equal opportunity to contribute.

    Only rows with a non-null 'state' column are eligible.
    Returns the sampled DataFrame (index reset).
    """
    rng = np.random.default_rng(random_state)

    eligible = df[df['state'].notna() & (df['state'] != '')].copy()

    if eligible.empty:
        return eligible

    if len(eligible) <= n:
        return eligible.reset_index(drop=True)

    # Stage 1: assign a quota to each month
    eligible['_ym'] = pd.to_datetime(
        eligible['collection_date'], errors='coerce'
    ).dt.to_period('M').astype(str)

    months = sorted(eligible['_ym'].dropna().unique())
    n_months = len(months)

    # start with a floor-equal quota for every month
    base_quota = n // n_months
    remainder  = n % n_months
    quotas = {m: base_quota for m in months}

    # distribute the remainder one slot at a time to the months
    # with the most samples (most capacity to absorb it)
    month_sizes = eligible.groupby('_ym').size()
    for m in month_sizes.sort_values(ascending=False).index[:remainder]:
        quotas[m] += 1

    # cap each quota at the actual number of samples in that month
    # and carry over any surplus to months still under-sampled
    capped = True
    while capped:
        capped = False
        surplus = 0
        has_capacity = []
        for m in months:
            available = int(month_sizes.get(m, 0))
            if quotas[m] > available:
                surplus += quotas[m] - available
                quotas[m] = available
                capped = True
            elif quotas[m] < available:
                has_capacity.append(m)
        # redistribute surplus evenly across months with remaining capacity
        if surplus > 0 and has_capacity:
            extra = surplus // len(has_capacity)
            leftover = surplus % len(has_capacity)
            for m in has_capacity:
                quotas[m] += extra
            for m in has_capacity[:leftover]:
                quotas[m] += 1

    # Stage 2: within each month, sample by state population weight
    pieces = []
    for m, quota in quotas.items():
        if quota == 0:
            continue
        pool = eligible[eligible['_ym'] == m].copy()
        pop_weights = pool['state'].map(
            lambda s: STATE_POPULATION.get(s, 1)
        ).astype(float).values
        pop_weights = pop_weights + 1e-6
        pop_weights = pop_weights / pop_weights.sum()

        n_draw = min(quota, len(pool))
        chosen_idx = rng.choice(
            len(pool),
            size=n_draw,
            replace=False,
            p=pop_weights,
        )
        chosen = pool.iloc[chosen_idx]
        pieces.append(chosen)

    sampled = pd.concat(pieces, ignore_index=True)
    sampled = sampled.drop(columns=['_ym'])
    return sampled.reset_index(drop=True)

# no of samples to retain per ST after downsampling
N_SUBSAMPLE = 100

assembly_index = {}
for root, dirs, files in os.walk('./assemblies'):
    for file in files:
        if file.endswith('.fa.gz'):
            assembly_index[file] = os.path.join(root, file)

for meta_file in glob.glob('./tables/*metadata.csv'):
    targets = set()
    df = pd.read_csv(meta_file)
    df['collection_date'] = df['collection_date'].apply(normalize_collection_date)
    df[['year', 'month', 'day']] = (
        df['collection_date']
        .str.split('-', expand=True)
        .reindex(columns=[0, 1, 2])
    )
    has_month = (
        df[df['month'].notnull() & df['country'].notnull()]
        .drop_duplicates(subset='sample_accession')
    )
    targets.update(has_month['sample_accession'].tolist())

    has_month = has_month.copy()
    has_month['state'] = has_month['country'].apply(
    lambda s: parse_usa_state(s, _unmatched_log=UNMATCHED_USA)
)
    has_month['hhs_region'] = has_month['state'].apply(
        lambda s: STATE_TO_HHS.get(s) if s else None
    )
    if 'instrument_platform' in has_month.columns:
        has_month['has_long_read'] = has_month['instrument_platform'].apply(has_long_read)
    else:
        has_month['has_long_read'] = False

    rows = []
    for target in targets:
        matches = [path for fname, path in assembly_index.items() if target in fname]
        for path in matches:
            rows.append({'sample_accession': target, 'path': path})

    if rows:
        out_df = pd.DataFrame(rows)
        basename = os.path.basename(meta_file)
        prefix = basename.replace('_metadata.csv', '')
        out_path = f'./tables/{prefix}_isolates.list'
        out_df.to_csv(out_path, sep='\t', index=False, header=False)
        print(f"Wrote {len(out_df)} rows to {out_path}")

        # treetime write
        date_rows = (
            has_month[has_month['sample_accession'].isin(targets)]
            [['sample_accession', 'collection_date']]
            .copy()
        )
        date_rows['treetime_date'] = date_rows['collection_date'].apply(
            format_date_treetime
        )

        dates_path = f'./tables/{prefix}_dates.txt'
        with open(dates_path, 'w') as f:
            f.write("name,date\n")
            for _, row in date_rows.iterrows():
                f.write(f"{row['sample_accession']},{row['treetime_date']}\n")
        print(f"Wrote {len(date_rows)} dates to {dates_path}")

        meta_out_path = f'./tables/{prefix}_metadata_with_hhs.csv'
        has_month.to_csv(meta_out_path, index=False)
        print(f"Wrote metadata with HHS regions to {meta_out_path}")

    else:
        print(f"No matches found for {meta_file}")

os.makedirs('./gubbins', exist_ok=True)

for ref_file in glob.glob('./reference_genomes/*_reference.fa*'):
    st_name = os.path.basename(ref_file).split('_reference.fa')[0]
    if os.path.exists(os.path.join(f'./gubbins', '{st_name}.final_tree.tre')):
        print(f"Skipping {st_name}, results already exist.")
        continue

    # load metadata and isolate list for this ST
    meta_path    = f'./tables/{st_name}_metadata.csv'
    isolates_path = f'./tables/{st_name}_isolates.list'

    if not os.path.exists(meta_path) or not os.path.exists(isolates_path):
        print(f"  Missing metadata or isolate list for {st_name}, skipping.")
        continue

    st_meta = pd.read_csv(meta_path)
    st_meta['collection_date'] = st_meta['collection_date'].apply(normalize_collection_date)
    st_meta[['year', 'month', 'day']] = (
        st_meta['collection_date']
        .str.split('-', expand=True)
        .reindex(columns=[0, 1, 2])
    )
    st_meta = st_meta[
        st_meta['month'].notnull() & st_meta['country'].notnull()
    ].drop_duplicates(subset='sample_accession').copy()
    st_meta['state'] = st_meta['country'].apply(
    lambda s: parse_usa_state(s, _unmatched_log=UNMATCHED_USA)
)
    if 'instrument_platform' in st_meta.columns:
        st_meta['has_long_read'] = st_meta['instrument_platform'].apply(has_long_read)
    else:
        st_meta['has_long_read'] = False

    # only usa first
    st_meta_usa = st_meta[st_meta['state'].notna()].copy()

    if st_meta_usa.empty:
        print(f"  No USA samples with state info for {st_name}, skipping.")
        continue

    # weighted downsample
    sampled = weighted_subsample(st_meta_usa, n=N_SUBSAMPLE)
    n_long_read = sampled['has_long_read'].sum() if 'has_long_read' in sampled.columns else 0
    print(
        f"  {st_name}: {len(st_meta_usa)} USA samples with state → "
        f"downsampled to {len(sampled)} "
        f"({n_long_read} with long-read data)"
    )

    # write filtered isolate list
    isolates_df = pd.read_csv(
        isolates_path, sep='\t', header=None,
        names=['sample_accession', 'path']
    )
    sampled_isolates = isolates_df[
        isolates_df['sample_accession'].isin(sampled['sample_accession'])
    ]

    downsampled_list_path = f'./tables/{st_name}_isolates_downsampled.list'
    sampled_isolates.to_csv(
        downsampled_list_path, sep='\t', index=False, header=False
    )
    print(f"  Wrote {len(sampled_isolates)} rows to {downsampled_list_path}")

    # run SKA alignment on the downsampled set
    print(f"  Running Gubbins for {st_name}...")
    os.system(
        f"generate_ska_alignment.py "
        f"--reference {ref_file} "
        f"--input {downsampled_list_path} "
        f"--out {os.path.join('./gubbins', st_name)}"
    )

    os.system(
        f"run_gubbins.py "
        f"--prefix ./gubbins/{st_name} "
        f"--first-tree-builder iqtree-fast "
        f"--first-model GTR "
        f"--tree-builder raxmlng "
        f"--model GTR "
        f"--date ./tables/{st_name}_dates.txt "
        f"--iterations 10 "
        f"--converge-method recombination "
        f"--min-snps 2 "
        f"--threads 8 "
        f"./gubbins/{st_name}"
    )

print("Finished. Files saved to: gubbins/")