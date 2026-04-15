import pandas as pd
import requests
from io import StringIO
import os
import numpy as np
import sqlite3
import matplotlib.pyplot as plt
import seaborn as sns

mlst_calls = pd.read_csv('allMlstCalls.tsv.gz', sep='\t', compression='gzip', header=None, low_memory=False)

# want to keep only STs with at least 10,000 samples
counts = mlst_calls[2].value_counts().to_dict()
sts_with_10k = {st: count for st, count in counts.items() if count >= 10000}

st_groups = {}
for st_num, _ in sts_with_10k.items():
    st_groups[f"st{st_num}"] = (
        mlst_calls[mlst_calls[2].astype(str) == str(st_num)][0]
        .str.split('.').str[0]
        .tolist()
    )

# want to fetch collection dates and country for all samples in these STs, and also instrument platform and library strategy (to find long-read data)
con = sqlite3.connect("atb.metadata.202408.sqlite")

os.makedirs("tables", exist_ok=True)

for name, st in st_groups.items():
    placeholders = ",".join("?" * len(st))
    query = f"""
        SELECT 
            sample_accession,
            run_accession,
            collection_date,
            collection_date_start,
            collection_date_end,
            country,
            strain,
            isolation_source,
            instrument_platform,
            instrument_model,
            library_strategy
        FROM ena_20240801
        WHERE sample_accession IN ({placeholders})
        AND sample_accession IS NOT NULL
    """
    meta = pd.read_sql(query, con, params=st)
    meta.to_csv(f"tables/{name}_metadata.csv", index=False)
    print(f"Saved to tables/{name}_metadata.csv")

con.close()

dfs = []
for name, st in st_groups.items():
    st_dates = pd.read_csv(f'tables/{name}_metadata.csv')
    st_dates[['year', 'month', 'day']] = st_dates['collection_date'].str.split('-', expand=True).reindex(columns=[0,1,2])

    # need samples with granular temporal data (at least month-level) and country info to do spatiotemporal analyses
    st_with_months = st_dates[st_dates['month'].notnull() & st_dates['country'].notnull()]
    unique_with_months = st_with_months.drop_duplicates(subset='sample_accession')

    # helpful for american epidemiology, and also to see if we can get more granular location info (state-level) for some samples
    american = unique_with_months[unique_with_months['country'].str.startswith('USA')]
    american_with_state = american[american['country'].str.split(':').str.len() > 1]

    # how many long-read samples available
    longread_samples = (
        st_dates[st_dates['instrument_platform'].isin(['OXFORD_NANOPORE', 'PACBIO_SMRT'])]
        ['sample_accession'].unique()
    )
    with_month_and_longread = unique_with_months[
        unique_with_months['sample_accession'].isin(longread_samples)
    ]

    dfs.append({
        'name': name,
        'unique_samples_with_month': unique_with_months['sample_accession'].nunique(),
        'american_with_month_and_state': american_with_state['sample_accession'].nunique(),
        'with_month_and_longread': with_month_and_longread['sample_accession'].nunique(),
    })

results_df = pd.DataFrame(dfs)
print(results_df)

os.makedirs("figures", exist_ok=True)

cols = ['unique_samples_with_month', 'american_with_month_and_state', 'with_month_and_longread']
labels = ['With month info', 'With month + US state', 'With month + long-read']

x = np.arange(len(results_df))
width = 0.25

fig, ax = plt.subplots(figsize=(5, 4))

for i, (col, label) in enumerate(zip(cols, labels)):
    ax.bar(x + i * width, results_df[col], width=width, label=label)

ax.set_xticks(x + width)
ax.set_xticklabels(results_df['name'])
ax.set_ylabel("Unique samples")
ax.set_title("Sample metadata availability by ST")
ax.legend()
plt.tight_layout()
plt.savefig("figures/st_metadata_summary.png", dpi=300)

print("\nFinished. Files saved to: figures/st_metadata_summary.png")