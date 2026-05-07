import pandas as pd
import subprocess
import requests
import sqlite3
import os
import glob

OUT_DIR = "assemblies"
ATB_FILE_INDEX = "https://osf.io/download/r6gcp"
os.makedirs(OUT_DIR, exist_ok=True)

# all relevant samples with month-level collection date and country info, across the 5 most common sts
target_samples = set()
for meta_file in glob.glob('./tables/*metadata.csv'):
    df = pd.read_csv(meta_file)
    df[['year', 'month', 'day']] = df['collection_date'].str.split('-', expand=True).reindex(columns=[0,1,2])
    has_month = df[df['month'].notnull() & df['country'].notnull()].drop_duplicates(subset='sample_accession')
    target_samples.update(has_month['sample_accession'].tolist())

print(f"Total unique samples with month info: {len(target_samples)}")

# download atb index file
print("Downloading ATB file index...")
atb_index = pd.read_csv(ATB_FILE_INDEX, sep='\t')
print(f"Index loaded: {len(atb_index)} files")

# filter to assembly tarballs only
asm_files = atb_index[
    atb_index['project'].str.contains('Assembly', na=False) &
    atb_index['filename'].str.endswith('.tar.xz')
].copy()
print(f"Assembly tarballs found: {len(asm_files)}")
print(asm_files[['project', 'filename']].head(10).to_string())

checkpoint_file = "downloaded_samples.txt"
already_done = set()
if os.path.exists(checkpoint_file):
    with open(checkpoint_file) as f:
        already_done = set(line.strip() for line in f)

remaining = target_samples - already_done
print(f"Samples still to download: {len(remaining)}")

# check inside each tarball for samples, download, then checkpoint
for _, row in asm_files.iterrows():
    if not remaining:
        print("All target samples downloaded!")
        break

    tarball_name = row['filename']
    tarball_url  = row['url']
    tarball_path = os.path.join(OUT_DIR, os.path.basename(tarball_name))

    batch_dir = os.path.basename(tarball_name).replace('.tar.xz', '')

    expected = {f"{batch_dir}/{s}.fa": s for s in remaining}

    print(f"\nDownloading {tarball_name} ({row['size(MB)']:.0f} MB)...")
    r = requests.get(tarball_url, stream=True)
    with open(tarball_path, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)

    result = subprocess.run(['tar', '-tf', tarball_path], capture_output=True, text=True)
    tarball_contents = set(result.stdout.strip().split('\n'))

    to_extract = [f for f in expected if f in tarball_contents]

    if not to_extract:
        print(f"  No target samples in this tarball, skipping.")
        os.remove(tarball_path)
        continue

    print(f"  Extracting {len(to_extract)} samples...")
    subprocess.run(['tar', '-xf', tarball_path, '-C', OUT_DIR] + to_extract)

    for f in to_extract:
        fa_path = os.path.join(OUT_DIR, f)
        if os.path.exists(fa_path):
            subprocess.run(['gzip', fa_path])
            sample_id = expected[f]
            with open(checkpoint_file, 'a') as cp:
                cp.write(sample_id + '\n')
            remaining.discard(sample_id)

    os.remove(tarball_path)
    print(f"  Done. Remaining samples: {len(remaining)}")

print(f"\nFinished. Files saved to: {OUT_DIR}/")
if remaining:
    print(f"WARNING: {len(remaining)} samples were not found in any tarball:")
    print(list(remaining)[:10])