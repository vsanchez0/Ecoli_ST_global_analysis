import pandas as pd
import os
import glob

os.makedirs("rgi_results", exist_ok=True)

for file in glob.glob('./assemblies/*/*'):
    if file.endswith('.fa.gz'):
        sample_id = file.split('/')[-1].split('.')[0]
        if os.path.exists(f"rgi_results/{sample_id}.txt"):
            print(f"Skipping {sample_id}, results already exist.")
            continue
        os.system(f"rgi main -i {file} -o rgi_results/{sample_id} -n 8 --clean")

dfs = []
for result in glob.glob('./rgi_results/*.txt'):
    df = pd.read_csv(result, sep='\t')
    df['IsolateID'] = os.path.basename(result).split('.')[0]
    dfs.append(df)

rgi_df = pd.concat(dfs, ignore_index=True)
rgi_df_filtered = rgi_df.loc[rgi_df['Best_Identities'] >= 85]
rgi_df.to_csv("tables/combined_rgi_results.tsv", sep='\t', index=False)
rgi_df_filtered.to_csv("tables/combined_rgi_results_filtered.tsv", sep='\t', index=False)

print("\nFinished. Files saved to: tables/combined_rgi_results.tsv")