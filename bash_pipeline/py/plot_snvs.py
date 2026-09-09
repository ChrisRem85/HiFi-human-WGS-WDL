#!/usr/bin/env python3
"""Port of bcftools_stats_roh_small_variants's SNV distribution plot (bcftools.wdl)."""
import sys
import pandas as pd
import seaborn as sns
import matplotlib
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt

input_tsv, title, out_png = sys.argv[1], sys.argv[2], sys.argv[3]

df = pd.read_csv(input_tsv, sep='\t')
df[['REF', 'ALT']] = df['type'].str.split('>', expand=True)
df = pd.pivot(df, index='ALT', columns='REF', values='count')
sns.set_style('dark')
mask = np.identity(df.shape[0], dtype=bool)
fig, ax = plt.subplots(figsize=(8, 6))
sns.heatmap(df, mask=mask, cmap='coolwarm', annot=True, fmt=',', annot_kws=dict(fontsize='large'), linewidth=.5, ax=ax)
plt.xlabel('REF base', fontsize='large')
plt.ylabel('ALT base', fontsize='large')
plt.xticks(fontsize='large', rotation=0)
plt.yticks(fontsize='large', rotation=0)
plt.title(f'{title}\nDeepVariant SNV distribution', fontsize='large')
fig.tight_layout()
plt.savefig(out_png)
plt.close()
