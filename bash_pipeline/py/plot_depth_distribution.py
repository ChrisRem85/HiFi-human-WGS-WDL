#!/usr/bin/env python3
"""Port of mosdepth task's inline depth-distribution plot (mosdepth.wdl)."""
import sys
import pandas as pd
import seaborn as sns
import matplotlib
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt

regions_bed_gz, title, out_png = sys.argv[1], sys.argv[2], sys.argv[3]

sns.set_theme(style='darkgrid')
df = pd.read_csv(
    regions_bed_gz,
    sep='\t',
    names=['chr', 'start', 'end', 'depth'],
    usecols=['depth'],
    dtype={'depth': 'float32'},
    compression='gzip',
)
xmax = int(2 * df[df['depth'] > 0]['depth'].mode().values[0])  # 2x non-zero mode
df = df[(df['depth'] <= xmax)]
fig, axs = plt.subplots(2, 1, figsize=(8, 6))
sns.histplot(df, x='depth', bins=xmax - 1, stat='proportion', ax=axs[0])
axs[0].set_xlim(0, xmax)
axs[0].set_xlabel('')
axs[0].set_xticklabels([])
sns.ecdfplot(df, x='depth', complementary=True, ax=axs[1])
axs[1].set_xlim(0, xmax)
axs[1].set_yticks(np.arange(0, 1.1, 0.1))
axs[1].set_ylabel('Proportion ≥ depth')
plt.suptitle(f'{title}\nAligned depth distribution')
fig.tight_layout()
plt.savefig(out_png)
plt.close()
