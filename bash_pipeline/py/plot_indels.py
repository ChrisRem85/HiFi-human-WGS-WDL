#!/usr/bin/env python3
"""Port of bcftools_stats_roh_small_variants's indel distribution plot (bcftools.wdl)."""
import sys
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

input_tsv, title, out_png = sys.argv[1], sys.argv[2], sys.argv[3]
df = pd.read_csv(input_tsv, sep='\t')


def size_filter(df, col, min_len, max_len):
    return df[(df[col].abs() >= min_len) & (df[col].abs() < max_len)]


def plot_hist(ax, df, min_len, max_len, logy=False, xlabel=True):
    g = sns.histplot(size_filter(df, 'length', min_len, max_len), x='length', weights='count', binwidth=1,
                      binrange=(-max_len - 0.5, max_len + 0.5), ax=ax)
    g.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x/1000)}k' if x >= 1000 else f'{int(x)}'))
    g.set_xlim(-max_len, max_len)
    g.set_title(f'±[{min_len},{max_len}) bp')
    if not xlabel:
        g.set_xlabel('')
        g.tick_params(labelbottom=False)
    if logy:
        g.set_yscale('log')
        g.set_ylabel('log10(Count)')


sns.set_style('darkgrid')
fig, axs = plt.subplots(2, 1, figsize=(8, 6))
plot_hist(axs[0], df, 1, 50, xlabel=False)
plot_hist(axs[1], df, 1, 50, logy=True)
plt.suptitle(f'{title}\nDeepVariant indel distribution, ±[1,50) bp')
fig.tight_layout()
plt.savefig(out_png)
plt.close()
