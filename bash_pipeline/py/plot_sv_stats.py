#!/usr/bin/env python3
"""Port of sv_stats's SV size-distribution plot (bcftools.wdl)."""
import sys
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import seaborn as sns
from matplotlib.ticker import MaxNLocator, FuncFormatter

len_ins, len_del, len_dup, len_inv, min_length, title, out_png = sys.argv[1:8]
min_length = int(min_length)

sns.set_style('darkgrid')


def fmt_bp(x, pos):
    sign = '-' if x < 0 else ''
    ax_val = abs(x)
    if ax_val >= 1000000:
        return f'{sign}{ax_val/1000000:g}M'
    if ax_val >= 1000:
        return f'{sign}{ax_val/1000:g}k'
    return f'{sign}{ax_val:g}'


def read_lengths(path):
    with open(path) as f:
        values = []
        for line in f.read().split():
            for tok in line.split(','):
                try:
                    values.append(int(tok))
                except ValueError:
                    pass
        return values


N_BINS_TOTAL = 50
CAP = 2000000


def tier_edges(lo, hi):
    w = hi // N_BINS_TOTAL
    return list(range(lo, hi, w)) + [hi]


def plot_del_tier(ax, dele, lo, hi, ylabel=False, legend=False):
    cap = hi is None
    hi_eff = CAP if cap else hi
    edges = tier_edges(lo, hi_eff)
    vals = [min(-x, CAP) for x in dele if -x >= lo] if cap else [-x for x in dele if lo <= -x < hi]
    ax.hist(vals, bins=edges, color='C1', label='DEL')
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel('SV length (bp)')
    ax.xaxis.set_major_formatter(FuncFormatter(fmt_bp))
    if ylabel:
        ax.set_ylabel('Count')
    if legend:
        ax.legend(fontsize='small')
    ax.set_title(f'>={lo:,} bp (capped {CAP:,})' if cap else f'[{lo:,}, {hi:,}) bp', fontsize='small')


def plot_insdup_tier(ax, ins, dup, lo, hi, legend=False):
    cap = hi is None
    hi_eff = CAP if cap else hi
    edges = tier_edges(lo, hi_eff)
    if cap:
        ins_v = [min(x, CAP) for x in ins if x >= lo]
        dup_v = [min(x, CAP) for x in dup if x >= lo]
    else:
        ins_v = [x for x in ins if lo <= x < hi]
        dup_v = [x for x in dup if lo <= x < hi]
    ax.hist([ins_v, dup_v], bins=edges, stacked=True, color=['C0', 'C2'], label=['INS', 'DUP'])
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel('SV length (bp)')
    ax.xaxis.set_major_formatter(FuncFormatter(fmt_bp))
    if legend:
        ax.legend(fontsize='small')
    ax.set_title(f'>={lo:,} bp (capped {CAP:,})' if cap else f'[{lo:,}, {hi:,}) bp', fontsize='small')


def plot_inv_tier(ax, inv, lo, hi, legend=False):
    cap = hi is None
    hi_eff = CAP if cap else hi
    edges = tier_edges(lo, hi_eff)
    vals = [min(x, CAP) for x in inv if x >= lo] if cap else [x for x in inv if lo <= x < hi]
    ax.hist(vals, bins=edges, color='C3', label='INV')
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel('SV length (bp)')
    ax.xaxis.set_major_formatter(FuncFormatter(fmt_bp))
    if legend:
        ax.legend(fontsize='small')
    ax.set_title(f'>={lo:,} bp (capped {CAP:,})' if cap else f'[{lo:,}, {hi:,}) bp', fontsize='small')


ins = read_lengths(len_ins)
dele = read_lengths(len_del)
# SVLEN sign convention varies by caller; force negative here so DEL always plots left of 0.
dele = [-abs(x) for x in dele]
dup = read_lengths(len_dup)
inv = read_lengths(len_inv)

fig, axd = plt.subplot_mosaic(
    [['d1', 'd2', 'd3', 'd4'],
     ['i1', 'i2', 'i3', 'i4'],
     ['v1', 'v2', 'v3', 'v4']],
    figsize=(13, 12),
)
plot_del_tier(axd['d1'], dele, min_length, 1000, ylabel=True, legend=True)
plot_del_tier(axd['d2'], dele, 1000, 10000)
plot_del_tier(axd['d3'], dele, 10000, 100000)
plot_del_tier(axd['d4'], dele, 100000, None)
plot_insdup_tier(axd['i1'], ins, dup, min_length, 1000, legend=True)
plot_insdup_tier(axd['i2'], ins, dup, 1000, 10000)
plot_insdup_tier(axd['i3'], ins, dup, 10000, 100000)
plot_insdup_tier(axd['i4'], ins, dup, 100000, None)
plot_inv_tier(axd['v1'], inv, 0, 1000, legend=True)
plot_inv_tier(axd['v2'], inv, 1000, 10000)
plot_inv_tier(axd['v3'], inv, 10000, 100000)
plot_inv_tier(axd['v4'], inv, 100000, None)
for d_key, i_key in [('d1', 'i1'), ('d2', 'i2'), ('d3', 'i3'), ('d4', 'i4')]:
    shared_max = max(axd[d_key].get_ylim()[1], axd[i_key].get_ylim()[1])
    axd[d_key].set_ylim(0, shared_max)
    axd[i_key].set_ylim(0, shared_max)
fig.suptitle(f'{title}\nStructural variant size distribution')
fig.tight_layout()
plt.savefig(out_png)
plt.close()
