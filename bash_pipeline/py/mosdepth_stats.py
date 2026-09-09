#!/usr/bin/env python3
"""Port of mosdepth task's mean-depth / sex-inference logic (mosdepth.wdl)."""
import sys
import pandas as pd

summary_txt, max_norm_female_chrY_depth, mean_depth_out, inferred_sex_out = sys.argv[1:5]

df = pd.read_csv(summary_txt, sep='\t')

try:
    mean_depth = df[df['chrom'] == 'total']['mean'].values[0]
except IndexError:
    mean_depth = 0.0
with open(mean_depth_out, 'w') as f:
    f.write(str(mean_depth))

chrA_depth = df[df['chrom'].str.match(r'^(chr)?\d{1,2}$')]['mean'].mean()
chrY_depth = df[df['chrom'].str.match(r'^(chr)?Y$')]['mean'].mean()
if chrA_depth == 0 or pd.isna(chrA_depth):
    inferred_sex = ''
else:
    inferred_sex = 'MALE' if chrY_depth / chrA_depth > float(max_norm_female_chrY_depth) else 'FEMALE'
with open(inferred_sex_out, 'w') as f:
    f.write(inferred_sex)
