#!/usr/bin/env python3
"""Port of hiphase's phase-stat extraction (hiphase.wdl)."""
import sys
import pandas as pd

df = pd.read_csv(sys.argv[1], sep='\t')
print(df[df['chromosome'] == 'all'][sys.argv[2]].values[0])
