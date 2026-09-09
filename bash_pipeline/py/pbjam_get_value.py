#!/usr/bin/env python3
"""Port of pbjam_bam_stats's summary-stat extraction (pbjam.wdl)."""
import json
import decimal
import sys

data = json.load(open(sys.argv[1]), parse_float=decimal.Decimal)
print(data['combined_stats']['summary_stats'][sys.argv[2]])
