#!/usr/bin/env python3
"""Port of bcftools_stats_roh_small_variants's ROH-to-BED conversion (bcftools.wdl)."""
import sys

roh_out, min_length, min_qual, out_bed = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), sys.argv[4]

with open(roh_out) as f, open(out_bed, 'w') as out:
    out.write("#chr\tstart\tend\tqual\n")
    for line in f:
        if line.startswith("RG"):
            # RG [2]Sample [3]Chromosome [4]Start [5]End [6]Length (bp) [7]Number of markers [8]Quality
            _, _, chrom, start, end, length, _, score = line.strip().split('\t')
            if int(length) >= min_length and float(score) >= min_qual:
                out.write('\t'.join([chrom, start, end, score]) + '\n')
