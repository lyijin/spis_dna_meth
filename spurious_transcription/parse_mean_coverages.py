#!/usr/bin/env python3

"""
> parse_mean_coverages.py <

Calculate per-gene, per-exon mean coverages in the form
    gene \t exon_1 \t exon_2 \t ...

All numbers sum up to 1 million, so that genes can be compared between
individual files.
"""
import argparse
import csv
import re
import statistics

import natural_sort

parser = argparse.ArgumentParser(description="""
Calculate per-gene, per-exon mean coverages in the form
gene \t exon_1 \t exon_2 \t ...""")

parser.add_argument('input_tsv', metavar='input_tsv',
                    type=argparse.FileType('r'),
                    help='File with per-exon coverages.')

args = parser.parse_args()

coverages = {}

# read input file
tsv_reader = csv.reader(args.input_tsv, delimiter='\t')
next(tsv_reader)
for row in tsv_reader:
    if not row: continue
    
    gene_name = re.match('Spis\d+', row[0]).group(0)
    exon = int(row[2])
    cov = int(row[9])
    
    if gene_name not in coverages:
        coverages[gene_name] = {}
    
    if exon not in coverages[gene_name]:
        coverages[gene_name][exon] = []
    
    coverages[gene_name][exon].append(cov)

# calculate "RPKM" scaling factor
total_coverage = 0
for g in coverages:
    for e in coverages[g]:
        # take the mean of list of coverage values, then store it back
        coverages[g][e] = statistics.mean(coverages[g][e])
        total_coverage += coverages[g][e]

scaling_factor = 1000000 / total_coverage

# print results out
for g in natural_sort.natural_sort(coverages):
    num_exons = max(coverages[g])
    temp = []
    for e in range(1, num_exons + 1):
        if e in coverages[g]:
            temp.append(round(coverages[g][e] * scaling_factor, 3))
        else:
            temp.append(0)
    
    print (g, *temp, sep='\t')
