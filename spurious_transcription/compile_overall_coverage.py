#!/usr/bin/env python3

"""
> compile_overall_coverage.py <

Given the per-position coverages from the desired input files, assign equal
weightage to each file, and compile overall trends from all datasets.

Implements several filters:
I.   Gene has to have >= NUM_EXONS exon
II.  Genes must be expressed in all tested exons.
III. Gene has to have mean RPKM > 0.5 overall
IV.  Gene has to pass all criteria above in all files.
"""
import argparse
import csv
import glob
import statistics

import numpy as np

import natural_sort

parser = argparse.ArgumentParser(description="""
Given the per-position coverages from the desired input files, assign equal
weightage to each file, and compile overall trends from all datasets.""")

parser.add_argument('input_tsvs', metavar='input_tsvs',
                    type=argparse.FileType('r'), nargs='+',
                    help='individual files to be compiled.')
args = parser.parse_args()

NUM_EXONS = 6
NUM_FILES = len(args.input_tsvs)

norm_cov = {}

# read input files
for f in args.input_tsvs:
    tsv_reader = csv.reader(f, delimiter='\t')

    for row in tsv_reader:
        if not row: continue
    
        gene_name = row[0]
        exon_covs = np.array([float(x) for x in row[1:]], dtype=np.float64)
        
        # abandon ship if gene has fewer than NUM_EXONS (criteria I), or 
        # not expressed in all tested exons (criteria II), or not having
        # sufficient overall coverage (criteria III)
        if len(exon_covs) < NUM_EXONS: continue
        if not all(exon_covs[:NUM_EXONS]): continue
        if np.mean(exon_covs) < 0.5: continue
        
        # normalise against coverage of first exon (fold change)
        exon_covs /= exon_covs[0]
        
        # calculate ln (FC)
        exon_covs = np.log(exon_covs)
        
        # ... but change first exon back to 1 to act as a counter, as a check
        # for whether all files pass the criteria
        exon_covs[0] = 1
        
        # add the coverages to norm_cov
        if gene_name not in norm_cov:
            norm_cov[gene_name] = exon_covs
        else:
            # add normalised coverages across files
            # code looks complicated to deal with arrays of different lengths
            # (which arises because, say, file 1 had data from exons 1-13 but
            # file 2 had data from exons 1-12)
            try:
                norm_cov[gene_name] += exon_covs
            except:
                temp = norm_cov[gene_name].copy()
                if len(temp) > len(exon_covs):
                    temp[:len(exon_covs)] += exon_covs
                    norm_cov[gene_name] = temp
                else:
                    exon_covs[:len(temp)] += temp
                    norm_cov[gene_name] = exon_covs

for g in natural_sort.natural_sort(norm_cov):
    # criteria IV can be detected by looking at normalised coverage of exon 1
    # if n == # of files, it means it passed I, II, III in all files.
    if norm_cov[g][0] < NUM_FILES: continue
    
    # yay, if gene reaches this part, it passed all criteria
    norm_cov[g][0] = 0
    
    # print stuff out
    print (g, *[round(x/NUM_FILES, 3) for x in norm_cov[g]], sep='\t')
