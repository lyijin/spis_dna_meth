#!/usr/bin/env python3

"""
> compile_covs_for_glm.py <

Script compiles all relevant cov files into a file, which can then be fed
into the R GLM script.
"""
import csv
import glob

treatments = ['pH_7.2', 'pH_7.4', 'pH_7.8', 'pH_8.1']

for t in treatments:
    for f in glob.glob('../*{}*.final.annot.cov'.format(t)):
        tsv_reader = csv.reader(open(f), delimiter='\t')
        
        for row in tsv_reader:
            if not row: continue
            
            # filter for genic positions; ignore intergenic ones and "no_info"
            if row[6] != 'intergenic' and row[6] != 'no_info':
                if row[10][:4] == 'Exon': 
                    parsed_output = '\t'.join([row[0], row[1], row[4], row[5],
                                               row[6], t])
                    print (parsed_output)
