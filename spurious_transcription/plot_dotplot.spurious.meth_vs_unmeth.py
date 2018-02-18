#!/usr/bin/env python3

"""
> plot_dotplot.spurious.meth_vs_unmeth.py <

Takes in output file from compile_overall_coverage.py, and plots the relative
expression of exon 1 vs. subsequent exons, up to NUM_EXONS.

Trends of unmethylated genes and highly methylated genes (median meth > 80%)
are plotted, to see whether methylation helps suppress spurious transcription.
"""
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats
import seaborn as sns

NUM_EXONS = 6

# read in median methylation levels
median_meth = pd.read_table('../bias_density_medians/compiled_median_meth.tsv')

# name standardisation: changes SpisGeneX.Novel --> SpisX
median_meth['gene'] = median_meth['gene'].apply(
    lambda x: x.replace('Gene', '').split('.')[0])
median_meth = median_meth.set_index('gene')

# calculate mean of median meths across 12 samples
median_meth['mean'] = median_meth.mean(axis=1)

# get lists of meth genes
unmeth_genes = median_meth[median_meth['mean'] == 0].index.tolist()
meth_genes = median_meth[median_meth['mean'] > 0].index.tolist()
highly_meth_genes = median_meth[median_meth['mean'] > 80].index.tolist()

# read in mean expression coverages
mean_cov = pd.read_table('mean_cov.all.compiled.tsv',
                         usecols=range(NUM_EXONS+1),
                         header=None, 
                         names=['gene'] + \
                               [str(x) for x in range(1, NUM_EXONS+1)],
                         engine='python')

# remove genes with < NUM_EXONS exons
mean_cov = mean_cov.dropna()

# colours for plot: create empty dict
colours = {}

unmeth_mean_cov = mean_cov[mean_cov['gene'].isin(unmeth_genes)]
n = len(unmeth_mean_cov)
unmeth_mean_cov = unmeth_mean_cov.drop(['gene'], axis=1)
um_label = 'Unmethylated (n={:,})'.format(n)
unmeth_mean_cov = pd.melt(unmeth_mean_cov, value_name='relative_expr', var_name='exon')
unmeth_mean_cov['meth_status'] = um_label
colours[um_label] = '#fee0d2'

meth_mean_cov = mean_cov[mean_cov['gene'].isin(meth_genes)]
n = len(meth_mean_cov)
meth_mean_cov = meth_mean_cov.drop(['gene'], axis=1)
m_label = 'Methylated (n={:,})'.format(n)
meth_mean_cov = pd.melt(meth_mean_cov, value_name='relative_expr', var_name='exon')
meth_mean_cov['meth_status'] = m_label
colours[m_label] = '#fc9272'

highly_meth_mean_cov = mean_cov[mean_cov['gene'].isin(highly_meth_genes)]
n = len(highly_meth_mean_cov)
highly_meth_mean_cov = highly_meth_mean_cov.drop(['gene'], axis=1)
hm_label = 'Highly methylated (n={:,})'.format(n)
highly_meth_mean_cov = pd.melt(highly_meth_mean_cov, value_name='relative_expr', var_name='exon')
highly_meth_mean_cov['meth_status'] = hm_label
colours[hm_label] = '#de2d26'

mean_cov = pd.concat(
    [unmeth_mean_cov, meth_mean_cov, highly_meth_mean_cov], ignore_index=True)

# seaborn
sns.set_style('white')
sns.set_style('ticks')
fig, ax= plt.subplots(figsize=(6, 4))

sns.pointplot(x='exon', y='relative_expr', hue='meth_status', label='',
              data=mean_cov, palette=colours,
              ci=68, n_boot=10000, dodge=True)

ax.set_xlabel('Exon')
ax.set_ylabel('ln (relative expression vs. exon 1)')
ax.legend(loc='lower right')
sns.despine(offset=5, trim=True)

# save figure
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'spurious.meth_vs_unmeth.pdf'
fig.savefig(output_filename, bbox_inches='tight')

# perform statistical tests (t-tests)
for n in range(1, NUM_EXONS+1):
    um = unmeth_mean_cov[unmeth_mean_cov['exon'] == str(n)]['relative_expr']
    m = meth_mean_cov[meth_mean_cov['exon'] == str(n)]['relative_expr']
    hm = highly_meth_mean_cov[highly_meth_mean_cov['exon'] == str(n)]['relative_expr']
    
    print ('Exon', n)
    print ('t test p for um vs. m: ', scipy.stats.ttest_ind(um, m))
    print ('t test p for um vs. hm: ', scipy.stats.ttest_ind(um, hm))
    print ()
