#!/usr/bin/env python3

"""
> plot_median_meth_timeseries.py <

Plots a time series chart for genes that positively and negatively regulate 
JNK and MAPK.
"""
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# read in delta median methylation from the xlsx worksheet
data = pd.read_excel('spis.bias_density_methylation.xlsx', 
                     skiprows=3, parse_cols='A,K:N,P')

# subselect genes that are differentially methylated at pH 7.2
goi = []
for line in open('reinv.7.2v8.1.5cpg.0.05.txt'):
    if not line: continue
    
    goi.append(line.strip())
    
dmg = data[data['Gene'].isin(goi)]

# subselect based on a GO term of interest
# 'GO:0043508': -ve regulation of JUN kinase
# 'GO:0043507': +ve regulation of JUN kinase
# 'GO:0043407': -ve regulation of MAP kinase
# 'GO:0043406': +ve regulation of MAP kinase
go_annot = {'GO:0043407': 'negative regulation of MAP kinase',
            'GO:0043406': 'positive regulation of MAP kinase',
            'GO:0043508': 'negative regulation of JUN kinase',
            'GO:0043507': 'positive regulation of JUN kinase',
}

# create four sub-dataframes for these GO terms
sub_df = {}
for g in go_annot:
    temp = dmg[dmg['GO terms'].str.contains(g, na=False)]
    temp = temp.filter(regex='^pH 7', axis=1)
    
    sub_df[g] = pd.DataFrame({
        'delta_meth': list(temp['pH 7.2']) + list(temp['pH 7.4']) + \
                      list(temp['pH 7.8']),
        'pH': [7.2] * len(temp['pH 7.2']) + \
              [7.4] * len(temp['pH 7.4']) + \
              [7.8] * len(temp['pH 7.8']),
        'subject': [x for x in range(len(temp['pH 7.2']))] + \
                   [x for x in range(len(temp['pH 7.4']))] + 
                   [x for x in range(len(temp['pH 7.8']))]
        })
    sub_df[g]['GO annotation'] = g + ': ' + go_annot[g]

df1 = pd.concat([sub_df['GO:0043508'], sub_df['GO:0043507']], ignore_index=True)
df2 = pd.concat([sub_df['GO:0043407'], sub_df['GO:0043406']], ignore_index=True)

# seaborn stuff
sns.set_style('white')
#sns.set_style('ticks')
sns.set_context('paper')

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4), sharey=True)

ts1 = sns.tsplot(data=df1, time='pH', value='delta_meth', unit='subject', 
                 condition='GO annotation', color=['#ef8a62', '#67a9cf'], ax=ax1)
ts2 = sns.tsplot(data=df2, time='pH', value='delta_meth', unit='subject', 
                 condition='GO annotation', color=['#ef8a62', '#67a9cf'], ax=ax2)

sns.despine(bottom=True)
plt.xlabel('pH')
ax1.set_ylabel(r'$\Delta$ median methylation level (%)')
ax2.set_ylabel('')

ax1.axhline(y=0, color='k')
ax2.axhline(y=0, color='k')

# redraws legend, which acts to remove the title within the legend
ax1.legend(loc='upper right', ncol=1)
ax2.legend(loc='upper right', ncol=1)

# save figure
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'median_meth_timeseries.pdf'
fig.savefig(output_filename, bbox_inches='tight')
