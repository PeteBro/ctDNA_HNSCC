# import and define variables + functions

import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import lq_biopsy as lq

subjects = [110, 133, 134, 141, 142, 150, 155, 165, 171]

genes = ['AJUBA',
 'CASP8',
 'CDKN2A',
 'CUL3',
 'FAT1',
 'FBXW7',
 'HLA.A',
 'KMT2D',
 'NFE2L2',
 'NOTCH1',
 'NSD1',
 'PIK3CA',
 'PIK3R1',
 'PTEN',
 'RB1',
 'TGFBR2',
 'TP53',
 'TRAF3',
 'HRAS']

def classify_loc(array):
    loc = []
    for c,m in array:
        if bool(c) and bool(m):
            loc.append('both')
        elif bool(c) and not(bool(m)):
            loc.append('core')
        elif not(bool(c)) and bool(m):
            loc.append('margin')

    return np.array(loc)



def get_props(array):

    core = len(np.where(array == 'core')[0]) / len(array)
    margin = len(np.where(array == 'margin')[0]) / len(array)
    both = len(np.where(array == 'both')[0]) / len(array)
    rate = core+margin

    return {'core' : core, 'margin': margin, 'both' : both, 'rate' : rate}



intg = lambda x: bool(int(x))
intv = np.vectorize(intg)

# load tumour data

data = {subject: pd.read_csv(f'outputs/tissue_data/{subject}_filtered.tsv', sep = '\t') for subject in subjects}

# get variant locations

propdict = {}
per_gene = {subject : {} for subject in subjects}
total_het = {subject : {} for subject in subjects}
for subject in subjects:
    y = np.zeros([2,len(data[subject])])
    y[0] = intv(np.array([i.split(':') for i in data[subject][f'{subject}_Core']]).T[5])
    y[1] = intv(np.array([i.split(':') for i in data[subject][f'{subject}_Margin']]).T[5])
    locgen = np.array(data[subject]['GENE'])
    f = classify_loc(y.T)
    for gene in locgen:
        idx = np.where(locgen==gene)[0]
        per_gene[subject][gene] = f[idx]

# get proportional heterogeneity

    propdict[subject] = {gene : get_props(per_gene[subject][gene]) for gene in per_gene[subject]}
    x = np.mean([propdict[subject][i]['core'] for i in propdict[subject]])
    y = np.mean([propdict[subject][i]['margin'] for i in propdict[subject]])
    z = np.mean([propdict[subject][i]['both'] for i in propdict[subject]])
    l = np.mean([propdict[subject][i]['rate'] for i in propdict[subject]])
    total_het[subject]['core'] = x
    total_het[subject]['margin'] = y
    total_het[subject]['both'] = z
    total_het[subject]['rate'] = l

df = pd.DataFrame(total_het).T
fig, ax = plt.subplots(figsize=(12,8))
ax.set_ylim(0,1)
ax.bar([str(s) for s in subjects], np.array(df['core']),  label='Tumour Core Only')
ax.bar([str(s) for s in subjects], np.array(df['margin']), bottom=(np.array(df['core'])), label ='Tumour Margin Only')
ax.bar([str(s) for s in subjects], np.array(df['both']), bottom=(np.array(df['margin'])+np.array(df['core'])), label='Tumour Core and Margin')
plt.title('Average heterogeneity variants in tumour tissue samples')
fig.legend()
plt.savefig('outputs/tissue_data/mean_plot.svg', format='svg')
plt.savefig('outputs/tissue_data/mean_plot')
pd.DataFrame(total_het).to_csv('outputs/tissue_data/mean_plot.tsv', sep='\t')
het_measures = {region : {subject : {} for subject in subjects} for region in ['core', 'margin', 'both', 'rate']}
for region in ['core', 'margin', 'both', 'rate']:
    for subject in subjects:
        for gene in genes:
            if gene in propdict[subject]:
                het_measures[region][subject][gene] = propdict[subject][gene][region]
            else:
                het_measures[region][subject][gene] = np.NaN

titles = {'core' : 'Proportion of variants observed exclusively in tumour centre',
          'margin' : 'Proportion of variants observed exclusively in tumour peripheral tissue',
          'both' : 'Proportion of variants observed in both peripheral and central tumour tissue',
          'rate' : 'Rate of intratumoral heterogeneity for variants'}

for reg in ['core', 'margin', 'both', 'rate']:
   plt.ioff()
   plt.figure(figsize=(15,12))
   df = pd.DataFrame(het_measures[reg])
   sns.heatmap(df, annot = True, cmap = 'hot_r', vmin=0, vmax=1)
   plt.title(titles[reg])
   plt.savefig(f'outputs/tissue_data/{reg}.svg', format='svg')
   plt.savefig(f'outputs/tissue_data/{reg}')
   plt.clf()

# repeat for HNSCC associated variants

data = {subject: pd.read_csv(f'outputs/tissue_data/{subject}_filtered_HNSCC.tsv', sep = '\t') for subject in subjects}

propdict = {}
per_gene = {subject : {} for subject in subjects}
total_het = {subject : {} for subject in subjects}
for subject in subjects:
    y = np.zeros([2,len(data[subject])])
    y[0] = intv(np.array([i.split(':') for i in data[subject][f'{subject}_Core']]).T[5])
    y[1] = intv(np.array([i.split(':') for i in data[subject][f'{subject}_Margin']]).T[5])
    locgen = np.array(data[subject]['GENE'])
    f = classify_loc(y.T)
    for gene in locgen:
        idx = np.where(locgen==gene)[0]
        per_gene[subject][gene] = f[idx]
    
    propdict[subject] = {gene : get_props(per_gene[subject][gene]) for gene in per_gene[subject]}
    x = np.mean([propdict[subject][i]['core'] for i in propdict[subject]])
    y = np.mean([propdict[subject][i]['margin'] for i in propdict[subject]])
    z = np.mean([propdict[subject][i]['both'] for i in propdict[subject]])
    l = np.mean([propdict[subject][i]['rate'] for i in propdict[subject]])
    total_het[subject]['core'] = x
    total_het[subject]['margin'] = y
    total_het[subject]['both'] = z
    total_het[subject]['rate'] = l

df = pd.DataFrame(total_het).T
fig, ax = plt.subplots(figsize=(12,8))
ax.set_ylim(0,1)
ax.bar([str(s) for s in subjects], np.array(df['core']),  label='Tumour Core Only')
ax.bar([str(s) for s in subjects], np.array(df['margin']), bottom=(np.array(df['core'])), label ='Tumour Margin Only')
ax.bar([str(s) for s in subjects], np.array(df['both']), bottom=(np.array(df['margin'])+np.array(df['core'])), label='Tumour Core and Margin')
plt.title('Average heterogeneity HNSCC associated variants in tumour tissue samples')
fig.legend()
plt.savefig('outputs/tissue_data/mean_plot_HNSCC.svg', format='svg')
plt.savefig('outputs/tissue_data/mean_plot_HNSCC')
pd.DataFrame(total_het).to_csv('outputs/tissue_data/mean_plot.tsv', sep='\t')
het_measures = {region : {subject : {} for subject in subjects} for region in ['core', 'margin', 'both', 'rate']}
for region in ['core', 'margin', 'both', 'rate']:
    for subject in subjects:
        for gene in genes:
            if gene in propdict[subject]:
                het_measures[region][subject][gene] = propdict[subject][gene][region]
            else:
                het_measures[region][subject][gene] = np.NaN


titles = {'core' : 'Proportion of HNSCC associated variants observed exclusively in tumour centre',
          'margin' : 'Proportion of HNSCC associated variants observed exclusively in tumour peripheral tissue',
          'both' : 'Proportion of HNSCC associated variants observed in both peripheral and central tumour tissue',
          'rate' : 'Rate of intratumoral heterogeneity for HNSCC associated variants'}

for reg in ['core', 'margin', 'both', 'rate']:
   plt.ioff()
   plt.figure(figsize=(15,12))
   df = pd.DataFrame(het_measures[reg])
   sns.heatmap(df, annot = True, cmap = 'hot_r', vmin=0, vmax=1)
   plt.title(titles[reg])
   plt.savefig(f'outputs/tissue_data/{reg}_HNSCC.svg', format='svg')
   plt.savefig(f'outputs/tissue_data/{reg}_HNSCC')
   plt.clf()
