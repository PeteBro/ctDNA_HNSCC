# import and define variables + functions

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

for subject in subjects:
    cov = pd.DataFrame(lq.get_coverage(data[subject]))
    cov.to_csv(f'outputs/tissue_data/{subject}_coverage.tsv', sep = '\t', index = False)

# find where variants occur in the tissue samples on a per-gene basis

locs = {subject : {} for subject in subjects}
for subject in subjects:
    for gene in np.unique(data[subject]['GENE']):
        tempt = data[subject].iloc[np.where(data[subject]['GENE'] == f'{gene}')[0]]
        y = np.zeros([2,len(tempt)])
        y[0] = intv(np.array([i.split(':') for i in tempt[f'{subject}_Core']]).T[5])
        y[1] = intv(np.array([i.split(':') for i in tempt[f'{subject}_Margin']]).T[5])
        y = y.T
        locs[subject][gene] = classify_loc(y)

# calculate proportional occurence

props = {subject : {gene : get_props(locs[subject][gene]) for gene in locs[subject]} for subject in subjects}

het_measures = {region : {subject : {} for subject in subjects} for region in ['core', 'margin', 'both', 'rate']}
for region in ['core', 'margin', 'both', 'rate']:
    for subject in subjects:
        for gene in genes:
            if gene in props[subject]:
                het_measures[region][subject][gene] = props[subject][gene][region]
            else:
                het_measures[region][subject][gene] = np.NaN

# plot

titles = {'core' : 'Proportion of all variants observed exclusively in tumour centre',
          'margin' : 'Proportion of all variants observed exclusively in tumour peripheral tissue',
          'both' : 'Proportion of all variants observed in both peripheral and central tumour tissue',
          'rate' : 'Rate of intratumoral heterogeneity for all variants'}

for reg in ['core', 'margin', 'both', 'rate']:
   plt.ioff()
   plt.figure(figsize=(15,12))
   df = pd.DataFrame(het_measures[reg])
   sns.heatmap(df, annot = True, cmap = 'hot_r', vmin=0, vmax=1)
   plt.title(titles[reg])
   plt.savefig(f'outputs/tissue_data/{reg}.svg', format='svg')
   plt.savefig(f'outputs/tissue_data/{reg}')
   plt.clf()

# do the same for all variants regardless of gene

locs = {subject : {} for subject in subjects}
for subject in subjects:
        y = np.zeros([2,len(data[subject])])
        y[0] = intv(np.array([i.split(':') for i in data[subject][f'{subject}_Core']]).T[5])
        y[1] = intv(np.array([i.split(':') for i in data[subject][f'{subject}_Margin']]).T[5])
        y = y.T
        locs[subject] = classify_loc(y)

het_all = {subject : get_props(locs[subject]) for subject in subjects}

subjects = list(het_all.keys())
df = pd.DataFrame(het_all).T
fig, ax = plt.subplots(figsize=(12,8))
ax.set_ylim(0,1)
ax.bar([str(s) for s in subjects], np.array(df['core']),  label='Tumour Core Only')
ax.bar([str(s) for s in subjects], np.array(df['margin']), bottom=(np.array(df['core'])), label ='Tumour Margin Only')
ax.bar([str(s) for s in subjects], np.array(df['both']), bottom=(np.array(df['margin'])+np.array(df['core'])), label='Tumour Core and Margin')
plt.title('Average heterogeneity of all variants in tumour tissue samples')
fig.legend()
plt.savefig('outputs/tissue_data/mean_plot.svg', format='svg')
plt.savefig('outputs/tissue_data/mean_plot')


# repeat, this time examining only COSMIC annotated variants with known roles in HNSCC

data = {subject: pd.read_csv(f'outputs/tissue_data/{subject}_filtered_HNSCC.tsv', sep = '\t') for subject in subjects}
for subject in subjects:
    cov = pd.DataFrame(lq.get_coverage(data[subject]))
    cov.to_csv(f'outputs/tissue_data/{subject}_coverage_HNSCC.tsv', sep = '\t', index = False)

locs = {subject : {} for subject in subjects}
for subject in subjects:
    for gene in np.unique(data[subject]['GENE']):
        tempt = data[subject].iloc[np.where(data[subject]['GENE'] == f'{gene}')[0]]
        y = np.zeros([2,len(tempt)])
        y[0] = intv(np.array([i.split(':') for i in tempt[f'{subject}_Core']]).T[5])
        y[1] = intv(np.array([i.split(':') for i in tempt[f'{subject}_Margin']]).T[5])
        y = y.T
        locs[subject][gene] = classify_loc(y)

props = {subject : {gene : get_props(locs[subject][gene]) for gene in locs[subject]} for subject in subjects}

het_measures = {region : {subject : {} for subject in subjects} for region in ['core', 'margin', 'both', 'rate']}
for region in ['core', 'margin', 'both', 'rate']:
    for subject in subjects:
        for gene in genes:
            if gene in props[subject]:
                het_measures[region][subject][gene] = props[subject][gene][region]
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

locs = {subject : {} for subject in subjects}
for subject in subjects:
        y = np.zeros([2,len(data[subject])])
        y[0] = intv(np.array([i.split(':') for i in data[subject][f'{subject}_Core']]).T[5])
        y[1] = intv(np.array([i.split(':') for i in data[subject][f'{subject}_Margin']]).T[5])
        y = y.T
        locs[subject] = classify_loc(y)

het_all = {subject : get_props(locs[subject]) for subject in subjects}

subjects = list(het_all.keys())
df = pd.DataFrame(het_all).T
fig, ax = plt.subplots(figsize=(12,8))
ax.set_ylim(0,1)
ax.bar([str(s) for s in subjects], np.array(df['core']),  label='Tumour Core Only')
ax.bar([str(s) for s in subjects], np.array(df['margin']), bottom=(np.array(df['core'])), label ='Tumour Margin Only')
ax.bar([str(s) for s in subjects], np.array(df['both']), bottom=(np.array(df['margin'])+np.array(df['core'])), label='Tumour Core and Margin')
plt.title('Average heterogeneity of HNSCC associated variants in tumour tissue samples')
fig.legend()
plt.savefig('outputs/tissue_data/mean_plot_HNSCC.svg', format='svg')
plt.savefig('outputs/tissue_data/mean_plot_HNSCC')
