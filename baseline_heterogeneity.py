# import and define variables + functions

import numpy as np
import pandas as pd
from math import log10, floor
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import seaborn as sns
import lq_biopsy as lq

subjects = [110, 133, 134, 141, 142, 150, 155, 165, 171]

genes = ['CASP8',
 'CDKN2A',
 'FAT1',
 'FBXW7',
 'KMT2D',
 'NOTCH1',
 'NSD1',
 'PIK3CA',
 'TP53']

def classify_loc(array):
    loc = []
    for c,m,a in array:
        if bool(c) and bool(m) and bool(a):
            loc.append('all')
        elif bool(c) and not(bool(m)) and bool(a):
            loc.append('core_and_blood')
        elif not(bool(c)) and bool(m) and bool(a):
            loc.append('margin_and_blood')
        elif not(bool(c)) and not(bool(m)) and bool(a):
            loc.append('blood')
        elif not(bool(c)) and not(bool(m)) and not(bool(a)):
            loc.append('none')

    return np.array(loc)



def gen_classify_loc(array):
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

    all = len(np.where(array == 'all')[0]) / len(array)
    blood = len(np.where(array == 'blood')[0]) / len(array)
    core_and_blood = len(np.where(array == 'core_and_blood')[0]) / len(array)
    margin_and_blood = len(np.where(array == 'margin_and_blood')[0]) / len(array)

    return {'all' : all, 'blood' : blood, 'core_and_blood' : core_and_blood, 'margin_and_blood' : margin_and_blood}



def gen_get_props(array):

    if len(array) == 0:
        core, margin, both, rate = np.NaN, np.NaN, np.NaN, np.NaN

    else:
        core = len(np.where(array == 'core')[0]) / len(array)
        margin = len(np.where(array == 'margin')[0]) / len(array)
        both = len(np.where(array == 'both')[0]) / len(array)
        rate = core+margin

    return {'core' : core, 'margin': margin, 'both' : both, 'rate' : rate}



intg = lambda x: bool(int(x))
intv = np.vectorize(intg)

# load tumour data

data = {subject: pd.read_csv(f'outputs/baseline_only/{subject}_filtered.tsv', sep = '\t') for subject in subjects}
for subject in subjects:
    cov = pd.DataFrame(lq.get_coverage(data[subject]))
    cov.to_csv(f'outputs/baseline_only/{subject}_coverage.tsv', sep = '\t', index = False)

# find where variants occur in the samples

locs = {subject : {} for subject in subjects}
for subject in subjects:
    y = np.zeros([3,len(data[subject])])
    y[0] = intv(np.array([i.split(':') for i in data[subject][f'{subject}_Core']]).T[5])
    y[1] = intv(np.array([i.split(':') for i in data[subject][f'{subject}_Margin']]).T[5])
    y[2] = intv(np.array([i.split(':') for i in data[subject][f'{subject}A']]).T[5])
    in_blood = list(np.where(y[2])[0])
    z = np.zeros([3,len(in_blood)])
    z[0] = y[0][in_blood]
    z[1] = y[1][in_blood]
    z[2] = y[2][in_blood]
    locs[subject] = classify_loc(z.T)

    # plot dispersion with venn diagram

    total = len(y[0])
    plt.figure(figsize=(10,10))
    plt.ioff()
    venn3([set(list(np.where(y[0])[0])), set(list(np.where(y[1])[0])), set(list(np.where(y[2])[0]))], ('Tumour Core', 'Tumour Margin', 'ctDNA'))
    plt.title(f'Patient {subject}')
    plt.savefig(f'outputs/baseline_only/sample_overlap_{subject}')
    plt.savefig(f'outputs/baseline_only/sample_overlap_{subject}.svg', format='svg')
    plt.clf()

# plot dispersion with composite bar

het_all = {subject : get_props(locs[subject]) for subject in subjects}

subjects = list(het_all.keys())
df = pd.DataFrame(het_all).T
fig, ax = plt.subplots(figsize=(12,8))
ax.set_ylim(0,1)
ax.bar([str(s) for s in subjects], np.array(df['blood']), label='ctDNA Only')
ax.bar([str(s) for s in subjects], np.array(df['core_and_blood']), bottom=np.array(df['blood']),  label='Tumor Core Only')
ax.bar([str(s) for s in subjects], np.array(df['margin_and_blood']), bottom=(np.array(df['core_and_blood'])+np.array(df['blood'])), label ='Tumor Margin Only')
ax.bar([str(s) for s in subjects], np.array(df['all']), bottom=(np.array(df['margin_and_blood'])+np.array(df['core_and_blood'])+np.array(df['blood'])), label='Tumor Core and Margin')
plt.title('Proportion of variants in ctDNA also observed in solid tumour samples')
fig.legend()
plt.savefig('outputs/baseline_only/mean_plot.svg', format='svg')
plt.savefig('outputs/baseline_only/mean_plot')

# find and plot where variants are in samples on a per gene basis

locs = {subject : {} for subject in subjects}
for subject in subjects:
    for gene in np.unique(data[subject]['GENE']):
        tempt = data[subject].iloc[np.where(data[subject]['GENE'] == f'{gene}')[0]]
        y = np.zeros([3,len(tempt)])
        y[0] = intv(np.array([i.split(':') for i in tempt[f'{subject}_Core']]).T[5])
        y[1] = intv(np.array([i.split(':') for i in tempt[f'{subject}_Margin']]).T[5])
        y[2] = intv(np.array([i.split(':') for i in tempt[f'{subject}A']]).T[5])

        present = list(set.union(set(np.where(y[0])[0]).intersection(set(np.where(y[2])[0])),
        set(np.where(y[1])[0]).intersection(set(np.where(y[2])[0]))))

        x = np.zeros([2,len(present)])
        x[0] = y[0][present]
        x[1] = y[1][present]

        locs[subject][gene] = gen_classify_loc(x.T)

for i in locs:
    for j in locs[i].values():
        print(len(j))

props = {subject : {gene : gen_get_props(locs[subject][gene]) for gene in locs[subject]} for subject in subjects}

het_measures = {region : {subject : {} for subject in subjects} for region in ['core', 'margin', 'both', 'rate']}
for region in ['core', 'margin', 'both', 'rate']:
    for subject in subjects:
        for gene in genes:
            if gene in props[subject]:
                het_measures[region][subject][gene] = props[subject][gene][region]
            else:
                het_measures[region][subject][gene] = np.NaN

titles = {'core' : 'Proportion of variants present in ctDNA also observed in tumour centre',
          'margin' : 'Proportion of variants present in ctDNA also observed in peripheral tumour tissue',
          'both' : 'Proportion of variants present in ctDNA also observed in both peripheral and central tumour tissue',
          'rate' : 'Rate of intratumoral heterogeneity for variants present in both ctDNA and tumour tissue'}

for reg in ['core', 'margin', 'both', 'rate']:
   plt.ioff()
   plt.figure(figsize=(15,12))
   df = pd.DataFrame(het_measures[reg])
   sns.heatmap(df, annot = True, cmap = 'hot_r', vmin=0, vmax=1)
   plt.title(titles[reg])
   df.to_csv(f'outputs/baseline_only/{reg}.tsv', sep = '\t')
   plt.savefig(f'outputs/baseline_only/{reg}.svg', format='svg')
   plt.savefig(f'outputs/baseline_only/{reg}')
   plt.clf()

# repeat with only HNSCC associated variants

data = {subject: pd.read_csv(f'outputs/baseline_only/{subject}_filtered_HNSCC.tsv', sep = '\t') for subject in subjects}
for subject in subjects:
    cov = pd.DataFrame(lq.get_coverage(data[subject]))
    cov.to_csv(f'outputs/baseline_only/{subject}_coverage_HNSCC.tsv', sep = '\t', index = False)

locs = {subject : {} for subject in subjects}
for subject in subjects:
    y = np.zeros([3,len(data[subject])])
    y[0] = intv(np.array([i.split(':') for i in data[subject][f'{subject}_Core']]).T[5])
    y[1] = intv(np.array([i.split(':') for i in data[subject][f'{subject}_Margin']]).T[5])
    y[2] = intv(np.array([i.split(':') for i in data[subject][f'{subject}A']]).T[5])
    in_blood = list(np.where(y[2])[0])
    z = np.zeros([3,len(in_blood)])
    z[0] = y[0][in_blood]
    z[1] = y[1][in_blood]
    z[2] = y[2][in_blood]
    locs[subject] = classify_loc(z.T)
    total = len(y[0])
    plt.figure(figsize=(10,10))
    plt.ioff()
    venn3([set(list(np.where(y[0])[0])), set(list(np.where(y[1])[0])), set(list(np.where(y[2])[0]))], ('Tumour Core', 'Tumour Margin', 'ctDNA'))
    plt.title(f'HNSCC associated variants for Patient {subject}')
    plt.savefig(f'outputs/baseline_only/sample_overlap_{subject}_HNSCC')
    plt.savefig(f'outputs/baseline_only/sample_overlap_{subject}_HNSCC.svg', format='svg')
    plt.clf()

het_all = {subject : get_props(locs[subject]) for subject in subjects}

subjects = list(het_all.keys())
df = pd.DataFrame(het_all).T
fig, ax = plt.subplots(figsize=(12,8))
ax.set_ylim(0,1)
ax.bar([str(s) for s in subjects], np.array(df['blood']), label='ctDNA Only')
ax.bar([str(s) for s in subjects], np.array(df['core_and_blood']), bottom=np.array(df['blood']),  label='Tumor Core Only')
ax.bar([str(s) for s in subjects], np.array(df['margin_and_blood']), bottom=(np.array(df['core_and_blood'])+np.array(df['blood'])), label ='Tumor Margin Only')
ax.bar([str(s) for s in subjects], np.array(df['all']), bottom=(np.array(df['margin_and_blood'])+np.array(df['core_and_blood'])+np.array(df['blood'])), label='Tumor Core and Margin')
plt.title('Proportion of HNSCC associated variants in ctDNA also observed in solid tumour samples')
fig.legend()
plt.savefig('outputs/baseline_only/mean_plot_HNSCC.svg', format='svg')
plt.savefig('outputs/baseline_only/mean_plot_HNSCC')

locs = {subject : {} for subject in subjects}
for subject in subjects:
    for gene in np.unique(data[subject]['GENE']):
        tempt = data[subject].iloc[np.where(data[subject]['GENE'] == f'{gene}')[0]]
        y = np.zeros([3,len(tempt)])
        y[0] = intv(np.array([i.split(':') for i in tempt[f'{subject}_Core']]).T[5])
        y[1] = intv(np.array([i.split(':') for i in tempt[f'{subject}_Margin']]).T[5])
        y[2] = intv(np.array([i.split(':') for i in tempt[f'{subject}A']]).T[5])

        present = list(set.union(set(np.where(y[0])[0]).intersection(set(np.where(y[2])[0])),
        set(np.where(y[1])[0]).intersection(set(np.where(y[2])[0]))))

        x = np.zeros([2,len(present)])
        x[0] = y[0][present]
        x[1] = y[1][present]

        locs[subject][gene] = gen_classify_loc(x.T)

for i in locs:
    for j in locs[i].values():
        print(len(j))

props = {subject : {gene : gen_get_props(locs[subject][gene]) for gene in locs[subject]} for subject in subjects}

het_measures = {region : {subject : {} for subject in subjects} for region in ['core', 'margin', 'both', 'rate']}
for region in ['core', 'margin', 'both', 'rate']:
    for subject in subjects:
        for gene in genes:
            if gene in props[subject]:
                het_measures[region][subject][gene] = props[subject][gene][region]
            else:
                het_measures[region][subject][gene] = np.NaN

titles = {'core' : 'Proportion of HNSCC associated variants present in ctDNA also observed in tumour centre',
          'margin' : 'Proportion of HNSCC associated variants present in ctDNA also observed in peripheral tumour tissue',
          'both' : 'Proportion of HNSCC associated variants present in ctDNA also observed in both peripheral and central tumour tissue',
          'rate' : 'Rate of intratumoral heterogeneity for HNSCC associated variants present in both ctDNA and tumour tissue'}

for reg in ['core', 'margin', 'both', 'rate']:
   plt.ioff()
   plt.figure(figsize=(15,12))
   df = pd.DataFrame(het_measures[reg])
   sns.heatmap(df, annot = True, cmap = 'hot_r', vmin=0, vmax=1)
   plt.title(titles[reg])
   df.to_csv(f'outputs/baseline_only/{reg}_HNSCC.tsv', sep = '\t')
   plt.savefig(f'outputs/baseline_only/{reg}_HNSCC.svg', format='svg')
   plt.savefig(f'outputs/baseline_only/{reg}_HNSCC')
   plt.clf()