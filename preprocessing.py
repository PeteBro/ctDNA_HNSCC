# import + define variables

import os
import numpy as np
import lq_biopsy as lq

path = os.path.join('lq_biopsy', 'data', 'allsamples.vcf.hg38_multianno.vcf')
data = lq.read_vcf(path)
genelib = lq.get_gene_library()

# read in vcf file
# split subjects into dataframes

with open('lq_biopsy/utilities/bamlist.txt', 'r') as order:
    subject_order = [l.split('/')[-1].split('.')[0] for l in order]
    subject_order = {f'Sample{subject_order.index(n)+1}' : n for n in subject_order}
    data = data.rename(columns = subject_order)

# annotate gene column

data = lq.get_gene_col(data, genelib)

# retrieve relevant sample information for tissue data, baseline only data, and recurrence profile data

subjects = [110, 133, 134, 141, 142, 150, 155, 165, 171]
subject_data = {sub : data[[*list(data.columns[:10]), *list(data.columns[np.where(np.char.find(list(data.columns), str(sub)) == 0)[0]])]] for sub in subjects}
solid_data = {sub : subject_data[sub][list(subject_data[sub].columns[:13])] for sub in subjects}
baseline_data = {sub : subject_data[sub][list(subject_data[sub].columns[:14])] for sub in subjects}
recurrence_data = {sub : subject_data[sub].copy() for sub in subjects}

# save data

for subject in subjects:
    solid_data[subject].to_csv(f'outputs/tissue_data/{subject}.tsv', sep='\t', index=False)
    baseline_data[subject].to_csv(f'outputs/baseline_only/{subject}.tsv', sep='\t', index=False)
    recurrence_data[subject].to_csv(f'outputs/recurrence_profiles/{subject}.tsv', sep='\t', index=False)