# import + define variables

import pandas as pd
import lq_biopsy as lq

subjects = [110, 133, 134, 141, 142, 150, 155, 165, 171]

# load data

data = {subject : pd.read_csv(f'outputs/baseline_only/{subject}.tsv', sep='\t') for subject in subjects}

# drop variants with NaN values
for subject in subjects:
    data[subject] = lq.drop_nans(data[subject])

# get filtering formattable data

data_tree = {subject : lq.get_filtering_dict(data[subject]) for subject in subjects}

for subject in subjects:
    # drop variants with coverage < 10 in any sample
    data_tree[subject] = lq.filter_coverage(data_tree[subject], 10, [f'{subject}_Buffy'])
    # drop variants with no variant supporting reads in any sample
    data_tree[subject] = lq.filter_supporting_reads(data_tree[subject], 1, [f'{subject}_Buffy'])
    # find variants with any presence in the control "Buffy" sample
    data_tree[subject] = lq.identify_somatic_candidates(data_tree[subject], f'{subject}_Buffy')
    # recover variants with > 20 coverage and < 5% VAF
    # if these variants have >= 4 variant supporting reads in tumourogenic samples then restore them
    # additionaly, restore any variants with > 1 variant supporting reads in tumourogenic samples which are recognized to be associated with HNSCC in the COSMIC database
    data_tree[subject] = lq.recover_somatic_candidates(data_tree[subject], 20, 5, 4, f'{subject}_Buffy', id_true=True)
    # repeat for variants with < 20 coverage and < 10% VAF
    data_tree[subject] = lq.recover_somatic_candidates(data_tree[subject], 1, 10, 10, f'{subject}_Buffy', max_coverage=19, id_true=True)
    # remove any remaining variants recognized as somatic
    data_tree[subject] = lq.remove_somatic_candidates(data_tree[subject])
    # remove any variants with a VAF < 5% in all samples
    data_tree[subject] = lq.filter_VAF(data_tree[subject], 0.1, [f'{subject}_Buffy', 'somatic_candidate'])
    #data_tree[subject] = lq.filter_VAF_venns(data_tree[subject], 5, {'entry' : f'{subject}A', 'min' : 5}, [f'{subject}_Buffy', 'somatic_candidate'])
    # fit filtering parameters to data
    tempdata = lq.fit_filter(data[subject], data_tree[subject])
    # save filtered data
    tempdata.to_csv(f'outputs/baseline_only/{subject}_filtered.tsv', sep='\t', index=False)
    # keep only HNSCC associated variants recognized by COSMIC
    data_tree[subject] = lq.filter_annotation(data_tree[subject])
    # save HNSCC filtered variants
    tempdata = lq.fit_filter(data[subject], data_tree[subject])
    tempdata.to_csv(f'outputs/baseline_only/{subject}_filtered_HNSCC.tsv', sep='\t', index=False)
