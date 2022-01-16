__all__ = ['read_vcf',
'get_gene_library',
'get_gene_col',
'get_coverage',
'drop_nans',
'get_filtering_dict',
'filter_coverage',
'filter_supporting_reads',
'filter_VAF',
'filter_annotation',
'identify_somatic_candidates',
'recover_somatic_candidates',
'remove_somatic_candidates',
'fit_filter']

import io
import numpy as np
import pandas as pd

def read_vcf(path, sep = '\t'):  

    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]

    wc = len(lines)
    data = pd.read_csv(io.StringIO(''.join(lines)),
                     dtype={'#CHROM': str,
                            'POS': int,
                            'ID': str,
                            'REF': str,
                            'ALT': str,
                            'QUAL': str,
                            'FILTER': str,
                            'INFO': str},
                            sep = sep).rename(columns={'#CHROM': 'CHROM'})
    print(f'{wc} lines loaded')
    return data



def get_gene_library(path='lq_biopsy/utilities/coords.bed', delimiter='\t'):

    with open(path, 'r') as genes:
        lines = [l[:-1] for l in genes]
        lines = [l.split(delimiter) for l in lines]
        exon_dict = {chr : {} for chr in np.unique(np.array(lines).T[0])}

    for exon in lines:
        if exon[3] not in exon_dict[exon[0]]:
            exon_dict[exon[0]][exon[3]] = []
        exon_dict[exon[0]][exon[3]].append(range(int(exon[1]), int(exon[2])+1))
    
    return exon_dict



def get_gene_col(data, library, colname = 'GENE'):
    chroms = list(data['CHROM'])
    pos = list(data['POS'])
    keys = []
    for i in range(len(pos)):
        det=0
        ik = 0
        while det == 0:
            key = list(library[chroms[i]].keys())[ik]
            det = 1 in [pos[i] in reg for reg in library[chroms[i]][key]]
            if det:
                keys.append(key)
            ik += 1
    data.insert(2,colname,keys,True)

    return data



def get_coverage(data):

    sep = lambda c: [l.split(':') for l in c]
    coverage_dict = {}
    for gene in np.unique(data['GENE']):
        tempt = data.iloc[np.where(data['GENE'] == f'{gene}')[0]]
        coverage_dict[gene] = {}
        for column in tempt.columns[10:]:
            matrix = np.stack(sep(tempt[column])).T
            coverage_dict[gene][column] = np.mean(np.array([int(val) for val in matrix[2]]))

    return coverage_dict




def drop_nans(data, nan_format = './.:.:', verbose = True):

    get_nan = lambda c: np.where(np.char.find(c, nan_format) != -1)[0]

    nans = []
    for column in list(data.columns[10:]):
        to_check = list(data[column])
        nan_idcs = get_nan(to_check)
        nans.append(list(nan_idcs))

    nans = list(set.union(*[set(na) for na in nans]))

    data = data.drop(nans)
    data = data.reset_index(drop = True)

    if verbose:
        print(f'{len(nans)} variants removed due to missingness')

    return data



def get_filtering_dict(data):

    sep = lambda c: [l.split(':') for l in c]
    filter_dict = {'pass' : np.array([True]*len(data))}

    cosmic_info = data['INFO']
    cosmic_info = [np.array(l.split(';'))[np.char.find(l.split(';'), 'cosmic') != -1] for l in cosmic_info]
    id_true = np.where([np.any((np.char.find(l, 'upper_aerodigestive_tract') != -1) == True) for l in cosmic_info])[0]
    filter_dict['annotated'] = np.array([False]*len(filter_dict['pass']))
    filter_dict['annotated'][id_true] = True

    for column in list(data.columns[10:]):
        matrix = np.stack(sep(data[column])).T

        filter_dict[column] = {'coverage' : np.array([int(val) for val in matrix[2]]),
                            'variant_supporting_reads' : np.array([int(val) for val in matrix[5]]),
                            'VAF' : np.array([float(val.strip('%')) for val in matrix[6]])}

    return filter_dict



def filter_coverage(filter_dict, min_coverage, ignore = None, fail_conditional = 'or', verbose = True):

    fail = []
    for key in [i for i in filter_dict.keys() if not(i in ['pass', 'annotated', *[i for i in ignore]])]:

        fail.append(list(np.where(filter_dict[key]['coverage'] < min_coverage)[0]))

    if fail_conditional == 'and':
        fail = list(set.intersection(*[set(f) for f in fail]))
    elif fail_conditional == 'or':
        fail = list(set.union(*[set(f) for f in fail]))

    filter_dict['pass'][fail] = False

    if verbose:
        print(f'{len(fail)} variants removed due to poor coverage')

    return filter_dict



def filter_supporting_reads(filter_dict, min_supporting_reads, ignore = None, fail_conditional = 'and', verbose = True):

    fail = []
    for key in [i for i in filter_dict.keys() if not(i in ['pass', 'annotated', *[i for i in ignore]])]:

        fail.append(list(np.where(filter_dict[key]['variant_supporting_reads'] < min_supporting_reads)[0]))

    if fail_conditional == 'and':
        fail = list(set.intersection(*[set(f) for f in fail]))
    elif fail_conditional == 'or':
        fail = list(set.union(*[set(f) for f in fail]))

    filter_dict['pass'][fail] = False

    if verbose:
        print(f'{len(fail)} variants removed due to lack of supporting read depth')

    return filter_dict



def filter_VAF(filter_dict, min_freq, ignore = None, fail_conditional = 'and', verbose = True):
    
    fail = []
    for key in [i for i in filter_dict.keys() if not(i in ['pass', 'annotated', *[i for i in ignore]])]:

        fail.append(list(np.where(filter_dict[key]['VAF'] < min_freq)[0]))

    if fail_conditional == 'and':
        fail = list(set.intersection(*[set(f) for f in fail]))
    elif fail_conditional == 'or':
        fail = list(set.union(*[set(f) for f in fail]))

    filter_dict['pass'][fail] = False

    if verbose:
        print(f'{len(fail)} variants removed due to low VAF')

    return filter_dict



#def filter_VAF_temp(filter_dict, min_freq, exceptions, ignore = None, fail_conditional = 'and', verbose = True):
#    
#    fail = []
#    for key in [i for i in filter_dict.keys() if not(i in ['pass', 'annotated', *[i for i in ignore]])]:
#        if key != exceptions['entry']:
#
#            fail.append(list(np.where(filter_dict[key]['VAF'] < min_freq)[0]))
#
#    fail.append(list(np.where(filter_dict[exceptions['entry']]['VAF'] < exceptions['min'])[0]))
#
#    if fail_conditional == 'and':
#        fail = list(set.intersection(*[set(f) for f in fail]))
#    elif fail_conditional == 'or':
#        fail = list(set.union(*[set(f) for f in fail]))
#
#    filter_dict['pass'][fail] = False
#
#    if verbose:
#        print(f'{len(fail)} variants removed due to low VAF')
#
#    return filter_dict



def filter_annotation(filter_dict, verbose = True):

    fail = np.where(~filter_dict['annotated'])[0]
    filter_dict['pass'][fail] = False

    if verbose:
        print(f'Identified {len(np.where(fail)[0])} unannotated variants')

    return filter_dict




def identify_somatic_candidates(filter_dict, control, verbose = True):

    candidates = list(np.where(filter_dict[control]['VAF'] > 0)[0])

    filter_dict['somatic_candidate'] = np.array([False]*len(filter_dict[control]['VAF']))
    filter_dict['somatic_candidate'][candidates] = True

    if verbose:
        print(f'Identified {len(candidates)} somatic candidates')

    return filter_dict



def recover_somatic_candidates(filter_dict, min_coverage, max_frequency, min_supporting_reads, control, max_coverage = None, id_true = False, verbose = True):

    if max_coverage:
        passes_coverage = list(set.intersection(set(list(np.where(min_coverage <= filter_dict[control]['coverage'])[0])), set(list(np.where(filter_dict[control]['coverage'] <= max_coverage)[0]))))
    else:
        passes_coverage = list(np.where(min_coverage <= filter_dict[control]['coverage'])[0])
    
    passes_frenquency = (list(np.where(filter_dict[control]['VAF'] <= max_frequency)[0]))

    passes_supporting_reads = []
    passes_id = []
    for key in [i for i in filter_dict.keys() if not(i in ['pass', 'annotated', 'somatic_candidate', control])]:
        passes_supporting_reads.append(list(np.where(filter_dict[key]['variant_supporting_reads'] >= min_supporting_reads)[0]))

        if id_true:
            pid = []
            pid.append(list(np.where(filter_dict['annotated'])[0]))
            pid.append(list(np.where(filter_dict[key]['variant_supporting_reads'] > 0)[0]))
            pid = list(set.intersection(*[set(p) for p in pid]))
            passes_id.append(pid)


    passes_supporting_reads = set.union(*[set(p) for p in passes_supporting_reads])
    passes_id = set.union(*[set(p) for p in passes_id])

    pass1 = set.intersection(set(passes_coverage), set(passes_frenquency), passes_supporting_reads)
    pass2 = set.intersection(set(passes_coverage), set(passes_frenquency), passes_id)

    recover = list(set.union(pass1, pass2))

    filter_dict['somatic_candidate'][recover] = False

    if verbose:
        print(f'Recovered {len(recover)} variants')

    return filter_dict



def remove_somatic_candidates(filter_dict):

    filt = list(np.where(filter_dict['somatic_candidate'])[0])
    filter_dict['pass'][filt] = False

    return filter_dict



def fit_filter(data, filter_dict):

    filt = filter_dict['pass']
    data = data[filt]
    data = data.drop(data.columns[np.where(['Buffy' in s for s in data.columns])][0], axis = 1)

    return data