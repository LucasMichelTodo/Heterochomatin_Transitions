#### Calculate correlations between ChIP Coverage ####

import os
import pybedtools as pb
import subprocess as sp
import pandas as pd
from itertools import combinations
from scipy.stats.stats import pearsonr

#### FUNCTIONS ####

def get_bdg_cov(bdg_file, cov_col):
    bed = pb.BedTool(bdg_file)
    ## Get coverage (from cov_col)
    cov = [float(feat.fields[cov_col]) for feat in bed]
    return(cov)

def make_union_bed(indir, bed1, bed2, out_fld):
    if not out_fld.endswith('/'):
        out_fld = out_fld+'/'
    b1 = indir+bed1
    b2 = indir+bed2
    id1 = bed1.rsplit('.', 1)[0]
    id2 = bed2.rsplit('.', 1)[0]
    outfile = out_fld+f'union_{id1}_{id2}.bdg'
    cmd = ['bedtools', 'unionbedg', '-i', b1, b2, '>', outfile]
    sp.call(' '.join(cmd), shell = True)

#### CALLS ####

wd = '/mnt/Disc4T/Projects/PhD_Project/Paper/Paper_Analysis/ChIP_Seq_Correlations/'
os.chdir(wd)

bdg_fld = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/RPKMs_NormInput_noDup_bs10_smth_200_pseudo_10/'

me_fls = [
    '1.2B_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
    '10G_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
    'A7K9_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
    'B11_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
    'E5K9_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg'
    ]

ac_fls = [
    '1.2B_ac_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
    '10G_ac_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
    'B11_ac_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
]

## We first create "union beds" so that coverage is calculated for the same exact bins in both samples.
union_dir = './Union_Beds_NormInput/'
os.makedirs(union_dir, exist_ok=True)

for pair in combinations(me_fls, 2):
    print(pair)
    make_union_bed(bdg_fld, pair[0], pair[1], union_dir)

for pair in combinations(ac_fls, 2):
    print(pair)
    make_union_bed(bdg_fld, pair[0], pair[1], union_dir)

union_beds_me = sorted([b for b in os.listdir(union_dir) if '_me_' in b])
union_beds_ac = sorted([b for b in os.listdir(union_dir) if '_ac_' in b])

## Now we can calculate correlations for pairs of samples.

## Me
corrs = {}
for u_bed in union_beds_me:
    name = u_bed.replace('union_', '').replace('_me_sort_q5_RPKMs_normInput', '').replace('.bdg', '')
    cov1 = get_bdg_cov(union_dir+u_bed, 3)
    cov2 = get_bdg_cov(union_dir+u_bed, 4)
    corr = pearsonr(cov1,cov2)
    corrs[name] = corr

## Create Correlation Mtx
corr_with = [
    '1.2B',
    '10G',
    'A7',
    'E5',
    'B11',
]

cols_dict = {}
for track in corr_with:
    cols_dict[track] = {key:None for key in corr_with}
    for k in cols_dict[track].keys():
        if k == track:
            cols_dict[track][k] = 1
        else:
            for kk, v in corrs.items():
                if (track in kk and k in kk):
                    cols_dict[track][k] = v[0]

corr_df = pd.DataFrame.from_dict(cols_dict)
corr_df.to_csv('norm_by_input_corrs.csv')

## Ac
corrs = {}
for u_bed in union_beds_ac:
    name = u_bed.replace('union_', '').replace('_me_sort_q5_RPKMs_normInput', '').replace('.bdg', '')
    cov1 = get_bdg_cov(union_dir+u_bed, 3)
    cov2 = get_bdg_cov(union_dir+u_bed, 4)
    corr = pearsonr(cov1,cov2)
    corrs[name] = corr

## Create Correlation Mtx
corr_with = [
    '1.2B',
    '10G',
    'B11',
]

cols_dict = {}
for track in corr_with:
    cols_dict[track] = {key:None for key in corr_with}
    for k in cols_dict[track].keys():
        if k == track:
            cols_dict[track][k] = 1
        else:
            for kk, v in corrs.items():
                if (track in kk and k in kk):
                    cols_dict[track][k] = v[0]

corr_df = pd.DataFrame.from_dict(cols_dict)
corr_df.to_csv('norm_by_input_ac_corrs.csv')


## Met vs Ac
bdg_fld = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/RPKMs_NormInput_noDup_bs10_smth_200_pseudo_10/'

me_fls = [
    '1.2B_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
    '10G_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
    'B11_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
    ]

ac_fls = [
    '1.2B_ac_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
    '10G_ac_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
    'B11_ac_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
]

union_dir = './Union_Beds_NormInput/'
os.makedirs(union_dir, exist_ok=True)

for pair in zip(me_fls, ac_fls):
    print(pair)
    make_union_bed(bdg_fld, pair[0], pair[1], union_dir)

union_beds_me_ac = sorted([b for b in os.listdir(union_dir) if '_me_' in b and '_ac_' in b])

## Me
corrs = {}
for u_bed in union_beds_me_ac:
    name = u_bed.replace('union_', '').replace('_me_sort_q5_RPKMs_normInput', '').replace('.bdg', '')
    cov1 = get_bdg_cov(union_dir+u_bed, 3)
    cov2 = get_bdg_cov(union_dir+u_bed, 4)
    corr = pearsonr(cov1,cov2)
    corrs[name] = corr

pd.DataFrame(corrs).T.to_csv('norm_by_input_me_ac_corrs.csv')
