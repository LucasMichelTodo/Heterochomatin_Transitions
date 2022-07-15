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

bdg_fld = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/RPKMs_noDup_bs10_smth_200/'
me_fls = [
    '1.2B_me_sort_q5_noDup_rpkm_normInput_bs10_smth200.bdg',
    '10G_me_sort_q5_noDup_rpkm_normInput_bs10_smth200.bdg',
    'A7K9_me_sort_q5_noDup_rpkm_normInput_bs10_smth200.bdg',
    'B11_me_sort_q5_noDup_rpkm_normInput_bs10_smth200.bdg',
    'E5K9_me_sort_q5_noDup_rpkm_normInput_bs10_smth200.bdg',
]
ac_fls = [
    '1.2B_ac_sort_q5_noDup_rpkm_normInput_bs10_smth200.bdg',
    '10G_ac_sort_q5_noDup_rpkm_normInput_bs10_smth200.bdg',
    'B11_ac_sort_q5_noDup_rpkm_normInput_bs10_smth200.bdg',
]
in_files = [
    '1.2B_in_sort_q5_noDup_rpkm_normInput_bs10_smth200.bdg',
    '10G_in_sort_q5_noDup_rpkm_normInput_bs10_smth200.bdg',
    'A7K9_in_sort_q5_noDup_rpkm_normInput_bs10_smth200.bdg',
    'B11_in_sort_q5_noDup_rpkm_normInput_bs10_smth200.bdg',
    'E5K9_in_sort_q5_noDup_rpkm_normInput_bs10_smth200.bdg',
]

## We first create "union beds" so that coverage is calculated for the same exact bins in both samples.

union_dir = './Union_Beds/'
os.makedirs(union_dir, exist_ok=True)

for pair in combinations(me_fls+ac_fls+in_files, 2):
    print(pair)
    make_union_bed(bdg_fld, pair[0], pair[1], union_dir)

union_beds = sorted(os.listdir(union_dir))

## Now we can calculate correlations for pairs of samples.

corrs = {}
for u_bed in union_beds:
    name = u_bed.replace('union_', '').replace('.bdg', '')
    print(name)
    cov1 = get_bdg_cov(union_dir+u_bed, 3)
    cov2 = get_bdg_cov(union_dir+u_bed, 4)
    corr = pearsonr(cov1,cov2)
    corrs[name] = corr

## Create Correlation Mtx
corr_with = [
    '1.2B_ac', '1.2B_me', '1.2B_in',
    '10G_ac', '10G_me', '10G_in',
    'A7K9_me', 'A7K9_in',
    'E5K9_me', 'E5K9_in',
    'B11_ac', 'B11_me', 'B11_in'
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
corr_df.to_csv('met_ac_in_corrs.csv')
