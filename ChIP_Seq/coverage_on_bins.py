#### Calculate coverage over regions defined in a BED file ####

import pybedtools as pb
import os

wd = '/mnt/Disc4T/Projects/PhD_Project/Paper/Paper_Analysis/Data_Files/Binned_Beds/'
os.chdir(wd)

genes_bed_fnames = [f for f in os.listdir()]
genes_bed_fnames

## Me coverage
for gb in genes_bed_fnames:
    genes_bed = pb.BedTool(gb)

    cov_path = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/RPKMs_NormInput_noDup_bs10_smth_200_pseudo_10/'

    cov_tracks = [
        'E5K9_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
        'A7K9_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
        '1.2B_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
        '10G_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
        'B11_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg'
    ]

    for track in cov_tracks:
        coverage = pb.BedTool(cov_path+track)
        name = track.split('_')[0]
        print(f'Joining {name} coverage...')
        cov = genes_bed.sort().map(coverage, c = 4, o='mean')
        outname = gb.replace('.bed', '_coverage_')
        cov.saveas(f'../Coverages/{outname}{name}.bed')

## Ac coverage
genes_bed_fnames  = [
    'binned_1000tss_0orf_allowoverlaps_False.bed',
    'binned_0tss_500orf_allowoverlaps_True.bed',
    'binned_3prime1000_allowoverlaps_False.bed'
]

for gb in genes_bed_fnames:
    genes_bed = pb.BedTool(gb)

    cov_path = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/RPKMs_NormInput_noDup_bs10_smth_200_pseudo_10/'
    cov_tracks = [
        '1.2B_ac_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
        '10G_ac_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
        'B11_ac_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg'
    ]

    for track in cov_tracks:
        coverage = pb.BedTool(cov_path+track)
        name = track.split('_')[0]
        print(f'Joining {name} coverage...')
        cov = genes_bed.sort().map(coverage, c = 4, o='mean')
        outname = gb.replace('.bed', '_coverage_acetylation_')
        cov.saveas(f'../Coverages/{outname}{name}.bed')
