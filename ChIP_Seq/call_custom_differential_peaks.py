from Custom_Differential_Peak_Calling.dif_peak_calling import *
from itertools import combinations

## Call differential on all samples
## Load MACS2 peaks file and normalized by input coverage

wd = '/mnt/Disc4T/Projects/Miniprojects/Custom_Differential_Peak_Calling/'
os.chdir(wd)

datadir = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/'

macs2_dir = datadir+'/Peak_Calling_MACS2/'
cov_dir = datadir+'/RPKMs_NormInput_noDup_bs10_smth_200_pseudo_10/'

macs2_fls = sorted([f for f in os.listdir(macs2_dir) if f.endswith('narrowPeak') and '_me_' in f])
cov_fls = sorted([f for f in os.listdir(cov_dir) if f.endswith('_pseudo10.bdg') and '_me_' in f])

macs2_fls = [
    '1.2B_me_Macspeaks_peaks.narrowPeak',
    '10G_me_Macspeaks_peaks.narrowPeak',
    'A7K9_me_Macspeaks_peaks.narrowPeak',
    'E5K9_me_Macspeaks_peaks.narrowPeak',
    'B11_me_Macspeaks_peaks.narrowPeak',
]

cov_fls = [
    '1.2B_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
    '10G_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
    'A7K9_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
    'E5K9_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
    'B11_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
]

genome = './Pf3D7.genome'
winsize = 100
stepsize = 100
minprobdif = 0.3
mergedist = 500
minlen = 1000
outfld = f'./DifPeaks_W{winsize}_S{stepsize}_PD{minprobdif}_Mg{mergedist}_Ml{minlen}/'

for c in zip(combinations(macs2_fls, 2), combinations(cov_fls, 2)):
    print(c[0])
    print(c[1])

    peakfile1 = macs2_dir+c[0][0]
    peakfile2 = macs2_dir+c[0][1]
    covfile1 = cov_dir+c[1][0]
    covfile2 = cov_dir+c[1][1]
    prefix1 = c[0][0].split('_')[0]
    prefix2 = c[0][1].split('_')[0]

    get_differential_peaks(
        peakfile1, peakfile2,
        covfile1, covfile2,
        prefix1, prefix2,
        genome, winsize, stepsize,
        minprobdif, mergedist, minlen,
        outfld
    )

#### Annotate differential peaks ####

import pybedtools as pb
import pandas as pd
import os

## Function

def annotate_bed(peaks_bed, ref_bed, ncols, gid_col):

    ref = pb.BedTool(ref_bed)
    ref = ref.sort()

    bed = pb.BedTool(peaks_bed)
    anot = bed.intersect(ref, wao=True)

    parsed_anot = []
    for interval in anot:

        originalfields = interval.fields[0:ncols]

        if interval.fields[ncols] == '.':
            gid = "intergenic"
        else:
            gid = interval.fields[ncols+gid_col-1]

        parsed_anot.append(originalfields + [gid])

    df = pd.DataFrame(parsed_anot)
    outfile = peaks_bed.replace(".bed", "_10005p_500orf_gene_crossed.tsv")
    df.to_csv(outfile, sep="\t", header=False, index=False)

## Calls

wd = '/mnt/Disc4T/Projects/PhD_Project/Paper/Paper_Analysis/Data_Files/'
os.chdir(wd)

gff = './Binned_Beds/binned_1000fp_500orf.bed'
peaksdir = './DifPeaks_W100_S100_PD0.3_Mg500_Ml1000/'

ncols = 3
gid_col = 4

filtered_peaks = [f for f in os.listdir(peaksdir) if f.endswith('.bed')]
for f in filtered_peaks:
    annotate_bed(
        peaks_bed = peaksdir+f,
        ref_bed = gff,
        ncols = ncols,
        gid_col = gid_col
    )
