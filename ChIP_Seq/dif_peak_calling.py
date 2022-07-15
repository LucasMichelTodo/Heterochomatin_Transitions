import time
import os
import pybedtools as pb
import numpy as np
import subprocess as sp
from tqdm import tqdm
from scipy.stats import skewnorm
import matplotlib.pyplot as plt
import pandas as pd
from multiprocessing import Pool
from itertools import repeat

## Functions

def create_window_ref(window_size, genome_file, step):
    bed = pb.BedTool()
    win_ref = bed.window_maker(w=window_size, g = genome_file, s = step)
    return(win_ref)

def plot_skewnorm(cov_vector, params, plotname):

    # Plot the histogram.
    plt.hist(cov_vector, bins=25, density=True, alpha=0.6, color='g')

    # Plot the PDF.
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = skewnorm.pdf(x, *params)
    plt.plot(x, p, 'k', linewidth=2)
    plt.savefig(plotname)
    plt.clf()

def get_probs(coverage, params):
    return(skewnorm.cdf(coverage, *params))

def compare_cov(feat, skewnorm_params1, skewnorm_params2):
    c1 = float(feat.fields[3])
    c2 = float(feat.fields[4])
    prob1 = skewnorm.cdf(c1, *skewnorm_params1)
    prob2 = skewnorm.cdf(c2, *skewnorm_params2)
    probdif = prob1 - prob2
    return(probdif)

def get_differential_peaks(
        peakfile1, peakfile2,
        covfile1, covfile2,
        prefix1, prefix2,
        genome_file,
        window_size, stepsize,
        minprobdif,
        mergedist,
        minlen,
        outfld):

    start = time.time()

    ## Create folders for output
    os.makedirs(outfld, exist_ok = True)
    if not outfld.endswith('/'): outfld = outfld+'/'
    os.makedirs(outfld+'Data/', exist_ok = True)
    os.makedirs(outfld+'Data/Common_Coverage/', exist_ok = True)
    os.makedirs(outfld+'Data/Common_Peaks/', exist_ok = True)
    os.makedirs(outfld+'Data/Plots/', exist_ok = True)
    os.makedirs(outfld+'Data/Window_Coverage/', exist_ok = True)
    os.makedirs(outfld+'Data/PreFilter_Difpeaks/', exist_ok = True)


    ## Generate skewnorm distribution model for coverage in peaks
    ## Peaks coverage
    print('Fitting distribution...')

    pf1 = pb.BedTool(peakfile1)
    pf2 = pb.BedTool(peakfile2)
    cf1 = pb.BedTool(covfile1)
    cf2 = pb.BedTool(covfile2)

    cov_peaks1 = np.array([float(f.fields[10]) for f in pf1.map(cf1, c = 4, o='mean')])
    cov_peaks2 = np.array([float(f.fields[10]) for f in pf2.map(cf2, c = 4, o='mean')])

    skewnorm_params1 = skewnorm.fit(cov_peaks1)
    skewnorm_params2 = skewnorm.fit(cov_peaks2)

    suffix = '_peak_coverage_fit.png'
    plot_skewnorm(cov_peaks1, skewnorm_params1, f'{outfld}Data/Plots/{prefix1}{suffix}')
    plot_skewnorm(cov_peaks2, skewnorm_params2, f'{outfld}Data/Plots/{prefix2}{suffix}')

    ## Create window reference
    print('Creating windowed coverage...')
    suffix = '_window_coverage.bdg'
    win_ref = create_window_ref(window_size, genome_file, stepsize)
    win_cov1 = win_ref.map(cf1, c=4, o='mean').saveas(f'{outfld}Data/Window_Coverage/{prefix1}{suffix}')
    win_cov2 = win_ref.map(cf2, c=4, o='mean').saveas(f'{outfld}Data/Window_Coverage/{prefix2}{suffix}')

    ## Create join coverage bdg
    print('Creating join coverage bed...')
    outfile = outfld+f'Data/Common_Coverage/{prefix1}_{prefix2}_common_coverage.bdg'
    cmd = ['bedtools', 'unionbedg', '-i',
           f'{outfld}Data/Window_Coverage/{prefix1}{suffix}',
           f'{outfld}Data/Window_Coverage/{prefix2}{suffix}',
           '>', outfile]
    sp.call(' '.join(cmd), shell = True)
    union_bed = pb.BedTool(outfile)

    ## Create common peaks bed
    print('Creating common peaks bed...')

    files_to_cross = []
    str_beds = peakfile1 + ' ' + peakfile2
    outfile = f'{outfld}Data/Common_Peaks/{prefix1}_{prefix2}_common_peaks.bed'
    cmd = 'awk \'{print}\' '+f'{str_beds} > {outfile}'
    sp.call(cmd, shell = True)
    common_peaks = pb.BedTool(outfile).sort()

    ## Subset union_bed to common_peaks
    common_cov_common_peaks = union_bed.intersect(common_peaks)

    ## Call differential peaks
    print('Calling differential peaks (this might take a while)...')

    c1s = [float(x.fields[3]) for x in common_cov_common_peaks]
    c2s = [float(x.fields[4]) for x in common_cov_common_peaks]

    with Pool(8) as p:
        prob1s = p.starmap(get_probs, zip(c1s, repeat(skewnorm_params1)))
    with Pool(8) as p:
        prob2s = p.starmap(get_probs, zip(c2s, repeat(skewnorm_params2)))

    prob1s = np.array(prob1s)
    prob2s = np.array(prob2s)

    probdifs = prob1s - prob2s
    peaks1over2 = probdifs > minprobdif
    peaks2over1 = -probdifs > minprobdif
    chroms = [x.chrom for x in common_cov_common_peaks]
    starts = [x.start for x in common_cov_common_peaks]
    stops = [x.stop for x in common_cov_common_peaks]

    pd.DataFrame(common_cov_common_peaks)
    df = pd.read_table(
        common_cov_common_peaks.fn,
        names=['#chrom', 'start', 'stop', 'cov1', 'cov2']
               )
    rawout1 = f'{outfld}Data/PreFilter_Difpeaks/difpeaks_{prefix1}over{prefix2}_w{window_size}_s{stepsize}_pd{minprobdif}.bed'
    rawout2 = f'{outfld}Data/PreFilter_Difpeaks/difpeaks_{prefix2}over{prefix1}_w{window_size}_s{stepsize}_pd{minprobdif}.bed'

    df[peaks1over2].to_csv(rawout1, sep='\t', index=False)
    df[peaks2over1].to_csv(rawout2, sep='\t', index=False)

    rawbed1 = pb.BedTool(rawout1)
    rawbed2 = pb.BedTool(rawout2)

    ## Merge peaks and filter by length
    out1 = f'{outfld}difpeaks_{prefix1}over{prefix2}_w{window_size}_s{stepsize}_pd{minprobdif}_mg{mergedist}_ml{minlen}.bed'
    out2 = f'{outfld}difpeaks_{prefix2}over{prefix1}_w{window_size}_s{stepsize}_pd{minprobdif}_mg{mergedist}_ml{minlen}.bed'
    out1 = rawbed1.sort().merge(d=mergedist).filter(lambda f: f.stop - f.start > minlen).saveas(out1)
    out2 = rawbed2.sort().merge(d=mergedist).filter(lambda f: f.stop - f.start > minlen).saveas(out2)

    end = time.time()
    print('Finished! Elapsed time:')
    print(end - start)

