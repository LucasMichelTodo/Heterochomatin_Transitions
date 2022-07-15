#### Functions ####
import subprocess
import os
import itertools
import pybedtools as pb
import subprocess as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import pandas as pd
import itertools

## Functions
def get_RPKMs(bam, bs, smooth, norm, outfld):
    
    outfld = outfld if outfld.endswith('/') else outfld + '/'
    name = bam.rsplit('/', 1)[1]
    outname = outfld+name.replace('.bam', f'_rpkm_raw_bs{bs}_smth{smooth}.bdg')

    cmd = (
        f'bamCoverage -b {bam} '
        '--outFileFormat bedgraph '
        f'--normalizeUsing {norm} '
        '-p 8 '
        f'-bs {bs} '
        f'--smoothLength {smooth} '
        f'-o {outname}'
    )
    print(cmd)
    subprocess.call(cmd, shell=True)


def get_RPKMs_normInput(
        bam_IP,
        bam_in,
        bs = 50,
        smooth = 150,
        norm = 'RPKM',
        of = 'bedgraph',
        outfld = './',
        pseudo = 1,
        num_process = 8,
):
    outfld = outfld if outfld.endswith('/') else outfld + '/'
    name = bam_IP.rsplit('/', 1)[1]
    outname = name.replace('.bam', f'_rpkm_normInput_bs{bs}_smth{smooth}_pseudo{pseudo}.bdg')
    cmd = (
        f'bamCompare -b1 {bam_IP} -b2 {bam_in} '
        '--outFileFormat bedgraph '
        '--scaleFactorsMethod None '
        f'--normalizeUsing {norm} '
        f'-p {num_process} '
        f'-bs {bs} '
        f'--smoothLength {smooth} '
        f'--pseudocount {pseudo} '
        f'-o {outfld+outname} '
        f'-of {of}'
    )

    print(cmd)
    subprocess.call(cmd, shell=True)

def params_to_name(params_dict):
    """
    Take a dicttionary with params and convert it to a string to add to name files.
    """
    suffix = ''
    for k,v in params_dict.items():
        k_str = k.replace('-', '').strip()
        v_str = v.strip()
        suffix += f'_{k_str}{v_str}'
    return(suffix)

def get_coverage(
        bam_IP,
        outfld = './',
        params = {
            '-bs':'50',
            '--smoothLength':'100',
            '--normalizeUsing':'RPKM',
            '-p':'8',
            '-of':'bedgraph',
        },
):
    """
    Wrapper for deepTools bamCoverage command. 
    If you want to run with custom params,
    copy the default params dict and then modify and add/remove
    the ones you want to change.
    """
    name = bam_IP.rsplit('/', 1)[1].rsplit('.', 1)[0]
    suffix = params_to_name(params)
    ext = '.bdg' if params['-of'] == 'bedgraph' else '.bw'
    outfile = outfld+name+suffix+ext

    cmd = f'bamCoverage -b {bam_IP} -o {outfile}'
    for k,v in params.items():
        cmd += ' '+k+' '+v

    print(cmd)
    subprocess.call(cmd, shell=True)


def get_coverage_normInput(
        bam_IP,
        bam_in,
        outfld = './',
        params = {
            '-bs':'50',
            '--smoothLength':'100',
            '--scaleFactorsMethod':'None',
            '--normalizeUsing':'RPKM',
            '-p':'8',
            '-of':'bedgraph',
        },
):
    """
    Wrapper for deepTools bamCompare command. 
    If you want to run with custom params,
    copy the default params dict and then modify and add/remove
    the ones you want to change.
    """
    name = bam_IP.rsplit('/', 1)[1].rsplit('.', 1)[0]
    suffix = params_to_name(params)
    ext = '.bdg' if params['-of'] == 'bedgraph' else '.bw'
    outfile = outfld+name+suffix+ext

    cmd = f'bamCompare -b1 {bam_IP} -b2 {bam_in} -o {outfile}'
    for k,v in params.items():
        cmd += ' '+k+' '+v

    print(cmd)
    subprocess.call(cmd, shell=True)

def bigWigCompare(
        bw_IP,
        bw_in,
        outfld = './',
        params = {
            '-bs':'50',
            '-p':'8',
            '-of':'bedgraph',
        },
):
    """
    Wrapper for deepTools bigwigCompare command. 
    If you want to run with custom params,
    copy the default params dict and then modify and add/remove
    the ones you want to change.
    """
    name = bw_IP.rsplit('/', 1)[1].rsplit('.', 1)[0]
    suffix = params_to_name(params)
    ext = '.bdg' if params['-of'] == 'bedgraph' else '.bw'
    outfile = outfld+name+suffix+ext

    cmd = f'bigwigCompare -b1 {bw_IP} -b2 {bw_in} -o {outfile}'
    for k,v in params.items():
        cmd += ' '+k+' '+v
        
    print('\n--------------\n')
    print(cmd)
    print('\n--------------\n')
    subprocess.call(cmd, shell=True)
    
def subtract_bams(bam1, bam2, bs, smooth, outname):

    cmd = (f'bamCompare -b1 {bam1} -b2 {bam2} '
           '--operation subtract '
           f'-bs {bs} '
           '-p 8 '
           '--outFileFormat bedgraph '
           '--scaleFactorsMethod None '
           '--normalizeUsing RPKM '
           f'--smoothLength {smooth} '
           f'-o {outname}')
    
    sp.call(cmd, shell=True)


def getName(filename):
    if '/' in filename:
        filename = filename.rsplit("/", 1)[1] #Remove path
    out = filename.rsplit(".", 1)[0] #Remove extension
    return(out)


def getDepth(peaksfile):
    dlist = []
    with open(peaksfile) as f:
        for line in f:
            if line.startswith("# total fragments "):
                d = line.split(":")[1].strip()
                dlist.append(int(d))
    depth = round(min(dlist) / 1000000, 2)
    return(str(depth))


def macs2callpeak(t, c, params):

    # Make sure we are using apropiate MACS2 version (2.1.2)
    print("You are using MACS version:")
    cmd = "macs2 --version"
    subprocess.call(cmd, shell=True)
    print("\n")
    n = t.rsplit('/', 1)[1].replace('.bam', '')
    cmd = ("macs2 callpeak -t {} -c {} -n {} ") .format(t, c, n) + params

    subprocess.call(cmd, shell=True)


def macs2DifPeaks(indir, t1, c1, t2, c2, g, l, outdir, c='3'):

    # Make sure we are using apropiate MACS2 version (2.1.2)
    print("You are using MACS version:")
    cmd = "macs2 --version"
    subprocess.call(cmd, shell=True)
    print("\n")

    t1pile = indir + getName(t1) + "_treat_pileup.bdg"
    c1pile = indir + getName(t1) + "_control_lambda.bdg"

    t2pile = indir + getName(t2) + "_treat_pileup.bdg"
    c2pile = indir + getName(t2) + "_control_lambda.bdg"

    peaks1 = indir + getName(t1) + "_peaks.xls"
    peaks2 = indir + getName(t2) + "_peaks.xls"

    d1 = getDepth(peaks1)
    d2 = getDepth(peaks2)

    prefix = getName(t1)+"_vs_"+getName(t2)+"_g{}_l{}_c{}" .format(g,l,c)

    cmd = ("macs2 bdgdiff "
           "--t1 {} --c1 {} "
           "--t2 {} --c2 {} "
           "--d1 {} --d2 {} "
           "--outdir {} "
           "--o-prefix {} "
           "-g {} -l {} --cutoff {}") .format(t1pile, c1pile,
                                              t2pile, c2pile,
                                              d1, d2,
                                              outdir,
                                              prefix,
                                              g, l, c)

    print(cmd)
    subprocess.call(cmd, shell=True)

def manormDifPeaks(peaks1, peaks2, reads1, reads2, params):

    #Check which MaNorm version we are using:
    print("\n")
    print("You are using MANorm version:")
    cmd = "manorm --version"
    subprocess.call(cmd, shell=True)
    print("\n")

    cmd = ("manorm --p1 {} --p2 {} --r1 {} --r2 {} ") .format(peaks1, peaks2,
                                                              reads1, reads2)
    cmd = cmd + params
    subprocess.call(cmd, shell=True)

def create_log2_tracks(treat, ctl):
    outfile = treat.replace('.bam', '_log2.bdg')
    cmd = 'bamCompare -b1 {} -b2 {} -o {} -of bedgraph -p8' .format(treat, ctl, outfile)
    subprocess.call(cmd, shell=True)


class ChIP_seq():
    def __init__(self, ID, ip=None, bf=None, s=None, c=None, h=None, com=None):
        self.ID = ID
        self.ip_type = ip
        self.bamfile = bf
        self.strain = s
        self.condition = c
        self.hpi = h
        self.conmments = com

    
def pairwise(iterable):
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return zip(a, a)

