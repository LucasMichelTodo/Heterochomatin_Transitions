import os
from chip_seq_processing import *

wd = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/'
os.chdir(wd)

bamdir = './Bams/'
bams = [b for b in os.listdir(bamdir) if b.endswith('_sort_q5.bam')]
IPs = [f for f in bams if '_me_' in f or '_ac_' in f]
IPs.sort()
inputs = [f for f in bams if '_in_' in f]
inputs.sort()
me_files = [f for f in bams if '_me_' in f]
me_files.sort()
ac_files = [f for f in bams if '_ac_' in f]
ac_files.sort()


#### Get RPKMs ####


wd = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/'
os.chdir(wd)

bamdir = './Bams/'
outdir = './RPKMs_noDup_bs10_smth_200/'

bams = [b for b in os.listdir(bamdir) if b.endswith('_noDup.bam')]

bs = 10
smooth = 200
norm = 'RPKM'

for bam in bams:
    out = outdir+bam.replace('.bam', '_RPKMs.bdg')
    get_RPKMs(bamdir+bam, bs, smooth, norm, outdir)

#### Get normalized by input RPKMs ####

outdir = './RPKMs_NormInput_noDup_bs10_smth_200_pseudo_10/'
IPs = [b for b in os.listdir(bamdir) if b.endswith('_noDup.bam') and '_in_' not in b]
inputs = [b for b in os.listdir(bamdir) if b.endswith('_noDup.bam') and '_in_' in b]

bs = 10
smooth = 200
norm = 'RPKM'
of = 'bedgraph'
outfld = outdir
pseudo = 10
num_process = 8

for ip in IPs:
    prefix = ip.split('_')[0]
    inpt = [f for f in inputs if prefix in f][0]
    print(ip, inpt)

    get_RPKMs_normInput(
        bamdir+ip,
        bamdir+inpt,
        bs,
        smooth,
        norm,
        of,
        outfld,
        pseudo
    )
