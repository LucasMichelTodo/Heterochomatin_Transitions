#### Clean Reads with BBDUK ####

import subprocess as sp
import os

## Functions

def call_BBDUK(in1, in2, out1, out2, outm, ref, params):
    cmd = ("bbduk.sh in={} in2={} "
           "out={} out2={} outm={} "
           "ref={}") .format(in1, in2, out1, out2, outm, ref)

    cmd = cmd+" "+params
    sp.call(cmd, shell=True)

## Calls

params = "ktrim=r k=22 mink=6 overwrite=t"
ref = "/home/lucas/Programs/bbmap/resources/adapters.fa"

root_path = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Raw_Data/"
read1s = []
read2s = []
for path, subdirs, files in os.walk(root_path):
    for f in files:
        if all(x in f for x in ["read1", ".fastq"]):
            read1s.append(os.path.join(path, f))
        elif all(x in f for x in ["read2", ".fastq"]):
            read2s.append(os.path.join(path, f))
        else:
            print(f)


read1s = sorted(read1s)
read2s = sorted(read2s)
outpath = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Clean_Reads/"

for pair in zip(read1s, read2s):
    in1, in2 = pair[0], pair[1]
    out1 = outpath + pair[0].rsplit("/", 1)[1].replace(".fastq", "_clean.fastq")
    out2 = outpath + pair[1].rsplit("/", 1)[1].replace(".fastq", "_clean.fastq")
    outm = outpath + pair[0].rsplit("/", 1)[1].replace(".fastq", "_badreads.fastq")
    call_BBDUK(in1, in2, out1, out2, outm, ref, params)

#### Align with Bowtie2 ####

## Funtions

def call_Bowtie2(in1, in2, out, params):
    cmd = "bowtie2 -1 {} -2 {} -S {}" .format(in1, in2, out)
    cmd = cmd+" "+params
    print(cmd)
    sp.call(cmd, shell=True)

## Calls
params = ("-p 4 --very-sensitive --local "
          "-5 4 -3 4 -I 50 -X 200 "
          "-x /home/lucas/Programs/bowtie2-2.3.0-legacy/Pf3D7")

inpath = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/Clean_Reads/"
files = [f for f in os.listdir(inpath)if f.endswith("_clean.fastq.gz")]

read1s = sorted([f for f in files if "read1" in f])
read2s = sorted([f for f in files if "read2" in f])

outpath = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/Current/"

for pair in zip(read1s, read2s):
    name = pair[0].split("_")[0]

    in1, in2 = inpath+pair[0], inpath+pair[1]
    out = outpath+name+".sam"
    call_Bowtie2(in1, in2, out, params)

#### Convert to BAM, sort and index ####

import subprocess as sp
import sys
from tqdm import tqdm

## Functions

def from_sam_to_bam(samfile):
    name = samfile.rsplit(".")[0]
    cmd = "samtools view -bS {} > {}" .format(samfile, name+".bam")
    sp.call(cmd, shell=True)

    ### Erase SAM after creating BAM
    # cmd = "rm {}" .format(samfile)
    # sp.call(cmd, shell=True)

    cmd = "samtools sort {} > {}" .format(name+".bam", name+"_sort.bam")
    sp.call(cmd, shell=True)

    ### Erase bam after creating sortedBAM
    cmd = "rm {}" .format(name+".bam")
    sp.call(cmd, shell=True)

    cmd = "samtools index {} > {}" .format(name+"_sort.bam", name+"_sort.bam.bai")
    sp.call(cmd, shell=True)

    ## Filter only >=q5 reads
    cmd = "samtools view -b -q 5 {} > {}" .format(name+"_sort.bam", name+"_q5_sort.bam")
    sp.call(cmd, shell=True)

## Calls

#indir = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/Bams/"
indir = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/New_Bams/"
samfiles = [f for f in os.listdir(indir) if f.endswith(".sam")]

for f in tqdm(samfiles):
    sam = indir+f
    from_sam_to_bam(sam)

#### Remove duplicates using GATK MarkDuplicates ####

import os
import subprocess as sp

def remove_duplicates(indir, outdir, bam):

    gatk_path = '/home/lucas/Programs/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar'
    i = indir+bam
    o = outdir+bam.replace(".bam", "_noDup.bam")
    m = outdir+bam.replace(".bam", "_metrics.txt")

    cmd = (f"java -jar {gatk_path} MarkDuplicates "
           "REMOVE_DUPLICATES=true I={} O={} M={}") .format(i, o, m)

    sp.call(cmd, shell=True)


indir = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/Bams/'
bams = sorted([b for b in os.listdir(indir) if b.endswith("sort_q5.bam")])
outdir = indir

bams = ['E5HA_ac_renamed_sort_q5.bam']


for bam in bams:
    remove_duplicates(indir, outdir, bam)
