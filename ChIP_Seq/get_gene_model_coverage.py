#### Create binned bed file for gene model ####

import os
import pybedtools as py
import numpy as np

def bin_region(start, stop, nbins=10):

    cuts = np.linspace(start, stop, num=(nbins+1), dtype = "int")
    bins = []

    for i in range(0,len(cuts)-1):
        bins.append((cuts[i], cuts[i+1]))

    return(bins)

def genome_to_dict(genome_file):

    chr_sizes = {}

    with open(genome_file) as infile:
        for line in infile:
            vals = line.strip().split()
            chr_sizes[vals[0]] = int(vals[1])

    return(chr_sizes)

def getGeneId(gene_bed):
    gid = gene_bed.fields[8].split(";")[0].replace("ID=", "")
    return(gid)

def elongate_and_bin_GFF(gff, genome, nbins=5):

    # Ensure output is overwritten
    prefix = "/mnt/Disc4T/Projects/PhD_Project/Data/"
    name = gff.rsplit( ".", 1 )[0].rsplit("/", 1)[1]
    sufix = "_bin"+str(nbins)+"_2prevGenes"
    outname = "".join([prefix, name, sufix])+".bed"

    print("Ouput will be written in:\n %s" % outname)

    fl = open(outname, "w+")
    fl.close()

    ## Load genome
    chr_sizes = genome_to_dict(genome)

    ## Load GFF for annotation
    ref = py.BedTool(gff)
    ref = ref.filter(lambda x: x[2] == "gene")
    ref = ref.sort()

    ## Create a variable for number of gene-bis
    gbins = nbins*3

    for gene in ref:

        ## Get gene ID
        gid = getGeneId(gene)
        chrom = gene.chrom
        start = int(gene.start)+1 #Compensate a base (maybe it comes from gff/bed)
        stop = int(gene.stop)

        ## Create a mini-bed for the gene and look for
        ## closest 2 genes before and after it (but not overlapping)
        linebed = "\t".join([gene.chrom, str(gene.start), str(gene.stop)])
        gene_bed = py.BedTool(linebed, from_string=True)

        pregenes = gene_bed.closest(ref, D = "ref", id = True, io = True, k = 2)
        postgenes = gene_bed.closest(ref, D = "ref", iu = True, io = True, k = 2)

        ## Resolve cases in which a gene has none or just one gene after/before it.
        ## In case there is no gene before/after:
        ## Set pre region to gene-start -2000 or start of chrom.
        ## Set post region to gene-stop + 2000 or end of chrom.

        if len(pregenes) == 1:

            pregene = pregenes[0]

            if pregene.fields[3] == ".": #No gene before
                pre_stop = max([start-2000, 0])
                pre_start = pre_stop
                prepre_start, prepre_stop = 0, 0

            else: #Just 1 gene before
                pre_start, pre_stop = pregene.fields[6], pregene.fields[7]
                prepre_start, prepre_stop = 0, 0

        else:

            pre_start, prepre_start = pregenes[0].fields[6], pregenes[1].fields[6]
            pre_stop, prepre_stop = pregenes[0].fields[7], pregenes[1].fields[7]

        if len(postgenes) == 1:

            postgene = postgenes[0]

            if postgene.fields[3] == ".": #No gene after
                post_start = min([stop+2000, chr_sizes[chrom]])
                post_stop = post_start
                postpost_start, postpost_stop = 0, 0

            else: #Just 1 gene after
                post_start, post_stop = postgene.fields[6], postgene.fields[7]
                postpost_start, postpost_stop = 0, 0

        else:

            post_stop, postpost_stop = postgenes[0].fields[7], postgenes[1].fields[7]
            post_start, postpost_start = postgenes[0].fields[6], postgenes[1].fields[6]


        ## Create bins

        prepre_gen_cov = bin_region(int(prepre_start), int(prepre_stop), nbins=2)
        pre_gen_cov = bin_region(int(pre_start), int(pre_stop), nbins=2)

        postpost_gen_cov = bin_region(int(postpost_start), int(postpost_stop), nbins=2)
        post_gen_cov = bin_region(int(post_start), int(post_stop), nbins=2)

        pre = bin_region(int(pre_stop), start, nbins=nbins)
        body = bin_region(start, stop, nbins=nbins)
        post = bin_region(stop, int(post_start), nbins=nbins)

        ## Print output taking into accound strandness
        ## We cannot sort here (reverse if strand is "-")
        ## because we will tabix afterwards.
        ## We create a column for sorting afterwards.

        output = []

        regions = [
            prepre_gen_cov,
            pre_gen_cov,
            pre,
            body,
            post,
            post_gen_cov,
            postpost_gen_cov,
        ]

        if gene.strand == "-":
            order = range(gbins+8, 0, -1)
        else:
            order = range(1, gbins+9)

        i = 0
        for reg in regions:
            for interval in reg:
                output.append([chrom,
                               str(interval[0]),
                               str(interval[1]),
                               gid, str(order[i])])
                i += 1

        for line in output:
            with open(outname, "a+") as outfile:
                outfile.write("\t".join(line)+"\n")


gff = "/mnt/Disc4T/Projects/PhD_Project/Data/PlasmoDB-46_Pfalciparum3D7_withGDV1_ncRNAs.gff"
gen = "/home/lucas/ISGlobal/Gen_Referencies/Pf3D7.genome"

elongate_and_bin_GFF(gff, gen)

#### BGZIP and TABIX coverage files ####

import os
import subprocess as sp

wd = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/RPKMs_NormInput_noDup_bs10_smth_200_pseudo_10'
os.chdir(wd)

# cov_files = [
#     'E5K9_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
#     'A7K9_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
#     '1.2B_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
#     '10G_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
#     'B11_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg'
# ]

cov_files = [
    '1.2B_ac_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
    '10G_ac_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
    'B11_ac_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg'
]


for f in cov_files:
    outf = f.replace('.bdg', '.bdg.gz')
    cmd = f'bgzip -c -f {f} > {outf}'
    print(cmd)
    sp.call(cmd, shell=True)

gz_files = [f for f in os.listdir() if f.endswith('.bdg.gz') and '_ac_' in f]

for f in gz_files:
    cmd = f'tabix -p bed {f}'
    print(cmd)
    sp.call(cmd, shell=True)

#### Cross Coverage data with binned BED file ####

import numpy as np
import pandas as pd
import os
import pybedtools as pb
from collections import defaultdict

wd = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/RPKMs_NormInput_noDup_bs10_smth_200_pseudo_10/'
os.chdir(wd)

indir = './'
bin_bed = pb.BedTool("/mnt/Disc4T/Projects/PhD_Project/Paper/Paper_Analysis/Data_Files/PlasmoDB-46_Pfalciparum3D7_withGDV1_ncRNAs_bin5_2prevGenes.bed")

#cov_files = [f for f in os.listdir(indir) if f.endswith(".bdg.gz")]
cov_files = [f for f in os.listdir(indir) if f.endswith(".bdg.gz") and '_ac_' in f]
cov_files

for cov in cov_files:

    #cov = cov_files[0]
    flstr = indir+cov
    cov = pb.BedTool(flstr)
    print(flstr, "Converted to bed!")
    outfile = flstr.replace(".bdg.gz", "_5binned_cov_2prevpost.csv")
    #logfile = outfile.replace(".bed", ".log")
    genevals = defaultdict(list)

    for interval in bin_bed:

        #interval = bin_bed[0]
        gene = interval.name
        pos = interval.score

        # Not all of them have a match!
        try:
            match = cov.tabix_intervals(interval)
            val = np.mean([float(x.fields[3]) for x in match])
            genevals[gene].append((val, pos))

        except:
            pass

    # Rearrange values deppending on strandness (we have to "flip" genes on "-" strand)
    sorted_genevals = {}
    for gene, val in genevals.items():
        svals = sorted(val, key = lambda x:int(x[1]))
        vals = [x[0] for x in svals]
        sorted_genevals[gene] = vals

    # Write output
    df = pd.DataFrame.from_dict(sorted_genevals, orient='index')
    df.to_csv(outfile)
    print("Done with file: {}" .format(flstr))
