#### Create "only_genes_gff" ####

import pybedtools as pb
import subprocess as sp
import os

wd = '/mnt/Disc4T/Projects/PhD_Project/Paper/Paper_Analysis/Data_Files/'
os.chdir(wd)

# Create genome dict
genome={}
genome_file = './PlasmoDB-46_Pfalciparum3D7_Genome.fasta' #later versions have 5'UTRs and 3'UTRs
with open(genome_file, 'r+') as infile:
    for line in infile:
        if line.startswith('>'):
            linelist = line.strip().split(' | ')
            chrom = linelist[0].replace('>', '')
            seize = linelist[3].replace('length=', '')
            genome[chrom] = (0, int(seize))

# Import GFF
gff_file = './PlasmoDB-46_Pfalciparum3D7_withGDV1_ncRNAs.gff'
gff = pb.BedTool(gff_file)

## Filter genes only
sel = ['gene']
gene_gff = gff.filter(lambda x: x.fields[2] in sel)

## Discard Apicoplast
gene_gff = gene_gff.filter(lambda x: x.chrom != 'Pf3D7_API_v3')
gene_gff.saveas(gff_file.replace('.gff', '_only_genes.gff'))

## Sort Gene-GFF
unsorted = gff_file.replace('.gff', '_only_genes.gff')
sorted = gff_file.replace('.gff', '_only_genes_sorted.gff')
cmd = f'python3 ./gff_sorter.py {unsorted} > {sorted}'
sp.call(cmd, shell=True)

#### Find genes that overlap neighboring genes ####

import pybedtools as pb
import os

wd = '/mnt/Disc4T/Projects/PhD_Project/Paper/Paper_Analysis/Data_Files/'
os.chdir(wd)

gff = pb.BedTool('./PlasmoDB-46_Pfalciparum3D7_withGDV1_ncRNAs_only_genes_sorted.gff')

exceptions = set()
for idx, gene in enumerate(gff):
    if idx > 0 and idx < len(gff)-1:
        gid = gene.fields[8].split(';')[0].replace('ID=', '')

        if gene.chrom == gff[idx-1].chrom:
            if gene.start < gff[idx-1].stop or gene.stop <= gff[idx-1].stop:
                exceptions.add(gid)
                print(gid)

        if gene.chrom == gff[idx+1].chrom:
            if gene.start >= gff[idx+1].start or gene.stop > gff[idx+1].start:
                exceptions.add(gid)
                print(gid)

with open('gene_exceptions.txt', 'w+') as outfile:
    for gene in exceptions:
        outfile.write(gene+'\n')

#### Get binned BED for desired genomic regions ####

import pybedtools as pb
import os

wd = '/mnt/Disc4T/Projects/PhD_Project/Paper/Paper_Analysis/Data_Files/'
os.chdir(wd)

def get_5prime_ORF_3prime(gff, genome_file, excep_file, fiveP=1000, threeP=1000, orf=[0,500], allowoverlap = False):

    if not allowoverlap:
    ## Get a list of genes that fall inside another gene
        with open(excep_file, 'r+') as f:
            exceptions = [g.strip() for g in f.readlines()]

    ## Create dict with lenght of each chromosome
    genome={}
    with open(genome_file, 'r+') as infile:
        for line in infile:
            if line.startswith('>'):
                linelist = line.strip().split(' | ')
                chrom = linelist[0].replace('>', '')
                seize = linelist[3].replace('length=', '')
                genome[chrom] = (0, int(seize))

    ## Set some variables and start!
    ref = pb.BedTool(gff)
    current_chrom = ''
    ngenes = len(ref)
    str_bed = ''
    first_in_chrom = False
    last_in_chrom = False

    for idx, gene in enumerate(ref):

        ## Get gene id
        gid = gene.fields[8].split(';')[0].replace('ID=', '')

        ## Check Orientation:
        strand = gene.fields[6]

        ## Check if first/last in chromosome
        chrom = gene.chrom

        if current_chrom != chrom:
            first_in_chrom = True

        if idx == ngenes-1:
            ## First check if we are in the last gene!
            last_in_chrom = True
        else:
            if ref[idx+1].chrom != chrom:
                last_in_chrom = True

        ## Set new start5, stop5 and star3, stop3 depending on strand:

        if strand == '+':

            prestart = gene.start-fiveP
            prestop = gene.start
            poststart = gene.stop
            poststop = gene.stop+threeP

        else:
            prestart = gene.start-threeP
            prestop = gene.start
            poststart = gene.stop
            poststop = gene.stop+fiveP

        ## Set ORF start-stop
        if strand == '+':
            orf_start = gene.start + orf[0]
            if orf_start > gene.stop: orf_start = orf_stop -1
            orf_stop = gene.start + orf[1]
            if orf_stop > gene.stop : orf_stop = gene.stop
        else:
            orf_start = gene.stop - orf[1]
            if orf_start < gene.start: orf_start = gene.start
            orf_stop = gene.stop - orf[0]
            if orf_stop < gene.start: orf_stop = gene.start + 1

        if not allowoverlap:
        ## Check overlapp previous gene if +strand or next gene if -strand
        ## Except for genes in exception list

            if gid not in exceptions:
                if first_in_chrom:
                    pass
                else:
                    if prestart < ref[idx-1].stop:
                        prestart = ref[idx-1].stop

                if last_in_chrom:
                    pass
                else:
                    if poststop > ref[idx+1].start:
                        poststop = ref[idx+1].start

        ## Check we dont go < 0
        if prestart < 0: prestart = 0
        if orf_start < 0: orf_start = 0

        ## Check we don't go > chrom length
        if poststop > genome[chrom][1]: poststop = genome[chrom][1]
        if orf_stop > genome[chrom][1]: orf_stop = genome[chrom][1]

        ## Check start always start < stop
        if prestart >= prestop:
            prestop = prestart+1
            print(f'Pre region error! In gene :{gid}')
            print(prestart, prestop)
        if orf_start >= orf_stop:
            print(f'ORF region error! In gene :{gid}')
            print(orf_start, orf_stop)
            orf_stop = orf_start+1
        if poststart >= poststop:
            poststop = poststart+1
            print(f'Post region error! In gene :{gid}')
            print(poststart, poststop)


        ## Reset variables
        first_in_chrom = False
        last_in_chrom = False
        current_chrom = chrom

        ## Prepare output
        presuffix = '_5prime' if strand == '+' else '_3prime'
        postsuffix = '_3prime' if strand == '+' else '_5prime'

        preline = '\t'.join([gene.chrom,
                             str(prestart),
                             str(prestop),
                             gid+presuffix, '.', strand])

        ORFline = '\t'.join([gene.chrom,
                             str(orf_start),
                             str(orf_stop),
                             gid, '.', strand])

        postline = '\t'.join([gene.chrom,
                              str(poststart),
                              str(poststop),
                              gid+postsuffix,'.', strand])

        str_bed += ('\n'.join([preline, ORFline, postline])+'\n')

    ## Convert strig output into bedtools object and return
    out_bed = pb.BedTool(str_bed, from_string=True)
    #out_bed.saveas(f'binned_5prime{fiveP}_ORF_{orf[0]}_{orf[1]}_3prime{threeP}.bed')
    return(out_bed)

ref = './PlasmoDB-46_Pfalciparum3D7_withGDV1_ncRNAs_only_genes_sorted.gff'
excep_file = './gene_exceptions.txt'
genome_file = './PlasmoDB-46_Pfalciparum3D7_Genome.fasta'

fiveP = 1000
threeP = 1000
orf = [0,500]
allowoverlap = False

binned_bed = get_5prime_ORF_3prime(ref, genome_file, excep_file, fiveP, threeP, orf, allowoverlap)
binned_bed.saveas(f'./Binned_Beds/binned_5prime{fiveP}_ORF_{orf[0]}_{orf[1]}_3prime{threeP}_allowoverlap{allowoverlap}.bed')

## Keep only 3' coverage:
file_3_orf_5 = 'binned_5prime1000_ORF_0_500_3prime1000_allowoverlapFalse_new.bed'
outfile = 'binned_3prime1000_allowoverlapFalse_new.bed'

with open(file_3_orf_5, 'r+') as infile:
    with open(outfile, 'w+') as output:
        for line in infile:
            if '_3prime' in line:
                output.write(line.replace('_3prime', ''))
