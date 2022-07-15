#### Check for Del and Dupl ####

import os
import pybedtools as pb
import numpy as np
import subprocess as sp

##### Functions #####

def add_names(bed, featname):
    i = 1
    str_bed = ''
    for feat in bed:
        newline = [feat.chrom, feat.start, feat.stop, featname+'_'+str(i), feat.name]
        newline = [str(x) for x in newline]
        str_bed += '\t'.join(newline)+'\n'
        i += 1
    outbed = pb.BedTool(str_bed, from_string=True)
    return(outbed)

def set_score_x(feat, x):
        feat.score = x
        return(feat)

def get_dupl_del(bed, fact_up, fact_dw, score_col, mergelen, minlen):

    bedname = bed.rsplit('/', 1)[1]
    print(bedname)
    bed = pb.BedTool(bed)
    cov = [float(feat.fields[score_col]) for feat in bed]

    th_up = np.mean(cov)*fact_up
    th_dw = np.mean(cov)*fact_dw

    thbed_up = bed.filter(lambda x: float(x.name) >= th_up)
    thbed_dw = bed.filter(lambda x: float(x.name) <= th_dw)

    ## Cluster peaks together
    clu_bed_up = thbed_up.merge(d = mergelen, c = 4, o = 'mean')
    clu_bed_dw = thbed_dw.merge(d = mergelen, c = 4, o = 'mean')

    ## Filter peaks by length
    len_bed_up = clu_bed_up.filter(lambda x: float(x.stop) - float(x.start) >= minlen)
    len_bed_dw = clu_bed_dw.filter(lambda x: float(x.stop) - float(x.start) >= minlen)

    ## Add name (and move score to 5th column)

    blueprint = '_bymean_{}_fact_{}_minlen{}_mergelen_{}.bed'

    suffix_up =  blueprint .format('dupl', fact_up, minlen, mergelen)
    outname_up = ddfld+bedname.replace('.bdg', suffix_up)

    suffix_dw = blueprint .format('del', fact_dw, minlen, mergelen)
    outname_dw = ddfld+bedname.replace('.bdg', suffix_dw)

    dupl_bed = add_names(len_bed_up, 'duplication').each(set_score_x, 1).saveas(outname_up)
    del_bed = add_names(len_bed_dw, 'deletion').each(set_score_x, -1).saveas(outname_dw)

    ## Merge dupl and del
    strfact = '{}up_{}dw' .format(fact_up, fact_dw)
    suffix = blueprint .format('dupl_del', strfact, minlen, mergelen)
    outname = ddfld+bedname.replace('.bdg', suffix)
    cmd = f'cat {outname_up} {outname_dw} > {outname}'
    sp.call(cmd, shell = True)

    ## Sort and set score for all features to 1
    outbed = pb.BedTool(outname).sort().saveas(outname)


##### Calls #####

wd = '/mnt/Disc4T/Projects/PhD_Project/Paper/Paper_Analysis/'
os.chdir(wd)

ddfld = './Data_Files/Duplication_Deletion_Regions_Mean_Separate_DuplDel/'
os.makedirs(ddfld, exist_ok=True)

rpkms_fld = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/RPKMs_noDup_bs10_smth_200/'

in_files = [f for f in os.listdir(rpkms_fld) if '_in_' in f and f.endswith('.bdg')]

## Params

fact_up = 1.75
fact_dw = 0.1
score_col = 3
mergelen = 200
minlen = 500

for f in in_files:
    get_dupl_del(rpkms_fld+f, fact_up, fact_dw, score_col, mergelen, minlen)

#### Check dupl/del present in ALL strains ####

## Deletions present in all strains likely represent repetitive regions poorly mapped,
## by the aligner and instead of real deletions.

ddfld = './Data_Files/Duplication_Deletion_Regions_Mean_Separate_DuplDel/'

os.listdir(ddfld)

prefixes = ('1.2B', '10G', 'A7K9', 'E5K9', 'B11')
suffix = '_bymean_dupl_del_fact_1.75up_0.1dw_minlen500_mergelen_200.bed'
files_to_cross = [f for f in os.listdir(ddfld) if f.endswith(suffix) and f.startswith(prefixes)]

##### Join all dupl_del files and select regions dupl/del in ALL samples

str_beds = ' '.join([ddfld+f for f in files_to_cross])
cmd = 'awk \'{print}\' '+f'{str_beds} > {ddfld}allstrains_supl_del.bed'
sp.call(cmd, shell = True)

join_bed = pb.BedTool(ddfld+'allstrains_supl_del.bed').sort()

result = join_bed.genome_coverage(bg=True, g='./Data_Files/Pf3D7.genome')
result.saveas(ddfld+'final_allsrtains_dupl_del.bdg')
filter_bed = result.filter(lambda x: int(x.name) >= 5)
filter_bed.saveas(ddfld+'final_allstrains_dupl_del_>5.bdg')

##### Cross each dupl/del file with the ALLstraisn dupl/del file

## Functions

def check_overlapp_perc(feat, bed, perc_th, exclude = False):
    """
    Check wether "feat" overlaps any feature in "bed" and if the overlapp spans
    >= "perc_th" % of "feat" (in length) If exclude = True keep peaks that don't
    overlapp another.
    """
    is_match = False
    for interval in bed:
        if interval.chrom == feat.chrom:
            # Check overlapp
            if feat.start <= interval.stop and interval.start <= feat.stop:
                #Check percentage overlapp
                bigger_start = max([feat.start, interval.start])
                smaller_end = min([feat.stop, interval.stop])
                perc_overlapp = (smaller_end - bigger_start/len(feat))*100

                if perc_overlapp >= perc_th: is_match = True
            else:
                pass
        else:
            pass

    if exclude: is_match = not is_match
    return(is_match)

## Calls

files_to_cross
perc = 80

## Retain only NON-overlapping peaks

for bed in files_to_cross:

    ref = 'final_allstrains_dupl_del_>5.bdg'
    refbed = pb.BedTool(ddfld+ref)

    outname = bed.replace('.bed', '_filtered.bed')
    bed1 = pb.BedTool(ddfld+bed)
    outbed = bed1.filter(lambda b: check_overlapp_perc(b, refbed, perc, exclude=True))
    outbed.saveas(ddfld+outname)
    print(outname)

#### Cross dupl/del with genes ####

import os
import pybedtools as pb

wd = '/mnt/Disc4T/Projects/PhD_Project/Paper/Paper_Analysis/'
os.chdir(wd)

ddfld = './Data_Files/Duplication_Deletion_Regions_Mean_Separate_DuplDel/'
dupl_del_files = [f for f in os.listdir(ddfld) if f.endswith('_filtered.bed')]

ref_dir = './Data_Files/PlasmoDB-52_Pfalciparum3D7.gff'
gff = pb.BedTool(ref_dir)

types = set([entry.fields[2] for entry in gff])
gene_types = [
    'ncRNA_gene',
    'protein_coding_gene',
    'pseudogene'
]

outdir = ddfld+'/Crossed_with_genes/'
os.makedirs(outdir, exist_ok=True)

for bed_f in dupl_del_files:

    gene_gff = gff.filter(lambda x: x.fields[2] in gene_types)
    print(bed_f)
    outfile = outdir+bed_f.replace('.bed', '_genes.tsv')
    dd_bed = pb.BedTool(ddfld+bed_f)
    cross = dd_bed.intersect(gene_gff, wao = True)

    with open(outfile, 'w+') as out_file:
        for x in cross:
            #print(x)
            if x.fields[5] != '.':
                #print(x)
                attrs_field = x.fields[13].split(';')
                attrs_dict = {x.split('=')[0]:x.split('=')[1] for x in attrs_field}
                out_line = [
                    attrs_dict['ID'],
                    attrs_dict.get('Name', ''),
                    attrs_dict.get('description', '')
                ]
                #print(out_line)
                out_file.write('\t'.join(out_line)+'\n')
