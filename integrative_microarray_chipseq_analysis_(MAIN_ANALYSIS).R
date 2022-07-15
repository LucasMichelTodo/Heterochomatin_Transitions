#### Imports and Dirs ####

library(ggplot2)
library(tidyverse)
library(readxl)
library(reshape2)
require(gridExtra)
library(readxl)
library(ggh4x)
library(ggrepel)
library(tsne)
library(scales)
library(viridis)
library(cluster)
library(NbClust)
library(factoextra)
library(class)
library(eulerr)
library(ggpubr)
library("colorspace")

wd <- '/mnt/Disc4T/Projects/PhD_Project/Paper/Paper_Analysis/'
setwd(wd)

## Data Dirs
microarrays_dir <- '../../Microarrays/New_Old_separate_approach/'
sicer_dir <- '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/DiffPeaks_SICER/PhD_Project_Strains/Overlapped_with_MACS2_Peaks/'
dupl_del_dir <- './Data_Files/Duplication_Deletion_Regions_Mean_Separate_DuplDel/Crossed_with_genes/'
coverage_dir <- './Data_Files/Coverages/'
genemodel_dir <- '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/RPKMs_NormInput_noDup_bs10_smth_200_pseudo_10/'
telomeres_dir <- './Data_Files/Telomeres/Coverages/'

## Output Dirs
tables_dir <- './Output_Tables/'
plots_dir <- './Plots/'


## Create Output Dirs
dir.create(tables_dir)
dir.create(paste0(tables_dir, '/Sicer_no_trans'))

dir.create(plots_dir)
dir.create(paste0(plots_dir, '/Coverage_PCAs'))
dir.create(paste0(plots_dir, '/Correlations_byStrain'))
dir.create(paste0(plots_dir, '/Trans_Het_byStrain'))
dir.create(paste0(plots_dir, 'Transcription_Donuts_ByStrain/'))
dir.create(paste0(plots_dir, 'MaxFC_Cor_Plots/'))
dir.create(paste0(plots_dir, 'Correlations_Intervals/'))
dir.create(paste0(plots_dir, 'Gene_Model/'))
dir.create(paste0(plots_dir, 'Gene_Model/PCAs/'))
dir.create(paste0(plots_dir, 'Gene_Model/Family_Heatmaps/'))
dir.create(paste0(plots_dir, 'Gene_Model/States_OnOff/'))
dir.create(paste0(plots_dir, 'Met_Ac/'))
dir.create(paste0(plots_dir, 'Met_Ac/PCAs/'))
dir.create(paste0(plots_dir, 'Met_Ac/Scatterplots/'))
dir.create(paste0(plots_dir, 'Met_Ac/Boxplots/'))

info_df <- read_tsv('./Data_Files/PlasmoDB-52_Pfalciparum3D7_parsed_annotation.tsv')

## Flag tRNAs
info_df <- info_df %>%
  mutate(Is_tRNA = grepl('tRNA', Annot, fixed = T) & Type == 'ncRNA_gene')

#### Load gene families data ####

gene_fam <- read_excel('./Data_Files/Supplementary_table_2_CVG_list_161120_ap.xlsx', sheet = 2)
gene_fam <- gene_fam %>%
  rename(Gene_id = `Gene ID`,
         Gene_name = `Gene Name or Symbol`,
         SubFamily = `Family Detail`) %>%
  mutate(SubFamily = case_when(SubFamily == 'var pseudo,truncated or -like.' ~ 'var-like',
                               TRUE ~ SubFamily)) %>%
  mutate(Gene_name = ifelse(Gene_name == 'N/A', NA, Gene_name)) %>%
  select(Gene_id, Gene_name, Family, SubFamily)

gene_fam %>%
  filter(Family == 'OTHER') %>%
  mutate(NewFam = case_when(is.na(SubFamily) ~ Gene_name,
                            !is.na(SubFamily) ~ SubFamily))

bigfams <- c(
  'VAR',
  'FIKK',
  'HYP',
  'PHIST',
  'RIFIN',
  'STEVOR',
  'PFMC-2TM'
)

info_df <- info_df %>%
  left_join(gene_fam) %>%
  mutate(Name = ifelse(is.na(Name), Gene_name, Name)) %>%
  mutate(Name = ifelse(is.na(Name) & Family != 'OTHER', Family, Name)) %>%
  mutate(Name = ifelse(Gene_id == 'PF3D7_0935390', 'GDV1as', Name)) %>%
  mutate(Label = ifelse(is.na(Name), Gene_id, paste(Gene_id, Name, sep = ': '))) %>%
  mutate(Family_Grouped = case_when(
           Family %in% bigfams ~ Family,
           !is.na(Family) & !Family %in% bigfams ~ 'Other CVGs',
           Gene_id == 'PF3D7_0935390' ~ 'Other CVGs',
           is.na(Family) ~ 'Not CVGs'))

table(info_df$Name == info_df$Gene_name)

info_df %>%
  filter(Name != Gene_name) %>%
  select(Gene_id, Name, Gene_name) %>%
  print(n = 40)

#### Load variant genes data ####
cvgs <- read.csv2(paste0(microarrays_dir, 'New_Arrays/Files/taula_CVG_final.csv'), stringsAsFactors = F)
cvgs <- cvgs %>%
  select(Gene_id = Gene.ID, Variant = Final.Customized) %>%
  mutate(Variant = ifelse(Variant == 'YES', TRUE, FALSE))


varlist <- filter(cvgs, Variant) %>% select(Gene_id) %>% pull()
varlist <- c(varlist, 'PF3D7_0935390')

info_df <- info_df %>%
  mutate(Variant = Gene_id %in% varlist)

## Probably Variant
## Some genes are from variant families but we have them as non-variant
info_df %>%
  filter(!Variant & !is.na(Family)) %>%
  select(Gene_id, Annot, Family) %>%
  print(n = 100)

## Check no variant gene has NA Family
info_df <- info_df %>%
  mutate(Family = ifelse(Gene_id == 'PF3D7_0935390', 'GDV1as', Family))

info_df %>%
  filter(Variant & is.na(Family))

gams <- read_csv(paste0(microarrays_dir, 'Oriol_Suplementary/gam_table.csv')) %>%
  select(Gene_id, Gam_specific)

info_df <- info_df %>%
  left_join(gams)

## Load DuplDel Data
suffix <- '_in_sort_q5_noDup_rpkm_normInput_bs10_smth200_bymean_dupl_del_fact_1.75up_0.1dw_minlen500_mergelen_200_filtered_genes.tsv'
prefix <- c('1.2B', '10G', 'A7K9', 'E5K9', 'B11')

get_dupl_del <- function(prefix){
  fname <- paste0(prefix, suffix)
  df <- read_tsv(paste0(dupl_del_dir, fname), col_names = F) %>%
    rename(Gene_id = X1, Name = X2, Annot = X3)
}

dupl_del <- lapply(prefix, get_dupl_del)
names(dupl_del) <- prefix

info_df <- info_df %>%
  mutate(DuplDel_12B = Gene_id %in% dupl_del$`1.2B`$Gene_id) %>%
  mutate(DuplDel_10G = Gene_id %in% dupl_del$`10G`$Gene_id) %>%
  mutate(DuplDel_A7 = Gene_id %in% dupl_del$A7K9$Gene_id) %>%
  mutate(DuplDel_E5 = Gene_id %in% dupl_del$E5K9$Gene_id) %>%
  mutate(DuplDel_B11 = Gene_id %in% dupl_del$B11$Gene_id)

#### Load Differential Peaks ####

difpeaks_dir <- './Data_Files/DifPeaks_W100_S100_PD0.3_Mg500_Ml1000/'
remove_transcript <- function(geneid){sub('\\.\\d$', '', geneid)}

difpeak_fls <- list.files(difpeaks_dir, pattern = 'gene_crossed.tsv')

get_peaks <- function(difpeaks_file){
  dp_colnames <- c('Chrom', 'Start', 'Stop', 'Gene_id')
  read_tsv(difpeaks_file, col_names = dp_colnames) %>%
    filter(Gene_id != 'intergenic') %>%
    select(Gene_id) %>%
    pull() %>%
    unique()
}

difpeaks <- lapply(difpeak_fls, function(x) get_peaks(paste0(difpeaks_dir, x)))
difpeaks_name <- function(x){
  paste(str_split(x, '_')[[1]][2], collapse = '_')
}
names(difpeaks) <- sapply(difpeak_fls, difpeaks_name)
difpeaks

difpeaks_df <- info_df %>%
  select(Gene_id) %>%
  mutate(peak_12B_10G = ifelse(Gene_id %in% difpeaks$`1.2Bover10G`, T, F)) %>%
  mutate(peak_12B_A7 = ifelse(Gene_id %in% difpeaks$`1.2BoverA7K9`, T, F)) %>%
  mutate(peak_12B_E5 = ifelse(Gene_id %in% difpeaks$`1.2BoverE5K9`, T, F)) %>%
  mutate(peak_12B_B11 = ifelse(Gene_id %in% difpeaks$`1.2BoverB11`, T, F)) %>%

  mutate(peak_10G_12B = ifelse(Gene_id %in% difpeaks$`10Gover1.2B`, T, F)) %>%
  mutate(peak_10G_A7 = ifelse(Gene_id %in% difpeaks$`10GoverA7K9`, T, F)) %>%
  mutate(peak_10G_E5 = ifelse(Gene_id %in% difpeaks$`10GoverE5K9`, T, F)) %>%
  mutate(peak_10G_B11 = ifelse(Gene_id %in% difpeaks$`10GoverB11`, T, F)) %>%

  mutate(peak_A7_12B = ifelse(Gene_id %in% difpeaks$`A7K9over1.2B`, T, F)) %>%
  mutate(peak_A7_10G = ifelse(Gene_id %in% difpeaks$`A7K9over10G`, T, F)) %>%
  mutate(peak_A7_E5 = ifelse(Gene_id %in% difpeaks$`A7K9overE5K9`, T, F)) %>%
  mutate(peak_A7_B11 = ifelse(Gene_id %in% difpeaks$`A7K9overB11`, T, F)) %>%

  mutate(peak_E5_12B = ifelse(Gene_id %in% difpeaks$`E5K9over1.2B`, T, F)) %>%
  mutate(peak_E5_10G = ifelse(Gene_id %in% difpeaks$`E5K9over10G`, T, F)) %>%
  mutate(peak_E5_A7 = ifelse(Gene_id %in% difpeaks$`E5K9overA7K9`, T, F)) %>%
  mutate(peak_E5_B11 = ifelse(Gene_id %in% difpeaks$`E5K9overB11`, T, F)) %>%

  mutate(peak_B11_12B = ifelse(Gene_id %in% difpeaks$`B11over1.2B`, T, F)) %>%
  mutate(peak_B11_10G = ifelse(Gene_id %in% difpeaks$`B11over10G`, T, F)) %>%
  mutate(peak_B11_A7 = ifelse(Gene_id %in% difpeaks$`B11overA7K9`, T, F)) %>%
  mutate(peak_B11_E5 = ifelse(Gene_id %in% difpeaks$`B11overE5K9`, T, F))

## Check Peaks that appear in both strains at the same time
sub_combs <- list(
  c('12B', '10G'),
  c('12B', 'A7'),
  c('12B', 'E5'),
  c('12B', 'B11'),
  c('10G', 'A7'),
  c('10G', 'E5'),
  c('10G', 'B11'),
  c('A7', 'E5'),
  c('A7', 'B11'),
  c('E5', 'B11')
)

for (i in 1:length(sub_combs)){
  print(sub_combs[i])
  difpeaks_df %>%
    select(Gene_id, contains(sub_combs[[i]][1]) & contains(sub_combs[[i]][2])) %>%
    set_names(c('Gene_id', 'Col1', 'Col2')) %>%
    filter(Col1 & Col2) %>%
    print()
  print('--------')
  print('')
  print('')
}

## Manually correct them
difpeaks_df[difpeaks_df$Gene_id == 'PF3D7_0936700',]$peak_B11_12B <- F
difpeaks_df[difpeaks_df$Gene_id == 'PF3D7_0936700',]$peak_B11_10G <- F
difpeaks_df[difpeaks_df$Gene_id == 'PF3D7_0401500',]$peak_B11_12B <- F

info_df %>%
  count(Family)

novar_family <- info_df %>%
  filter((!Variant & !is.na(Family)) | (Variant & is.na(Family)))

write_csv(novar_family, paste0(tables_dir, 'noVariant_with_family.csv'))
write_tsv(info_df, paste0(tables_dir, 'info_df.tsv'))

#### Load transcription data ####
trans_df_old <- read_csv(paste0(microarrays_dir, 'Old_Arrays/R_results_OldArrays_Variantome/final_summary_table.csv'))
trans_df_new <- read_csv(paste0(microarrays_dir, 'New_Arrays/R_results_NewArray/final_summary_table.csv'))

trans_df <- full_join(trans_df_old, trans_df_new) %>%
  select(-Name, -Annot, -GamGene)

## Subset trans_df to genes that appear on info_df
trans_df %>%
  filter(!Gene_id %in% info_df$Gene_id) %>%
  select(Gene_id) %>%
  write_csv(paste0(tables_dir, 'gene_ids_not_in_info_df.csv'))

trans_df <- trans_df %>%
  filter(Gene_id %in% info_df$Gene_id)

trans_df <- trans_df %>%
  mutate(`E5-B11_MaxVal` = -`B11-E5_MaxVal`) %>%
  mutate(`E5-B11_MaxTime` = `B11-E5_MaxTime`) %>%
  select(-`B11-E5_MaxVal`, -`B11-E5_MaxTime`)

## Load Areas Data
arees_old <- read_csv(paste0(microarrays_dir, 'Old_Arrays/R_results_OldArrays_Variantome/areaDiferences_geneLevel.csv')) %>%
  rename(Gene_id = ...1)

arees_new <- read_csv(paste0(microarrays_dir, 'New_Arrays/R_results_NewArray/areaDiferences_geneLevel.csv')) %>%
  rename(Gene_id = ...1)

arees_df <- full_join(arees_old, arees_new)

## Load Red Percentile Data
red_percent_old <- read_csv(paste0(microarrays_dir, 'Old_Arrays/R_results_OldArrays_Variantome/red_percentiles.csv'))
red_percent_new <- read_csv(paste0(microarrays_dir, 'New_Arrays/R_results_NewArray/red_percentiles.csv'))

red_df <- full_join(red_percent_old, red_percent_new, by = 'Gene_id')

## Load Max-Time data
old_maxtime <- read_csv(paste0(microarrays_dir, 'Old_Arrays/R_results_OldArrays_Variantome/old_arrays_maxtime.csv'))

new_maxtime <- read_csv(paste0(microarrays_dir, 'New_Arrays/R_results_NewArray/new_arrays_maxtime.csv'))

maxtime_df <- full_join(old_maxtime, new_maxtime, by = 'Gene_id', suffix = c('_Old', '_New'))

old_breaks <- read_csv(paste0(microarrays_dir, 'Old_Arrays/R_results_OldArrays_Variantome/old_area_breaks.csv'))

new_breaks <- read_csv(paste0(microarrays_dir, 'New_Arrays/R_results_NewArray/new_area_breaks.csv'))

breaks_df <- bind_cols(old_breaks, new_breaks)
colnames(breaks_df) <- c('Old_Area_Breaks', 'New_Area_Breaks')

#### Load Filtered FC Tables ####
old_arrays <- paste0(microarrays_dir, 'Old_Arrays/R_results_OldArrays_Variantome/')

filtered_12B_10G <- read_tsv(paste0(old_arrays, '12B_10G_final_df.tsv')) %>%
  filter(PassRed & PassDuplDel)
filtered_12B_3D7B <- read_tsv(paste0(old_arrays, '12B_3D7B_final_df.tsv')) %>%
  filter(PassRed & PassDuplDel)
filtered_10G_3D7B <- read_tsv(paste0(old_arrays, '10G_3D7B_final_df.tsv')) %>%
  filter(PassRed & PassDuplDel)

new_arrays <- paste0(microarrays_dir, 'New_Arrays/R_results_NewArray/')
filtered_A7_B11 <- read_tsv(paste0(new_arrays, 'A7_B11_final_df.tsv')) %>%
  filter(PassRed & PassDuplDel)
filtered_A7_E5 <- read_tsv(paste0(new_arrays, 'A7_E5_final_df.tsv')) %>%
  filter(PassRed & PassDuplDel)
filtered_B11_E5 <- read_tsv(paste0(new_arrays, 'B11_E5_final_df.tsv')) %>%
  filter(PassRed & PassDuplDel)

## Swich B11vsE5 for E5vsB11
filtered_E5_B11 <- filtered_B11_E5 %>%
  rename(`E5-B11_MaxVal` = `B11-E5_MaxVal`, `E5-B11_MaxTime` = `B11-E5_MaxTime`) %>%
  mutate(`E5-B11_MaxVal` = -`E5-B11_MaxVal`)


## Create filtered lists for each comparison, we need to add MaxTime and tRNA filters
difs_12B_10G <- filtered_12B_10G %>%
  filter(abs(`12B-10G_MaxVal`) > 2) %>%
  left_join(info_df %>% select(Gene_id, Is_tRNA)) %>%
  filter(PassAll) %>%
  filter(!Is_tRNA) %>%
  select(Gene_id) %>%
  pull()

difs_12B_3D7B <- filtered_12B_3D7B %>%
  filter(abs(`12B-3D7B_MaxVal`) > 2) %>%
  left_join(info_df %>% select(Gene_id, Is_tRNA)) %>%
  filter(PassAll) %>%
  filter(!Is_tRNA) %>%
  select(Gene_id) %>%
  pull()

difs_10G_3D7B <- filtered_10G_3D7B %>%
  filter(abs(`10G-3D7B_MaxVal`) > 2) %>%
  left_join(info_df %>% select(Gene_id, Is_tRNA)) %>%
  filter(PassAll) %>%
  filter(!Is_tRNA) %>%
  select(Gene_id) %>%
  pull()

difs_A7_E5 <- filtered_A7_E5 %>%
  filter(abs(`A7-E5_MaxVal`) > 2) %>%
  left_join(info_df %>% select(Gene_id, Is_tRNA)) %>%
  filter(PassAll) %>%
  filter(!Is_tRNA) %>%
  select(Gene_id) %>%
  pull()

difs_A7_B11 <- filtered_A7_B11 %>%
  filter(abs(`A7-B11_MaxVal`) > 2) %>%
  left_join(info_df %>% select(Gene_id, Is_tRNA)) %>%
  filter(PassAll) %>%
  filter(!Is_tRNA) %>%
  select(Gene_id) %>%
  pull()

difs_E5_B11 <- filtered_E5_B11 %>%
  filter(abs(`E5-B11_MaxVal`) > 2) %>%
  left_join(info_df %>% select(Gene_id, Is_tRNA)) %>%
  filter(PassAll) %>%
  filter(!Is_tRNA) %>%
  select(Gene_id) %>%
  pull()


all_difs <- unique(c(
  difs_12B_10G,
  difs_12B_3D7B,
  difs_10G_3D7B,
  difs_A7_E5,
  difs_A7_B11,
  difs_E5_B11
))

length(all_difs)

## Select MaxDif for each gene (once filtered by redfilter and dupl_del)
## It takes long to run so we read premade table, uncoment if need to rerun

maxDif_df <- tibble()
for (gid in trans_df$Gene_id){

  ## Create a list with each contrast difference df
  difs <- list(
    filtered_12B_10G,
    filtered_12B_3D7B,
    filtered_10G_3D7B,
    filtered_A7_B11,
    filtered_A7_E5,
    filtered_E5_B11
  )

  ## Filter by Gene_id
  dif_dfs <- lapply(difs, function(x) x %>% filter(Gene_id == gid))

  ## Function to get MaxVal of each df and join them in a vector
  get_MaxVal <- function(x){
    maxVal <- x %>%
      select(contains('MaxVal')) %>%
      pull()
    if (identical(maxVal, numeric(0))) {maxVal <- NA}
    return(maxVal)
  }

  maxVect <- sapply(dif_dfs, get_MaxVal)
  max_idx <- which.max(abs(maxVect))

  ## Handle genes that don't pass filters (set to NA)
  if (identical(maxVect, rep(NA, 6))) {
    out_row <- tibble(
      Gene_id = gid,
      Max_aMAFC = NA,
      Max_Time = NA,
      On_trans = NA,
      Off_trans = NA,
      PassMaxtime = FALSE
    )
    ## Handle rest of genes
  } else {
    ## Get on/off strain names
    maxDif <- dif_dfs[[max_idx]] %>% select(contains('_MaxVal'))
    maxVal <- maxDif %>% pull()

    maxStrains <- colnames(maxDif)
    strains <- strsplit(maxStrains, split = '_', fixed = T)[[1]][1]
    strains <- strsplit(strains, split = '-', fixed = T)[[1]]
    if (maxVal >= 0){
      On <- strains[1]
      Off <- strains[2]
    } else {
      On <- strains[2]
      Off <- strains[1]
    }

    ## Collect MaxDif row from the appropiate df
    out_row <- dif_dfs[[max_idx]] %>%
      select(Gene_id, contains('_MaxVal'), contains('_MaxTime'), PassMaxtime) %>%
      setNames(c('Gene_id', 'Max_aMAFC', 'Max_Time', 'PassMaxtime')) %>%
      mutate(
        On_trans = On,
        Off_trans = Off
      )

  }
  maxDif_df <- bind_rows(maxDif_df, out_row)
}

maxDif_df <- maxDif_df %>%
  arrange(-abs(Max_aMAFC)) %>%
  left_join(trans_df, by = 'Gene_id')

maxDif_df %>%
  write_tsv(paste0(tables_dir, 'max_differences_df.tsv'))

maxDif_df <- read_tsv(paste0(tables_dir, 'max_differences_df.tsv'))
maxDif_df_old <- read_tsv(paste0(tables_dir, 'max_differences_df_old.tsv'))

maxDif_df %>%
  print(width = 400)

maxDif_df %>%
  filter(abs(Max_aMAFC) > 2)

pre3D7_substitution <- maxDif_df
new_array_areas <- read_csv(paste0(microarrays_dir, 'New_Arrays/R_results_NewArray/area_geneLevel.csv')) %>%
  select(-Name, -Annot, -Variant)

max_trans_newarray <- maxDif_df %>%
  left_join(new_array_areas)

get_new_onoff <- function(gid) {

  #gid <- 'PF3D7_0601200'
  on <- max_trans_newarray %>%
    filter(Gene_id == gid) %>%
    select(On_trans) %>%
    pull()

  off <- max_trans_newarray %>%
    filter(Gene_id == gid) %>%
    select(Off_trans) %>%
    pull()

  time <- max_trans_newarray %>%
    filter(Gene_id == gid) %>%
    select(Max_Time) %>%
    pull()

  onoff <- TRUE
  if (is.na(on) | is.na(off)){
    onoff <- FALSE
  } else {
    if (on == '3D7B') {func <- which.max}
    if (off == '3D7B') {func <- which.min}
    if (on != '3D7B' & off != '3D7B') {onoff <- FALSE}
  }

  if (is.na(time)){
    out <- NA
  } else if (!onoff){
    out <- NA
  } else {
    strains <- c('A7', 'B11', 'E5')

    vect <- max_trans_newarray %>%
      filter(Gene_id == gid) %>%
      select(contains(time))

    ## Filter by red percent if 3D7B is 'On' strain
    red_pcnts <- red_df %>%
      filter(Gene_id == gid) %>%
      select(A7, B11, E5)

    red_mask <- red_pcnts > 15
    if (on != '3D7B'){red_mask <- c(TRUE, TRUE, TRUE)}

    ## Filter by dupl/del
    dupl_del_mask <- c(
      !gid %in% dupl_del$A7K9$Gene_id,
      !gid %in% dupl_del$B11$Gene_id,
      !gid %in% dupl_del$E5K9$Gene_id
    )

    ## Apply both filters
    whole_mask <- red_mask & dupl_del_mask

    vect <- vect[whole_mask]
    strains <- strains[whole_mask]

    ## Final output
    ifelse(all(is.na(vect)) | !any(whole_mask), out <- NA, out <- strains[func(vect)])
  }
  return(out)
}

newonoffs <- sapply(max_trans_newarray$Gene_id, get_new_onoff)
maxDif_df['New_OnOffs'] <- newonoffs
max_trans_newarray['New_OnOffs'] <- newonoffs

max_trans_newarray %>%
  select(Gene_id, On_trans, Off_trans, New_OnOffs) %>%
  print(n = 50)

maxDif_df <- maxDif_df %>%
  mutate(On_trans = ifelse(On_trans == '3D7B',
                           New_OnOffs,
                           On_trans)) %>%
  mutate(Off_trans = ifelse(Off_trans == '3D7B',
                           New_OnOffs,
                           Off_trans)) %>%
  mutate(Is_3D7B = !is.na(New_OnOffs)) %>%
  select(-New_OnOffs)

maxDif_df %>%
  select(-contains('_MaxVal'), -contains('_MaxTime'))

## Add tRNA info
maxDif_df <- maxDif_df %>%
  left_join(info_df %>% select(Gene_id, Is_tRNA), by = 'Gene_id')

maxDif_df %>%
  filter(Is_tRNA)
table(maxDif_df$Is_tRNA)

colnames(maxDif_df)

finalFC_df <- maxDif_df %>%
  filter(!is.na(On_trans) & !is.na(Off_trans)) %>%
  filter(!Is_tRNA) %>%
  filter(PassMaxtime)

## Missing Genes
msg <- all_difs[!all_difs %in% finalFC_df$Gene_id]

maxDif_df %>%
  filter(Gene_id %in% msg) %>%
  print(width = 999)

info_df %>%
  filter(Gene_id %in% msg) %>%
  print(width = 999)

write_csv(finalFC_df, paste0(tables_dir, 'final_filtered_transcription_differences.csv'))

colnames(info_df)

finalFC_df %>%
  filter(PassMaxtime) %>%
  filter(!Is_tRNA) %>%
  left_join(info_df, by = 'Gene_id') %>%
  filter(!Variant) %>%
  select(Gene_id, Name, Annot, Family, Gam_specific)

finalFC_df %>%
  left_join(info_df) %>%
  select(-PassMaxtime, -Is_3D7B, -Gene_name) %>%
  write_tsv(paste0(tables_dir, 'final_filtered_mAFC4_annotated.tsv'))

#### Load heterochromatin data 1000bp+500CDS ####

cov_12b <- read_tsv(paste0(
  coverage_dir,
  'binned_1000tss_500orf_coverage_1.2B.bed'
), col_names = F)
cov_10g <- read_tsv(paste0(
  coverage_dir,
  'binned_1000tss_500orf_coverage_10G.bed'
), col_names = F)
cov_a7 <- read_tsv(paste0(
  coverage_dir,
  'binned_1000tss_500orf_coverage_A7K9.bed'
), col_names = F)
cov_e5 <- read_tsv(paste0(
  coverage_dir,
  'binned_1000tss_500orf_coverage_E5K9.bed'
), col_names = F)
cov_b11 <- read_tsv(paste0(
  coverage_dir,
  'binned_1000tss_500orf_coverage_B11.bed'
), col_names = F)

join_cols = c('X1', 'X2', 'X3', 'X4')

het_df <- cov_12b %>%
  full_join(cov_10g, by = join_cols) %>%
  full_join(cov_a7, by = join_cols) %>%
  full_join(cov_e5, by = join_cols) %>%
  full_join(cov_b11, by = join_cols)

colnames(het_df) <- c('Chrom', 'Start', 'Stop', 'Gene_id',
                      'Het_12B', 'Het_10G', 'Het_A7', 'Het_E5', 'Het_B11')

het_df <- het_df %>% select(Gene_id, contains('Het'), everything())

hdf <- het_df %>% select(Gene_id, contains('Het'))
cnames <- str_replace(colnames(hdf), 'Het', 'Cov_1000fp_500orf')
colnames(hdf) <- cnames

#### Load heterochromatin data 5pORF3p ####
## Load 3prime, ORF, 5prime Coverage

load3ORF5 <- function(cov_5ORF3_file){
  cov_5ORF3 <- read_tsv(cov_5ORF3_file, col_names = F)
  cov_5ORF3 <- cov_5ORF3 %>%
    setNames(c('Chrom', 'Start', 'Stop', 'Gene_id', 'Intensity', 'Strand', 'Cov')) %>%
    select(-Intensity)
  prime5 <- cov_5ORF3 %>% filter(grepl('5prime', fixed = T, Gene_id))
  ORF <- cov_5ORF3 %>% filter(!grepl('5prime', fixed = T, Gene_id) & !grepl('3prime', fixed = T, Gene_id))
  prime3 <- cov_5ORF3 %>% filter(grepl('3prime', fixed = T, Gene_id))
  ORF['Cov_5prime'] <- prime5$Cov
  ORF['Cov_3prime'] <- prime3$Cov
  ORF <- ORF %>% rename(Cov_ORF = Cov)
  return(ORF)
}

e5_5ORF3 <- load3ORF5(paste0(coverage_dir, 'binned_5prime1000_ORF_3prime1000_coverage_E5K9.bed'))
a7_5ORF3 <- load3ORF5(paste0(coverage_dir, 'binned_5prime1000_ORF_3prime1000_coverage_A7K9.bed'))
b11_5ORF3 <- load3ORF5(paste0(coverage_dir, 'binned_5prime1000_ORF_3prime1000_coverage_B11.bed'))
x12b_5ORF3 <- load3ORF5(paste0(coverage_dir, 'binned_5prime1000_ORF_3prime1000_coverage_1.2B.bed'))
x10g_5ORF3 <- load3ORF5(paste0(coverage_dir, 'binned_5prime1000_ORF_3prime1000_coverage_10G.bed'))

join_cols = c('Chrom', 'Start', 'Stop', 'Gene_id', 'Strand')

cov_5ORF3_df <- x12b_5ORF3 %>%
  full_join(x10g_5ORF3, by = join_cols, suffix = c('_12B', '_10G')) %>%
  full_join(a7_5ORF3, by = join_cols, suffix = c('', '_A7')) %>%
  full_join(e5_5ORF3, by = join_cols, suffix = c('', '_E5')) %>%
  full_join(b11_5ORF3, by = join_cols, suffix = c('', '_B11')) %>%
  rename(Cov_ORF_A7 = Cov_ORF, Cov_5prime_A7 = Cov_5prime, Cov_3prime_A7 = Cov_3prime) %>%
  select(all_of(join_cols), contains('5prime'), contains('ORF'), contains('3prime'))

cov_5orf3 <- cov_5ORF3_df %>% select(Gene_id, contains('Cov'))
cnames <- str_replace(colnames(cov_5orf3), '5prime', '1000fp') %>%
  str_replace('ORF', 'allorf') %>%
  str_replace('3prime', '1000tp')
cov_5orf3 <- cov_5orf3 %>% set_names(cnames)

#### Load ORF abs-bp Coverages ####

read_abs_orf_cov <- function(strain, len){
  col_names <- c('Chrom', 'Start', 'Stop', 'Gene_id', 'Cov')
  path <- paste0(
    coverage_dir,
    'binned_0tss_', len,
    'orf_allowoverlaps_True_coverage_',
    strain, '.bed'
  )
  df <- read_tsv(path, col_names = col_names) %>%
    select(Gene_id, Cov)
  return(df)
}

cov_500_12B <- read_abs_orf_cov('1.2B', '500')
cov_500_10G <- read_abs_orf_cov('10G', '500')
cov_500_A7 <- read_abs_orf_cov('A7K9', '500')
cov_500_E5 <- read_abs_orf_cov('E5K9', '500')
cov_500_B11 <- read_abs_orf_cov('B11', '500')

cov_1000_12B <- read_abs_orf_cov('1.2B', '1000')
cov_1000_10G <- read_abs_orf_cov('10G', '1000')
cov_1000_A7 <- read_abs_orf_cov('A7K9', '1000')
cov_1000_E5 <- read_abs_orf_cov('E5K9', '1000')
cov_1000_B11 <- read_abs_orf_cov('B11', '1000')

cov_1500_12B <- read_abs_orf_cov('1.2B', '1500')
cov_1500_10G <- read_abs_orf_cov('10G', '1500')
cov_1500_A7 <- read_abs_orf_cov('A7K9', '1500')
cov_1500_E5 <- read_abs_orf_cov('E5K9', '1500')
cov_1500_B11 <- read_abs_orf_cov('B11', '1500')

cov_orf_abs_df <- cov_500_12B %>%
  full_join(cov_500_10G, by = 'Gene_id') %>%
  full_join(cov_500_A7, by = 'Gene_id') %>%
  full_join(cov_500_E5, by = 'Gene_id') %>%
  full_join(cov_500_B11, by = 'Gene_id') %>%

  full_join(cov_1000_12B, by = 'Gene_id') %>%
  full_join(cov_1000_10G, by = 'Gene_id') %>%
  full_join(cov_1000_A7, by = 'Gene_id') %>%
  full_join(cov_1000_E5, by = 'Gene_id') %>%
  full_join(cov_1000_B11, by = 'Gene_id') %>%

  full_join(cov_1500_12B, by = 'Gene_id') %>%
  full_join(cov_1500_10G, by = 'Gene_id') %>%
  full_join(cov_1500_A7, by = 'Gene_id') %>%
  full_join(cov_1500_E5, by = 'Gene_id') %>%
  full_join(cov_1500_B11, by = 'Gene_id') %>%

  set_names('Gene_id',
            'Cov_500orf_12B',
            'Cov_500orf_10G',
            'Cov_500orf_A7',
            'Cov_500orf_E5',
            'Cov_500orf_B11',
            'Cov_1000orf_12B',
            'Cov_1000orf_10G',
            'Cov_1000orf_A7',
            'Cov_1000orf_E5',
            'Cov_1000orf_B11',
            'Cov_1500orf_12B',
            'Cov_1500orf_10G',
            'Cov_1500orf_A7',
            'Cov_1500orf_E5',
            'Cov_1500orf_B11'
            )

cov_orf_abs_df %>%
  filter(!complete.cases(.))


cov_orf_abs_df %>%
  select(contains('12B'))

cov_orf_abs_df %>%
  select(contains('10G'))

#### Load plasmoDB TSS coverages ####
load_pDB <- function(cov_pDB_file){

  col_names <- c('Chrom', 'Start', 'Stop', 'Gene_id', 'Type', 'Cov')
  cov_pDB <- read_tsv(cov_pDB_file, col_names = col_names) %>%
    mutate(Gene_id = str_replace(Gene_id, '.*(PF3D7_\\d{7}).*', '\\1'))

  prime5 <- cov_pDB %>%
    filter(Type == 'five_prime_UTR') %>%
    rename(p5utr_Cov = Cov) %>%
    select(Gene_id, p5utr_Cov)

  prime3 <- cov_pDB %>%
    filter(Type == 'three_prime_UTR') %>%
    rename(p3utr_Cov = Cov) %>%
    select(Gene_id, p3utr_Cov)

  outdf <- prime5 %>%
    full_join(prime3, by = 'Gene_id') %>%
    group_by(Gene_id) %>%
    summarize(p5utr_Cov = mean(p5utr_Cov), p3utr_Cov = mean(p3utr_Cov))

  return(outdf)
}

cov_plasmoDB_12B <- load_pDB(paste0(coverage_dir, 'plasmoDB_UTRs_sorted_coverage_1.2B.bed'))
cov_plasmoDB_10G <- load_pDB(paste0(coverage_dir, 'plasmoDB_UTRs_sorted_coverage_10G.bed'))
cov_plasmoDB_A7 <- load_pDB(paste0(coverage_dir, 'plasmoDB_UTRs_sorted_coverage_A7K9.bed'))
cov_plasmoDB_E5 <- load_pDB(paste0(coverage_dir, 'plasmoDB_UTRs_sorted_coverage_E5K9.bed'))
cov_plasmoDB_B11 <- load_pDB(paste0(coverage_dir, 'plasmoDB_UTRs_sorted_coverage_B11.bed'))

cov_pDB_df <- cov_plasmoDB_12B %>%
  full_join(cov_plasmoDB_10G, by = 'Gene_id', suffix = c('_12B', '_10G')) %>%
  full_join(cov_plasmoDB_A7, by = 'Gene_id', suffix = c('', '_A7')) %>%
  full_join(cov_plasmoDB_E5, by = 'Gene_id', suffix = c('', '_E5')) %>%
  full_join(cov_plasmoDB_B11, by = 'Gene_id', suffix = c('', '_B11')) %>%
  rename(p5utr_Cov_A7 = p5utr_Cov, p3utr_Cov_A7 = p3utr_Cov)


load_pDB_TSS <- function(cov_pDB_file){
  col_names <- c('Chrom', 'Start', 'Stop', 'Gene_id', 'TSS_Cov')
  cov_pDB <- read_tsv(cov_pDB_file, col_names = col_names) %>%
    mutate(Gene_id = str_replace(Gene_id, '.*(PF3D7_\\d{7}).*', '\\1')) %>%
    select(Gene_id, TSS_Cov) %>%
    group_by(Gene_id) %>%
    summarize(TSS_Cov = mean(TSS_Cov))

  return(cov_pDB)
}

cov_plasmoDB_TSS_12B <- load_pDB_TSS(paste0(coverage_dir, 'plasmoDB_TSSs_sorted_coverage_1.2B.bed'))
cov_plasmoDB_TSS_10G <- load_pDB_TSS(paste0(coverage_dir, 'plasmoDB_TSSs_sorted_coverage_10G.bed'))
cov_plasmoDB_TSS_A7 <- load_pDB_TSS(paste0(coverage_dir, 'plasmoDB_TSSs_sorted_coverage_A7K9.bed'))
cov_plasmoDB_TSS_E5 <- load_pDB_TSS(paste0(coverage_dir, 'plasmoDB_TSSs_sorted_coverage_E5K9.bed'))
cov_plasmoDB_TSS_B11 <- load_pDB_TSS(paste0(coverage_dir, 'plasmoDB_TSSs_sorted_coverage_B11.bed'))

cov_pDB_TSS_df <- cov_plasmoDB_TSS_12B %>%
  full_join(cov_plasmoDB_TSS_10G, by = 'Gene_id', suffix = c('_12B', '_10G')) %>%
  full_join(cov_plasmoDB_TSS_A7, by = 'Gene_id', suffix = c('', '_A7')) %>%
  full_join(cov_plasmoDB_TSS_E5, by = 'Gene_id', suffix = c('', '_E5')) %>%
  full_join(cov_plasmoDB_TSS_B11, by = 'Gene_id', suffix = c('', '_B11')) %>%
  rename(TSS_Cov_A7 = TSS_Cov)

cov_pDB_df <- cov_pDB_df %>%
  full_join(cov_pDB_TSS_df)

#### Load Many-Bin coverages ####

get_bins <- function(dfin, bin_txt, outcol){
  dfout <- dfin %>% filter(grepl(bin_txt, fixed = T, Gene_id)) %>%
    mutate(Gene_id = str_replace(Gene_id, bin_txt, '')) %>%
    rename(!!outcol := Cov) %>%
    select(Gene_id, Strand, !!outcol)
}

get_manybins <- function(strain){

  fpath <- paste0(
    coverage_dir,
    'binned_5prime500_1000_1500_2000_ORF_3prime500_1000_1500_allowoverlapTrue_coverage_',
    strain, '.bed'
  )

  cov_5orf3_manybins <- read_tsv(fpath, col_names = F)%>%
    setNames(c('Chrom', 'Start', 'Stop', 'Gene_id', 'Intensity', 'Strand', 'Cov')) %>%
    select(-Intensity)

  cov_2000fp_manybins <- get_bins(cov_5orf3_manybins, '_5prime_2000', 'Cov_2000fp_manybins')
  cov_1500fp_manybins <- get_bins(cov_5orf3_manybins, '_5prime_1500', 'Cov_1500fp_manybins')
  cov_1000fp_manybins <- get_bins(cov_5orf3_manybins, '_5prime_1000', 'Cov_1000fp_manybins')
  cov_500fp_manybins <- get_bins(cov_5orf3_manybins, '_5prime_500', 'Cov_500fp_manybins')
  cov_1qorf_manybins <- get_bins(cov_5orf3_manybins, '_q1', 'Cov_1qorf_manybins')
  cov_2qorf_manybins <- get_bins(cov_5orf3_manybins, '_q2', 'Cov_2qorf_manybins')
  cov_3qorf_manybins <- get_bins(cov_5orf3_manybins, '_q3', 'Cov_3qorf_manybins')
  cov_4qorf_manybins <- get_bins(cov_5orf3_manybins, '_q4', 'Cov_4qorf_manybins')
  cov_1500tp_manybins <- get_bins(cov_5orf3_manybins, '_3prime_1500', 'Cov_1500tp_manybins')
  cov_1000tp_manybins <- get_bins(cov_5orf3_manybins, '_3prime_1000', 'Cov_1000tp_manybins')
  cov_500tp_manybins <- get_bins(cov_5orf3_manybins, '_3prime_500', 'Cov_500tp_manybins')

  join_cols <- c('Strand', 'Gene_id')
  df <- cov_2000fp_manybins %>%
    full_join(cov_1500fp_manybins, by = join_cols) %>%
    full_join(cov_1000fp_manybins, by = join_cols) %>%
    full_join(cov_500fp_manybins, by = join_cols) %>%
    full_join(cov_1qorf_manybins, by = join_cols) %>%
    full_join(cov_2qorf_manybins, by = join_cols) %>%
    full_join(cov_3qorf_manybins, by = join_cols) %>%
    full_join(cov_4qorf_manybins, by = join_cols) %>%
    full_join(cov_1500tp_manybins, by = join_cols) %>%
    full_join(cov_1000tp_manybins, by = join_cols) %>%
    full_join(cov_500tp_manybins, by = join_cols)

  return(df)
}

cov_manybins_12B <- get_manybins('1.2B')
cov_manybins_10G <- get_manybins('10G')
cov_manybins_A7 <- get_manybins('A7K9')
cov_manybins_E5 <- get_manybins('E5K9')
cov_manybins_B11 <- get_manybins('B11')


## cov_manybins_B11 %>%
##   filter(Gene_id == 'PF3D7_1130800') %>%
##   print(width = 800)


join_cols <- c('Gene_id', 'Strand')
cov_manybins <- cov_manybins_12B %>%
  full_join(cov_manybins_10G, by = join_cols, suffix = c('_12B', '_10G')) %>%
  full_join(cov_manybins_A7, by = join_cols, suffix = c('', '_A7')) %>%
  full_join(cov_manybins_E5, by = join_cols, suffix = c('', '_E5')) %>%
  full_join(cov_manybins_B11, by = join_cols, suffix = c('', '_B11')) %>%
  rename(Cov_2000fp_manybins_A7 = Cov_2000fp_manybins,
         Cov_1500fp_manybins_A7 = Cov_1500fp_manybins,
         Cov_1000fp_manybins_A7 = Cov_1000fp_manybins,
         Cov_500fp_manybins_A7 = Cov_500fp_manybins,
         Cov_1qorf_manybins_A7 = Cov_1qorf_manybins,
         Cov_2qorf_manybins_A7 = Cov_2qorf_manybins,
         Cov_3qorf_manybins_A7 = Cov_3qorf_manybins,
         Cov_4qorf_manybins_A7 = Cov_4qorf_manybins,
         Cov_1500tp_manybins_A7 = Cov_1500tp_manybins,
         Cov_1000tp_manybins_A7 = Cov_1000tp_manybins,
         Cov_500tp_manybins_A7 = Cov_500tp_manybins,
         ) %>%
  select(-Strand)

## Convert to numeric
cov_manybins <- cov_manybins %>%
  mutate(across(contains('Cov'), as.numeric))

cov_manybins %>%
  filter(!complete.cases(.))

#### Load Abs ORF Coverages ####
load_abs_ORF <- function(cov_abs_orf_file){
  cov_abs_ORF <- read_tsv(cov_abs_orf_file, col_names = F)
  cov_abs_ORF <- cov_abs_ORF %>%
    setNames(c('Chrom', 'Start', 'Stop', 'Gene_id', 'Intensity', 'Strand', 'Cov')) %>%
    select(-Intensity)
  ORF <- cov_abs_ORF %>%
    filter(!grepl('5prime', fixed = T, Gene_id) & !grepl('3prime', fixed = T, Gene_id))
  ORF <- ORF %>%
    rename(Cov_ORF = Cov) %>%
    mutate(Cov_ORF = as.double(Cov_ORF)) %>%
    select(Gene_id, Cov_ORF)
  return(ORF)
}

abs_ORF_500_12B <- load_abs_ORF(paste0(
  coverage_dir,
  'binned_5prime1000_ORF_0_500_3prime1000_allowoverlapTrue_coverage_1.2B.bed'))
abs_ORF_500_10G <- load_abs_ORF(paste0(
  coverage_dir,
  'binned_5prime1000_ORF_0_500_3prime1000_allowoverlapTrue_coverage_10G.bed'))
abs_ORF_500_A7 <- load_abs_ORF(paste0(
  coverage_dir,
  'binned_5prime1000_ORF_0_500_3prime1000_allowoverlapTrue_coverage_A7K9.bed'))
abs_ORF_500_E5 <- load_abs_ORF(paste0(
  coverage_dir,
  'binned_5prime1000_ORF_0_500_3prime1000_allowoverlapTrue_coverage_E5K9.bed'))
abs_ORF_500_B11 <- load_abs_ORF(paste0(
  coverage_dir,
  'binned_5prime1000_ORF_0_500_3prime1000_allowoverlapTrue_coverage_B11.bed'))

abs_ORF_1000_12B <- load_abs_ORF(paste0(
  coverage_dir,
  'binned_5prime1000_ORF_500_1000_3prime1000_allowoverlapTrue_coverage_1.2B.bed'))
abs_ORF_1000_10G <- load_abs_ORF(paste0(
  coverage_dir,
  'binned_5prime1000_ORF_500_1000_3prime1000_allowoverlapTrue_coverage_10G.bed'))
abs_ORF_1000_A7 <- load_abs_ORF(paste0(
  coverage_dir,
  'binned_5prime1000_ORF_500_1000_3prime1000_allowoverlapTrue_coverage_A7K9.bed'))
abs_ORF_1000_E5 <- load_abs_ORF(paste0(
  coverage_dir,
  'binned_5prime1000_ORF_500_1000_3prime1000_allowoverlapTrue_coverage_E5K9.bed'))
abs_ORF_1000_B11 <- load_abs_ORF(paste0(
  coverage_dir,
  'binned_5prime1000_ORF_500_1000_3prime1000_allowoverlapTrue_coverage_B11.bed'))

abs_ORF_1500_12B <- load_abs_ORF(paste0(
  coverage_dir,
  'binned_5prime1000_ORF_1000_1500_3prime1000_allowoverlapTrue_coverage_1.2B.bed'))
abs_ORF_1500_10G <- load_abs_ORF(paste0(
  coverage_dir,
  'binned_5prime1000_ORF_1000_1500_3prime1000_allowoverlapTrue_coverage_10G.bed'))
abs_ORF_1500_A7 <- load_abs_ORF(paste0(
  coverage_dir,
  'binned_5prime1000_ORF_1000_1500_3prime1000_allowoverlapTrue_coverage_A7K9.bed'))
abs_ORF_1500_E5 <- load_abs_ORF(paste0(
  coverage_dir,
  'binned_5prime1000_ORF_1000_1500_3prime1000_allowoverlapTrue_coverage_E5K9.bed'))
abs_ORF_1500_B11 <- load_abs_ORF(paste0(
  coverage_dir,
  'binned_5prime1000_ORF_1000_1500_3prime1000_allowoverlapTrue_coverage_B11.bed'))

join_cols = c('Gene_id')

cov_abs_ORF <- abs_ORF_500_12B %>%
  full_join(abs_ORF_500_10G, by = join_cols, suffix = c('_0_500_12B', '_0_500_10G')) %>%
  full_join(abs_ORF_500_A7, by = join_cols, suffix = c('', '_0_500_A7')) %>%
  full_join(abs_ORF_500_E5, by = join_cols, suffix = c('', '_0_500_E5')) %>%
  full_join(abs_ORF_500_B11, by = join_cols, suffix = c('', '_0_500_B11')) %>%

  full_join(abs_ORF_1000_12B, by = join_cols, suffix = c('', '_500_1000_12B')) %>%
  full_join(abs_ORF_1000_10G, by = join_cols, suffix = c('', '_500_1000_10G')) %>%
  full_join(abs_ORF_1000_A7, by = join_cols, suffix = c('', '_500_1000_A7')) %>%
  full_join(abs_ORF_1000_E5, by = join_cols, suffix = c('', '_500_1000_E5')) %>%
  full_join(abs_ORF_1000_B11, by = join_cols, suffix = c('', '_500_1000_B11')) %>%

  full_join(abs_ORF_1500_12B, by = join_cols, suffix = c('', '_1000_1500_12B')) %>%
  full_join(abs_ORF_1500_10G, by = join_cols, suffix = c('', '_1000_1500_10G')) %>%
  full_join(abs_ORF_1500_A7, by = join_cols, suffix = c('', '_1000_1500_A7')) %>%
  full_join(abs_ORF_1500_E5, by = join_cols, suffix = c('', '_1000_1500_E5')) %>%
  full_join(abs_ORF_1500_B11, by = join_cols, suffix = c('', '_1000_1500_B11')) %>%

  rename(Cov_ORF_0_500_A7 = Cov_ORF) %>%
  select(Gene_id, contains('Cov_'))

#### Join all Coverages ####
cor_cov_df <- cov_orf_abs_df %>%
  full_join(hdf, by = 'Gene_id') %>%
  full_join(cov_5orf3, by = 'Gene_id') %>%
  full_join(cov_pDB_df, by = 'Gene_id') %>%
  full_join(cov_manybins, by = 'Gene_id') %>%
  full_join(cov_abs_ORF, by = 'Gene_id')

positive_cov_df <- cor_cov_df %>%
  mutate(across(-Gene_id, ~ ifelse(.x > 0, .x, 0)))

full_df <- inner_join(cor_cov_df, trans_df)  %>%
  left_join(info_df)

#### Create 'fixed' color palette ####
col_factor <- factor(unique(info_df$Family_Grouped),
                     levels=c(
                       'Not CVGs',
                       'Other CVGs',
                       'FIKK',
                       'HYP',
                       'PHIST',
                       'STEVOR',
                       'RIFIN',
                       'VAR',
                       'PFMC-2TM'
                     ))

my_colors <- c('gray', scales::brewer_pal(palette = 'Set3', type = 'qual')(9))
names(my_colors) <- levels(col_factor)
my_scale <- scale_fill_manual(name = "Family_Grouped", values = my_colors)
my_scale2 <- scale_color_manual(name = "Family_Grouped", values = my_colors)

#### Coverage PCAs ####
colnames(full_df %>% select(contains('Cov')))

str_selector <- 'allorf'

cov_df <- full_df %>%
  select(Gene_id, contains(str_selector)) %>%
  filter(complete.cases(.))

noNA_df <- cov_df %>%
  select(contains('Cov'))

pca <- prcomp(t(noNA_df))
cmp1 <- format(round(summary(pca)$importance[2,1]*100, 2), nsmall = 2)
cmp2 <- format(round(summary(pca)$importance[2,2]*100, 2), nsmall = 2)

df_pca <- as.data.frame(pca$x)
df_pca$Strain <- sapply(colnames(cov_df %>% select(contains('Cov'))),
       function(x) gsub('Cov_1000fp_500orf_', '', x))

p <- ggplot(df_pca, aes(x=PC1,y=PC2, col = Strain, group = Strain))
p <- p + geom_point(aes(), size = 3)
p <- p + scale_x_continuous(name=paste0("PC1: ", cmp1, "%"))
p <- p + scale_y_continuous(name=paste0("PC2: ", cmp2, "%"))
p <- p + theme_classic()
p <- p + theme(text = element_text(size=20))
p <- p + geom_label_repel(aes(label = Strain))
p <- p + theme(legend.position = "none")
#p

ggsave(p, filename = paste0(plots_dir, 'Coverage_PCAs/', "coverage_", str_selector, '_PCA.png'), device = "png")

#### FilteredFC On/Off Donut Plots ####
on_off_donut <- function(col){
  data <- finalFC_df %>%
    count(get(col)) %>%
    set_names('category', 'count')

  ## Compute percentages
  data$fraction = data$count / sum(data$count)

  ## Compute the cumulative percentages (top of each rectangle)
  data$ymax = cumsum(data$fraction)

  ## Compute the bottom of each rectangle
  data$ymin = c(0, head(data$ymax, n=-1))

  ## Compute label position
  data$labelPosition <- (data$ymax + data$ymin) / 2

  ## Compute a good label
  data$label <- paste0(data$category, "\n value: ", data$count)

  ## Make the plot
  p <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category))
  p <- p +  geom_rect(color="white")
  ##p <- p + geom_label( x=4.2, aes(y=labelPosition, label=label), size=3)
  p <- p + coord_polar(theta="y") # Try to remove that to understand how the chart is built initially
  p <- p + xlim(c(2, 4)) # Try to remove that to see how to make a pie chart
  p <- p + theme_void()
  p <- p + scale_fill_viridis(discrete=TRUE)
  p
}


ggsave(paste0(plots_dir, 'on_difgenes_donut.pdf'), on_off_donut('On_trans'), device = 'pdf')
ggsave(paste0(plots_dir, 'off_difgenes_donut.pdf'), on_off_donut('Off_trans'), device = 'pdf')

#### Correlations By Strain ####

get_cor_mtx <- function(df, contrast, trans_filter, difpeaks_filter){

  ##contrast <- c('12B', '10G')
  ##df <- full_df

  trans_col <- paste0(contrast[1], '-', contrast[2], '_MaxVal')
  difpeaks_col <- paste0('Difpeaks_', contrast[1], '_', contrast[2])


  ## df %>%
  ##   select(contains('Cov_ORF_'))

  df_dif <- df %>%
  ## Absolute bp ORF
    mutate(Cov_500orf_Dif = get(paste0('Cov_500orf_', contrast[1])) - get(paste0('Cov_500orf_', contrast[2]))) %>%
    mutate(Cov_1000orf_Dif = get(paste0('Cov_1000orf_', contrast[1])) - get(paste0('Cov_1000orf_', contrast[2]))) %>%
    mutate(Cov_1500orf_Dif = get(paste0('Cov_1500orf_', contrast[1])) - get(paste0('Cov_1500orf_', contrast[2]))) %>%

  ## Absolute bp ORF, by interval
    mutate(Cov_0_500orf_Dif = get(paste0('Cov_ORF_0_500_', contrast[1])) - get(paste0('Cov_ORF_0_500_', contrast[2]))) %>%
    mutate(Cov_500_1000orf_Dif = get(paste0('Cov_ORF_500_1000_', contrast[1])) - get(paste0('Cov_ORF_500_1000_', contrast[2]))) %>%
    mutate(Cov_1000_1500orf_Dif = get(paste0('Cov_ORF_1000_1500_', contrast[1])) - get(paste0('Cov_ORF_1000_1500_', contrast[2]))) %>%

  ## 1000bp 5prime + 500bp ORF
    mutate(Cov_1000fp_500orf_Dif = get(paste0('Cov_1000fp_500orf_', contrast[1])) - get(paste0('Cov_1000fp_500orf_', contrast[2]))) %>%

  ## 1000bp 5prime/ ORF / 1000bp 3prime
    mutate(Cov_1000fp_Dif = get(paste0('Cov_1000fp_', contrast[1])) - get(paste0('Cov_1000fp_', contrast[2]))) %>%
    mutate(Cov_allorf_Dif = get(paste0('Cov_allorf_', contrast[1])) - get(paste0('Cov_allorf_', contrast[2]))) %>%
    mutate(Cov_1000tp_Dif = get(paste0('Cov_1000tp_', contrast[1])) - get(paste0('Cov_1000tp_', contrast[2]))) %>%

  ## PlsmoDB 5prime and TSS
    mutate(Cov_plasmoDB_5p_Dif = get(paste0('p5utr_Cov_', contrast[1])) - get(paste0('p5utr_Cov_', contrast[2]))) %>%
    mutate(Cov_plasmoDB_TSS_Dif = get(paste0('TSS_Cov_', contrast[1])) - get(paste0('TSS_Cov_', contrast[2]))) %>%
    mutate(Cov_plasmoDB_3p_Dif = get(paste0('p3utr_Cov_', contrast[1])) - get(paste0('p3utr_Cov_', contrast[2]))) %>%

  ## Manybins (2000, 1500, 1000, 500, 1q, 2q, 3q, 4q, 500, 1000, 1500)
    mutate(Cov_2000fp_manybins_Dif = get(paste0('Cov_2000fp_manybins_', contrast[1])) - get(paste0('Cov_2000fp_manybins_', contrast[2]))) %>%
    mutate(Cov_1500fp_manybins_Dif = get(paste0('Cov_1500fp_manybins_', contrast[1]))  - get(paste0('Cov_1500fp_manybins_', contrast[2]))) %>%
    mutate(Cov_1000fp_manybins_Dif = get(paste0('Cov_1000fp_manybins_', contrast[1])) - get(paste0('Cov_1000fp_manybins_', contrast[2]))) %>%
    mutate(Cov_500fp_manybins_Dif = get(paste0('Cov_500fp_manybins_', contrast[1]))- get(paste0('Cov_500fp_manybins_', contrast[2]))) %>%

    mutate(Cov_1qorf_manybins_Dif = get(paste0('Cov_1qorf_manybins_', contrast[1])) - get(paste0('Cov_1qorf_manybins_', contrast[2]))) %>%
    mutate(Cov_2qorf_manybins_Dif = get(paste0('Cov_2qorf_manybins_', contrast[1]))  - get(paste0('Cov_2qorf_manybins_', contrast[2]))) %>%
    mutate(Cov_3qorf_manybins_Dif = get(paste0('Cov_3qorf_manybins_', contrast[1])) - get(paste0('Cov_3qorf_manybins_', contrast[2]))) %>%
    mutate(Cov_4qorf_manybins_Dif = get(paste0('Cov_4qorf_manybins_', contrast[1]))- get(paste0('Cov_4qorf_manybins_', contrast[2]))) %>%

    mutate(Cov_1500tp_manybins_Dif = get(paste0('Cov_1500tp_manybins_', contrast[1]))  - get(paste0('Cov_1500tp_manybins_', contrast[2]))) %>%
    mutate(Cov_1000tp_manybins_Dif = get(paste0('Cov_1000tp_manybins_', contrast[1])) - get(paste0('Cov_1000tp_manybins_', contrast[2]))) %>%
    mutate(Cov_500tp_manybins_Dif = get(paste0('Cov_500tp_manybins_', contrast[1]))- get(paste0('Cov_500tp_manybins_', contrast[2]))) %>%

    filter(abs(get(trans_col)) > trans_filter) %>%
    select(Gene_id, matches(trans_col), contains('_Dif'))

  ## df_dif %>%
  ##   select(Gene_id, Cov_0_500orf_Dif, Cov_500_1000orf_Dif, Cov_1000_1500orf_Dif) %>%
  ##   summary()

  ## table(df_dif$Cov_500_1000orf_Dif)

  if (difpeaks_filter){df_dif <- df_dif %>% filter(get(difpeaks_col))}

  cormtx <- cor(df_dif %>% select(-Gene_id) %>% drop_na(), method = 'pearson')
  print(paste0('Contrast ', contrast[1], ' vs ', contrast[2]))
  print(paste0('Number of genes: ', dim(df_dif %>% select(-Gene_id) %>% drop_na())[1]))
  print(as.data.frame(cormtx[,1]))
  print('+++++++++++++++++++++++++++++')

  outdf <- as.data.frame(cormtx)
  outdf['Corr_with'] <- rownames(outdf)
  outdf <- outdf %>% select(Corr_with, everything())
  outname = paste0('corr_ALL_', contrast[1], '_', contrast[2],
                   '_transth_', trans_filter, '_difpeaks_', difpeaks_filter, '.csv')
  write_csv(outdf, file = outname)
}

contrasts <- list(
  c('12B', '10G'),
  c('A7', 'E5'),
  c('A7', 'B11'),
  c('E5', 'B11')
)

trans_filter <- 0
difpeaks_filter <- F

for (contrast in contrasts) {get_cor_mtx(full_df, contrast, trans_filter, difpeaks_filter)}

cor_12B_10G <- get_cor_mtx(full_df, c('12B', '10G'), trans_filter, difpeaks_filter)
cor_A7_E5 <- get_cor_mtx(full_df, c('A7', 'E5'), trans_filter, difpeaks_filter)
cor_A7_B11 <- get_cor_mtx(full_df, c('A7', 'B11'), trans_filter, difpeaks_filter)
cor_B11_E5 <- get_cor_mtx(full_df, c('E5', 'B11'), trans_filter, difpeaks_filter)

all_Corr <- bind_cols(Corr_With = cor_12B_10G[-1,1],
                      Cor_12B_10G = cor_12B_10G[-1,2],
                      Cor_A7_E5 = cor_A7_E5[-1,2],
                      Cor_A7_B11 = cor_A7_B11[-1,2],
                      Cor_B11_E5 = cor_B11_E5[-1,2])

write_csv(all_Corr, file = paste0(tables_dir, 'correlations_transth_', trans_filter, '_allowoverlaps.csv'))

#### Correlation By Strain Plots ####

myCorPlot <- function(df, region, s1, s2, trans_th){

  df <- plot_5ORF3_df

  trans <- paste0(s1,'-',s2,'_MaxVal')
  difpeak <- paste0('Difpeaks_', s1, '_', s2)
  title_str <- paste0(s1, ' vs ', s2)
  outfile <- paste0('./Plots/corplot_', s1, '_', s2, '_', region, '.png')

  plot_df <- df %>%
    mutate(Dif_het = get(paste0('Cov_', region, '_', s1)) - get(paste0('Cov_', region, '_', s2))) %>%
    select(Gene_id, matches(trans), Dif_het, matches(difpeak)) %>%
    setNames(c('Gene_id', 'Trans', 'Het', 'Difpeak')) %>%
    filter(abs(Trans) > trans_th)

  p <- ggplot(plot_df, aes(x=Het,
                           y=Trans)) +
    scale_alpha_discrete(range = c(0.2, 1)) +
    ggtitle(title_str) +
    ylab('aMAFC') +
    xlab(paste0('Heterochromatin difference at ', region)) +
    geom_point()
  ggsave(outfile, p, device = 'png')
  print(p)
}

plot_5ORF3_df <- inner_join(cov_5ORF3_df, trans_df)  %>%
  left_join(info_df)

regions <- c('5prime', 'ORF', '3prime')

contrasts <- list(
  c('12B', '10G'),
  c('A7', 'E5'),
  c('A7', 'B11'),
  c('E5', 'B11')
)

## ## Call corr plots for 5p/ORF/3p coverage
## for (region in regions){
##   for (contrast in contrasts){
##     s1 <- contrast[1]
##     s2 <- contrast[2]
##     trans_th <- 0
##     myCorPlot(plot_5ORF3_df, region, s1, s2, trans_th)
##   }
## }

region <- 'ORF'
s1 <- 'A7'
s2 <- 'B11'
trans_th <- 1
myCorPlot(plot_5ORF3_df, region, s1, s2, trans_th)

#### Correlation By Strain Plots ####

myCorPlot_all <- function(df, region, s1, s2, trans_th){

  trans <- paste0(s1,'-',s2,'_MaxVal')
  difpeak <- paste0('peakFC_', s1, '_', s2)
  title_str <- paste0(s1, ' vs ', s2)
  outfile <- paste0(
    plots_dir,
    'Correlations_byStrain/corplot_',
    s1, '_', s2, '_', region,
    'transth_', trans_th,
    '.png'
  )

  plot_df <- df %>%
    mutate(Dif_het = get(paste0(region, s1)) - get(paste0(region, s2))) %>%
    select(Gene_id, matches(trans), Dif_het, matches(difpeak)) %>%
    setNames(c('Gene_id', 'Trans', 'Het', 'Difpeak')) %>%
    mutate(Difpeak = !is.na(Difpeak)) %>%
    filter(abs(Trans) > trans_th)

  p <- ggplot(plot_df, aes(x=Het,
                      y=Trans,
                      color = Difpeak,
                      alpha = Difpeak)) +
    scale_alpha_discrete(range = c(0.2, 1)) +
    ggtitle(title_str) +
    ylab('aMAFC') +
    xlab(paste0('Heterochromatin difference at ', region)) +
    geom_point()
  ggsave(outfile, p, device = 'png')
  print(p)
}

## Join transcription, coverage and difpeaks
trans_het_difpeaks_df <- full_df %>%
  left_join(sicer_df, by = 'Gene_id')


## Make plots for all strains and coverages
cnms <- str_replace(colnames(full_df %>% select(contains('Cov') & contains('12B'))), '12B', '')

contrasts <- list(
  c('12B', '10G'),
  c('A7', 'E5'),
  c('A7', 'B11'),
  c('E5', 'B11')
)

## Plot correlations for all Coverages
for (c in cnms){
  for (contrast in contrasts){
    s1 <- contrast[1]
    s2 <- contrast[2]
    trans_th <- 1
    myCorPlot_all(trans_het_difpeaks_df, c, s1, s2, trans_th)
  }
}


region <- 'Cov_allorf_'
s1 <- 'A7'
s2 <- 'E5'
trans_th <- 1
myCorPlot_all(trans_het_difpeaks_df, region, s1, s2, trans_th)

#### Heatmap Funtions ####

myHeatmap <- function(df, family_facet){
  p <- ggplot(df, aes(x = variable, y = Label, fill = value)) +
    geom_tile(colour="snow3") +
    #geom_tile() +
    theme(
      text=element_text(size=24),

      legend.position='bottom',
      legend.title = element_blank(),

      panel.background=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),

      axis.title = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )
  if (family_facet){p <- p+facet_grid(Family ~., scales = "free_y", space = "free")}
  return(p)
}
hetHeatmap <- function(df, family_facet){
  p <- myHeatmap(df, family_facet)
  p <- p + scale_fill_gradient(low = "white",
                               high = "orange",
                               na.value="grey",
                               limits = c(0,4),
                               oob=squish) +
    scale_y_discrete(position = "right") +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),

      panel.border=element_blank(),
      panel.grid.major=element_blank(),

      strip.background = element_blank(),
      strip.text.x = element_blank(),
      strip.text.y = element_blank(),
    )
  return(p)
}
hetDifHeatmap <- function(df, family_facet){
  p <- myHeatmap(df, family_facet)
  p <- p + scale_fill_gradient2(low = "chartreuse3",
                                mid = "white",
                                high = "darkred",
                                na.value="grey",
                                limits = c(-4, 4),
                                oob=squish) +
    scale_y_discrete(position = "left") +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),

      panel.border=element_blank(),
      panel.grid.major=element_blank(),

      strip.background = element_blank(),
      strip.text.x = element_blank(),
      strip.text.y = element_blank(),
    )
  return(p)
}
transHeatmap <- function(df, family_facet){
  p <- myHeatmap(df, family_facet)
  p <- p + scale_fill_gradient2(
             low = "blue",
             mid = "black",
             high = "yellow",
             na.value="grey",
             limits = c(-4, 4),
             oob=squish
           )
  p <- p + scale_y_discrete(position = "left")
  p <- p + theme(
             axis.text.y = element_blank(),
             axis.ticks.y = element_blank(),

             panel.border=element_blank(),
             panel.grid.major=element_blank(),

             strip.background = element_blank(),
             strip.text.x = element_blank(),
             strip.text.y = element_blank(),
             )
  return(p)
}

family_heatmap <- function(mdf){
  p <- myHeatmap(mdf, family_facet)
  ##p <- p + geom_text(aes(label=Label))
  p <- p + my_scale
  #p <- p + scale_y_discrete(limits=(rev(levels(mdf$Label))))
  p <- p + theme(
             axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks.y = element_blank(),

             panel.border=element_blank(),
             panel.grid.major=element_blank(),

             strip.background = element_blank(),
             strip.text.x = element_blank(),
             strip.text.y = element_blank()
           )
  return(p)
}

info_heatmap <- function(m_df, family_facet){
  p <- myHeatmap(m_df, family_facet)
  p <- p + scale_fill_manual(
             values = c('white', 'black')
           )
  p <- p + theme(
             axis.text.y = element_blank(),
             axis.ticks.y = element_blank(),

             panel.border=element_blank(),
             panel.grid.major=element_blank(),

             strip.background = element_blank(),
             strip.text.x = element_blank(),
             strip.text.y = element_blank(),
             )
  return(p)
}

#### Transcription and Hetterochrom. By Strain Heatmaps ####

all_transcov_df <- finalFC_df %>%
  left_join(positive_cov_df, by = 'Gene_id') %>%
  left_join(info_df, by = 'Gene_id')

finalHeatmap <- function(
                         contrast,
                         het_col,
                         abs_trans_filter,
                         abs_het_filter,
                         family,
                         family_facet
                         ) {

  ## contrast <- c('A7', 'E5')
  ## het_col <- 'Cov_1000fp_500orf'
  ## abs_trans_filter <- 2
  ## abs_het_filter <- 0
  ## family_facet <- F
  ## family <- NA

  trans_col <- paste0(contrast[1], '-', contrast[2], '_MaxVal')
  het_col_1 <- paste0(het_col, '_', contrast[1])
  het_col_2 <- paste0(het_col, '_', contrast[2])
  difgenes <- paste0('difs_', contrast[1], '_', contrast[2])

  subset_all <- all_transcov_df %>%
    mutate(Het_dif = get(het_col_1) - get(het_col_2)) %>%
    ## filter(abs(get(trans_col)) > abs_trans_filter &
    ##        (abs(Het_dif) > abs_het_filter | is.na(Het_dif))) %>%
    filter(Gene_id %in% get(difgenes) &
           (abs(Het_dif) >= abs_het_filter | is.na(Het_dif))) %>%
    select(Gene_id,
           matches(trans_col),
           matches(het_col_1),
           matches(het_col_2),
           Het_dif,
           Label, Name, Annot, Family, SubFamily, Gam_specific
           )

  subset_all %>% filter(Gene_id == 'PF3D7_0813300') %>%
    print(width = 999)

  ## Subset by family is needed
  if (!is.na(family)){
    subset_all <- subset_all %>%
      filter(Family == family)
  }

  ## Ordering
  if (sorting == 'clust') {
    mtx <- subset_all %>%
      select(where(is.numeric))

    ## Make hierarquical Clustering
    dmtx <- dist(scale(mtx), method = "euclidean")
    cl <- hclust(dmtx, method = 'average')
    subset_all$Label <- factor(subset_all$Label, levels = subset_all$Label[cl$order])

  } else if (sorting == 'trans') {
    subset_all <- subset_all %>%
      arrange(desc(get(trans_col))) %>%
      mutate(Label = factor(Label, levels = Label))

  } else if (sorting == 'het') {
    subset_all <- subset_all %>%
      arrange(desc(Het_dif)) %>%
      mutate(Label = factor(Label, levels = Label))
  }

  df_name <- paste0(
    plots_dir,
    'Trans_Het_byStrain/',
    contrast[1], '_',
    contrast[2], '_',
    het_col, '_',
    'transth_', as.character(abs_trans_filter),
    '_new.csv'
  )
  sorted_df <- subset_all %>% arrange(fct_rev(Label))
  write_csv(sorted_df, df_name)

  subset_trans <- subset_all %>%
    select(Gene_id, matches(trans_col), Label, Name, Annot, Family, SubFamily, Gam_specific)

  subset_het <- subset_all %>%
    select(Gene_id, matches(het_col_1), matches(het_col_2), Label, Name, Annot, Family, SubFamily, Gam_specific)

  subset_het %>%
    print(width = 200)

  subset_hetDif <- subset_all %>%
    select(Gene_id, Het_dif, Label, Name, Annot, Family, SubFamily, Gam_specific)

  subset_info <- subset_all %>%
    select(Gene_id, Gam_specific, Label, Name, Annot, Family, SubFamily)

  melt_vars <- c('Gene_id', 'Label', 'Name', 'Annot', 'Family', 'SubFamily', 'Gam_specific')

  mdf_het <- melt(subset_het, id.vars = melt_vars)
  mdf_trans <- melt(subset_trans, id.vars = melt_vars)
  mdf_hetDif <- melt(subset_hetDif, id.vars = melt_vars)
  mdf_info <- melt(subset_info, id.vars = c('Gene_id', 'Label', 'Name', 'Annot', 'Family', 'SubFamily'))

  het_plot <- hetHeatmap(mdf_het, family_facet) +
    theme(
      axis.text.y = element_text(),
      axis.ticks.y = element_line(),
      axis.text.x = element_blank()
    )
  trans_plot <- transHeatmap(mdf_trans, family_facet)
  hetDif_plot <- hetDifHeatmap(mdf_hetDif, family_facet)
  info_plot <- info_heatmap(mdf_info, family_facet)

  all_plot <- grid.arrange(
    trans_plot, hetDif_plot, het_plot, #info_plot,
    nrow = 1, widths = c(1,1,3))

  outname <- paste0(
    plots_dir,
    'Trans_Het_byStrain/',
    contrast[1], '_',
    contrast[2], '_',
    het_col, '_',
    'transth_', as.character(abs_trans_filter),
    '_new.pdf'
  )
  print(outname)
  ggsave(outname, all_plot, device = 'pdf')
}


colnames(cor_cov_df)

contrast <- c('12B', '10G')
het_col <- 'Cov_1000fp_500orf'
abs_trans_filter <- 2
abs_het_filter <- 0
family_facet <- F
family <- NA

## trans, het or clust
sorting <- 'trans'
plt <- finalHeatmap(contrast, het_col, abs_trans_filter, abs_het_filter, family, family_facet)
het_trans_12b_10g <- read_csv(paste0(
  plots_dir, 'Trans_Het_byStrain/',
  '12B_10G_Cov_1000fp_500orf_transth_2.csv'
))

cor(het_trans_12b_10g$`12B-10G_MaxVal`, het_trans_12b_10g$Het_dif)
cor.test(het_trans_12b_10g$`12B-10G_MaxVal`, het_trans_12b_10g$Het_dif)

contrasts <- list(
  c('12B', '10G'),
  c('A7', 'E5'),
  c('A7', 'B11'),
  c('E5', 'B11')
)

for (contrast in contrasts){
  het_col <- 'Cov_1000fp_500orf'
  abs_trans_filter <- 2
  abs_het_filter <- 0
  family_facet <- F
  family <- NA
  finalHeatmap(contrast, het_col, abs_trans_filter, abs_het_filter, family, family_facet)
}


het_trans_e5_B11 <- read_csv(paste0(
  plots_dir, 'Trans_Het_byStrain/',
  'E5_B11_Cov_1000fp_500orf_transth_2_new.csv'
))

difs_E5_B11
get('difs_E5_B11')

het_trans_e5_B11$Gene_id[!het_trans_e5_B11$Gene_id %in% difs_E5_B11]
difs_E5_B11[!difs_E5_B11 %in% het_trans_e5_B11$Gene_id]


finalFC_df %>%
  filter(Gene_id == 'PF3D7_0813300') %>%
  print(width = 999)

filtered_E5_B11 %>%
  filter(Gene_id == 'PF3D7_0813300') %>%
  print(width = 999)

all_transcov_df %>%
  filter(Gene_id == 'PF3D7_0813300') %>%
  print(width = 999) %>%
  select(Gene_id, Cov_1000fp_500orf_E5, Cov_1000fp_500orf_B11) %>%
  print(width = 999)

#### 1.2B vs 10G SICER Peaks Heatmap ####

## Load Data

read_cov <- function(strain){
  cov_5p <- read_tsv(
    paste0(coverage_dir,
           'binned_1000tss_0orf_allowoverlaps_False_coverage_',
           strain,
           '.bed'
           ),
    col_names = F
  ) %>%
    select(X4, X5) %>%
    set_names(c('Gene_id', 'Cov_5p'))

  cov_ORF1 <- read_tsv(
    paste0(coverage_dir,
           'binned_CDS_first_half_coverage_',
           strain,
           '.bed'
           ),
    col_names = F
  ) %>%
    select(X4, X5) %>%
    set_names(c('Gene_id', 'Cov_ORF1'))

  cov_ORF2 <- read_tsv(
    paste0(coverage_dir,
           'binned_CDS_second_half_coverage_',
           strain,
           '.bed'
           ),
    col_names = F
  ) %>%
    select(X4, X5) %>%
    set_names(c('Gene_id', 'Cov_ORF2'))

  cov_3p <- read_tsv(
    paste0(coverage_dir,
           'binned_3prime1000_allowoverlaps_False_coverage_',
           strain,
           '.bed'
           ),
    col_names = F
  ) %>%
    select(X4, X7) %>%
    set_names(c('Gene_id', 'Cov_3p'))

  sicer_cov_ <- cov_5p %>%
    left_join(cov_ORF1) %>%
    left_join(cov_ORF2) %>%
    left_join(cov_3p)
}

make_difpeaks_hetdif_df <- function(strain1, strain2){

  undotted1 <- gsub('.', '', strain1, fixed = T)
  undotted1 <- gsub('K9', '', undotted1, fixed = T)

  undotted2 <- gsub('.', '', strain2, fixed = T)
  undotted2 <- gsub('K9', '', undotted2, fixed = T)

  cov_1 <- read_cov(strain1)
  cov_2 <- read_cov(strain2)

  subtraction_cov <- tibble(cov_1[,-1] - cov_2[,-1]) %>%
    mutate(Gene_id = cov_1$Gene_id) %>%
    select(Gene_id, everything())

  x <- difpeaks_df %>%
    select(Gene_id, contains(undotted1) & contains(undotted2))

  difpeaks <- x %>%
    set_names('Gene_id', 'dp1', 'dp2') %>%
    filter(dp1 | dp2) %>%
    pull(Gene_id)

  heat_df <- subtraction_cov %>%
    filter(Gene_id %in% difpeaks) %>%
    left_join(trans_df %>%
              select(
                Gene_id,
                contains(undotted1) & contains(undotted2) & contains('MaxVal')
              ) %>%
              set_names(c('Gene_id', 'Trans'),
                        ), by = 'Gene_id') %>%
    left_join(info_df %>% select(
                            Gene_id,
                            Label,
                            contains('DuplDel') & contains(undotted1),
                            contains('DuplDel') & contains(undotted2),
                            Gam_specific,
                            Annot
                          ) %>%
              set_names(c('Gene_id', 'Label',
                          'DuplDel_1', 'DuplDel_2',
                          'Gam_Specific', 'Annot'),
                        ), by = 'Gene_id') %>%
    select(Gene_id, Trans, everything()) %>%
    arrange(-Trans) %>%
    mutate(Label = factor(Label, levels = Label))

  heat_df %>%
    select(Gene_id, Trans, Annot) %>%
    filter(is.na(Trans)) %>%
    write_csv(paste0(
      tables_dir, 'figure3_',
      undotted1, '_', undotted2,
      '_notrans.csv'))

  sicer_heat_df <- heat_df %>%
    filter(!is.na(Trans)) %>%
    filter(!DuplDel_1) %>%
    filter(!DuplDel_2) %>%
    select(-DuplDel_1, -DuplDel_2)

  return(sicer_heat_df)
}

make_difpeaks_heat <- function(df, strain1, strain2){

  undotted1 <- gsub('.', '', strain1, fixed = T)
  undotted1 <- gsub('K9', '', undotted1, fixed = T)

  undotted2 <- gsub('.', '', strain2, fixed = T)
  undotted2 <- gsub('K9', '', undotted2, fixed = T)

  mtx <- df %>%
    select(Cov_5p, Cov_ORF1)

  ## Make hierarquical Clustering
  dmtx <- dist(mtx, method = "euclidean")
  cl <- hclust(dmtx, method = 'complete')
  df$Label <- factor(df$Label, levels = df$Label[cl$order])

  sorted_df <- df %>% arrange(fct_rev(Label))
  outname <- paste0(
    plots_dir, 'Trans_Het_byStrain/',
    'difpeaks_',
    undotted1, '_', undotted2,
    '_heatmap_nonas_hc_5pORF1.csv'
  )
  write_csv(sorted_df, outname)

  x <- df %>%
    pivot_longer(c(-Gene_id, -Label, -Annot)) %>%
    mutate(name = factor(name, levels = unique(name))) %>%
    filter(grepl('Cov', name))

  y <- df %>%
    pivot_longer(c(-Gene_id, -Label, -Annot)) %>%
    mutate(name = factor(name, levels = unique(name))) %>%
    filter(grepl('Trans', name))

  z <- df %>%
    pivot_longer(c(-Gene_id, -Label, -Annot)) %>%
    mutate(name = factor(name, levels = unique(name))) %>%
    filter(grepl('Gam', name))

  p1 <- ggplot(x, aes(x = name, y = Label, fill = value))
  p1 <- p1 + geom_tile(colour="snow3")
  p1 <- p1 + scale_fill_gradient2(
               low = "chartreuse3",
               mid = "white",
               high = "darkred",
               midpoint = 0,
               limits = c(-4, 4),
               oob=squish
             )
  p1 <- p1 + theme(
               ##axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
               panel.grid.major=element_blank(),
               panel.background=element_blank(),
               legend.position = 'bottom'
             )
  p1 <- p1 + labs(y = '', x = '')

  p2 <- ggplot(y, aes(x = name, y = Label, fill = value))
  p2 <- p2 + geom_tile(colour="snow3")
  p2 <- p2 + scale_fill_gradient2(
               low = "blue",
               mid = "black",
               high = "Yellow",
               midpoint = 0,
               limits = c(-4, 4),
               oob=squish
             )
  p2 <- p2 + theme(
               ##axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
               axis.title.y = element_blank(),
               panel.grid.major=element_blank(),
               panel.background=element_blank(),
               legend.position = 'bottom'
             )
  p2 <- p2 + labs(y = '', x = '')

  p3 <- ggplot(z, aes(x = name, y = Label, fill = value))
  p3 <- p3 + geom_tile(colour="snow3")
  p3 <- p3 + scale_fill_gradient(
               low = "white",
               high = "Black",
               )
  p3 <- p3 + scale_y_discrete(position = "right")
  p3 <- p3 + theme(
               ##axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               panel.grid.major=element_blank(),
               panel.background=element_blank(),
               legend.position = 'bottom'
             )
  p3 <- p3 + labs(y = '', x = '')

  difpeaks_plot <- ggarrange(p1, p2, p3, ncol = 3, nrow = 1, widths = c(4, 1, 3))
  difpeaks_plot
  ggsave(paste0(
    plots_dir, 'Trans_Het_byStrain/',
    'difpeaks_',
    undotted1, '_', undotted2,
    '_heatmap_nonas_hc_5pORF1.pdf'),
    difpeaks_plot,
    device = 'pdf',
    ##width = 30, units = 'cm'
    )
}

## Run all heatmaps

contrasts <- list(
  c('1.2B', '10G'),
  c('A7K9', 'E5K9'),
  c('A7K9', 'B11'),
  c('E5K9', 'B11')
)

for (contrast in contrasts){
  print(contrast)
  df <- make_difpeaks_hetdif_df(contrast[1], contrast[2])
  make_difpeaks_heat(df, contrast[1], contrast[2])
}

#### Transcription Donut Plots ####
## Donut Plot

contrasts_donuts <- function(contrast, het_col, trans_th, df){

  trans_col <- paste0(contrast[1], '-', contrast[2], '_MaxVal')
  het_col_1 <- paste0(het_col, '_', contrast[1])
  het_col_2 <- paste0(het_col, '_', contrast[2])

  infocols <- colnames(info_df)

  x <- df %>%
    select(Gene_id,
           Family_Grouped,
           matches(trans_col),
           matches(het_col_1),
           matches(het_col_2),
           one_of(infocols)
           ) %>%
    filter(abs(get(trans_col)) > trans_th)
  data <- as.data.frame(table(x$Family_Grouped))
  colnames(data) <- c('category', 'count')

  ## Compute percentages
  data$fraction = data$count / sum(data$count)

  ## Compute the cumulative percentages (top of each rectangle)
  data$ymax = cumsum(data$fraction)

  ## Compute the bottom of each rectangle
  data$ymin = c(0, head(data$ymax, n=-1))

  ## Compute label position
  data$labelPosition <- (data$ymax + data$ymin) / 2

  ## Compute a good label
  data$label <- paste0(data$category, "\n value: ", data$count)

  ## Make the plot
  p <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect(color="white") +
    ##geom_label( x=4.2, aes(y=labelPosition, label=label), size=3) +
    coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
    xlim(c(2, 4)) +# Try to remove that to see how to make a pie chart
    theme_void() + my_scale#scale_fill_viridis(discrete=TRUE)

  print(p)
  outname <- paste0(
    plots_dir,
    'Transcription_Donuts_ByStrain/',
    contrast[1], '_',
    contrast[2], '_',
    'families',
    '.svg'
  )

  ggsave(outname, p, device = 'svg')
}

contrasts <- list(
  c('12B', '10G'),
  c('A7', 'E5'),
  c('A7', 'B11'),
  c('E5', 'B11')
  )

het_col <- 'Cov_1000fp_500orf'
contrast <-   c('12B', '10G')
trans_th <- 2
df <- trans_het_difpeaks_df

for(contrast in contrasts){
  print(contrast)
  contrasts_donuts(contrast, het_col, trans_th, df)
}

#### Get coverage for ON/OFF genes with sign ####
## Select transcription cols
## Create empty DF
hetcols <- colnames(cor_cov_df %>% select(contains('12B')))
hetcols <- str_replace(hetcols, '_12B', '')

col_names <- c('Gene_id',
               paste0(hetcols, '_On'),
               paste0(hetcols, '_Off'),
               paste0(hetcols, '_Dif'))

cn <- setNames(rep('', length(col_names)), col_names)
df <- bind_rows(cn)[0,]

## Traverse max trans. df and get coverage cols

## which(tdf$Gene_id == 'PF3D7_0100300')
## i <- 63
## tdf[63,]

for (i in 1:dim(finalFC_df)[1]){

  gid <- as.character(finalFC_df$Gene_id[i])
  on <- finalFC_df$On_trans[i]
  off <- finalFC_df$Off_trans[i]

  ## Select contrast in which difference sign must be switched
  contrast <- paste0(on, '-', off)
  positive_contrasts <- c(
    '12B-10G', '12B-A7', '12B-E5', '12B-B11',
    '10G-A7', '10G-E5', '10G-B11',
    'A7-E5', 'A7-B11',
    'E5-B11'
  )

  ifelse(
    contrast %in% positive_contrasts,
    neg_dif <- FALSE,
    neg_dif <- TRUE
    )

  ## Main Loop
  if (gid %in% cor_cov_df$Gene_id & !is.na(on) & !is.na(off)){

    onvect <- cor_cov_df %>%
      filter(Gene_id == gid) %>%
      select(contains(on))

    offvect <- cor_cov_df %>%
      filter(Gene_id == gid) %>%
      select(contains(off))

    diffvect <- onvect-offvect
    if (neg_dif) {diffvect <- -diffvect}

    row <- c(gid,
             unlist(onvect, use.names = F),
             unlist(offvect, use.names = F),
             unlist(diffvect, use.names = F))
    row <- setNames(row, col_names)
    df <- df %>% add_row(bind_rows(row))
  }
}

signed_maxtrans_cov <- df %>%
  mutate(across(-Gene_id,  as.numeric)) %>%
  full_join(finalFC_df, by = 'Gene_id') %>%
  left_join(info_df)

#### Get coverage for ON/OFF genes unsigned ####
## Select transcription cols
## Create empty DF
hetcols <- colnames(cor_cov_df %>% select(contains('12B')))
hetcols <- str_replace(hetcols, '_12B', '')

col_names <- c('Gene_id',
               paste0(hetcols, '_On'),
               paste0(hetcols, '_Off'),
               paste0(hetcols, '_Dif'))

cn <- setNames(rep('', length(col_names)), col_names)
df <- bind_rows(cn)[0,]

## Traverse max trans. df and get coverage cols

## which(finalFC_df$Gene_id == 'PF3D7_0100300')
## i <- 63
## finalFC_df[63,]

for (i in 1:dim(finalFC_df)[1]){

  gid <- as.character(finalFC_df$Gene_id[i])
  on <- finalFC_df$On_trans[i]
  off <- finalFC_df$Off_trans[i]

  ## Select contrast in which difference sign must be switched
  ## contrast <- paste0(on, '-', off)
  ## positive_contrasts <- c(
  ##   '12B-10G', '12B-A7', '12B-E5', '12B-B11',
  ##   '10G-A7', '10G-E5', '10G-B11',
  ##   'A7-E5', 'A7-B11',
  ##   'B11-E5'
  ## )

  ## ifelse(
  ##   contrast %in% positive_contrasts,
  ##   neg_dif <- FALSE,
  ##   neg_dif <- TRUE
  ##   )

  ## Main Loop
  if (gid %in% positive_cov_df$Gene_id & !is.na(on) & !is.na(off)){

    onvect <- positive_cov_df %>%
      filter(Gene_id == gid) %>%
      select(contains(on))

    offvect <- positive_cov_df %>%
      filter(Gene_id == gid) %>%
      select(contains(off))

    diffvect <- onvect-offvect
    #if (neg_dif) {diffvect <- -diffvect}

    row <- c(gid,
             unlist(onvect, use.names = F),
             unlist(offvect, use.names = F),
             unlist(diffvect, use.names = F))
    row <- setNames(row, col_names)
    df <- df %>% add_row(bind_rows(row))
  }
}

unsigned_maxtrans_cov <- df %>%
  mutate(across(-Gene_id,  as.numeric)) %>%
  full_join(finalFC_df, by = 'Gene_id') %>%
  mutate(Max_aMAFC = abs(Max_aMAFC)) %>%
  left_join(info_df)

genes_nocov <- unsigned_maxtrans_cov %>%
     select(Gene_id, contains('Cov_1000fp_500orf')) %>%
     filter(!complete.cases(.)) %>%
     select(Gene_id) %>% pull()

unsigned_maxtrans_cov %>%
  select(Gene_id, contains('Cov_1000fp_500orf')) %>%
  filter(!complete.cases(.)) %>%
  print(width = 400)

#### Correlations On/Off genes ####
## Get Correlations

cor_signed_maxtrans <- signed_maxtrans_cov %>%
  filter(!Is_3D7B)

## Subset column by column
clnms <- colnames(cor_signed_maxtrans %>% select(contains('_Dif')))
for (c in clnms) {
  print(c)
  cor_genes <- cor_signed_maxtrans %>%
    select(Max_aMAFC, matches(c)) %>%
    filter(abs(Max_aMAFC) > 1.5) %>%
    drop_na()

  cormtx <- cor(cor_genes, method = 'pearson')
  print(paste0('Number of genes: ', dim(cor_genes)[1]))
  print(as.data.frame(cormtx[,1]))
  print('------------------')
}

## All toghether

trans_th <- 2

cor_genes <- cor_signed_maxtrans %>%
  select(Max_aMAFC, contains('_Dif')) %>%
  filter(abs(Max_aMAFC) > trans_th) %>%
  drop_na()

cor_genes

cormtx <- cor(cor_genes, method = 'pearson')
print(paste0('Number of genes: ', dim(cor_genes)[1]))
print(as.data.frame(cormtx[,1]))

outdf <- as.data.frame(cormtx)
outdf['Corr_with'] <- rownames(outdf)
outdf <- outdf %>% select(Corr_with, everything())
outname = paste0(tables_dir, 'corr_MaxFC_', trans_th, '_allowoverlaps.csv')
write_csv(outdf, outname)

## By subsets

get_corr <- function(df){
  cor_genes <- df %>%
    select(-Gene_id) %>%
    filter(abs(Max_aMAFC) > trans_th) %>%
    drop_na()

  cormtx <- cor(cor_genes, method = 'pearson')
  print(paste0('Number of genes: ', dim(cor_genes)[1]))

  outdf <- as.data.frame(cormtx)
  outdf['Corr_with'] <- rownames(outdf)
  outdf <- outdf %>% select(Corr_with, everything())
  outdf <- as_tibble(outdf[,c(1,2)])
  print(outdf)
}

trans_th <- 2

## Abs bins
abs_bins <- cor_signed_maxtrans %>%
  select(Gene_id, Max_aMAFC,
         Cov_1500fp_manybins_Dif,
         Cov_1000fp_manybins_Dif,
         Cov_500fp_manybins_Dif,
         contains('ORF', ignore.case = F) & contains('Dif'),
         Cov_500tp_manybins_Dif,
         Cov_1000tp_manybins_Dif
         )

abs_bins_cor <- get_corr(abs_bins)

## Rel bins

rel_bins <- cor_signed_maxtrans %>%
  select(Gene_id, Max_aMAFC,
         p5utr_Cov_Dif,
         contains('qorf') & contains('Dif'),
         p3utr_Cov_Dif
         )

rel_bins_cor <- get_corr(rel_bins)

## TSS

tss_bins <- cor_signed_maxtrans %>%
  select(Gene_id, Max_aMAFC,
         TSS_Cov_Dif,
         Cov_500fp_manybins_Dif,
         p5utr_Cov_Dif
         )

tss_bins_cor <- get_corr(tss_bins)


## Plots

df <- cor_signed_maxtrans %>% filter(abs(Max_aMAFC) > trans_th)

for (c in clnms){
  p <- ggplot(df,
              aes(x = Max_aMAFC, y = get(c), color = Gam_specific)) +
    geom_point() +
    geom_label_repel(data=subset(df, get(c) > 1),
                     aes(label = Gene_id),
                     box.padding   = 0.35,
                     point.padding = 0.5,
                     segment.color = 'grey50')
  ggsave(paste0(plots_dir, 'MaxFC_Cor_Plots/cor_signed_maxFC_', c, '_allowoverlaps.png'), p, device = 'png')
}

#### Make correlation by bin plots STEP ####
corr_data <- as_tibble(outdf[,1:2])

corr_data %>%
  arrange(Max_aMAFC) %>%
  print(n = 40)

bins <- c(
  'Cov_2000fp_manybins_Dif',
  'Cov_1500fp_manybins_Dif',
  'Cov_1000fp_manybins_Dif',
  'Cov_500fp_manybins_Dif',
  'Cov_500orf_Dif',
  'Cov_1000orf_Dif',
  'Cov_1500orf_Dif',
  'Cov_500tp_manybins_Dif',
  'Cov_1000tp_manybins_Dif',
  'Cov_1500tp_manybins_Dif'
  )


## Geom_step approach

pos <- c(-2000, -1500, -1000, -500, 0, 500, 1000, 2500, 3000, 3500)

cor_bins <- corr_data %>%
  filter(Corr_with %in% bins) %>%
  arrange(factor(Corr_with, levels = bins)) %>%
  mutate(Pos = pos)

cor_bins

cor_bins2 <- cor_bins %>%
  add_row(Corr_with = 'start', Max_aMAFC = 0, Pos = -2500) %>%
  add_row(Corr_with = 'end_of_gene', Max_aMAFC = 0, Pos = 1500) %>%
  add_row(Corr_with = 'end', Max_aMAFC = 0, Pos = 4000) %>%
  add_row(Corr_with = 'endtail', Max_aMAFC = 0, Pos = 4500) %>%
  arrange(Pos)

xbreaks <- c(
  '',
  'ATG-2000',
  'ATG-1500',
  'ATG-1000',
  'ATG-500',
  'ATG',
  'ATG+500',
  'ATG+1000',
  'ATG+1500',
  'End',
  'End+500',
  'End+1000',
  'End+1500',
  ''
  )

p <- ggplot(cor_bins2, aes(x = Pos, y = -Max_aMAFC, group = 1)) +
  geom_step() +
  ylim(0, 1) +
  ylab('- Pearson Corr') + xlab('Region') +
  theme_minimal() +
  scale_x_continuous(breaks = cor_bins2$Pos, labels = xbreaks) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'gray') +
  geom_vline(xintercept = 2500, linetype = 'dashed', color = 'gray')

p
ggsave(paste0(plots_dir, 'Correlations_Intervals/corr_genemodel_steps_allowoverlaps.png'), p, device = 'png')

#### Make correlation by bin plots BAR ####
## Geom Bar approach

cor_bins_plot <- function(df, labs, regions){
  corplot_df <- df %>%
    filter(Corr_with != 'Max_aMAFC') %>%
    mutate(Corr_with = factor(Corr_with, levels = Corr_with)) %>%
    arrange(Corr_with) %>%
    mutate(Labs = factor(labs, levels = labs)) %>%
    mutate(Region = regions)

  p <- ggplot(corplot_df, aes(x = Labs, y = -Max_aMAFC, color = Region)) +
    geom_boxplot(key_glyph = "point") +
    ylim(0, 1) +
    ylab('- Pearson Corr\n') + xlab('\nRegion') +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      text = element_text(size=30),
            legend.position = "none"
          )

  ggsave(
    paste0(plots_dir, 'Correlations_Intervals/', deparse(substitute(df)), '_allowoverlaps.pdf'),
    p, device = 'pdf')
  p
}


## Abs

labs = c(
  '-1500 to -1000',
  '-1000 to -500',
  '-500 to ATG',
  'ATG to 500',
  '500 to 1000',
  '1000 to 1500',
  'End to +500',
  '+500 to +1000'
)

regions = c(rep('5prime', 3), rep('ORF', 3), rep('3prime', 2))

cor_bins_plot(abs_bins_cor, labs, regions)

## Rel

labs = c(
  '5\'UTR',
  'ORF 1/4',
  'ORF 2/4',
  'ORF 3/4',
  'ORF 4/4',
  '3\'UTR'
)

regions = c(rep('5prime', 1), rep('ORF', 4), rep('3prime', 1))

cor_bins_plot(rel_bins_cor, labs, regions)

## TSS

labs = c(
  'TSS',
  '-500 to ATG',
  '5\'UTR'
)

regions = c(rep('5prime', 3), rep('ORF', 0), rep('3prime', 0))

cor_bins_plot(tss_bins_cor, labs, regions) + geom_boxplot(color = 'green')

#### Make correlation On/Off genes Heatmap ####
trans_th <- 2

cor_signed_maxtrans

cor_genes <- cor_signed_maxtrans %>%
  select(Max_aMAFC, contains('_Dif')) %>%
  filter(abs(Max_aMAFC) > trans_th) %>%
  drop_na()

cor_genes

cormtx <- cor(cor_genes, method = 'pearson')
print(paste0('Number of genes: ', dim(cor_genes)[1]))
print(as.data.frame(cormtx[,1]))

outdf <- as.data.frame(cormtx)
outdf['Corr_with'] <- rownames(outdf)
outdf <- outdf %>% select(Corr_with, everything()) %>% as_tibble()

plotdf <- outdf %>% select(Corr_with, Max_aMAFC) %>% arrange(Max_aMAFC)
mplotdf <- melt(plotdf)
mplotdf$Corr_with <- factor(mplotdf$Corr_with, levels = mplotdf$Corr_with)

p <- ggplot(mplotdf, aes(x = variable, y = Corr_with)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label=value)) +
  scale_fill_gradient2(low = "chartreuse3",
                       mid = "white",
                       high = "darkred",
                       na.value="grey90")

#p
ggsave(paste0(plots_dir, 'Correlations_Intervals/ranked_intervals_heatmap.png'))

#### Gene-Model Analysis: Load Data ####
## Coverage

gm_12B <- read_csv(paste0(
  genemodel_dir,
  '1.2B_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10_5binned_cov_2prevpost.csv'
))
gm_10G <- read_csv(paste0(
  genemodel_dir,
  '10G_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10_5binned_cov_2prevpost.csv'
))
gm_A7 <- read_csv(paste0(
  genemodel_dir,
  'A7K9_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10_5binned_cov_2prevpost.csv'
))
gm_E5 <- read_csv(paste0(
  genemodel_dir,
  'E5K9_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10_5binned_cov_2prevpost.csv'
))
gm_B11 <- read_csv(paste0(
  genemodel_dir,
  'B11_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10_5binned_cov_2prevpost.csv'
))


gm_strains <- list(
  df12B=gm_12B,
  df10G=gm_10G,
  dfA7=gm_A7,
  dfE5=gm_E5,
  dfB11=gm_B11
)

## Change dfs in a list
gm_strains <- lapply(gm_strains, function(df) {
    colnames(df)[1] <- "Gene_id"
    df
})

## Convert list into individual objects again
list2env(gm_strains, envir=.GlobalEnv)

##Create numeric non-na mtxs
nona.mtxs <- lapply(gm_strains, function(df) {
    mtx <- as.matrix(df[complete.cases(df),-1])
    rownames(mtx) <- df[complete.cases(df),1] %>% pull()
    mtx
})

nona_gm_strains <- lapply(gm_strains, function(tibble) {
  tibble %>%
    filter(complete.cases(tibble))
  tibble
})


names(nona.mtxs) <- c("mtx12B", "mtx10G", 'mtxA7', 'mtxE5', 'mtxB11')
names(nona_gm_strains) <- c('nona_12B', 'nona_10G', 'nona_A7', 'nona_E5', 'nona_B11')

## Create side-by-side matrices
#mtx_12b10g <- cbind(nona.mtxs$mtx12B, nona.mtxs$mtx10G)

mtx_all <- cbind(
  nona.mtxs$mtx12B,
  nona.mtxs$mtx10G,
  nona.mtxs$mtxA7,
  nona.mtxs$mtxE5,
  nona.mtxs$mtxB11
)

colnames(mtx_all) <- 1:dim(mtx_all)[2]

bins <- paste0('bin', c(1:23))
all_cols <- c(
  paste0(bins, '_12B'),
  paste0(bins, '_10G'),
  paste0(bins, '_A7'),
  paste0(bins, '_E5'),
  paste0(bins, '_B11')
)

gmodel_all <- gm_strains$df12B %>%
  full_join(gm_strains$df10G, by = 'Gene_id') %>%
  full_join(gm_strains$dfA7, by = 'Gene_id') %>%
  full_join(gm_strains$dfE5, by = 'Gene_id') %>%
  full_join(gm_strains$dfB11, by = 'Gene_id') %>%
  setNames(c('Gene_id', all_cols))

nona_gm_all <- gmodel_all %>%
  filter(complete.cases(gmodel_all))

#### Gene-Model Analysis: PCAs ####
## Filtered Differential Genes
gm_trans <- finalFC_df %>%
  left_join(gmodel_all) %>%
  left_join(info_df) %>%
  filter(across(contains('bin'), complete.cases))

gm_trans_mtx <- gm_trans %>%
  select(contains('bin'))

gm_trans_pca <- prcomp(gm_trans_mtx)
gm_trans_df <- as_tibble(gm_trans_pca$x[, c(1,2)])
gm_trans <- bind_cols(gm_trans_df,  gm_trans)

#### PCA plots ####
## All-genes, transcription filter
p <- ggplot(gm_trans, aes(x=PC1, y=PC2, color = Family_Grouped))
p <- p + geom_point() + my_scale2
##p
ggsave(paste0(plots_dir, 'Gene_Model/PCAs/pca_families.svg'), p, device = 'svg')

p <- ggplot(gm_trans, aes(x=PC1, y=PC2, color = Gam_specific))
p <- p + geom_point()
##p
ggsave(paste0(plots_dir, 'Gene_Model/PCAs/pca_gametocytes.svg'), p, device = 'svg')

p <- ggplot(gm_trans, aes(x=PC1, y=PC2, color = Max_aMAFC))
p <- p + geom_point()
p <- p + scale_color_gradient2(midpoint = 0, low="red", mid = "black", high="green")
##p
ggsave(paste0(plots_dir, 'Gene_Model/PCAs/pca_transcription.svg'), p, device = 'svg')

#### Gene-Model Analysis: Max. trans Dif. Approach ####
tdf <- maxDif_df %>%
  filter(!Is_tRNA) %>%
  select(Gene_id, Max_aMAFC, On_trans, Off_trans, Is_3D7B)

gmodel_all_pos <- gmodel_all %>%
  mutate(across(-Gene_id, ~ ifelse(.x > 0, .x, 0)))

## Create empty DF
hetcols <- gmodel_all %>%
  select(contains('12B') & contains('bin')) %>%
  colnames(.) %>%
  gsub('12B', '', ., fixed=TRUE)

col_names <- c('Gene_id',
               paste0(hetcols, 'On'),
               paste0(hetcols, 'Off'),
               paste0(hetcols, 'Dif'))

cn <- setNames(rep('', length(col_names)), col_names)
df <- bind_rows(cn)[0,]

## Traverse max trans. df and get coverage cols

for (i in 1:dim(tdf)[1]){

  gid <- as.character(tdf$Gene_id[i])
  on <- tdf$On_trans[i]
  off <- tdf$Off_trans[i]

  ## Select contrast in which difference sign must be switched
  ## contrast <- paste0(on, '-', off)
  ## positive_contrasts <- c(
  ##   '12B-10G', '12B-A7', '12B-E5', '12B-B11',
  ##   '10G-A7', '10G-E5', '10G-B11',
  ##   'A7-E5', 'A7-B11',
  ##   'B11-E5'
  ## )

  ## ifelse(
  ##   contrast %in% positive_contrasts,
  ##   neg_dif <- FALSE,
  ##   neg_dif <- TRUE
  ## )

  ### Main Loop
  if (gid %in% gmodel_all_pos$Gene_id & !is.na(on) & !is.na(off)){

    onvect <- gmodel_all_pos %>%
      filter(Gene_id == gid) %>%
      select(contains(on))

    offvect <- gmodel_all_pos %>%
      filter(Gene_id == gid) %>%
      select(contains(off))

    diffvect <- onvect-offvect
    #if (neg_dif) {diffvect <- -diffvect}

    row <- c(gid,
             unlist(onvect, use.names = F),
             unlist(offvect, use.names = F),
             unlist(onvect-offvect, use.names = F))
    row <- setNames(row, col_names)
    df <- df %>% add_row(bind_rows(row))
  }
}

df <- df %>%
  mutate(across(-Gene_id, as.numeric))

gm_pos_maxtrans <- df %>%
  full_join(tdf, by = 'Gene_id') %>%
  left_join(info_df) %>%
  mutate(Max_aMAFC = abs(Max_aMAFC))

#### Gene-Model Analysis: K-Means Heatmap Functions ####
my_info_heatmap <- function(m_df, family_facet){
  p <- myHeatmap(m_df, family_facet)
  p <- p + scale_fill_manual(
             values = c('white', 'black')
           )
  p <- p + theme(
             axis.text.y = element_blank(),
             axis.ticks.y = element_blank(),

             panel.border=element_blank(),
             panel.grid.major=element_blank(),

             strip.background = element_blank(),
             strip.text.x = element_blank(),
             strip.text.y = element_blank(),
             )
  return(p)
}

test_kmeans <- function(mtx){

  ## Set max k for k analysis
  maxk <- mtx %>% distinct() %>% nrow()
  maxk <- min(maxk, 20)

  ## Elbow method, keep reducing maxk if it fails
  while (maxk > 0) {
    test <- try(
      ## elbow plot
      elbow <- fviz_nbclust(mtx, kmeans, method = "wss", k.max = maxk) +
        labs(subtitle = "Elbow method")
    , silent = T)
    ## If kmax throws an error, reduce it by 1
    if (class(test)[1] != "gg") {
      maxk <- maxk -1
    } else {
      break
    }
  }

  ## Silhouette method (if kmax works for elbow it will work for silhouette)
  silhouette <- fviz_nbclust(mtx, kmeans, method = "silhouette", k.max = maxk) +
    labs(subtitle = "Silhouette method")

  result = list(elbow = elbow, silhouette = silhouette)
}

heatMap_allstrains_kmeans <- function(df, aFC_th, fam, fam_facet, nclu, tbar, fbar, labels){

  ## set.seed(123)
  ## df <- kmeans_df
  ## aFC_th <- 2
  ## fam <- NA
  ## fam_facet <- F
  ## nclu <- 7
  ## tbar <- T
  ## fbar <- T

  ## Returns a list object with list = (df, plot)

  subset_all <- df %>%
    filter(abs(Max_aMAFC) > aFC_th)

  subset_all %>%
    select(Gene_id, Max_aMAFC)

  if (!is.na(fam)){
    subset_all <- subset_all %>%
      filter(Family == fam)
  }

  ## subset_all <- subset_all %>%
  ##   filter(across(
  ##     contains('bin'),
  ##     complete.cases
  ##   ))

  mtx <- subset_all %>%
    select(contains('bin'),
           -contains('bin1_'),
           -contains('bin2_'),
           -contains('bin3_'),
           -contains('bin4_'),
           -contains('bin20_'),
           -contains('bin21_'),
           -contains('bin22_'),
           -contains('bin23_'),
           )
    #select(contains('On'), contains('Off'))
    #select(contains('Dif'))

  ## Make K-means Clustering

  km_fit <- kmeans(mtx, nclu, iter.max = 100000)
  cl <- km_fit$cluster
  subset_all['Cluster'] <- cl
  subset_all <- subset_all %>%
    mutate(Cluster = case_when(
             Cluster == 1 ~ 1,
             Cluster == 2 ~ 7,
             Cluster == 3 ~ 2,
             Cluster == 4 ~ 4,
             Cluster == 5 ~ 6,
             Cluster == 6 ~ 3,
             Cluster == 7 ~ 5,
             )) %>%
    arrange(Cluster) %>%
    mutate(Label = factor(Label, levels = rev(Label)))
  ## We use rev() because for some reason, order gets reversed inside facets.

  ## Test optimal k's
  k_tests <- test_kmeans(mtx)
  elbow <- k_tests$elbow
  silhouette <- k_tests$silhouette

  ## Subset DF for different heatmaps

  subset_trans <- subset_all %>%
    select(Gene_id, Max_aMAFC, Label, Name, Annot, Family, SubFamily, Cluster)

  subset_het_on <- subset_all %>%
    select(Gene_id, Label,
           contains('On') & contains('bin'),
           Name, Annot, Family, SubFamily, Cluster)

  subset_het_off <- subset_all %>%
    select(Gene_id, Label,
           contains('Off') & contains('bin'),
           Name, Annot, Family, SubFamily, Cluster)

  subset_hetDif <- subset_all %>%
    select(Gene_id, contains('_Dif'), Label, Name, Annot, Family, SubFamily, Cluster)

  subset_fambar <- subset_all %>%
    select(Gene_id, Label,
           Family_Grouped,
           Name, Annot, Family, SubFamily, Cluster) %>%
    #mutate(Label = factor(Label, levels = Label)) %>%
    replace_na(list(Family_Grouped = 'Non CVGs'))

  subset_info <- subset_all %>%
    select(
      Gene_id, Gam_specific, Label, Name, Annot, Family, SubFamily, Cluster
    )

  ## Melt Data-Frames

  melt_vars <- c('Gene_id', 'Label', 'Name', 'Annot', 'Family', 'SubFamily', 'Cluster')
  mdf_het_on <- melt(subset_het_on, id.vars = melt_vars)
  mdf_het_off <- melt(subset_het_off, id.vars = melt_vars)
  mdf_trans <- melt(subset_trans, id.vars = melt_vars)
  mdf_hetDif <- melt(subset_hetDif, id.vars = melt_vars)
  mdf_fambar <- melt(subset_fambar, id.vars = melt_vars)
  mdf_info <- melt(subset_info, id.vars = melt_vars)

  head(mdf_fambar)

  ## Create individual Heatmaps

  het_plot_on <- hetHeatmap(mdf_het_on, fam_facet) +
    facet_grid(Cluster ~., scales = "free_y", space = "free") +
    theme(axis.text.x = element_blank())

  het_plot_off <- hetHeatmap(mdf_het_off, fam_facet) +
    facet_grid(Cluster ~., scales = "free_y", space = "free")+
    theme(axis.text.x = element_blank())

  trans_plot <- transHeatmap(mdf_trans, fam_facet) +
    facet_grid(Cluster ~., scales = "free_y", space = "free") +
    theme(axis.text.x = element_blank()) +
    ## scale_fill_gradient(low = "white",
    ##                     high = "blue",
    ##                     na.value="gray",
    ##                     limits = c(0,NA))
    scale_fill_gradient2(
      low = "yellow",
      mid = 'white',
      high = "blue",
      na.value="black"
    )

  hetDif_plot <- hetDifHeatmap(mdf_hetDif, fam_facet) +
    facet_grid(Cluster ~., scales = "free_y", space = "free") +
    theme(axis.text.x = element_blank())

  fambar_plot <- family_heatmap(mdf_fambar) +
    facet_grid(Cluster ~., scales = "free_y", space = "free")

  info_plot <- my_info_heatmap(mdf_info, fam_facet) +
    facet_grid(Cluster ~., scales = "free_y", space = "free")

  ## Create fake Heatmap for Labels

  labels_df <- mdf_trans %>%
    select(Label, Family, Cluster) %>%
    mutate(value = 1, variable = 'fake_val')

  labels_plot <- hetHeatmap(labels_df, fam_facet) +
    facet_grid(Cluster ~., scales = "free_y", space = "free")  +
    theme(
      axis.text.y = element_text(),
      axis.ticks.y = element_line(),
      axis.text.x = element_blank()
    )

  ## Arrange plots (with or without transcription heatmap and fambar)

  pl_wd <- tibble(
    P_names = c('trans', 'het_dif', 'het_on', 'het_off', 'fambar', 'info_plot', 'labels'),
    Plots = list(
      trans_plot,
      hetDif_plot,
      het_plot_on,
      het_plot_off,
      fambar_plot,
      info_plot,
      labels_plot),
    widths = c(0.1,1,1,1,0.1,0.1,1)
  )

  if (!tbar){pl_wd <- pl_wd %>% filter(P_names != 'trans')}
  if (!fbar){pl_wd <- pl_wd %>% filter(P_names != 'fambar')}
  if (!labels){pl_wd <- pl_wd %>% filter(P_names != 'labels')}

  all_plot <- grid.arrange(grobs = pl_wd$Plots, nrow = 1, widths = pl_wd$widths)

  ## Create output
  result <- list(df = subset_all, plot = all_plot, elbow = elbow, silhouette = silhouette)
  return(result)
}


heatMap_families_kmeans <- function(df, aFC_th, fam, fam_facet, nclu, tbar, fbar, labels){

  ## df <- gm_pos_maxtrans %>%
  ##   mutate(Label <- paste(Gene_id, SubFamily, sep <- ': '))
  ## aFC_th <- 0
  ## fam <- 'ACS'
  ## fam_facet <- F
  ## nclu <- as.numeric(fam_k[2])
  ## tbar <- T
  ## fbar <- F
  ## labels <- F

  ## Returns a list object with list = (df, plot)

  subset_all <- df %>%
    filter(abs(Max_aMAFC) > aFC_th)

  subset_all %>%
    select(Gene_id, Max_aMAFC)

  if (!is.na(fam)){
    subset_all <- subset_all %>%
      filter(Family == fam)
  }

  subset_all <- subset_all %>%
    filter(across(
      contains('bin'),
      complete.cases
    ))

  mtx <- subset_all %>%
    select(contains('bin'),
           -contains('bin1_'),
           -contains('bin2_'),
           -contains('bin3_'),
           -contains('bin4_'),
           -contains('bin20_'),
           -contains('bin21_'),
           -contains('bin22_'),
           -contains('bin23_'),
           )
    #select(contains('On'), contains('Off'))
    #select(contains('Dif'))

  ## Make K-means Clustering

  km_fit <- kmeans(mtx, nclu, iter.max = 100000)
  cl <- km_fit$cluster
  subset_all['Cluster'] <- cl
  subset_all <- subset_all %>%
    mutate(Cluster = case_when(
             Cluster == 1 ~ 1,
             Cluster == 2 ~ 7,
             Cluster == 3 ~ 2,
             Cluster == 4 ~ 4,
             Cluster == 5 ~ 6,
             Cluster == 6 ~ 3,
             Cluster == 7 ~ 5,
             )) %>%
    arrange(Cluster) %>%
    mutate(Label = factor(Label, levels = rev(Label)))
  ## We use rev() because for some reason, order gets reversed inside facets.

  ## Test optimal k's
  k_tests <- test_kmeans(mtx)
  elbow <- k_tests$elbow
  silhouette <- k_tests$silhouette

  ## Subset DF for different heatmaps

  subset_trans <- subset_all %>%
    select(Gene_id, Max_aMAFC, Label, Name, Annot, Family, SubFamily, Cluster)

  subset_het_on <- subset_all %>%
    select(Gene_id, Label,
           contains('On') & contains('bin'),
           Name, Annot, Family, SubFamily, Cluster)

  subset_het_off <- subset_all %>%
    select(Gene_id, Label,
           contains('Off') & contains('bin'),
           Name, Annot, Family, SubFamily, Cluster)

  subset_hetDif <- subset_all %>%
    select(Gene_id, contains('_Dif'), Label, Name, Annot, Family, SubFamily, Cluster)

  subset_fambar <- subset_all %>%
    select(Gene_id, Label,
           Family_Grouped,
           Name, Annot, Family, SubFamily, Cluster) %>%
    #mutate(Label = factor(Label, levels = Label)) %>%
    replace_na(list(Family_Grouped = 'Non CVGs'))

  subset_info <- subset_all %>%
    select(
      Gene_id, Gam_specific, Label, Name, Annot, Family, SubFamily, Cluster
    )

  ## Melt Data-Frames

  melt_vars <- c('Gene_id', 'Label', 'Name', 'Annot', 'Family', 'SubFamily', 'Cluster')
  mdf_het_on <- melt(subset_het_on, id.vars = melt_vars)
  mdf_het_off <- melt(subset_het_off, id.vars = melt_vars)
  mdf_trans <- melt(subset_trans, id.vars = melt_vars)
  mdf_hetDif <- melt(subset_hetDif, id.vars = melt_vars)
  mdf_fambar <- melt(subset_fambar, id.vars = melt_vars)
  mdf_info <- melt(subset_info, id.vars = melt_vars)

  head(mdf_fambar)

  ## Create individual Heatmaps

  het_plot_on <- hetHeatmap(mdf_het_on, fam_facet) +
    facet_grid(Cluster ~., scales = "free_y", space = "free")
    #theme(axis.text.x = element_blank(), legend.position="none")

  het_plot_off <- hetHeatmap(mdf_het_off, fam_facet) +
    facet_grid(Cluster ~., scales = "free_y", space = "free")
    #theme(axis.text.x = element_blank(), legend.position="none")

  trans_plot <- transHeatmap(mdf_trans, fam_facet) +
    facet_grid(Cluster ~., scales = "free_y", space = "free")
    #theme(axis.text.x = element_blank(), legend.position="none")

  hetDif_plot <- hetDifHeatmap(mdf_hetDif, fam_facet) +
    facet_grid(Cluster ~., scales = "free_y", space = "free")
    #theme(axis.text.x = element_blank(), legend.position="none")

  fambar_plot <- family_heatmap(mdf_fambar) +
    facet_grid(Cluster ~., scales = "free_y", space = "free")
    #theme(legend.position="none")

  info_plot <- my_info_heatmap(mdf_info, fam_facet) +
    facet_grid(Cluster ~., scales = "free_y", space = "free")
    #theme(axis.text.x = element_blank(),legend.position="none")

  ## Create fake Heatmap for Labels

  labels_df <- mdf_trans %>%
    select(Label, Family, Cluster) %>%
    mutate(value = 1, variable = 'fake_val')

  labels_plot <- hetHeatmap(labels_df, fam_facet) +
    facet_grid(Cluster ~., scales = "free_y", space = "free") +
    theme(
      axis.text.y = element_text(),
      axis.ticks.y = element_line(),
      axis.text.x = element_blank()
    )
    #theme(legend.position="none")

  ## Arrange plots (with or without transcription heatmap and fambar)

  pl_wd <- tibble(
    P_names = c('trans', 'het_dif', 'het_on', 'het_off', 'fambar', 'info_plot', 'labels'),
    Plots = list(
      trans_plot,
      hetDif_plot,
      het_plot_on,
      het_plot_off,
      fambar_plot,
      info_plot,
      labels_plot),
    widths = c(0.1,1,1,1,0.1,0.1,1)
  )

  if (!tbar){pl_wd <- pl_wd %>% filter(P_names != 'trans')}
  if (!fbar){pl_wd <- pl_wd %>% filter(P_names != 'fambar')}
  if (!labels){pl_wd <- pl_wd %>% filter(P_names != 'labels')}

  all_plot <- grid.arrange(grobs = pl_wd$Plots, nrow = 1, widths = pl_wd$widths)

  ## Create output
  result <- list(df = subset_all, plot = all_plot, elbow = elbow, silhouette = silhouette)
  return(result)
}

#### Gene-Model Analysis: K-Means Clustering Heatmap ####
## Current Heatmap

finalFC_df %>%
  left_join(info_df, by = 'Gene_id') %>%
  count(Variant)

finalFC_df %>%
  filter(Off_trans == 'B11')

set.seed(123)

kmeans_df <- gm_pos_maxtrans %>%
  filter(Gene_id %in% finalFC_df$Gene_id)

nclu <- 7

x <- heatMap_allstrains_kmeans(
  df = kmeans_df,
  aFC_th = 2,
  fam = NA,
  fam_facet = F,
  nclu = nclu,
  tbar = T,
  fbar = T,
  labels = T
)

plot(x$plot)
x$df %>%
  select(Gene_id, Label, Cluster)

## Save resulting Heatmap and DF
ggsave(
  paste0(plots_dir, 'Gene_Model/heatmap_with_neighbors_Kmeans_filteredFC_reordered_new.pdf'),
  x$plot,
  device = 'pdf',
  height = 80, width = 60, units = 'cm'
)

all_gmodel_kmeans <- x$df

all_gmodel_kmeans %>%
  filter(Gene_id == 'PF3D7_0324800') %>%
  print(width = 999)


all_gmodel_kmeans %>%
  select(-contains('bin')) %>%
  write_tsv(paste0(
    plots_dir,
    'Gene_Model/heatmap_with_neighbors_Kmeans_filteredFC_reordered_table.csv'
  ))

## Check which genes are we loosing
plot_df <- read_tsv(paste0(
  plots_dir,
  'Gene_Model/heatmap_with_neighbors_Kmeans_filteredFC_reordered_table.csv'
))

plot_df %>%
  count(Cluster)

finalFC_df %>%
  filter(!Gene_id %in% plot_df$Gene_id) %>%
  left_join(info_df) %>%
  select(Gene_id, Is_3D7B, Is_tRNA, Annot, contains('DuplDel')) %>%
  left_join(gm_pos_maxtrans) %>%
  write_csv(paste0(tables_dir, 'gene_model_heatmap_lost_genes.csv'))



## Check output
x$elbow
x$silhouette
plot(x$plot)

#### Gene-Model Analysis: K-Means Heatmap by Families ####
## Plot by Family

set.seed(123)

fam_plot <- heatMap_families_kmeans(
  df = gm_pos_maxtrans %>%
    mutate(Label = paste(Gene_id, SubFamily, sep = ': ')),
  aFC_th = 0,
  fam = 'CLAG',
  fam_facet = F,
  nclu = 1,
  tbar = T,
  fbar = F,
  labels = T
)

fam_plot$elbow
fam_plot$silhouette
fam_plot$df
plot(fam_plot$plot)


ggsave(
  paste0(plots_dir, 'Gene_Model/Family_Heatmaps/CLAG', '_aMAFC_0', '_k', 1, '_withlegend.pdf'),
  fam_plot$plot, width = 60, height = 20,
  units = 'cm', device = 'pdf', limitsize = F
)

## Current famillies and k's
## fams_ks <- list(
##   c('ACS', 3),
##   c('CLAG', 2),
##   c('HYP', 3),
##   c('OTHER', 7),
##   c('PFMC-2TM', 3),
##   c('PHIST', 6),
##   c('RIFIN', 5),
##   c('STEVOR', 3),
##   c('VAR', 5)
## )

fams_ks <- list(
  c('ACS', 1),
  c('CLAG', 1),
  c('HYP', 1),
  c('OTHER', 1),
  c('PFMC-2TM', 1),
  c('PHIST', 1),
  c('RIFIN', 1),
  c('STEVOR', 1),
  c('VAR', 1)
)

for (fam_k in fams_ks){

  print(fam_k)

  nrows <- gm_pos_maxtrans %>%
    filter(Family == fam_k[1]) %>%
    count(Family) %>%
    pull(n)

  fam_plot <- heatMap_families_kmeans(
    df = gm_pos_maxtrans %>%
      mutate(Label = paste(Gene_id, SubFamily, sep = ': ')),
    aFC_th = 0,
    fam = fam_k[1],
    fam_facet = F,
    nclu = as.numeric(fam_k[2]),
    tbar = T,
    fbar = F,
    labels = F
  )

  ggsave(
    paste0(plots_dir, 'Gene_Model/Family_Heatmaps/', fam_k[1], '_aMAFC_0', '_k', fam_k[2], '_withlegend.pdf'),
    fam_plot$plot, width = 60, height = 0.7*nrows,
    units = 'cm', device = 'pdf', limitsize = F
  )
  ggsave(
    paste0(plots_dir, 'Gene_Model/Family_Heatmaps/', fam_k[1], '_elbow', '.pdf'),
    fam_plot$elbow, device = 'pdf'
  )
  ggsave(
    paste0(plots_dir, 'Gene_Model/Family_Heatmaps/', fam_k[1], '_silhouette', '.pdf'),
    fam_plot$silhouette, device = 'pdf'
  )
  write_csv(
    fam_plot$df %>% select(-contains('bin')),
    paste0(plots_dir, 'Gene_Model/Family_Heatmaps/', fam_k[1], '_table.csv')
    )
}

## Plot the "NA" Family

aMAFC_th <- 2
family_facet <- F
family <- NA
nclust <- 2
df <- gm_pos_maxtrans %>%
  filter(is.na(Family) & !Is_tRNA)
x <- heatMap_allstrains_trans(df, aMAFC_th, family, family_facet)

#### Gene-Model Analysis: Kmeans Cluster Tendency Plots ####
tendency_plot_loess <- function(df, cluster){

  max_y <- max(df %>% select(contains('_On') | contains('_Off')), na.rm = T)
  min_y <- min(df %>% select(contains('_On') | contains('_Off')), na.rm = T)

  plt_df <-  df %>%
    select(Gene_id, contains('On'), contains('Off'), Cluster) %>%
    filter(Cluster == cluster) %>%
    select(-Cluster)

  plot_mdf <- melt(plt_df)
  head(plot_mdf)

  plot_mdf['State'] <- sapply(plot_mdf$variable, function(x) str_split(x, '_')[[1]][2])
  plot_mdf <- plot_mdf %>%
    mutate(variable = as.numeric(sub("bin(\\d+).*", "\\1", variable)))

  head(plot_mdf)

  p <- ggplot(plot_mdf, aes(x = variable, y = value, group = State))
  p <- p + geom_smooth(aes(color = State), method = 'loess', level = 0.95)
  p <- p + scale_color_manual(values = c('red', 'green'))
  p <- p + ggtitle(paste0('Cluster ', cluster))
  p <- p + geom_vline(xintercept = 4.5, linetype="solid")
  p <- p + geom_vline(xintercept = 9.5, linetype="dotted")
  p <- p + geom_vline(xintercept = 14.5, linetype="dotted")
  p <- p + geom_vline(xintercept = 19.5, linetype="solid")
  p <- p + theme_classic()
  p <- p + theme(axis.title.x = element_blank())
  p <- p + ylab('H3K9me3 enrichment')
  p <- p + scale_x_continuous(
             breaks=c(2.5,4.5,7,12,17,19.5, 21.5),
             labels=c('Prev\nGenes', "ATG", "5'UTR", "CDS", "3'UTR", "END", 'Post\nGenes'),
             expand = c(0, 0)
           )
  ##p <- p + scale_y_continuous(limits = c(min_y-0.3, max_y+0.3), expand = c(0, 0))
  ##p <- p + scale_y_continuous(limits = c(-0.3, 3), expand = c(0, 0))
  p <- p + coord_cartesian(ylim = c(-0.3, 4))
  p
}

## tendency_plot_mean_sd <- function(df, cluster){

##   onoff_df <- df %>%
##     select(-contains('_Dif'))

##   mean_df <- onoff_df %>%
##     group_by(Cluster) %>%
##     summarize(across(contains('bin'), list(mean))) %>%
##     pivot_longer(!Cluster, values_to = 'mean')

##   sd_df <- onoff_df %>%
##     group_by(Cluster) %>%
##     summarize(across(contains('bin'), list(sd))) %>%
##     pivot_longer(!Cluster, values_to = 'sd')

##   plot_df <- mean_df %>%
##     full_join(sd_df, by = c('Cluster', 'name'))

##   max_y <- max(plot_df$mean + plot_df$sd)
##   min_y <- min(plot_df$mean - plot_df$sd)

##   plot_df <- plot_df %>%
##     filter(Cluster == cluster) %>%
##     select(-Cluster)

##   plot_df['State'] <- sapply(plot_df$name, function(x) str_split(x, '_')[[1]][2])
##   plot_mdf <- plot_df %>%
##     mutate(variable = as.numeric(sub("bin(\\d+).*", "\\1", name)))

##   #print(min_y)

##   p <- ggplot(plot_mdf, aes(x = variable, group = State))
##   p <- p + geom_line(aes(y = mean, color = State), size = 1)
##   p <- p + geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = State),
##                        alpha = .2)
##   p <- p + scale_color_manual(values = c('red', 'green'))
##   p <- p + ggtitle(paste0('Cluster ', cluster))
##   p <- p + geom_vline(xintercept = 4.5, linetype="solid")
##   p <- p + geom_vline(xintercept = 9.5, linetype="dotted")
##   p <- p + geom_vline(xintercept = 14.5, linetype="dotted")
##   p <- p + geom_vline(xintercept = 19.5, linetype="solid")
##   p <- p + theme_classic()
##   p <- p + theme(axis.title.x = element_blank())
##   p <- p + ylab('H3K9me3 enrichment')
##   p <- p + scale_x_continuous(
##              breaks=c(2.5,4.5,7,12,17,19.5, 21.5),
##              labels=c('Prev\nGenes', "ATG", "5'UTR", "CDS", "3'UTR", "END", 'Post\nGenes'),
##              expand = c(0, 0)
##            )
##   p <- p + scale_y_continuous(limits = c(min_y, max_y), expand = c(0, 0))
##   p
## }


## 1 cluster plots
tendency_plot_loess(all_gmodel_kmeans, 6)
##tendency_plot_mean_sd(all_gmodel_kmeans, 6)


## All cluster plots, overlapped
## K-means loess
plots <- lapply(1:nclu, function(x) tendency_plot_loess(all_gmodel_kmeans, x))
clusters_plot <- grid.arrange(grobs = plots, nrow = length(plots), ncol = 1)

ggsave(
  paste0(plots_dir, 'Gene_Model/clusters_plot_overlapped_with_neighbors_kmeans_new.pdf'),
  clusters_plot,
  height = 80, width = 15, units = 'cm',
  device = 'pdf'
)

## ## K-means means +- sd
## plots <- lapply(1:nclust, function(x) tendency_plot_mean_sd(all_gmodel_kmeans, x))
## clusters_plot <- grid.arrange(grobs = plots, nrow = length(plots), ncol = 1)

## ggsave(
##   paste0(plots_dir, 'Gene_Model/clusters_plot_overlapped_with_neighbors_kmeans_means_sd.pdf'),
##   clusters_plot,
##   height = 80, width = 15, units = 'cm',
##   device = 'pdf'
## )

#### Gene-Model Analysis: K-Means Box and whiskers Plots ####
box_df <- all_gmodel_kmeans %>%
  select(Gene_id, Max_aMAFC, Cluster)

p <- ggplot(all_gmodel_kmeans, aes(x = as.character(Cluster), y= Max_aMAFC))
p <- p + stat_boxplot(geom = "errorbar", width = 0.5) #add perpendicular whiskers
p <- p + geom_boxplot(outlier.colour="red")
p <- p + theme_classic()
p <- p + xlab('Cluster')
p <- p + ylab('Level of expression (aMAFC)')
p

ggsave(paste0(plots_dir, 'Gene_Model/cluster_boxplots.pdf'), p, device = 'pdf')

#### Gene-Model Analysis: K-Means Cluster Donuts ####
## Donut Plot

donut_plot <- function(df, cluster) {
    x <- df %>%
      filter(Cluster == cluster)

    data <- as.data.frame(table(x$Family_Grouped))
    colnames(data) <- c('category', 'count')

    count(info_df, Family)

    ## Compute percentages
    data$fraction = data$count / sum(data$count)

    ## Compute the cumulative percentages (top of each rectangle)
    data$ymax = cumsum(data$fraction)

    ## Compute the bottom of each rectangle
    data$ymin = c(0, head(data$ymax, n=-1))

    ## Compute label position
    data$labelPosition <- (data$ymax + data$ymin) / 2

    ## Compute a good label
    data$label <- paste0(data$category, "\n value: ", data$count)

    ## Make the plot
    p <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
      geom_rect(color="white") +
      #geom_label( x=4.2, aes(y=labelPosition, label=label), size=3) +
      coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
      xlim(c(2, 4)) +# Try to remove that to see how to make a pie chart
      theme_void() + my_scale

    outname <- paste0('./Plots/Donuts/gmodel_hetmap_cluster_',
                      as.character(cluster),
                      '.svg')
    #ggsave(outname, p, device = 'svg')
    p
}

count(all_gmodel_kmeans, Cluster)

p1 <- donut_plot(all_gmodel_kmeans, 1)
p2 <- donut_plot(all_gmodel_kmeans, 2)
p3 <- donut_plot(all_gmodel_kmeans, 3)
p4 <- donut_plot(all_gmodel_kmeans, 4)
p5 <- donut_plot(all_gmodel_kmeans, 5)
p6 <- donut_plot(all_gmodel_kmeans, 6)
p7 <- donut_plot(all_gmodel_kmeans, 7)

all_donuts <- grid.arrange(p1, p2, p3, p4, p5, p6, p7,  nrow = 7, ncol = 1)
ggsave(paste0(plots_dir, 'Gene_Model/gmodel_heatmap_alldonuts.svg'), all_donuts, device = 'svg')


## Plot all genes to get legend
## all_genes_donut <- info_df %>%
##   mutate(Family_Dount = case_when(
##            Family %in% col_factor ~ Family,
##            Family == 'OTHER' ~ 'Other CVGs',
##            is.na(Family) ~ 'Not CVGs'))


data <- as.data.frame(table(info_df$Family_Grouped))
colnames(data) <- c('category', 'count')

count(info_df, Family)

## Compute percentages
data$fraction = data$count / sum(data$count)

## Compute the cumulative percentages (top of each rectangle)
data$ymax = cumsum(data$fraction)

## Compute the bottom of each rectangle
data$ymin = c(0, head(data$ymax, n=-1))

## Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

## Compute a good label
data$label <- paste0(data$category, "\n value: ", data$count)

## Make the plot
p <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect(color="white") +
                                        #geom_label( x=4.2, aes(y=labelPosition, label=label), size=3) +
  coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
  xlim(c(2, 4)) +# Try to remove that to see how to make a pie chart
  theme_void() + my_scale

ggsave(paste0(plots_dir, 'families_donut_for_legend.svg'), p, device = 'svg')
##p

#### Gene State Tables Analysis ####

## Check for which genes we have a strain in an active state and a strain in a silenced state
## Load whole gene state tables
state_old_whole <- read_tsv('/mnt/Disc4T/Projects/Active_gene_list/Results_Tables/state_df_old_arrays_rna25_red_up50_red_dw25_reddw15_areaFC1_areaFCbig3_redpcntdif30_maxtime1_dupldel_filtered.csv') %>%
  select(-Variant, -Gene_name, -Annot)

state_new_whole <- read_tsv('/mnt/Disc4T/Projects/Active_gene_list/New_Arrays_Results/state_df_new_arrays__rna25_red_up50_red_dw25_reddw15_areaFC1_areaFCbig3_redpcntdif30_maxtime1_dupldel_filtered.csv') %>%
  select(
    -Variant, -Gene_name, -Annot,
    -Chrom, -Start, -Stop, -Type, -Strand, -Name,
    -Is_tRNA, -Family, -SubFamily, -Label,
    -Family_Grouped, -Gam_specific,
    -Difpeaks_12B_10G, -Difpeaks_A7_E5, -Difpeaks_A7_B11, -Difpeaks_E5_B11
  )


state_whole <- state_old_whole %>%
  full_join(state_new_whole, by = 'Gene_id', suffix = c('_old', '_new'))


## ## Different state in 3D7B and A7/E5/B11
myMax <- function(x){
  if (all(is.na(x))){
    out <- NA
  } else {
    out <- max(x, na.rm=T)
  }
  return(out)
}

myMin <- function(x){
  if (all(is.na(x))){
    out <- NA
  } else {
    out <- min(x, na.rm=T)
  }
  return(out)
}

## Load gene expression tables
genexp_old <- read_csv(paste0(microarrays_dir, 'Old_Arrays/R_results_OldArrays_Variantome/geneLevel_exp.csv')) %>%
  select(Gene_id, contains('tp'))

genexp_old <- genexp_old %>%
  mutate(MaxExp_12B = apply(select(., contains('12B')), 1, myMax)) %>%
  mutate(MaxExp_10G = apply(select(., contains('10G')), 1, myMax)) %>%
  mutate(MaxExp_3D7B = apply(select(., contains('3D7B')), 1, myMax)) %>%
  mutate(MaxExp_old = apply(select(., contains('MaxExp')), 1, myMax)) %>%
  mutate(MinExp_old = apply(select(., contains('MaxExp')), 1, myMin))

genexp_new <- read_csv(paste0(microarrays_dir, 'New_Arrays/R_results_NewArray/geneLevel_exp.csv')) %>%
  select(Gene_id, contains('tp'))

genexp_new <- genexp_new %>%
  mutate(MaxExp_A7 = apply(select(., contains('A7')), 1, myMax)) %>%
  mutate(MaxExp_E5 = apply(select(., contains('E5')), 1, myMax)) %>%
  mutate(MaxExp_B11 = apply(select(., contains('B11')), 1, myMax)) %>%
  mutate(MaxExp_new = apply(select(., contains('MaxExp')), 1, myMax)) %>%
  mutate(MinExp_new = apply(select(., contains('MaxExp')), 1, myMin))

genexp_df <- genexp_old %>% select(Gene_id, contains('Max'), contains('Min')) %>%
  full_join(genexp_new %>% select(Gene_id, contains('Max'), contains('Min')), by = 'Gene_id')

genexp_df

state_exp <- state_whole %>%
  left_join(genexp_df, by = 'Gene_id')

## Genes in which 3D7B is in a different state than A7/E5/B11
## Repressed in 3D7B but active in A7/E5/B11
bigLogDif_dw <- 3

## Filter
inactive_A7 <- state_exp %>%
  filter(category_3D7B == 'CVG_Silenced') %>%
  filter(category_A7 == 'CVG_Active') %>%
  filter((MaxExp_old - MaxExp_A7) > bigLogDif_dw) %>%
  select(Gene_id) %>%
  pull()

inactive_E5 <- state_exp %>%
  filter(category_3D7B == 'CVG_Silenced') %>%
  filter(category_E5 == 'CVG_Active') %>%
  filter((MaxExp_old - MaxExp_E5) > bigLogDif_dw) %>%
  select(Gene_id) %>%
  pull()

inactive_B11 <- state_exp %>%
  filter(category_3D7B == 'CVG_Silenced') %>%
  filter(category_B11 == 'CVG_Active') %>%
  filter((MaxExp_old - MaxExp_B11) > bigLogDif_dw) %>%
  select(Gene_id) %>%
  pull()

## Replace
state_exp <- state_exp %>%
  mutate(state_A7 = ifelse(
           Gene_id %in% inactive_A7,
           gsub('Active', 'Undetermined', state_A7),
           state_A7
         )) %>%
  mutate(category_A7 = ifelse(
           Gene_id %in% inactive_A7,
           gsub('Active', 'Undetermined', category_A7),
           category_A7
         ))

state_exp <- state_exp %>%
  mutate(state_E5 = ifelse(
           Gene_id %in% inactive_E5,
           gsub('Active', 'Undetermined', state_E5),
           state_E5
         )) %>%
  mutate(category_E5 = ifelse(
           Gene_id %in% inactive_E5,
           gsub('Active', 'Undetermined', category_E5),
           category_E5
         ))

state_exp <- state_exp %>%
  mutate(state_B11 = ifelse(
           Gene_id %in% inactive_B11,
           gsub('Active', 'Undetermined', state_B11),
           state_B11
         )) %>%
  mutate(category_B11 = ifelse(
           Gene_id %in% inactive_B11,
           gsub('Active', 'Undetermined', category_B11),
           category_B11
         ))


## Actius a 3D7B i inactius a A7/E5/B11
## Per soca: actius a 3D7B però inactiu (1.2B/10G) molt més baix que A7/E5/B11
## Treure només els 'positius'
## bigLogDif_up <- 3

## ## Filter
## active_A7 <- state_exp %>%
##   filter(category_3D7B == 'CVG_Active') %>%
##   filter(category_A7 == 'CVG_Repressed') %>%
##   filter((MaxExp_A7 - MinExp_old) > bigLogDif_up) %>%
##   select(Gene_id) %>%
##   pull()

## active_E5 <- state_exp %>%
##   filter(category_3D7B == 'CVG_Active') %>%
##   filter(category_E5 == 'CVG_Repressed') %>%
##   filter((MaxExp_E5 - MinExp_old) > bigLogDif_up) %>%
##   select(Gene_id) %>%
##   pull()

## active_B11 <- state_exp %>%
##   filter(category_3D7B == 'CVG_Active') %>%
##   filter(category_B11 == 'CVG_Repressed') %>%
##   filter((MaxExp_B11 - MinExp_old) > bigLogDif_up) %>%
##   select(Gene_id) %>%
##   pull()


## state_exp %>%
##   select(Gene_id, MaxExp_old, MinExp_old, MaxExp_A7) %>%
##   filter(Gene_id == 'PF3D7_0424900')


## state_exp


## Replace

## Replace
## state_exp <- state_exp %>%
##   mutate(state_A7 = ifelse(
##            Gene_id %in% active_A7,
##            gsub('Repressed', 'Undetermined', state_A7),
##            state_A7
##          )) %>%
##   mutate(category_A7 = ifelse(
##            Gene_id %in% active_A7,
##            gsub('Repressed', 'Undetermined', category_A7),
##            category_A7
##          ))

## state_exp <- state_exp %>%
##   mutate(state_E5 = ifelse(
##            Gene_id %in% active_E5,
##            gsub('Repressed', 'Undetermined', state_E5),
##            state_E5
##          )) %>%
##   mutate(category_E5 = ifelse(
##            Gene_id %in% active_E5,
##            gsub('Repressed', 'Undetermined', category_E5),
##            category_E5
##          ))

## state_exp <- state_exp %>%
##   mutate(state_B11 = ifelse(
##            Gene_id %in% active_B11,
##            gsub('Repressed', 'Undetermined', state_B11),
##            state_B11
##          )) %>%
##   mutate(category_B11 = ifelse(
##            Gene_id %in% active_B11,
##            gsub('Repressed', 'Undetermined', category_B11),
##            category_B11
##          ))

## Load gene state tables
my_difstate_filter <- function(state_vect){
  active <- any(grepl('Active', state_vect))
  silenced <- any(grepl('CVG_Silenced', state_vect)) | any(grepl('Inactive', state_vect))
  return(active & silenced)
}

state_exp %>%
  select(contains('category'))

colnames(state_exp)

state_final <- state_exp %>%
  select(-contains('3D7B')) %>%
  rename_with(~ gsub('category_', 'Gene_State_', .x)) %>%
  mutate(Gene_id = ifelse(Gene_id == 'PF3D7_0935400_as', 'PF3D7_0935390', Gene_id)) %>%
  mutate(DifState = apply(select(., contains('State')), 1, my_difstate_filter)) %>%
  filter(Gene_id %in% info_df$Gene_id) %>%
  left_join(info_df, by = 'Gene_id') %>%
  mutate(Gene_State_12B = case_when(
           is.na(Gene_State_12B) & !Variant ~ 'Undetermined',
           is.na(Gene_State_12B) & Variant ~ 'CVG_Undetermined',
           TRUE ~ Gene_State_12B
         )) %>%
  mutate(Gene_State_10G = case_when(
           is.na(Gene_State_10G) & !Variant ~ 'Undetermined',
           is.na(Gene_State_10G) & Variant ~ 'CVG_Undetermined',
           TRUE ~ Gene_State_10G
         )) %>%
  mutate(Gene_State_A7 = case_when(
           is.na(Gene_State_A7) & !Variant ~ 'Undetermined',
           is.na(Gene_State_A7) & Variant ~ 'CVG_Undetermined',
           TRUE ~ Gene_State_A7
         )) %>%
  mutate(Gene_State_E5 = case_when(
           is.na(Gene_State_E5) & !Variant ~ 'Undetermined',
           is.na(Gene_State_E5) & Variant ~ 'CVG_Undetermined',
           TRUE ~ Gene_State_E5
         )) %>%
  mutate(Gene_State_B11 = case_when(
           is.na(Gene_State_B11) & !Variant ~ 'Undetermined',
           is.na(Gene_State_B11) & Variant ~ 'CVG_Undetermined',
           TRUE ~ Gene_State_B11
         ))

state_final %>%
  write_tsv(paste0(tables_dir, 'join_gene_states.tsv'))

## Parse for supplementary
colnames(state_final)

supp_state <- state_final %>%
  select(Gene_id, contains('Gene_State'), Variant, Annot, everything()) %>%
  rename(CVG = Variant) %>%
  write_tsv(paste0(tables_dir, 'join_gene_states_for_supp.tsv'))

colnames(supp_state)


colnames(state_final)

stinfo_df <- state_final %>%
  left_join(info_df) %>%
  filter(DifState) %>%
  select(Gene_id, contains('State'), Variant, Annot) %>%
  print(width = 200)

stinfo_df %>%
  filter(DifState) %>%
  count(Variant)

vars <- stinfo_df %>%
  filter(DifState) %>%
  filter(Variant) %>%
  select(Gene_id, contains('State'), Annot) %>%
  print(width = 200)

novars <- stinfo_df %>%
  filter(DifState) %>%
  filter(!Variant) %>%
  select(Gene_id, contains('State'), Annot) %>%
  print(width = 200)

vh <- gm_pos_maxtrans %>%
  filter(Gene_id %in% vars$Gene_id)

nvh <- gm_pos_maxtrans %>%
  filter(Gene_id %in% novars$Gene_id)

x <- heatMap_allstrains_kmeans(
  df = vh,
  aFC_th = 0,
  fam = NA,
  fam_facet = F,
  nclu = 6,
  tbar = F,
  fbar = T,
  labels = T
)

x$silhouette

ggsave(
  paste0(plots_dir, 'Gene_Model/heatmap_dif_gene_states_variant.pdf'),
  x$plot,
  device = 'pdf',
  height = 80, width = 60, units = 'cm'
)

## Genereate empty df for each out df (on/off/diff)
hetcols <- gmodel_all %>%
  select(contains('12B') & contains('bin')) %>%
  colnames(.) %>%
  gsub('12B', '', ., fixed=TRUE)

col_names <- c(
  'Gene_id',
  'Strain',
  paste0(hetcols, 'Cov')
  )

cn <- setNames(rep('', length(col_names)), col_names)
df <- bind_rows(cn)[0,]

df_on <- df
df_off <- df
df_diff <- df

vars

for (i in 1:dim(vars)[1]){
  #i <- 1
  gid <- as.character(vars$Gene_id[i])
  states <- vars[i,] %>% select(-Gene_id)
  on_cols <- colnames(states[which(states == 'CVG_Active')])
  off_cols <- colnames(states[which(states == 'CVG_Silenced')])
  on_strains <- gsub('Gene_State_', '', on_cols)
  off_strains <- gsub('Gene_State_', '', off_cols)

  ## Get Gene-model bins data
  ## On_df
  for (strain in on_strains){

    vect <- gmodel_all %>%
      filter(Gene_id == gid) %>%
      select(contains(strain))

    row <- c(
      gid,
      strain,
      unlist(vect, use.names = F)
      )

    row <- setNames(row, col_names)
    df_on <- df_on %>% add_row(bind_rows(row))
  }
  ## Off_df
  for (strain in off_strains){

    vect <- gmodel_all %>%
      filter(Gene_id == gid) %>%
      select(contains(strain))

    row <- c(
      gid,
      strain,
      unlist(vect, use.names = F)
      )

    row <- setNames(row, col_names)
    df_off <- df_off %>% add_row(bind_rows(row))
  }
}
## Convert to numeric
df_on <- df_on %>% mutate(across(contains('bin'), as.numeric))
df_off <- df_off %>% mutate(across(contains('bin'), as.numeric))


## target <- 'PF3D7_1301500'

## df_off %>%
##   filter(Gene_id == target) %>%
##   print(width = 90000)

## vars %>%
##   filter(Gene_id == target) %>%
##   print(width = 90000)

## trans_df %>%
##   filter(Gene_id == target) %>%
##   select(Gene_id, contains('MaxVal')) %>%
##   print(width = 90000)


## Filter out genes with only one ON / OFF strain
repeated_on <- df_on %>%
  group_by(Gene_id) %>%
  summarise(N_strains = n()) %>%
  filter(N_strains > 1) %>%
  select(Gene_id) %>%
  pull()

repeated_off <- df_off %>%
  group_by(Gene_id) %>%
  summarise(N_strains = n()) %>%
  filter(N_strains > 1) %>%
  select(Gene_id) %>%
  pull()


## Calculate means and substract
on_means <- df_on %>%
  group_by(Gene_id) %>%
  summarise(
    across(contains('bin'), ~ mean(.x, na.rn = T))
  )

off_means <- df_off %>%
  group_by(Gene_id) %>%
  summarise(
    across(contains('bin'), ~ mean(.x, na.rn = T))
  )

onoff_means <- bind_cols(
  on_means %>% select(Gene_id),
  on_means %>% select(-Gene_id) - off_means %>% select(-Gene_id)
)

onmean_mtx <- df_on %>%
  left_join(on_means, by = 'Gene_id', suffix = c('', '_mean')) %>%
  select(contains('mean'))

on_mtx <- df_on %>%
  select(contains('bin'))

offmean_mtx <- df_off %>%
  left_join(off_means, by = 'Gene_id', suffix = c('', '_mean')) %>%
  select(contains('mean'))

off_mtx <- df_off %>%
  select(contains('bin'))

## Create mean subtracted DFs and remove rows with a single strain per gene
onmean_df <- bind_cols(df_on %>% select(-contains('bin')), on_mtx - onmean_mtx) %>%
  filter(Gene_id %in% repeated_on)

offmean_df <- bind_cols(df_off %>% select(-contains('bin')), off_mtx - offmean_mtx) %>%
  filter(Gene_id %in% repeated_off)

## Make Heatmaps

my_difmean_heat <- function(df){

  df_heat <- df %>%
    mutate(N_id = row_number()) %>%
    left_join(info_df) %>%
    mutate(N_id = paste0(N_id, ': ', Label)) %>%
    select(
      -contains('bin1_'),
      -contains('bin2_'),
      -contains('bin3_'),
      -contains('bin4_'),
      -contains('bin20_'),
      -contains('bin21_'),
      -contains('bin22_'),
      -contains('bin23_'),
      )

  bin_cols <- df_heat %>%
    select(contains('bin')) %>%
    colnames

  nids <- df_heat$N_id

  write_csv(
    df_heat,
    paste0(
      plots_dir,
      'Gene_Model/States_OnOff/',
      deparse(substitute(df)),
      '.csv'
    )
  )

  mdf <- df_heat %>%
    pivot_longer(contains('bin'), names_to = 'variable', values_to = 'value') %>%
    mutate(variable = factor(variable, levels = bin_cols)) %>%
    mutate(N_id = factor(N_id, levels = rev(nids)))

  p <- ggplot(mdf, aes(x = variable, y = N_id, fill = value))
  p <- p + geom_tile(colour="snow3")
  p <- p + theme(
             #axis.text.y = element_blank(),
             #axis.ticks.y = element_blank(),

             axis.title = element_blank(),
             axis.line.x = element_blank(),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             axis.title.x = element_blank(),


             plot.background=element_blank(),

             panel.border=element_blank(),
             panel.grid.major=element_blank(),
             panel.grid.minor=element_blank(),
             panel.background=element_blank(),

             strip.background = element_blank(),
             strip.text.x = element_blank(),
             strip.text.y = element_blank(),

             text=element_text(size=24),

             legend.position='bottom',
             legend.title = element_blank()
             )
  p <- p + scale_fill_gradient2(
             low = "chartreuse3",
             mid = "white",
             high = "darkred",
             na.value="grey90",
             limits = c(-3,3)
           )
  p
}


pon <- my_difmean_heat(onmean_df)
poff <- my_difmean_heat(offmean_df)
ponoff <- my_difmean_heat(onoff_means)

pon
poff
ponoff

ggsave(
  paste0(plots_dir, 'Gene_Model/States_OnOff/genes_on_meandif_withnegatives.pdf'),
  pon, device = 'pdf', height = 80, width = 60, units = 'cm'
)
ggsave(
  paste0(plots_dir, 'Gene_Model/States_OnOff/genes_off_meandif_withnegatives.pdf'),
  poff, device = 'pdf', height = 80, width = 60, units = 'cm'
)
ggsave(
  paste0(plots_dir, 'Gene_Model/States_OnOff/genes_on_minus_ff_withnegatives.pdf'),
  ponoff, device = 'pdf', height = 80, width = 60, units = 'cm'
)

my_difmean_heat(df_on) +
  scale_fill_gradient(
    low = "white",
    high = "darkred",
    na.value="grey90"
  )
my_difmean_heat(df_off) +
  scale_fill_gradient(
    low = "white",
    high = "darkred",
    na.value="grey90"
  )

## Make Correlation Boxplots

table(unique(df_on$Gene_id) == unique(df_off$Gene_id))

## In ON/OFF states, get correlation between gene-model coverage
## of the same gene in all the different strains that share a state
get_cors <- function(df) {
  cors_out <- c()
  for (gid in unique(df$Gene_id)) {
    mtx <- df %>%
      filter(Gene_id == gid) %>%
      select(contains('bin')) %>%
      select(
        -contains('bin1_'),
        -contains('bin2_'),
        -contains('bin3_'),
        -contains('bin4_'),
        -contains('bin20_'),
        -contains('bin21_'),
        -contains('bin22_'),
        -contains('bin23_'),
        ) %>%
      t()

    cors <-cor(mtx)
    cors[col(cors)==row(cors)] <- NA
    mean_cor <- mean(c(cors), na.rm = T)
    cors_out <- c(cors_out, mean_cor)
  }
  return(cors_out)
}

on_cors <- get_cors(df_on)
off_cors <- get_cors(df_off)

## For the same genes, get correlation between the gene-model mean coverage
## between active and silenced state

on_means_noneighbors <- on_means %>%
  select(contains('bin')) %>%
  select(
    -contains('bin1_'),
    -contains('bin2_'),
    -contains('bin3_'),
    -contains('bin4_'),
    -contains('bin20_'),
    -contains('bin21_'),
    -contains('bin22_'),
    -contains('bin23_'),
    )

off_means_noneighbors <- off_means %>%
  select(contains('bin')) %>%
  select(
    -contains('bin1_'),
    -contains('bin2_'),
    -contains('bin3_'),
    -contains('bin4_'),
    -contains('bin20_'),
    -contains('bin21_'),
    -contains('bin22_'),
    -contains('bin23_'),
    )


onoff_cors <- c()
for (i in 1:dim(on_means_noneighbors)[1]){
  x <- on_means_noneighbors[i,] %>% t()
  y <- off_means_noneighbors[i,] %>% t()
  c <- cor(x, y)
  onoff_cors <- c(onoff_cors, c)
}
onoff_cors

### Make graphic

mean(on_cors, na.rm = T)
mean(off_cors, na.rm = T)
mean(onoff_cors)

length(on_cors)
length(off_cors)
length(onoff_cors)

cor_df <- tibble(
  On_Cors = on_cors,
  Off_Cors = off_cors,
  OnOff_Cors = onoff_cors
)

mean(cor_df$On_Cors, na.rm = T)
mean(cor_df$Off_Cors, na.rm = T)
mean(cor_df$OnOff_Cors, na.rm = T)


cor_mdf <- melt(cor_df)

p <- ggplot(cor_mdf, aes(x = variable, y = value))
p <- p + geom_boxplot(outlier.shape = NA, aes(ymax = quantile(value, 0.95, na.rm = T),
                                              ymin = quantile(value, 0.05, na.rm = T)))
p <- p + geom_jitter()
p <- p + coord_cartesian(ylim = c(0.7, 1))
p
ggsave(
  paste0(plots_dir, 'similarity_on_off_boxplot.pdf'),
  device = 'pdf'
)

#### MaxTransDif filtered genes, same state genes vs difstate genes ####
set.seed(123)

table(finalFC_df$Is_3D7B)
table(finalFC_df$Is_tRNA)
table(finalFC_df$PassMaxtime)
colnames(finalFC_df)


variant <- info_df %>%
  filter(Variant) %>%
  select(Gene_id) %>%
  pull()

## Subset DFs of interest
difgenes_states <- state_final %>%
  select(Gene_id, contains('Gene_State')) %>%
  filter(Gene_id %in% finalFC_df$Gene_id) %>%
  filter(Gene_id %in% variant)

## Remove "neighboring" genes
het_bins <- gmodel_all_pos %>%
  select(
    -contains('bin1_'),
    -contains('bin2_'),
    -contains('bin3_'),
    -contains('bin4_'),
    -contains('bin20_'),
    -contains('bin21_'),
    -contains('bin22_'),
    -contains('bin23_'),
    )

## Create empty DF for result
col_names <- c(
  'Original_Row',
  'Gene_id',
  'Same_Strain_1',
  'Same_Strain_2',
  'Dif_Strain',
  paste0('SameState_', c(1:15)),
  paste0('DifState_', c(1:15))
)
cn <- setNames(rep('', length(col_names)), col_names)
df <- bind_rows(cn)[0,]

## Main "for loop"
for (i in 1:dim(difgenes_states)[1]){

  gid <- difgenes_states[i, 'Gene_id'] %>% pull()

  off <- difgenes_states[i,] %>%
    pivot_longer(-Gene_id) %>%
    filter(grepl('Silenced', value, fixed = T) | grepl('Inactive', value, fixed = T)) %>%
    select(name) %>%
    pull()

  on <- difgenes_states[i,] %>%
    pivot_longer(-Gene_id) %>%
    filter(grepl('Active', value, fixed = T)) %>%
    select(name) %>%
    pull()

  off <- gsub('Gene_State_', '', off)
  on <- gsub('Gene_State_', '', on)

  if (length(off) == 0) {off <- NA}
  if (length(on) == 0) {on <- NA}

  if (any(!is.na(on)) & any(!is.na(off)) & length(c(on, off)) >= 3){

    if (length(on) >= length(off)){
      nodif <- 'on'
    } else {
      nodif <- 'off'
    }

    if (nodif == 'on'){
      on_strains <- sample(on, 2)
      off_strain <- sample(off, 1)

      onvect1 <- het_bins %>%
        filter(Gene_id == gid) %>%
        select(contains(on_strains[1]))

      onvect2 <- het_bins %>%
        filter(Gene_id == gid) %>%
        select(contains(on_strains[2]))

      offvect <- het_bins %>%
        filter(Gene_id == gid) %>%
        select(contains(off_strain))

      nodif_vect <- onvect1 - onvect2
      difvect <- onvect1 - offvect
      same_strains <- on_strains
      difstrain <- off_strain
    } else {
      off_strains <- sample(off, 2)
      on_strain <- sample(on, 1)

      offvect1 <- het_bins %>%
        filter(Gene_id == gid) %>%
        select(contains(off_strains[1]))

      offvect2 <- het_bins %>%
        filter(Gene_id == gid) %>%
        select(contains(off_strains[2]))

      onvect <- het_bins %>%
        filter(Gene_id == gid) %>%
        select(contains(on_strain))

      nodif_vect <- offvect1 - offvect2
      difvect <- onvect - offvect1
      same_strains <- off_strains
      difstrain <- on_strain
    }

    row <- c(
      i,
      gid,
      same_strains[1],
      same_strains[2],
      difstrain,
      unlist(nodif_vect, use.names = F),
      unlist(difvect, use.names = F)
    )
    row <- setNames(row, col_names)
    df <- df %>% add_row(bind_rows(row))
  }
}


## Convert to numeric
same_state_df <- df %>%
  mutate(across(contains('State'), as.numeric)) %>%
  left_join(info_df)

same_state_df %>%
  select(-Original_Row) %>%
  left_join(difgenes_states) %>%
  select(-Is_tRNA) %>%
  select(Gene_id, contains('Gene_State'), everything()) %>%
  rename(CVG = Variant) %>%
  write_csv(paste0(tables_dir, 'same_sate_vs_different_state_new_new.csv'))

## Make Heatmap
ss <- same_state_df %>%
  select(
    Gene_id,
    contains('SameState'),
    contains('DifState'),
    Gam_specific,
  ) %>%
  pivot_longer(-Gene_id) %>%
  mutate(name = factor(name, levels = unique(name))) %>%
  mutate(Gene_id = factor(Gene_id, levels = rev(unique(Gene_id))))


x <- ss %>%
  filter(grepl('State', name)) %>%
  mutate(
    SameState = ifelse(grepl('SameState', name), 'Same State', 'Different State')
  ) %>%
  mutate(SameState = factor(SameState, levels=(c('Same State', 'Different State'))))

y <- ss %>%
  filter(name == 'Gam_specific')

p <- ggplot(x, aes(x = name, y = Gene_id, fill = value))
p <- p + geom_tile(colour="snow3")
p <- p + scale_fill_gradient2(
           low = "limegreen",
           mid = "white",
           high = "red",
           midpoint = 0
         )
p <- p + facet_grid(cols = vars(x$SameState), scales = "free", space="free_y")
p <- p + theme(
           axis.text.x = element_blank(),
           axis.ticks.x = element_blank(),
           axis.title.y = element_blank(),
           #axis.text.y = element_blank(),
           axis.ticks.y = element_blank(),
           panel.grid.major=element_blank(),
           panel.background=element_blank(),
           legend.position = 'bottom'
           )
p <- p + labs(y = '', x = '')

p2 <- ggplot(y, aes(x = name, y = Gene_id, fill = value))
p2 <- p2 + geom_tile(colour="snow3")
p2 <- p2 + scale_fill_gradient(
             low = "white",
             high = "Orange",
             )
p2 <- p2 + facet_grid(cols = vars(y$name), scales = "free", space="free_y")
p2 <- p2 + scale_y_discrete(position = "right")
p2 <- p2 + theme(
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             axis.title.y = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks.y = element_blank(),
             panel.grid.major=element_blank(),
             panel.background=element_blank(),
             legend.position = 'bottom'
           )
p2 <- p2 + labs(y = '', x = '')

same_state_heat <- ggarrange(p, p2, nrow = 1, widths = c(10,1))
same_state_heat

ggsave(
  paste0(plots_dir, 'same_dif_state_filteredFC_onlyvariant_on_minus_off.pdf'),
  same_state_heat, device = 'pdf'
)


## Make tendency Plots

same <- ss %>%
  filter(grepl('SameState', name)) %>%
  mutate(State = 'same') %>%
  mutate(value = abs(value))%>%
  ggplot(aes(x = name, y = value, group = State))

same <- same + geom_smooth(aes(color = State), method = 'loess', level = 0.95)
same <- same + theme_classic()
same <- same + coord_cartesian(ylim = c(0,2))
#same <- same + scale_y_continuous(limits = c(0,5), expand = c(0, 0))

dif <- ss %>%
  filter(grepl('DifState', name)) %>%
  mutate(State = 'dif') %>%
  mutate(value = abs(value)) %>%
  ggplot(aes(x = name, y = value, group = State))
dif <- dif + geom_smooth(aes(color = State), method = 'loess', level = 0.95)
dif <- dif + theme_classic()
dif <- dif + coord_cartesian(ylim = c(0,2))
#dif <- dif + scale_y_continuous(limits = c(0,5), expand = c(0, 0))
dif
same_dif_states_tendency <- ggarrange(same, dif, nrow = 1)

ggsave(
  paste0(plots_dir, 'same_dif_state_filteredFC_onlyvariant_tendencies_absval.pdf'),
  device = 'pdf'
  )

#### Metilation and Acetilation Analysis ####

## Load Me and Ac Data

read_ac <- function(strain) {
  strain_dotless = gsub('.', '', strain, fixed = T)

  cov_ac_5p <- read_tsv(paste0(
    coverage_dir,
    'binned_1000tss_0orf_allowoverlaps_False_coverage_acetylation_', strain, '.bed'
  ), col_names = F) %>%
    select(X4, X5) %>%
    set_names(c('Gene_id', 'Ac_5p')) %>%
    mutate(Strain = strain_dotless)

  cov_ac_orf <- read_tsv(paste0(
    coverage_dir,
    'binned_0tss_500orf_allowoverlaps_True_coverage_acetylation_', strain, '.bed'
  ), col_names = F) %>%
    select(X4, X5) %>%
    set_names(c('Gene_id', 'Ac_ORF')) %>%
    mutate(Strain = strain_dotless)

  cov_ac_3p <- read_tsv(paste0(
    coverage_dir,
    'binned_3prime1000_allowoverlaps_False_coverage_acetylation_', strain, '.bed'
  ), col_names = F) %>%
    select(X4, X7) %>%
    set_names(c('Gene_id', 'Ac_3p')) %>%
    mutate(Strain = strain_dotless)

  outdf <- cov_ac_5p %>%
    full_join(cov_ac_orf, join_cols = c('Gene_id', 'Strain')) %>%
    full_join(cov_ac_3p, join_cols = c('Gene_id', 'Strain')) %>%
    select(Gene_id, Strain, everything())

  return(outdf)
}

read_me <- function(strain) {
  strain_dotless = gsub('.', '', strain, fixed = T)

  cov_me_5p <- read_tsv(paste0(
    coverage_dir,
    'binned_1000tss_0orf_allowoverlaps_False_coverage_', strain, '.bed'
  ), col_names = F) %>%
    select(X4, X5) %>%
    set_names(c('Gene_id', 'Me_5p')) %>%
    mutate(Strain = strain_dotless)

  cov_me_orf <- read_tsv(paste0(
    coverage_dir,
    'binned_0tss_500orf_allowoverlaps_True_coverage_', strain, '.bed'
  ), col_names = F) %>%
    select(X4, X5) %>%
    set_names(c('Gene_id', 'Me_ORF')) %>%
    mutate(Strain = strain_dotless)

  cov_me_3p <- read_tsv(paste0(
    coverage_dir,
    'binned_3prime1000_allowoverlaps_False_coverage_', strain, '.bed'
  ), col_names = F) %>%
    select(X4, X7) %>%
    set_names(c('Gene_id', 'Me_3p')) %>%
    mutate(Strain = strain_dotless)

  outdf <- cov_me_5p %>%
    full_join(cov_me_orf, join_cols = c('Gene_id', 'Strain')) %>%
    full_join(cov_me_3p, join_cols = c('Gene_id', 'Strain')) %>%
    select(Gene_id, Strain, everything())

  return(outdf)
}

cov_12b_me <- read_me('1.2B')
cov_10g_me <- read_me('10G')
cov_b11_me <- read_me('B11')

cov_12b_ac <- read_ac('1.2B')
cov_10g_ac <- read_ac('10G')
cov_b11_ac <- read_ac('B11')

me_df <- bind_rows(cov_12b_me, cov_10g_me, cov_b11_me)
ac_df <- bind_rows(cov_12b_ac, cov_10g_ac, cov_b11_ac)
me_ac_df <- me_df %>%
  full_join(ac_df, by = c('Gene_id', 'Strain'))

## Load Gene States table (response/target variable)

labels_df <- state_exp %>%
  select(Gene_id, category_12B, category_10G, category_A7, category_E5, category_B11)

## Table for further analysis

labels_df <- labels_df %>%
  pivot_longer(!Gene_id) %>%
  rename(Label = value) %>%
  mutate(Strain = gsub('category_', '', name, fixed=T)) %>%
  select(-name)

## Join Both

final_df <- me_ac_df %>%
  right_join(labels_df, by=c('Gene_id'='Gene_id', 'Strain'='Strain'))

## Filter out unlabbeled genes and add alpha depending on "Label" (for plots)


final_df <- final_df %>%
  filter(Label != 'No_Category') %>%
  filter(Label != 'Not_Settable') %>%
  filter(Label != 'CVG_Semiactive') %>%
  filter(Label != 'CVG_Undetermined') %>%
  mutate(Alpha = case_when(
           Label == 'Active' ~ 0.1,
           Label == 'Inactive' ~ 0.5,
           TRUE ~ 1
         ))

final_df %>%
  count(Label)

write_tsv(final_df, paste0(tables_dir, 'met_ac_states_df.tsv'))

## Remove rows that represent a gene in the same "state" in different strains (keep first)

final_df <- final_df %>%
  mutate(Label = factor(Label, levels = c('Active', 'Inactive', 'CVG_Active', 'CVG_Silenced')))

unique_df <- final_df %>%
  group_by(Gene_id, Label) %>%
  filter(row_number()==1) %>%
  ungroup()

# Check it worked
unique_df %>%
  group_by(Gene_id, Label) %>%
  count() %>%
  filter(n > 1)

unique_df %>%
  count(Label)

write_tsv(unique_df, paste0(tables_dir, './met_ac_unique_state.tsv'))

## Remove API and mal_mito genes

unique_noApi <- unique_df %>%
  filter(!grepl('API', Gene_id)) %>%
  filter(!grepl('mal_', Gene_id))

write_csv(unique_noApi, paste0(tables_dir, 'unique_noApi.csv'))

## PCAs

set.seed(123)
my_colors = c('#ddb310', '#4053d3',  '#00b25d', '#b51d14')#,  '#fb49b0')
## blue, yellow, red, green, pink

make_pca <- function(df){
  nona_df <- df %>%
    filter(complete.cases(.))

  mtx <- nona_df %>% select(matches('Ac|Me'))
  mtx <- as.matrix(mtx)
  rownames(mtx) <- paste(nona_df$Gene_id, nona_df$Strain, sep = '_')
  pca <- prcomp(mtx)
  pca_df <- tibble(as.data.frame(pca$x)) %>%
    mutate(
      Gid = paste(nona_df$Gene_id, nona_df$Strain, sep = '_'),
      Label = nona_df$Label,
      Alpha = nona_df$Alpha
    ) %>%
    select(Gid, Label, everything())

  cmp1 <- format(round(summary(pca)$importance[2,1]*100, 2), nsmall = 2)
  cmp2 <- format(round(summary(pca)$importance[2,2]*100, 2), nsmall = 2)

  p <- ggplot(pca_df, aes(x=PC1,y=PC2, col = Label, group = Label))
  p <- p + geom_point(alpha = pca_df$Alpha)
  p <- p + labs(x=paste0("PC1: ", cmp1, "%"), y=paste0("PC2: ", cmp2, "%"))
  p <- p + scale_color_manual(values=my_colors)
  #p <- p + theme_minimal() + guides(alpha='none')
  p <- p + theme_bw() + guides(alpha='none')
  p <- p + theme(
             axis.text.x=element_blank(),
             axis.ticks.x=element_blank(),
             axis.text.y=element_blank(),
             axis.ticks.y=element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank()
           )
  return(list(df = pca_df, plot = p))
}

## PCA full dataset

full_pca <- make_pca(final_df)

full_pca$plot
ggsave(paste0(plots_dir, 'Met_Ac/PCAs/full_df_pca.pdf'), full_pca$plot)

full_pca$df %>%
  filter(PC2 < -10) %>%
  select(Gid) %>%
  print(n=200)

## PCA unique genes/states dataset

unique_pca <- make_pca(unique_df)
unique_pca$plot
ggsave(paste0(plots_dir, 'Met_Ac/PCAs/unique_df_pca.pdf'), full_pca$plot)

## PCA compact (5'/ORF/3') DFs
unique_noApi_pca <- unique_noApi %>%
  make_pca()

unique_noApi_pca$plot
ggsave(paste0(plots_dir, 'Met_Ac/PCAs/unique_noApi_df_pca.pdf'), full_pca$plot)

## Make Scatterplots

colnames(unique_df)

y_max <- unique_noApi %>%
  select(matches('Ac|Me')) %>%
  max(na.rm = T)

y_min <- unique_noApi %>%
  select(matches('Ac|Me')) %>%
  min(na.rm = T)


my_scatter <- function(df, ac_col, me_col, outname){
  dfname <- deparse(substitute(df))
  p <- ggplot(df, aes_string(x=ac_col, y=me_col, col = 'Label', alpha = 'Alpha'))
  p <- p + scale_color_manual(values=my_colors)
  #p <- p + theme_minimal() + guides(alpha='none')
  p <- p + theme_bw() + guides(alpha='none')
  p <- p + theme(
             #axis.title.x=element_blank(),
             #axis.text.x=element_blank(),
             #axis.ticks.x=element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank()
           )
  p <- p + geom_point()
  ## p <- p + ylim(y_min, y_max)
  ## p <- p + xlim(y_min, y_max)
  p <- p + coord_cartesian(ylim = c(y_min, y_max))
  ggsave(
    paste0(
      plots_dir,
      sprintf(
        'Met_Ac/Scatterplots/%s_%s_%s_scatter_new.pdf',
        dfname, ac_col, me_col
      )))
  return(p)
}

sc_5p <- my_scatter(unique_noApi, 'Ac_5p', 'Me_5p')
sc_orf <- my_scatter(unique_noApi, 'Ac_ORF', 'Me_ORF')
sc_3p <- my_scatter(unique_noApi, 'Ac_3p', 'Me_3p')
#my_scatter(unique_noApi, 'Ac_5p', 'Me_ORF')
sc_5p
ggarrange(sc_5p, sc_orf, sc_3p,
          labels = c("5p", "ORF", "3p"),
          common.legend = TRUE,
          ncol = 3, nrow = 1) %>%
  ggexport(filename = paste0(plots_dir, 'Met_Ac/Scatterplots/join_scatter_new.pdf'))

## Plot boxplots for each interval

y_max <- unique_noApi %>%
  select(matches('Ac|Me')) %>%
  max(na.rm = T)

y_min <- unique_noApi %>%
  select(matches('Ac|Me')) %>%
  min(na.rm = T)

my_comparisons <- list(
  c("Active", "Inactive"),
  c("Active", "CVG_Active"),
  c("Active", "CVG_Silenced"),
  c("Inactive", "CVG_Active"),
  c("Inactive", "CVG_Silenced"),
  c("CVG_Active", "CVG_Silenced")
)

my_boxplot <- function(df, column, yaxis){
  dfname <- deparse(substitute(df))
  p <- ggplot(unique_noApi, aes_string(x='Label', y=column, fill='Label'))
  p <- p + geom_boxplot(aes(ymax = quantile(eval(parse(text = column)), 0.95, na.rm = T),
                            ymin = quantile(eval(parse(text = column)), 0.05, na.rm = T)),
                        ##outlier.shape = NA
                        )
  p <- p + coord_cartesian(ylim = yaxis)
  p <- p + scale_fill_manual(values=my_colors)
  p <- p + theme_bw() + guides(alpha='none')
  p <- p + theme(
             axis.title.x=element_blank(),
             axis.text.x=element_blank(),
             axis.ticks.x=element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank()
           )
  ## p <- p + stat_compare_means(
  ##            label = 'p.format',
  ##            method = 't.test',
  ##            comparisons = my_comparisons
  ##          )

  ggsave(
    paste0(
      plots_dir,
      sprintf('Met_Ac/Boxplots/%s_%s_boxplot.pdf', dfname, column))
    )
  return(p)
}

bx_ac_5p <- my_boxplot(unique_noApi, 'Ac_5p', c(-4, 4))
bx_ac_ORF <- my_boxplot(unique_noApi, 'Ac_ORF', c(-4, 4))
bx_ac_3p <- my_boxplot(unique_noApi, 'Ac_3p', c(-4, 4))

fm = aov(lm(unique_noApi$Ac_ORF ~ unique_noApi$Label))
summary(fm)
intervals = TukeyHSD(fm)
intervals
plot(intervals)

bx_me_5p <- my_boxplot(unique_noApi, 'Me_5p', c(-5,5))
bx_me_ORF <- my_boxplot(unique_noApi, 'Me_ORF', c(-5,5))
bx_me_3p <- my_boxplot(unique_noApi, 'Me_3p', c(-5,5))


ggarrange(bx_ac_5p, bx_ac_ORF, bx_ac_3p,
          bx_me_5p, bx_me_ORF, bx_me_3p,
          labels = c("Ac_5p", "Ac_ORF", "Ac_3p", "Me_5p", "Me_ORF", "Me_3p"),
          common.legend = TRUE,
          ncol = 3, nrow = 2) %>%
  ggexport(filename = paste0(plots_dir, 'Met_Ac/Boxplots/join_boxplot.pdf'))

df <- unique_df %>%
  left_join(info_df, by='Gene_id') %>%
  select(-contains('Difpeaks'))

df %>%
  select(Gene_id, Ac_ORF, Me_ORF, Variant, Annot) %>%
  filter(Me_ORF > 0) %>%
  filter(!Variant) %>%
  arrange(-Me_ORF) %>%
  print(n = 100)

#### Telomere Analysis ####

## Coverage Analysis

read_telomeres <- function(strain){

  #strain <- '1.2B'
  suffix <- gsub('.', '', strain, fixed=T)
  suffix <- gsub('K9', '', suffix, fixed=T)

  df_ <- read_tsv(
    paste0(telomeres_dir, 'cov_telomeres_', strain, '.bed'),
    col_names = F
  ) %>%
    set_names(c('Chrom', 'Start', 'Stop', 'Cov'))

  tl_l <- df_ %>%
    filter(Start == 1) %>%
    select(Chrom, Cov) %>%
    rename(Cov_Left = Cov)

  tl_r <- df_ %>%
    filter(Start != 1) %>%
    select(Chrom, Cov) %>%
    rename(Cov_Right = Cov)

  tel_df <- tl_l %>%
    full_join(tl_r, by = 'Chrom') %>%
    set_names('Chrom', paste0('Cov_Left_', suffix), paste0('Cov_Right_', suffix))

  return(tel_df)
}

tel_heatmap <- function(df, side, meancenter, logit){

  #df <- tel_df

  tel_df <- df %>%
    select(Chrom, contains(side))

  ## Get row means and subtract from each val
  tel_df_meancentered <- tel_df %>%
    pivot_longer(cols=-Chrom) %>%
    group_by(Chrom)%>%
    mutate(MeanCov = mean(value)) %>%
    pivot_wider(names_from=name, values_from=value)  %>%
    ungroup %>%
    mutate(
      across(
        .cols = contains('Cov_'),
        .fns = function(x) x - MeanCov
      )
    ) %>%
    select(-MeanCov)


  ## Make Matrix for heatmap
  if (meancenter) {
    corr_mtx <- tel_df_meancentered %>%
      select(-Chrom) %>%
      as.matrix()
    rownames(corr_mtx) <- tel_df_meancentered$Chrom
  } else {
    corr_mtx <- tel_df %>%
      select(-Chrom) %>%
      as.matrix()
    rownames(corr_mtx) <- tel_df$Chrom
  }

  ## Log if necessary
  if (logit) {
    corr_mtx <- log2(corr_mtx)
  }

  ## Order Matrix (with clustering)
  ## set.seed(123)
  ## dmtx <- dist(t(corr_mtx), method = "euclidean")
  ## dendo <- hclust(dmtx)
  ## order <- hclust(dmtx)$order
  ## ordered_mtx <- corr_mtx[,order]
  ordered_mtx <- corr_mtx

  out_df <- as_tibble(ordered_mtx) %>%
    mutate(Chrom = rownames(ordered_mtx)) %>%
    select(Chrom, everything())

  ## Plot
  ordered_mtx_m <- melt(ordered_mtx)
  plt_df <- tibble(ordered_mtx_m) %>%
    mutate(Var1 = factor(Var1, levels = rev(tel_df$Chrom)))

  p <- ggplot(plt_df, aes(x=Var2, y=Var1, fill=value))
  p <- p + geom_tile()
  p <- p + theme(

             axis.title.x=element_blank(),
             axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
             axis.title.y=element_blank(),
             legend.position="top"
           )
  p <- p + scale_fill_gradient(
             low = "white",
             high = "orange",
             limits = c(min(corr_mtx), max(corr_mtx)),
             name = 'Coverage'
           )

  return(list(plot = p, df = out_df))
}

strains <- c('NF54', 'P63', '3D7imp', '1.2B', '10G', 'B11', 'A7K9', 'E5K9', 'E5HA')
tels <- lapply(strains, read_telomeres)

tel_df <- bind_cols(tels) %>%
  rename(Chrom = Chrom...1) %>%
  select(-contains('...')) %>%
  select(Chrom, contains('Left'), contains('Right'))

## Set deleted regions to NA ##
tel_df <- tel_df %>%
  mutate(Cov_Left_12B = case_when(
           Chrom == 'Pf3D7_02_v3' ~ NA_real_,
           Chrom == 'Pf3D7_03_v3' ~ NA_real_,
           Chrom == 'Pf3D7_08_v3' ~ NA_real_,
           TRUE ~ Cov_Left_12B
         )) %>%
  mutate(Cov_Left_10G = case_when(
           Chrom == 'Pf3D7_02_v3' ~ NA_real_,
           Chrom == 'Pf3D7_03_v3' ~ NA_real_,
           Chrom == 'Pf3D7_08_v3' ~ NA_real_,
           Chrom == 'Pf3D7_11_v3' ~ NA_real_,
           TRUE ~ Cov_Left_10G
         )) %>%
  mutate(Cov_Left_B11 = case_when(
           Chrom == 'Pf3D7_03_v3' ~ NA_real_,
           Chrom == 'Pf3D7_12_v3' ~ NA_real_,
           TRUE ~ Cov_Left_B11
         )) %>%
  mutate(Cov_Right_12B = case_when(
           Chrom == 'Pf3D7_04_v3' ~ NA_real_,
           TRUE ~ Cov_Right_12B
         )) %>%
  mutate(Cov_Right_10G = case_when(
           Chrom == 'Pf3D7_04_v3' ~ NA_real_,
           TRUE ~ Cov_Right_10G
         )) %>%
  mutate(Cov_Right_B11 = case_when(
           Chrom == 'Pf3D7_08_v3' ~ NA_real_,
           Chrom == 'Pf3D7_12_v3' ~ NA_real_,
           TRUE ~ Cov_Right_B11
         ))

tel_df[tel_df$Chrom == 'Pf3D7_05_v3', grepl('Left', colnames(tel_df))] <- NA_real_
tel_df[tel_df$Chrom == 'Pf3D7_13_v3', grepl('Left', colnames(tel_df))] <- NA_real_

tel_l <- tel_heatmap(
  df <- tel_df,
  side <- 'Left',
  meancentered <- F,
  logit <- F
)

tel_l$plot

tel_r <- tel_heatmap(
  df <- tel_df,
  side <- 'Right',
  meancentered <- F,
  logit <- F
)

tel_l$plot
plot(tel_l$dendo)
tel_l$df

tel_r_plot <- tel_r$plot +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

whole_p <- ggarrange(tel_l$plot, tel_r_plot, nrow = 1, ncol = 2)
ggsave(paste0(plots_dir, 'telomeres_heatmap_new.pdf'), whole_p, device = 'pdf')

## Deletions/Duplications Analysis

tel_dels <- read_csv(paste0(telomeres_dir, '../manual_deletions_table.csv'))

df <- tel_dels %>%
  select(Chrom, contains('left'))

tel_dupl_del_heatmap <- function(side){

  df <- tel_dels %>%
    select(Chrom, contains(side))

  ## Make Matrix for heatmap
  mtx_ <- df %>%
    select(-Chrom) %>%
    as.matrix()

  rownames(mtx_) <- df$Chrom

  ## Order Matrix (by clustering)
  ## set.seed(123)
  ## dmtx <- dist(t(mtx_), method = "euclidean")
  ## dendo <- hclust(dmtx)
  ##plot(dendo)
  ##order <- hclust(dmtx)$order
  ##ordered_mtx <- mtx_[,order]
  ordered_mtx <- mtx_


  ## Plot
  ordered_mtx_m <- melt(ordered_mtx)
  plt_df <- tibble(ordered_mtx_m) %>%
    mutate(Var1 = factor(Var1, levels = rev(df$Chrom)))

  p <- ggplot(plt_df, aes(x=Var2, y=Var1, fill=value))
  p <- p + geom_tile(colour="snow3", size=0.25)
  p <- p + theme(
             panel.background=element_blank(),
             ##panel.grid.minor=element_blank(),
             plot.background=element_blank(),
             axis.title.x=element_blank(),
             axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
             axis.title.y=element_blank(),
             legend.position="top"
           )
  p <- p + scale_fill_gradient2(
             mid = 0,
             name = 'Coverage'
           )
  return(p)
}


pl_l <- tel_dupl_del_heatmap('left')
pl_r <- tel_dupl_del_heatmap('right') +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

whole_p <- ggarrange(pl_l, pl_r, nrow = 1, ncol = 2)
whole_p
ggsave(paste0(plots_dir, 'telomeres_dupl_del_heatmap_new.pdf'), whole_p, device = 'pdf')

#### Save/load environtment ####
##save.image('paper_analysis_070722.RData')
##setwd('/mnt/Disc4T/Projects/PhD_Project/Paper/Paper_Analysis/')
##load('paper_analysis_070722.RData')
