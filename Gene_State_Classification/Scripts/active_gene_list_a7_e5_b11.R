#### Imports ####

library(readxl)
library(tidyverse)
library(eulerr)

#### Max Dif function ####

max_dif <- function(vect){
  mx <- max(vect, na.rm = T)
  mn <- min(vect, na.rm = T)
  if (is.infinite(mx) | is.infinite(mn)) {
    md <- NA
  } else {
    md <- mx - mn
  }
  return(md)
}

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

ave2higher <- function(vect){
  return(mean(sort(vect, decreasing = TRUE)[1:2]))
}

#### Folders ####

wd <- '/mnt/Disc4T/Projects/Active_gene_list/'
setwd(wd)
outdir = './New_Arrays_Results/'
dir.create(outdir)

#### Red Signal DF ####

new_fld <- '/mnt/Disc4T/Projects/PhD_Project/Microarrays/New_Old_separate_approach/New_Arrays/R_results_NewArray/'

## Add B11

B11_red <- read_csv(paste0(new_fld, 'geneLevel_redSignal_exp.csv')) %>%
  select(Gene_id, contains('B11')) %>%
  mutate(across(contains('B11'), function(x) 2**x)) %>%
  mutate(Red_B11 = apply(select(., contains('tp')), 1, ave2higher))

## Add A7

A7_red <- read_csv(paste0(new_fld, 'geneLevel_redSignal_exp.csv')) %>%
  select(Gene_id, contains('A7')) %>%
  mutate(across(contains('A7'), function(x) 2**x)) %>%
  mutate(Red_A7 = apply(select(., contains('tp')), 1, ave2higher))

## Add E5

E5_red <- read_csv(paste0(new_fld, 'geneLevel_redSignal_exp.csv')) %>%
  select(Gene_id, contains('E5')) %>%
  mutate(across(contains('E5'), function(x) 2**x)) %>%
  mutate(Red_E5 = apply(select(., contains('tp')), 1, ave2higher))

## Join

red_df <- A7_red %>%
  full_join(E5_red, by='Gene_id') %>%
  full_join(B11_red, by='Gene_id')

## Transform into percentiles

red_df <- red_df %>%
  mutate(Percent_A7 = (rank(Red_A7)/length(Red_A7))*100) %>%
  mutate(Percent_E5 = (rank(Red_E5)/length(Red_E5))*100) %>%
  mutate(Percent_B11 = (rank(Red_B11)/length(Red_B11))*100)

print(red_df) %>% print(width = 400)

## Add max percentile dif

red_df <- red_df %>%
  mutate(MaxRedPercentDif = apply(select(., contains('Percent_')), 1, max_dif)) %>%
  mutate(MeanRedPercent = apply(select(., contains('Percent_')), 1, mean)) %>%
  mutate(MaxRedPercent = apply(select(., contains('Percent_')), 1, myMax))

print(red_df, width = 200)
hist(red_df$MeanRedPercent)

write_tsv(red_df, paste0(outdir, 'red_df.tsv'))

#### Areas DF ####

B11_area <- read_csv(paste0(new_fld, 'area_geneLevel.csv')) %>%
  select(Gene_id, contains('B11')) %>%
  set_names(c('Gene_id', 'l_B11', 'r_B11', 'm_B11', 's_B11'))

A7_area <- read_csv(paste0(new_fld, 'area_geneLevel.csv')) %>%
  select(Gene_id, contains('A7')) %>%
  set_names(c('Gene_id', 'l_A7', 'r_A7', 'm_A7', 's_A7'))

E5_area <- read_csv(paste0(new_fld, 'area_geneLevel.csv')) %>%
  select(Gene_id, contains('E5')) %>%
  set_names(c('Gene_id', 'l_E5', 'r_E5', 'm_E5', 's_E5'))

area_df <- A7_area %>%
  left_join(E5_area, by='Gene_id') %>%
  left_join(B11_area, by='Gene_id')

## Get Max and Min
print(area_df, width = 200)

area_df <- area_df %>%
  mutate(MaxLeft = apply(select(., contains('l_')), 1, myMax)) %>%
  mutate(MinLeft = apply(select(., contains('l_')), 1, myMin)) %>%

  mutate(MaxRight = apply(select(., contains('r_')), 1, myMax)) %>%
  mutate(MinRight = apply(select(., contains('r_')), 1, myMin)) %>%

  mutate(MaxMid = apply(select(., contains('m_')), 1, myMax)) %>%
  mutate(MinMid = apply(select(., contains('m_')), 1, myMin)) %>%

  mutate(MaxSides = apply(select(., contains('s_')), 1, myMax)) %>%
  mutate(MinSides = apply(select(., contains('s_')), 1, myMin)) %>%

  mutate(DifLeft = MaxLeft - MinLeft) %>%
  mutate(DifRight = MaxRight - MinRight) %>%
  mutate(DifMid = MaxMid - MinMid) %>%
  mutate(DifSides = MaxSides - MinSides)

print(area_df, width = 200)

## Add max interval and difference

timeinterval <- 14.95

maxinterval <- area_df %>%
  select(Gene_id, contains('Dif')) %>%
  pivot_longer(-Gene_id, names_to = 'Interval', values_to = 'MaxDif') %>%
  group_by(Gene_id) %>%
  filter(rank(-MaxDif, ties.method = "first") == 1) %>%
  mutate(Interval = ifelse(is.na(MaxDif), 'No Data', Interval)) %>%
  mutate(Interval = case_when(Interval == 'DifLeft' ~ 'Left',
                              Interval == 'DifRight' ~ 'Right',
                              Interval == 'DifMid' ~ 'Mid',
                              Interval == 'DifSides' ~ 'Sides',
                              Interval == 'No Data' ~ 'No Data')) %>%
  mutate(areaFC = MaxDif/timeinterval)

maxinterval

area_df <- area_df %>%
  left_join(maxinterval, by = 'Gene_id')

print(area_df, width = 400)

## Select appropiate area for each gene and add max and min areas

area_df <- area_df %>%
  mutate(area_A7 = case_when(
           Interval == 'Left' ~ l_A7,
           Interval == 'Right' ~ r_A7,
           Interval == 'Mid' ~ m_A7,
           Interval == 'Sides' ~ s_A7,
           Interval == 'No Data' ~ NA_real_)) %>%
  mutate(area_E5 = case_when(
           Interval == 'Left' ~ l_E5,
           Interval == 'Right' ~ r_E5,
           Interval == 'Mid' ~ m_E5,
           Interval == 'Sides' ~ s_E5,
           Interval == 'No Data' ~ NA_real_)) %>%
  mutate(area_B11 = case_when(
           Interval == 'Left' ~ l_B11,
           Interval == 'Right' ~ r_B11,
           Interval == 'Mid' ~ m_B11,
           Interval == 'Sides' ~ s_B11,
           Interval == 'No Data' ~ NA_real_)) %>%
  mutate(MaxArea = apply(select(., contains('area_')), 1, myMax)) %>%
  mutate(MinArea = apply(select(., contains('area_')), 1, myMin))


print(area_df, width = 400)
table(duplicated(area_df$Gene_id))
write_tsv(area_df, paste0(outdir, 'area_df.tsv'))

#### Filter by Max-Time ####

breaks_df <- read_csv('/mnt/Disc4T/Projects/PhD_Project/Microarrays/New_Old_separate_approach/New_Arrays/R_results_NewArray/new_area_breaks.csv')

max_time <- read_csv('/mnt/Disc4T/Projects/PhD_Project/Microarrays/New_Old_separate_approach/New_Arrays/R_results_NewArray/new_arrays_maxtime.csv')

## Check which areas does maxtimepoint overlapp -> check if aAFC > th at this areas

point_overlap <- function(point, interval){
  point >= interval[1] & point <= interval[2]
}

pass_maxtime <- c()
gids <- c()
th_maxtime <- 1
for (gid in area_df$Gene_id){

  ## Create time-regions
  breaks <- breaks_df$Areas_Breaks
  left <- c(breaks[1], breaks[3])
  right <- c(breaks[3], breaks[5])
  mid <- c(breaks[2], breaks[4])
  sides_l <- c(breaks[1], breaks[2])
  sides_r <- c(breaks[4], breaks[5])

  ## Get maxtime
  if (gid %in% max_time$Gene_id){
    maxtime <- max_time %>%
      filter(Gene_id == gid) %>%
      pull()
  } else {
    maxtime <- NA
  }

  if (is.na(maxtime)){
    pass_maxtime <- c(pass_maxtime, NA)
    gids <- c(gids, gid)
  } else {
    ## Ensure maxtime is in the areas intervals
    if (maxtime < breaks[1]) {maxtime <- breaks[1]}
    if (maxtime > breaks[5]) {maxtime <- breaks[5]}

    ## Get overlappped areas
    areas <- list('left' = left, 'right' = right, 'mid' = mid,
                  'sides' = sides_l, 'sides' = sides_r)
    overlaps <- sapply(areas, function(x) point_overlap(maxtime, x))

    ## Get aAFC in overlapping areas by comparison
    maxes <- area_df %>%
      filter(Gene_id == gid) %>%
      select(contains(names(areas[overlaps])) & contains('Max')) %>%
      replace(is.na(.), 0) %>%
      as.numeric()

    mins <- area_df %>%
      filter(Gene_id == gid) %>%
      select(contains(names(areas[overlaps])) & contains('Min')) %>%
      replace(is.na(.), 0) %>%
      as.numeric()

    aFCs <- (maxes - mins)/timeinterval

    pass_maxtime <- c(pass_maxtime, any(abs(aFCs) > th_maxtime))
    gids <- c(gids, gid)
  }
}

maxtime_df <- tibble(
  Gene_id = gids,
  Pass_Maxtime = pass_maxtime
  )

maxtime_df %>%
  count(Pass_Maxtime)

area_df %>%
  filter(areaFC > 1) %>%
  filter(!Pass_Maxtime) %>%
  select(Gene_id, areaFC, Pass_Maxtime) %>%
  print(n = 999)

area_df %>%
  left_join(maxtime_df, by = 'Gene_id') %>%
  filter(areaFC < 1 & Pass_Maxtime)

area_df <- area_df %>%
  left_join(maxtime_df, by = 'Gene_id')

#### Load RNA-Seq Data ####

## Otto Data-Set

csv <- './Data/RNA_Seq_Percentiles/rnaseq_otto_normvals.csv'

otto <- read_delim(csv, delim=";") %>%
  select(Gene_id = `Gene ID`, contains('unique'))

otto <- otto %>%
  mutate(Max = apply(select(., contains('unique')), 1, myMax)) %>%
  mutate(Otto_Max_pcnt = (rank(Max)/length(Max))*100)

otto <- otto %>% select(Gene_id, Otto_Max_pcnt)
otto <- otto %>% group_by(Gene_id) %>% summarize_all(list(mean))

table(duplicated(otto$Gene_id))
otto %>% filter(duplicated(otto$Gene_id)) %>% print(width = 400)
otto %>% filter(Gene_id == 'PF3D7_0108400') %>% print(width = 400)

#hist(otto$Otto_Max_pcnt)

## Hoeijmakers Data-Set

csv <- './Data/RNA_Seq_Percentiles/rnaseq_hoeijmakers_normvals.csv'

hoeij <- read_delim(csv, delim=";") %>%
  select(Gene_id = `Gene ID`, contains('scaled'))

hoeij <- hoeij %>%
  mutate(Max = apply(select(., contains('scaled')), 1, myMax)) %>%
  mutate(Hoeij_Max_pcnt = (rank(Max)/length(Max))*100)

hoeij <- hoeij %>% select(Gene_id, Hoeij_Max_pcnt)
hoeij <- hoeij %>% group_by(Gene_id) %>% summarize_all(list(mean))
hoeij

#hist(hoeij$Hoeij_Max_pcnt)

## Toenhake Data-Sets

csv <- './Data/RNA_Seq_Percentiles/rnaseq_toen_normvals.csv'

toen <- read_delim(csv, delim=",") %>%
  select(Gene_id = `Gene ID`, contains('unique'))

toen <- toen %>%
  mutate(Max = apply(select(., contains('unique')), 1, myMax)) %>%
  mutate(Toen_Max_pcnt = (rank(Max)/length(Max))*100)

toen <- toen %>% select(Gene_id, Toen_Max_pcnt)
toen <- toen %>% group_by(Gene_id) %>% summarize_all(list(mean))

toen

#hist(toen$Toen_Max_pcnt)

## Bartfai Data-Sets

csv <- './Data/RNA_Seq_Percentiles/rnaseq_bartfai_normvals.csv'

bart <- read_delim(csv, delim=",") %>%
  select(Gene_id = `Gene ID`, contains('scaled'))

bart <- bart %>%
  mutate(Max = apply(select(., contains('scaled')), 1, myMax)) %>%
  mutate(Bart_Max_pcnt = (rank(Max)/length(Max))*100)

bart <- bart %>% select(Gene_id, Bart_Max_pcnt)
bart <- bart %>% group_by(Gene_id) %>% summarize_all(list(mean))

bart

#hist(bart$Bart_Max_pcnt)


## Join DF
rna_df <- otto %>%
  full_join(hoeij) %>%
  full_join(toen) %>%
  full_join(bart)

## Add mean and sd
rna_df <- rna_df %>%
  mutate(MeanPercent = apply(select(., -Gene_id), 1, mean)) %>%
  mutate(StdDevPercent = apply(select(., -Gene_id), 1, sd))


print(rna_df, width=200)

hist(rna_df$MeanPercent, breaks = 20)
hist(rna_df$StdDevPercent, breaks = 100)
table(duplicated(rna_df$Gene_id))

write_tsv(rna_df, paste0(outdir, 'rna_df.tsv'))

#### Load Annotation ####

annot_df <- read_csv('/mnt/Disc4T/Projects/PhD_Project/Binned_Coverage/info_df.csv') %>%
  mutate(Gene_id = ifelse(Gene_id == 'PF3D7_0935390', 'PF3D7_0935400_as', Gene_id))

write_tsv(annot_df, paste0(outdir, 'annot_df.tsv'))

## missing_ids %in% annot_df
## missing_ids

## annot_df %>%
##   rowwise() %>%
##   filter(grepl('gdv1', Annot)) %>%
##   select(Gene_id, Annot)

#### Create Join DF ####

print(red_df, width = 400)
print(area_df, width = 400)
print(rna_df, width = 400)

all_df <- select(red_df, Gene_id, contains('Percent'), MeanRedPercent) %>%
  full_join(select(area_df, Gene_id, Interval, Pass_Maxtime, contains('area')), by = 'Gene_id') %>%
  full_join(select(rna_df, Gene_id, MeanPercent), by = 'Gene_id')

print(all_df, width = 200)

## Here we create a dplyr function.
##To be able to use variables (for colnames) we needto use the special quote functions.
## Colnames to use inside functions must be "enquoted" before usage and preceded by !! when used.
## Colnames to assign must be "enquoted" first, preceded by !! and assigned by :=


## First create a col where we set categories for each gene according relative expression
## For each gene: gene-min----|---mid----|----gene-max

relexprs <- function(vect){
  if (any(is.na(vect))){
    return(NA)
  } else {
    labs = c('min', 'mid', 'max')
    lab <- cut(vect, 3, labels = labs)[1]
    return(as.character(lab))
  }
}

set_relexprs <- function(df, outcol, areacol){
  outcol <- enquo(outcol)
  areacol <- enquo(areacol)
  df %>%
    mutate(!! outcol := apply(select(., !! areacol, MaxArea, MinArea), 1, relexprs))
}

final_df <- all_df %>%
  set_relexprs(rel_A7, area_A7) %>%
  set_relexprs(rel_E5, area_E5) %>%
  set_relexprs(rel_B11, area_B11)

## Add annotation
final_df <- left_join(final_df, annot_df, by = 'Gene_id')
write_tsv(final_df, paste0(outdir, 'final_df.tsv'))

#### Save RData ####

##save.image('load_gene_state_new.RData')
##setwd('/mnt/Disc4T/Projects/Active_gene_list/')
##load('load_gene_state_new.RData')

#### Create Lists according to thresholds ####

print(final_df, width = 200)

th_rnapcnt <- 25
th_redpcnt_dw <- 25
th_redpcnt_up <- 50
#th_redrescue <- 40 # We use redpcnt_up for the same
th_reddown <- 15
th_areaFC <- 1
th_areaFC_redrescue <- 3
th_redpcntdif <- 30
th_maxtime <- 1

## Set areaFC of genes that don't pass MaxTime filter under the FC threshold.
final_df <- final_df %>%
  mutate(areaFC = ifelse(
           areaFC > th_areaFC & !Pass_Maxtime,
           th_areaFC - 0.1,
           areaFC
           ))

## Here we create a dplyr function.
##To be able to use variables (for colnames) we needto use the special quote functions.
## Colnames to use inside functions must be "enquoted" before usage and preceded by !! when used.
## Colnames to assign must be "enquoted" first, preceded by !! and assigned by :=

set_state <- function(df, statecol, redcol, relcol){

  statecol <- enquo(statecol)
  redcol <- enquo(redcol)
  relcol <- enquo(relcol)

  df <- df %>%
    mutate(!! statecol := case_when(

                ## No variant amb FC i RedPcntDif
                ## RNA-Seq & MeanRed
                !Variant &
                areaFC > th_areaFC &
                MeanPercent >= th_rnapcnt &
                MaxRedPercent >= th_reddown &
                MaxRedPercentDif >= th_redpcntdif &
                !! relcol == 'max' ~ 'Active_FC_Dif_max',

                !Variant &
                areaFC > th_areaFC &
                MeanPercent >= th_rnapcnt &
                MaxRedPercent >= th_reddown &
                MaxRedPercentDif >= th_redpcntdif &
                !! relcol == 'mid' ~ 'Active_FC_Dif_mid',

                !Variant &
                areaFC > th_areaFC &
                MeanPercent >= th_rnapcnt &
                MaxRedPercent >= th_reddown &
                MaxRedPercentDif >= th_redpcntdif &
                !! relcol == 'min' ~ 'Inactive_FC_Dif_min',

                ## No variant amb FC no RedPcntDif
                !Variant &
                areaFC > th_areaFC &
                MeanPercent >= th_rnapcnt &
                MaxRedPercent >= th_reddown &
                MaxRedPercentDif < th_redpcntdif ~ 'Active_FC',

                ## No variant amb FC < reddown
                !Variant &
                areaFC > th_areaFC &
                MeanPercent >= th_rnapcnt &
                MaxRedPercent < th_reddown ~ 'Inactive_FC_reddown',

                ## NO-RNA-Seq Rescued
                !Variant &
                areaFC > th_areaFC &
                MeanPercent < th_rnapcnt &
                MaxRedPercent >= th_redpcnt_up &
                MaxRedPercentDif >= th_redpcntdif &
                !! relcol == 'max' ~ 'Active_FC_Dif_rescued_max',

                !Variant &
                areaFC > th_areaFC &
                MeanPercent < th_rnapcnt &
                MaxRedPercent >= th_redpcnt_up &
                MaxRedPercentDif >= th_redpcntdif &
                !! relcol == 'mid' ~ 'Active_FC_Dif_rescued_mid',

                !Variant &
                areaFC > th_areaFC &
                MeanPercent < th_rnapcnt &
                MaxRedPercent >= th_redpcnt_up &
                MaxRedPercentDif >= th_redpcntdif &
                !! relcol == 'min' ~ 'Inactive_FC_Dif_rescued_min',

                ## FC NO RNA-Seq rescued, no redpcnt dif
                !Variant &
                areaFC > th_areaFC &
                MeanPercent < th_rnapcnt &
                MaxRedPercent >= th_redpcnt_up &
                MaxRedPercentDif < th_redpcntdif ~ 'Active_FC_rescued',

                ## FC  No RNA-Seq no rescue
                !Variant &
                areaFC > th_areaFC &
                MeanPercent < th_rnapcnt &
                MaxRedPercent < th_redpcnt_up ~ 'Inactive_FC',

                ## No Var, no FC, 3 categories
                !Variant &
                areaFC < th_areaFC &
                MeanPercent >= th_rnapcnt &
                MaxRedPercent >=th_redpcnt_up ~ 'Active',

                !Variant &
                areaFC < th_areaFC &
                MeanPercent >= th_rnapcnt &
                MaxRedPercent < th_redpcnt_up &
                MaxRedPercent >= th_redpcnt_dw ~ 'Undetermined',

                !Variant &
                areaFC < th_areaFC &
                MeanPercent >= th_rnapcnt &
                MaxRedPercent < th_redpcnt_dw ~ 'Inactive',

                ## !Variant &
                ## areaFC < th_areaFC &
                ## MeanPercent >= th_rnapcnt &
                ## MaxRedPercent < th_reddown ~ 'Inactive_reddown',

                !Variant &
                areaFC < th_areaFC &
                MeanPercent < th_rnapcnt &
                MaxRedPercent > th_redpcnt_up ~ 'Active_rescued',

                !Variant &
                areaFC < th_areaFC &
                MeanPercent < th_rnapcnt &
                MaxRedPercent < th_redpcnt_up ~ 'Inactive_rnaseq',

                ## Clonally Variant Genes
                ## With FC
                Variant &
                areaFC >= th_areaFC &
                MaxRedPercent >= th_reddown &
                !! relcol == 'max' ~ 'CVG_Active_FC',

                Variant &
                areaFC >= th_areaFC &
                MaxRedPercent >= th_reddown &
                !! relcol == 'mid' ~ 'CVG_Semiactive_FC',

                Variant &
                areaFC >= th_areaFC &
                MaxRedPercent >= th_reddown &
                !! relcol == 'min' ~ 'CVG_Silenced_FC',

                Variant &
                areaFC >= th_areaFC &
                MaxRedPercent < th_reddown ~ 'CVG_Silenced_FC_reddown',

                ## No FC
                Variant &
                areaFC < th_areaFC &
                MaxRedPercent >= th_redpcnt_up ~ 'CVG_Active_noFC',

                Variant &
                areaFC < th_areaFC &
                MaxRedPercent < th_redpcnt_up &
                MaxRedPercent >= th_redpcnt_dw  ~ 'CVG_Undetermined_noFC',

                Variant &
                areaFC < th_areaFC &
                MaxRedPercent < th_redpcnt_dw ~ 'CVG_Silenced_noFC',

                ## Exceptions
                ## CVGs with very high FC and very low Red Signal
                Variant &
                areaFC >= th_areaFC_redrescue &
                MaxRedPercent < th_reddown &
                !! relcol == 'max' ~ 'CVG_Undetermined_FC_reddown',

                Variant &
                areaFC >= th_areaFC_redrescue &
                MaxRedPercent < th_reddown &
                !! relcol == 'mid' ~ 'CVG_Silenced_FC_reddown',

                Variant &
                areaFC >= th_areaFC_redrescue &
                MaxRedPercent < th_reddown &
                !! relcol == 'min' ~ 'CVG_Silenced_FC_reddown',


                ###### Not settable
                ## No Var, FC NA
                !Variant &
                is.na(areaFC) &
                MeanPercent >= th_rnapcnt &
                MaxRedPercent >= th_redpcnt_dw ~ 'Active_FCna',

                !Variant &
                is.na(areaFC) &
                MeanPercent >= th_rnapcnt &
                MaxRedPercent < th_redpcnt_dw ~ 'Not_Settable',

                !Variant &
                is.na(areaFC) &
                MeanPercent >= th_rnapcnt &
                is.na(MaxRedPercent) ~ 'Active_FCna_Redna',

                !Variant &
                is.na(areaFC) &
                MeanPercent < th_rnapcnt &
                MaxRedPercent >= th_redpcnt_up ~ 'Not_Settable',

                !Variant &
                is.na(areaFC) &
                MeanPercent < th_rnapcnt &
                MaxRedPercent < th_redpcnt_up ~ 'Inactive_FCna',

                !Variant &
                is.na(areaFC) &
                MeanPercent < th_rnapcnt &
                is.na(MaxRedPercent) ~ 'Inactive_FCna_Redna',

                ## RNA-Seq NA

                !Variant &
                is.na(MeanPercent) &
                areaFC >= th_areaFC &
                MaxRedPercent >= th_redpcnt_up &
                !! relcol == 'max' ~ 'Active_FC_RNASeqNA_max',

                !Variant &
                is.na(MeanPercent) &
                areaFC >= th_areaFC &
                MaxRedPercent >= th_redpcnt_up &
                !! relcol == 'mid' ~ 'Active_FC_RNASeqNA_max',

                !Variant &
                is.na(MeanPercent) &
                areaFC >= th_areaFC &
                MaxRedPercent >= th_redpcnt_up &
                !! relcol == 'min' ~ 'Inactive_FC_RNASeqNA_min',

                !Variant &
                is.na(MeanPercent) &
                areaFC >= th_areaFC &
                MaxRedPercent < th_redpcnt_up &
                MaxRedPercent >= th_redpcnt_dw ~ 'Not_Settable',

                !Variant &
                is.na(MeanPercent) &
                areaFC >= th_areaFC &
                MaxRedPercent < th_redpcnt_dw ~ 'Inactive_FC_reddown',

                !Variant &
                is.na(MeanPercent) &
                areaFC >= th_areaFC &
                is.na(MaxRedPercent) ~ 'Not_Settable',

                !Variant &
                is.na(MeanPercent) &
                areaFC < th_areaFC &
                MaxRedPercent >= th_redpcnt_up ~ 'Active_RNASeqNA',

                !Variant &
                is.na(MeanPercent) &
                areaFC < th_areaFC &
                MaxRedPercent < th_redpcnt_up &
                MaxRedPercent >= th_redpcnt_dw ~ 'Undetermined_RNASeqNA',

                !Variant &
                is.na(MeanPercent) &
                areaFC < th_areaFC &
                MaxRedPercent < th_redpcnt_dw ~ 'Inactive_RNASeqNA_reddown',

                !Variant &
                is.na(MeanPercent) &
                areaFC < th_areaFC &
                is.na(MaxRedPercent) ~ 'Not_Settable',

                ## VARs
                Variant &
                is.na(areaFC) &
                MaxRedPercent >= th_redpcnt_up ~ 'CVG_Active_FCna',

                Variant &
                is.na(areaFC) &
                MaxRedPercent < th_redpcnt_up &
                MaxRedPercent >= th_redpcnt_dw ~ 'CVG_Undetermined_FCna',

                Variant &
                is.na(areaFC) &
                MaxRedPercent < th_redpcnt_dw ~ 'CVG_Silenced_FCna',

                Variant &
                is.na(areaFC) &
                is.na(MaxRedPercent) ~ 'CVG_Not_Settable',

                ## Double NAs
                is.na(areaFC) &
                is.na(MeanPercent) ~ 'CVG_Not_Settable',

                TRUE ~ 'Wrong!'))

  ## The 'TRUE ~ ...' handles rows that do not match any of previous patterns.
  ## Here we use it to make sure all rows are set (no "Wrong!" appearing)

  return(df)
}

set_category <- function(df, statecol, categorycol){

  statecol <- enquo(statecol)
  categorycol <- enquo(categorycol)

  df <- df %>%
    mutate(!! categorycol := case_when(
                startsWith(!! statecol, 'Active') ~ 'Active',
                startsWith(!! statecol, 'Inactive') ~ 'Inactive',
                startsWith(!! statecol, 'Undetermined') ~ 'Undetermined',
                startsWith(!! statecol, 'Not_Settable') ~ 'Undetermined',
                startsWith(!! statecol, 'CVG_Active') ~ 'CVG_Active',
                startsWith(!! statecol, 'CVG_Semiactive') ~ 'CVG_Undetermined',
                startsWith(!! statecol, 'CVG_Silenced') ~ 'CVG_Silenced',
                startsWith(!! statecol, 'CVG_Undetermined') ~ 'CVG_Undetermined',
                startsWith(!! statecol, 'CVG_Not_Settable') ~ 'CVG_Undetermined',
                TRUE ~ 'Undetermined'))
  return(df)
}

## We now set each gene to it's state

state_df <- final_df %>%
  set_state(state_A7, Percent_A7, rel_A7) %>%
  set_state(state_E5, Percent_E5, rel_E5) %>%
  set_state(state_B11, Percent_B11, rel_B11) %>%
  set_category(state_A7, category_A7) %>%
  set_category(state_E5, category_E5) %>%
  set_category(state_B11, category_B11)

print(state_df, width = 400)

table(state_df$state_A7)
table(state_df$category_A7)

missing_ids <- state_df %>%
  filter(
    state_A7 == 'Wrong!' |
    state_E5 == 'Wrong!' |
    state_B11 == 'Wrong!'
  ) %>%
  select(Gene_id) %>%
  pull()

missing_genes <- missing_ids[grepl('PF3D7', missing_ids)]

state_df <- state_df %>%
  filter(
    !Gene_id %in% missing_ids
  )

## Filter out genes with duplications/deletions

f_path <- '/mnt/Disc4T/Projects/PhD_Project/Paper/Paper_Analysis/Data_Files/Duplication_Deletion_Regions_Mean_Separate_DuplDel/Crossed_with_genes/'

file_list <- c(
  'A7K9_in_sort_q5_noDup_rpkm_normInput_bs10_smth200_bymean_dupl_del_fact_1.75up_0.1dw_minlen500_mergelen_200_filtered_genes.tsv',
  'E5K9_in_sort_q5_noDup_rpkm_normInput_bs10_smth200_bymean_dupl_del_fact_1.75up_0.1dw_minlen500_mergelen_200_filtered_genes.tsv',
  'B11_in_sort_q5_noDup_rpkm_normInput_bs10_smth200_bymean_dupl_del_fact_1.75up_0.1dw_minlen500_mergelen_200_filtered_genes.tsv'
)

dupl_dl_A7 <- read_tsv(paste0(f_path, file_list[1]), col_names = F) %>%
  select(X1) %>% pull()

dupl_dl_E5 <- read_tsv(paste0(f_path, file_list[2]), col_names = F) %>%
  select(X1) %>% pull()

dupl_dl_B11 <- read_tsv(paste0(f_path, file_list[3]), col_names = F) %>%
  select(X1) %>% pull()


## state_df <- state_df %>%
##   mutate(state_A7 = ifelse(
##            Gene_id %in% dupl_dl_A7,
##            'Not_Settable',
##            state_A7
##          )) %>%
##   mutate(category_A7 = ifelse(
##            Gene_id %in% dupl_dl_A7,
##            'Not_Settable',
##            category_A7
##          )) %>%
##   mutate(state_E5 = ifelse(
##            Gene_id %in% dupl_dl_E5,
##            'Not_Settable',
##            state_E5
##          )) %>%
##   mutate(category_E5 = ifelse(
##            Gene_id %in% dupl_dl_E5,
##            'Not_Settable',
##            category_E5
##          )) %>%
##   mutate(state_B11 = ifelse(
##            Gene_id %in% dupl_dl_B11,
##            'Not_Settable',
##            state_B11
##          )) %>%
##   mutate(category_B11 = ifelse(
##            Gene_id %in% dupl_dl_B11,
##            'Not_Settable',
##            category_B11
##          ))

state_df <- state_df %>%
  mutate(state_A7 = case_when(
           Gene_id %in% dupl_dl_A7 & !Variant ~ 'Undetermined_dupldel',
           Gene_id %in% dupl_dl_A7 & Variant ~ 'CVG_Undetermined_dupldel',
           TRUE ~ state_A7
         )) %>%
  mutate(category_A7 = case_when(
           Gene_id %in% dupl_dl_A7 & !Variant ~ 'Undetermined',
           Gene_id %in% dupl_dl_A7 & Variant ~ 'CVG_Undetermined',
           TRUE ~ category_A7
         )) %>%
  mutate(state_E5 = case_when(
           Gene_id %in% dupl_dl_E5 & !Variant ~ 'Undetermined_dupldel',
           Gene_id %in% dupl_dl_E5 & Variant ~ 'CVG_Undetermined_dupldel',
           TRUE ~ state_E5
         )) %>%
  mutate(category_E5 = case_when(
           Gene_id %in% dupl_dl_E5 & !Variant ~ 'Undetermined',
           Gene_id %in% dupl_dl_E5 & Variant ~ 'CVG_Undetermined',
           TRUE ~ category_E5
         )) %>%
  mutate(state_B11 = case_when(
           Gene_id %in% dupl_dl_B11 & !Variant ~ 'Undetermined_dupldel',
           Gene_id %in% dupl_dl_B11 & Variant ~ 'CVG_Undetermined_dupldel',
           TRUE ~ state_B11
         )) %>%
  mutate(category_B11 = case_when(
           Gene_id %in% dupl_dl_B11 & !Variant ~ 'Undetermined',
           Gene_id %in% dupl_dl_B11 & Variant ~ 'CVG_Undetermined',
           TRUE ~ category_B11
         ))



## Save results
outname <- sprintf(  'state_df_new_arrays__rna%s_red_up%s_red_dw%s_reddw%s_areaFC%s_areaFCbig%s_redpcntdif%s_maxtime%s_dupldel_filtered_new.csv',
  th_rnapcnt,
  th_redpcnt_up,
  th_redpcnt_dw,
  th_reddown,
  th_areaFC,
  th_areaFC_redrescue,
  th_redpcntdif,
  th_maxtime
)
write_tsv(state_df, paste0(outdir, outname))

#### Some checks ####

## We check no rows are set to "Wrong!"
state_df %>%
  filter(state_A7 == 'Not_Settable'|
         state_E5 == 'Not_Settable'|
         state_B11 == 'Not_Settable') %>%
  print(width = 400)

state_df %>%
  filter(category_A7 == 'No_Category' |
         category_A7 == 'No_Category' |
         category_B11 == 'No_Category') %>%
  print(width = 400)


## ## Create a table with number of each state per strain
## state_table <-  bind_rows(table(state_df$state_12B),
##                           table(state_df$state_10G),
##                           table(state_df$state_3D7B),
##                           table(state_df$state_B11)) %>%
##   replace_na(list(Var_Semiactive = 0)) %>%
##   mutate(Strain = c('12B', '10G', '3D7B', 'B11')) %>%
##   select(Strain, everything())

## ## Create a table with differences between 12B and 10G
## dif12B_10G <- state_df %>%
##   filter(state_12B != state_10G) %>%
##   select(Gene_id, contains('12B'), contains('10G'), Gene_name, Annot)

## ## Check Clags
## clags <- state_df %>%
##   filter(Gene_id == 'PF3D7_0302500' | Gene_id == 'PF3D7_0302200')

## write.csv(clags, './Results_Tables_All_Strains_withB11/clag_genes.csv')

## print(state_table, width = 200)
## summary(rna_df)
