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

#### Red Signal DF ####

## Read translation table
#map <- read.csv('./Data/oldnames_table.csv', stringsAsFactors = F)
map <- read_csv('/mnt/Disc4T/Projects/PhD_Project/Microarrays/New_Old_separate_approach/Old_Arrays/raw_gene_level_data.csv') %>%
  select(Old_id, Gene_id)

excl <- "./Data/3D7_Variantome_AllData_withGam.xls"

## Import Red Signal table
red <- read_excel(excl, sheet = 4)

colnames(red)[1] <- "Old_id"

red_df <- red %>%
  select(Old_id,
         Red_12B = `Aver.2Higher1.2B.`,
         Red_10G = `Aver.2Higher10G.`,
         Red_3D7B = `Aver.2Higher3D7-B.`) %>%
  left_join(map, by='Old_id')

## Manually add gene ids that do not appear on map (we blasted the probes to assign them)
red_df <- red_df %>%
  mutate(Gene_id = ifelse(Old_id == 'MAL8P1.310', 'PF3D7_0830200', Gene_id)) %>%
  mutate(Gene_id = ifelse(Old_id == 'PFI0905W', 'PF3D7_0918500', Gene_id)) %>%
  mutate(Gene_id = ifelse(Old_id == 'PFL1580W', 'PF3D7_1232700', Gene_id))

## Collapse gene ids that appear more than once (mean expression)
red_df <- red_df %>%
  select(-Old_id) %>%
  group_by(Gene_id) %>% summarize_all(list(mean))

## Filter out rows with untranslated gene ids (marked by '_oldname')
red_df <- red_df %>%
  filter(!grepl('_oldname', Gene_id))

## Transform into percentiles

red_df <- red_df %>%
  mutate(Percent_12B = (rank(Red_12B)/length(Red_12B))*100) %>%
  mutate(Percent_10G = (rank(Red_10G)/length(Red_10G))*100) %>%
  mutate(Percent_3D7B = (rank(Red_3D7B)/length(Red_3D7B))*100)

## Add max percentile dif

red_df <- red_df %>%
  mutate(MaxRedPercentDif= apply(select(., contains('Percent_')), 1, max_dif)) %>%
  mutate(MeanRedPercent = apply(select(., contains('Percent_')), 1, mean)) %>%
  mutate(MaxRedPercent = apply(select(., contains('Percent_')), 1, max))

print(red_df, width = 200)
hist(red_df$MeanRedPercent)

#### Areas DF ####

# Import Areas table

area <- read_excel(excl, sheet = 2)

colnames(area)[1] <- "Old_id"

area_df <- area %>%
  select(Old_id,
         l_12B = `left.1.2b`,
         r_12B = `right.1.2b`,
         m_12B = `mid.1.2b`,
         s_12B = `sides.1.2b`,
         l_10G = `left.10g`,
         r_10G = `right.10g`,
         m_10G = `mid.10g`,
         s_10G = `sides.10g`,
         l_3D7B = `left.3d7b`,
         r_3D7B = `right.3d7b`,
         m_3D7B = `mid.3d7b`,
         s_3D7B = `sides.3d7b`) %>%

  mutate_at(vars(-Old_id), as.numeric) %>%


  left_join(map, by='Old_id') %>%
  select(-Old_id) %>%
  group_by(Gene_id) %>% summarize_all(list(mean))


print(area_df, width = 200)

area_df <- area_df %>%
  mutate(MaxLeft = apply(select(., contains('l_')), 1, max)) %>%
  mutate(MinLeft = apply(select(., contains('l_')), 1, min)) %>%

  mutate(MaxRight = apply(select(., contains('r_')), 1, max)) %>%
  mutate(MinRight = apply(select(., contains('r_')), 1, min)) %>%

  mutate(MaxMid = apply(select(., contains('m_')), 1, max)) %>%
  mutate(MinMid = apply(select(., contains('m_')), 1, min)) %>%

  mutate(MaxSides = apply(select(., contains('s_')), 1, max)) %>%
  mutate(MinSides = apply(select(., contains('s_')), 1, min)) %>%

  mutate(DifLeft = MaxLeft - MinLeft) %>%
  mutate(DifRight = MaxRight - MinRight) %>%
  mutate(DifMid = MaxMid - MinMid) %>%
  mutate(DifSides = MaxSides - MinSides)

print(area_df, width = 200)

## Add max interval and difference

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
  mutate(areaFC = MaxDif/13.5)

maxinterval

area_df <- area_df %>%
  left_join(maxinterval, by = 'Gene_id')

print(area_df, width = 400)

## Select appropiate area for each gene and add max and min areas

area_df <- area_df %>%
  mutate(area_12B = case_when(
           Interval == 'Left' ~ l_12B,
           Interval == 'Right' ~ r_12B,
           Interval == 'Mid' ~ m_12B,
           Interval == 'Sides' ~ s_12B,
           Interval == 'No Data' ~ NA_real_)) %>%
  mutate(area_10G = case_when(
           Interval == 'Left' ~ l_10G,
           Interval == 'Right' ~ r_10G,
           Interval == 'Mid' ~ m_10G,
           Interval == 'Sides' ~ s_10G,
           Interval == 'No Data' ~ NA_real_)) %>%
  mutate(area_3D7B = case_when(
           Interval == 'Left' ~ l_3D7B,
           Interval == 'Right' ~ r_3D7B,
           Interval == 'Mid' ~ m_3D7B,
           Interval == 'Sides' ~ s_3D7B,
           Interval == 'No Data' ~ NA_real_)) %>%
  mutate(MaxArea = apply(select(., contains('area_')), 1, max)) %>%
  mutate(MinArea = apply(select(., contains('area_')), 1, min))


print(area_df, width = 400)
table(duplicated(area_df$Gene_id))

#### Filter by Max-Time ####

breaks_df <- read_csv('/mnt/Disc4T/Projects/PhD_Project/Microarrays/New_Old_separate_approach/Old_Arrays/R_results_OldArrays_Variantome/old_area_breaks.csv')

max_time <- read_csv('/mnt/Disc4T/Projects/PhD_Project/Microarrays/New_Old_separate_approach/Old_Arrays/R_results_OldArrays_Variantome/old_arrays_maxtime.csv')

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

    aFCs <- (maxes - mins)/13.5

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
  left_join(maxtime_df, by = 'Gene_id') %>%
  filter(areaFC > 1) %>%
  filter(!Pass_Maxtime) %>%
  select(Gene_id, areaFC, Pass_Maxtime)

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
  mutate(Max = apply(select(., contains('unique')), 1, max)) %>%
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
  mutate(Max = apply(select(., contains('scaled')), 1, max)) %>%
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
  mutate(Max = apply(select(., contains('unique')), 1, max)) %>%
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
  mutate(Max = apply(select(., contains('scaled')), 1, max)) %>%
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
write_csv(rna_df, './Results_Tables/rna_seq_percentiles.csv')


hist(rna_df$MeanPercent, breaks = 20)
hist(rna_df$StdDevPercent, breaks = 100)
table(duplicated(rna_df$Gene_id))

#### Load Annotation ####

annot_df <- read_delim('./Data/plasmoDB_geneAnnot.csv', delim = ';') %>%
  select(Gene_id = `Gene ID`,
         Gene_name = `Gene Name or Symbol`,
         Annot = `Product Description`) %>%
  distinct() # remove duplicated rows

print(annot_df, width=200)
table(duplicated(annot_df$Gene_id))

#### Create Join DF ####

print(red_df, width = 200)
print(area_df, width = 200)
print(rna_df, width = 200)

all_df <- select(red_df, Gene_id, contains('Percent'), MeanRedPercent) %>%
  full_join(select(area_df, Gene_id, Interval, Pass_Maxtime, contains('area')), by = 'Gene_id') %>%
  full_join(select(rna_df, Gene_id, MeanPercent), by = 'Gene_id')

## Add Vartiant Genes information

cvg <- read_excel("./Data/CVG_list_jan2020_final.xlsx", sheet = "Final")

final_df <- cvg %>%
  select("Gene_id" = `Gene ID`, "Variant" = `Final Customized`) %>%
  right_join(all_df, by = 'Gene_id') %>%
  mutate(Variant = recode(Variant, YES = TRUE, NO = FALSE, .missing = FALSE))

print(final_df, width = 200)

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

final_df <- final_df %>%
  set_relexprs(rel_12B, area_12B) %>%
  set_relexprs(rel_10G, area_10G) %>%
  set_relexprs(rel_3D7B, area_3D7B)

## Add annotation
final_df <- left_join(final_df, annot_df, by = 'Gene_id')

print(final_df, width = 200)

#### Save R Data ####

##save.image('load_gene_state_old.RData')
##setwd('/mnt/Disc4T/Projects/Active_gene_list/')
##load('load_gene_state_old.RData')

#### Create Lists according to thresholds ####
print(final_df, width = 200)

colnames(final_df)

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
  set_state(state_12B, Percent_12B, rel_12B) %>%
  set_state(state_10G, Percent_10G, rel_10G) %>%
  set_state(state_3D7B, Percent_3D7B, rel_3D7B) %>%
  set_category(state_12B, category_12B) %>%
  set_category(state_10G, category_10G) %>%
  set_category(state_3D7B, category_3D7B)



state_df %>%
  filter(Variant & category_3D7B == 'Undetermined') %>%
  select(Variant, contains('state'), contains('category'))

state_df %>%
  filter(grepl('Undetermined', state_12B)) %>%
  print(width = 400)

state_df %>%
  filter(category_12B == 'No_Category') %>%
  print(width = 400)


## Filter out deleted/duplicated genes

f_path <- '/mnt/Disc4T/Projects/PhD_Project/Paper/Paper_Analysis/Data_Files/Duplication_Deletion_Regions_Mean_Separate_DuplDel/Crossed_with_genes/'

file_list <- c( '1.2B_in_sort_q5_noDup_rpkm_normInput_bs10_smth200_bymean_dupl_del_fact_1.75up_0.1dw_minlen500_mergelen_200_filtered_genes.tsv',
               '10G_in_sort_q5_noDup_rpkm_normInput_bs10_smth200_bymean_dupl_del_fact_1.75up_0.1dw_minlen500_mergelen_200_filtered_genes.tsv'
)

dupl_dl_12B <- read_tsv(paste0(f_path, file_list[1]), col_names = F) %>%
  select(X1) %>% pull()

dupl_dl_10G <- read_tsv(paste0(f_path, file_list[2]), col_names = F) %>%
  select(X1) %>% pull()

## state_df <- state_df %>%
##   mutate(state_12B = ifelse(
##            Gene_id %in% dupl_dl_12B,
##            'Not_Settable',
##            state_12B
##          )) %>%
##   mutate(category_12B = ifelse(
##            Gene_id %in% dupl_dl_12B,
##            'Not_Settable',
##            category_12B
##          )) %>%
##   mutate(state_10G = ifelse(
##            Gene_id %in% dupl_dl_10G,
##            'Not_Settable',
##            state_10G
##          )) %>%
##   mutate(category_10G = ifelse(
##            Gene_id %in% dupl_dl_10G,
##            'Not_Settable',
##            category_10G
##          ))

state_df <- state_df %>%
  mutate(state_12B = case_when(
           Gene_id %in% dupl_dl_12B & !Variant ~ 'Undetermined_dupldel',
           Gene_id %in% dupl_dl_12B & Variant ~ 'CVG_Undetermined_dupldel',
           TRUE ~ state_12B
         )) %>%
  mutate(category_12B = case_when(
           Gene_id %in% dupl_dl_12B & !Variant ~ 'Undetermined',
           Gene_id %in% dupl_dl_12B & Variant ~ 'CVG_Undetermined',
           TRUE ~ category_12B
         )) %>%
  mutate(state_10G = case_when(
           Gene_id %in% dupl_dl_10G & !Variant ~ 'Undetermined_dupldel',
           Gene_id %in% dupl_dl_10G & Variant ~ 'CVG_Undetermined_dupldel',
           TRUE ~ state_10G
         )) %>%
  mutate(category_10G = case_when(
           Gene_id %in% dupl_dl_10G & !Variant ~ 'Undetermined',
           Gene_id %in% dupl_dl_10G & Variant ~ 'CVG_Undetermined',
           TRUE ~ category_10G
         ))


## Save results
outname <- sprintf(  'state_df_old_arrays_rna%s_red_up%s_red_dw%s_reddw%s_areaFC%s_areaFCbig%s_redpcntdif%s_maxtime%s_dupldel_filtered_new.csv',
  th_rnapcnt,
  th_redpcnt_up,
  th_redpcnt_dw,
  th_reddown,
  th_areaFC,
  th_areaFC_redrescue,
  th_redpcntdif,
  th_maxtime
)
write_tsv(state_df, paste0('./Results_Tables/', outname))

#### Some checks and miscelania ####
## We check no rows are set to "Wrong!"
state_df %>%
  filter(state_12B == 'Not_Settable' |
         state_10G == 'Not_Settable' |
         state_3D7B == 'Not_Settable') %>%
  print(width = 400)

state_df %>%
  filter(category_12B == 'No_Category' |
         category_10G == 'No_Category' |
         category_3D7B == 'No_Category') %>%
  print(width = 400)

## Create a table with number of each state per strain
state_table <-  bind_rows(table(state_df$state_12B),
                          table(state_df$state_10G),
                          table(state_df$state_3D7B)) %>%
  replace_na(list(Var_Semiactive = 0)) %>%
  mutate(Strain = c('12B', '10G', '3D7B')) %>%
  select(Strain, everything())

## Create a table with differences between 12B and 10G
dif12B_10G <- state_df %>%
  filter(state_12B != state_10G) %>%
  select(Gene_id, contains('12B'), contains('10G'), Gene_name, Annot)

## Check Clags
clags <- state_df %>%
  filter(Gene_id == 'PF3D7_0302500' | Gene_id == 'PF3D7_0302200')

write.csv(clags, './Results_Tables/clag_genes.csv')

print(state_table, width = 200)
summary(rna_df)

## Categories histograms

make_histogram <- function(df, column){
  col <- enquo(column)

  df %>%
    select(!!col) %>%
    count(!!col) %>%
    ggplot(aes(y=n, x=!!col, fill=!!col)) +
    geom_bar(stat='identity') +
    geom_text(aes(label=n), vjust=0) +
    theme(
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
}


hist12B <- make_histogram(state_df, category_12B)
hist10G <- make_histogram(state_df, category_10G)
hist3D7B <- make_histogram(state_df, category_3D7B)

ggsave(hist12B, filename = './Plots/histogram_12B.png', height=7, width=8)
ggsave(hist10G, filename = './Plots/histogram_10G.png', height=7, width=8)
ggsave(hist3D7B, filename = './Plots/histogram_3D7B.png', height=7, width=8)

excl <- "./Data/10Gvs1p2B.xlsx"

## Import 12B vs 10G differences by transcription table
trans_difs <- read_delim('./Data/transDif_12Bvs10G.csv', delim = ';') %>%
  rename(Old_id = `...1`) %>%
  left_join(map, by = 'Old_id')

print(trans_difs, width = 400)

trans_difs %>%
  select(Old_id, Gene_id)

dif12B_10G

table(trans_difs$Gene_id %in% dif12B_10G$Gene_id)
table(dif12B_10G$Gene_id %in% trans_difs$Gene_id)

allids <- unique(c(dif12B_10G$Gene_id, trans_difs$Gene_id))

compare_12Bvs10G <- state_df %>%
  filter(Gene_id %in% allids) %>%
  select(Gene_id,
         Variant, areaFC,
         MeanRedPercent,
         contains('12B'),
         contains('10G'),
         Gene_name,
         Annot) %>%
  mutate(TransDif = Gene_id %in% trans_difs$Gene_id) %>%
  mutate(Dif_state = category_12B != category_10G)

print(compare_12Bvs10G, width = 400)
write.csv(compare_12Bvs10G, './Results_Tables/gens_dif12B_10G.csv')

makeIntersects <- function(a,b,c){

  a_b <- intersect(a, b)
  a_c <- intersect(a, c)
  b_c <- intersect(b, c)
  a_b_c <- intersect(a_b, c)

  abc <- a_b_c
  ab <- a_b[!a_b %in% a_b_c]
  ac <- a_c[!a_c %in% a_b_c]
  bc <- b_c[!b_c %in% a_b_c]

  a <- a[!a %in% ab & !a %in% ac & !a %in% abc]
  b <- b[!b %in% ab & !b %in% bc & !b %in% abc]
  c <- c[!c %in% ac & !c %in% bc & !c %in% abc]

  return(list(a = a, b = b, c = c, ab = ab, bc = bc, ac = ac, abc = abc))

}

customEuler <- function(a,b,c, name){

  intersects <- makeIntersects(a,b,c)
  areas <- lapply(intersects, function(x) length(x))

  fit <- euler(c(A=areas$a, B=areas$b, C=areas$c,
                 "A&B"=areas$ab, "A&C"=areas$ac, "B&C"=areas$bc,
                 "A&B&C" = areas$abc))

  d <- plot(fit, fills = list(fill = c("#619CFF", "#F8766D", "#00BA38"), alpha = 0.5),
            edges = list(lwd = 0.1), quantities = list(quantities = T),
            labels = list(labels=c("12B", "10G", "3D7B")),
            main = name)

  ggsave(d, filename = paste0('./Plots/', "venn_", name, ".png"),
         device = "png", width = 10, height = 10, units = "cm")

  plot(d)
  print(fit)
}


## Var active

va_12B <- state_df %>%
  filter(category_12B == "CVG_Active") %>%
  select(Gene_id) %>%
  pull()

va_10G <- state_df %>%
  filter(category_10G == "CVG_Active") %>%
  select(Gene_id) %>%
  pull()

va_3D7B <- state_df %>%
  filter(category_3D7B == "CVG_Active") %>%
  select(Gene_id) %>%
  pull()


customEuler(va_12B, va_10G, va_3D7B, 'CVG_Active_Genes')

## Var inactive

vi_12B <- state_df %>%
  filter(category_12B == "CVG_Silenced") %>%
  select(Gene_id) %>%
  pull()

vi_10G <- state_df %>%
  filter(category_10G == "CVG_Silenced") %>%
  select(Gene_id) %>%
  pull()

vi_3D7B <- state_df %>%
  filter(category_3D7B == "CVG_Silenced") %>%
  select(Gene_id) %>%
  pull()


customEuler(vi_12B, vi_10G, vi_3D7B, 'CVG_Silenced_Genes')

## Same category
## We fuse the gene_id with it's state so that genes with same state will be exactly the same while genes with different state will be different.

print(state_df, width = 400)

genestate <- state_df %>%
  mutate(gs_12B = paste(Gene_id, category_12B, sep = "_")) %>%
  mutate(gs_10G = paste(Gene_id, category_10G, sep = "_")) %>%
  mutate(gs_3D7B = paste(Gene_id, category_3D7B, sep = "_")) %>%
  select(Gene_id,
         category_12B, category_10G, category_3D7B,
         gs_12B, gs_10G, gs_3D7B)

customEuler(genestate$gs_12B, genestate$gs_10G, genestate$gs_3D7B, 'Same_State')

x <- makeIntersects(genestate$gs_12B, genestate$gs_10G, genestate$gs_3D7B)

varstate <- state_df %>%
  filter(Variant) %>%
  mutate(gs_12B = paste(Gene_id, category_12B, sep = "_")) %>%
  mutate(gs_10G = paste(Gene_id, category_10G, sep = "_")) %>%
  mutate(gs_3D7B = paste(Gene_id, category_3D7B, sep = "_")) %>%
  select(Gene_id,
         category_12B, category_10G, category_3D7B,
         gs_12B, gs_10G, gs_3D7B)

customEuler(varstate$gs_12B, varstate$gs_10G, varstate$gs_3D7B, 'CVG_Same_State')
