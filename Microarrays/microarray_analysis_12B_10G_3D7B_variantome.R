#### Import libraries ####

print("Importing Libraries...")
list.of.packages <- c("reshape2", "ggfortify", "tidyverse", "RColorBrewer", "sp")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(reshape2)
library(ggfortify)
library(tidyverse)
library(RColorBrewer)
library(sp)
library(readxl)

if (!("Biobase" %in% installed.packages()[,"Package"])){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install()
}

library(Biobase)

#### Experiment Setup: MODIFY THIS PART #################################
##*********************************************************************##
##*********************************************************************##

wd <- ('/mnt/Disc4T/Projects/PhD_Project/Microarrays/New_Old_separate_approach/Old_Arrays/')
setwd(wd)
datadir <- "./Variantome_Original/"
outdir <- "./R_results_OldArrays_Variantome/"
annot <- read.csv("./Files/array_anotation.csv", sep="\t", header = F)
infiles <- "./Files/"

sample_names <- c(
  '12B_tp10', '12B_tp20', '12B_tp30', '12B_tp34', '12B_tp37', '12B_tp40', '12B_tp43',
  '10G_tp10', '10G_tp20', '10G_tp30', '10G_tp34', '10G_tp37', '10G_tp40', '10G_tp43',
  '3D7B_tp10', '3D7B_tp20', '3D7B_tp30', '3D7B_tp34', '3D7B_tp37', '3D7B_tp40', '3D7B_tp43'
)

nsamples <- length(sample_names)

times <- as.integer(sub("^.+_tp", "", sample_names))
types <- sub("_tp.+", "", sample_names)

## Load Annotation
print("Loading Array Annotation...")
gene_list <- readLines(paste0(infiles, "gene_list.txt"))
annot_df <- read_csv(paste0(infiles, 'info_df.csv')) %>%
  select(Gene_id, Name, Variant, Annot) %>%
  dplyr::filter(!is.na(Gene_id))

## Add Annotation
## gene_rosetta <- read_tsv(paste0(infiles, 'gene_ids_rosetta.txt'), col_names = F)
## colnames(gene_rosetta) <- c('Old_id', 'Gene_id')

gene_rosetta <- read_tsv('/mnt/Disc4T/Projects/Miniprojects/Gene_Rosetta/gene_rosetta.tsv')

## Set to NA rows with double ID (to be changed)

blasted <- read_csv('./Blast_Old_Primers/blasted_oligos.csv')
blasted <- blasted %>%
  mutate(Old_id = gsub('a(.*)_[0-9]', '\\1', blasted$Oligo_id))

unique_hit <- blasted %>%
  filter(Score >= 70 & Score_2 <= 70 ) %>%
  select(Old_id, Gene_id) %>%
  group_by(Old_id) %>%
  summarize(N_hits = n_distinct(Gene_id)) %>%
  filter(N_hits == 1) %>%
  select(Old_id) %>% pull()

unique_oligos <- blasted %>%
  filter(Old_id %in% unique_hit) %>%
  select(Old_id, Gene_id) %>%
  distinct()

gene_rosetta_blasted <- gene_rosetta %>%
  full_join(unique_oligos, by = 'Old_id', suffix = c('_aliases', '_blast')) %>%
  mutate(Gene_id = ifelse(is.na(Gene_id_blast), Gene_id_aliases, Gene_id_blast)) %>%
  select(Old_id, Gene_id, everything())

double_new_ids <- gene_rosetta_blasted %>%
   filter(grepl(',', Gene_id))

write_csv(double_new_ids, './genes_with_double_GeneID.csv')

gene_rosetta_blasted <- gene_rosetta_blasted %>%
  mutate(Gene_id = replace(Gene_id, grepl(',', Gene_id), NA))


write_csv(gene_rosetta_blasted, './final_rosseta_with_blasted_genes.csv')

## Load Gene-Level data

gene_level <- read_csv('./Variantome_Original/normalizedData_geneLevel.csv')
gene_level_noratio <- read_csv('./Variantome_Original/normalizedData_geneLevel_noRatio.csv')

gene_level %>%
  select(geneID) %>%
  write_csv('old_gene_ids.txt', , col_names = F)

gene_level <- gene_level %>%
  select(-contains('X3d7a'), -contains('w41')) %>%
  rename(Old_id = geneID, Old_name = Name) %>%
  mutate(Old_id = gsub('-a|-b|-a/b|-a/b/c', '', Old_id)) %>%
  left_join(gene_rosetta_blasted, by='Old_id') %>%
  select(Gene_id, everything())

## Add manually blasted ids
man_blast <- read_csv('./genes_manually_blasted.csv') %>%
  select(Old_id, Manual_Blast_Target) %>%
  filter(!is.na(Manual_Blast_Target))
gene_level[gene_level$Old_id %in% man_blast$Old_id,]$Gene_id <- man_blast$Manual_Blast_Target

gene_level <- gene_level %>% left_join(annot_df, by = 'Gene_id')
gene_level %>%
  left_join(annot_df, by = 'Gene_id') %>%
  write_csv('./raw_gene_level_data.csv')

gene_level %>%
  filter(is.na(Gene_id)) %>%
  write_csv('./raw_gene_level_noGene_id.csv')

## Remove genes without ids
gene_level <- gene_level %>%
  filter(!is.na(Gene_id))

## No Ratio

gene_level_noratio <- gene_level_noratio %>%
  select(-contains('X3d7a'), -contains('w41'))%>%
  rename(Old_id = geneID, Old_name = Name) %>%
  mutate(Old_id = gsub('-a|-b|-a/b|-a/b/c', '', Old_id)) %>%
  left_join(gene_rosetta_blasted, by='Old_id') %>%
  select(Gene_id, everything())

## Add manually blasted ids
gene_level_noratio[gene_level_noratio$Old_id %in% man_blast$Old_id,]$Gene_id <- man_blast$Manual_Blast_Target

gene_level_noratio <- gene_level_noratio %>% left_join(annot_df, by = 'Gene_id')
gene_level_noratio %>%
  left_join(annot_df, by = 'Gene_id') %>%
  write_csv('./raw_gene_level_noratio_data.csv')

gene_level_noratio %>%
  filter(is.na(Gene_id)) %>%
  write_csv('./raw_gene_level_noratio_noGene_id.csv')

## Remove genes without ids
gene_level_noratio <- gene_level_noratio %>%
  filter(!is.na(Gene_id))


## Summarize genes with duplicated new IDs

collapse_ids <- function(df) {
  non_num <- df %>%
    select(Gene_id, Name, Variant, Annot) %>%
    distinct()
  #print(dim(non_num))

  num <- df %>%
    select(Gene_id, contains('X')) %>%
    group_by(Gene_id) %>%
    summarise_each(funs(mean))

  #print(dim(num))
  df <- num %>%
    left_join(non_num, by = 'Gene_id') %>%
    filter(!is.na(Gene_id))


  return(df)
}


gene_level %>% group_by(Gene_id) %>% filter(n() > 1) %>%
  select(Gene_id, Name, Annot) %>%
  print(n = 99)

gene_level <- collapse_ids(gene_level)
gene_level_noratio <- collapse_ids(gene_level_noratio)

## Load Areas
areasDF <- read_csv('Variantome_Original/areas_subclons.csv')
areasDF <- areasDF %>%
  select(-contains('X3d7a'), -contains('w41')) %>%
  rename(Old_id = ...1, Old_name = Name) %>%
  mutate(Old_id = gsub('-a|-b|-a/b|-a/b/c', '', Old_id)) %>%
  left_join(gene_rosetta_blasted, by='Old_id') %>%
  select(Gene_id, everything())

## Add manually blasted ids
areasDF[areasDF$Old_id %in% man_blast$Old_id,]$Gene_id <- man_blast$Manual_Blast_Target

areasDF <- areasDF %>%
  filter(!is.na(Gene_id)) %>%
  select(Gene_id, contains('left'), contains('right'), contains('mid'), contains('sides'), -Old_id, -Old_name) %>%
  select(Gene_id, contains('1.2b'), contains('10g'), contains('3d7b')) %>%
  group_by(Gene_id) %>% filter(n() == 1) %>% ungroup()

colnames(areasDF) <- c("Gene_id",
                       "12B_Left", "12B_Right", "12B_Middle", "12B_Sides",
                       "10G_Left", "10G_Right", "10G_Middle", "10G_Sides",
                       "3D7B_Left","3D7B_Right", "3D7B_Middle", "3D7B_Sides")

areasDF <- as.data.frame(areasDF)
rownames(areasDF) <- areasDF$Gene_id
areasDF <-areasDF %>% select(-Gene_id)


##*********************************************************************##
##********************* END OF MODIFIABLE PART ************************##
##*********************************************************************##

#### Create folders for output ####

print("Creating folders for Output...")
dir.create(paste0(outdir))
dir.create(paste0(outdir, "/Plots"))
dir.create(paste0(outdir, "/Plots/Array_Plots"))
dir.create(paste0(outdir, "/Plots/MA_Plots"))
dir.create(paste0(outdir, "/Plots/Time_estim/"))
dir.create(paste0(outdir, "/Plots/Ratio"))
dir.create(paste0(outdir, "/Plots/Ratio/Gene_Level"))
dir.create(paste0(outdir, "/Plots/Ratio/Probe_Level"))
dir.create(paste0(outdir, "/Plots/Red_Signal"))
dir.create(paste0(outdir, "/Plots/Red_Signal/Gene_Level"))
dir.create(paste0(outdir, "/Plots/Red_Signal/Probe_Level"))
figPath <- paste0(outdir, "/Plots/")

#### Create Eset: xgene ####

exprsx <- as.matrix(gene_level %>% select(contains('X')))
colnames(exprsx) <- sample_names
rownames(exprsx) <- gene_level$Gene_id
fdata <- new("AnnotatedDataFrame", gene_level %>% select(Gene_id, Name, Variant, Annot))
rownames(fdata) <- gene_level$Gene_id
teor_time <- times
type <- types
pdata <- data.frame(type=type, teor_time=teor_time); rownames(pdata) <- sample_names
pdata <- new("AnnotatedDataFrame", pdata)
xgene <- new("ExpressionSet", exprs=exprsx, featureData=fdata, phenoData=pdata)
save(xgene, file=paste0(outdir, '/geneLevel.RData'))

exprsx <- as.matrix(gene_level_noratio %>% select(contains('X')) %>% select(contains('F635')))
rownames(exprsx) <- gene_level_noratio$Gene_id
colnames(exprsx) <- sample_names
xgene_red <- new("ExpressionSet", exprs=exprsx, featureData=fdata, phenoData=pdata)

write.csv(cbind(fdata@data, exprs(xgene)), paste0(outdir, "/geneLevel_exp.csv"), row.names=F)
write.csv(cbind(fdata@data, exprs(xgene_red)), paste0(outdir, "/geneLevel_redSignalexp.csv"), row.names=F)

#### Estimate times ####

ascendingTime <- function(x){
  current <- 0
  ncycle <- 0
  for (i in 1:length(x)){
    val  <- x[i]+(48*ncycle)
    if (val < current){
      current <- val
      val <- val+48
      ncycle <- ncycle+1
    }
    current <- val
    x[i] <- val
  }
  return(x)
}

estimatedTimes <- read_excel('Variantome_Original/Estimated_Times_3D7.xls')
estimatedTimes <- estimatedTimes %>%
  rename(Sample = `...1`, HPI = x) %>%
  filter(grepl('X1.2b', fixed = T, Sample) |
           grepl('X10g', fixed = T, Sample) |
           grepl('X3d7b', fixed = T, Sample)) %>%
  select(HPI) %>% pull()

#estimatedTimes <- getTimeEstimation(xgene,bozdechPath,LemieuxFunctionsPath,file.path(figPath),B=100)
estimatedTimes[estimatedTimes < 0] <- 0
hpi <- estimatedTimes

for (type in pData(xgene)$type){
  sel <- pData(xgene)$type == type
  typetime <- estimatedTimes[sel]
  time <- ascendingTime(typetime)
  estimatedTimes[sel] <- time
}

write.csv(estimatedTimes, paste0(outdir, "/Estimated_Times.csv"))
pData(xgene)$time <- estimatedTimes
pData(xgene_red)$time <- estimatedTimes
pData(xgene)$hpi <- hpi
pData(xgene_red)$hpi <- hpi

#### Areas: Functions ####

# imputePoint <- function(xs, ys, tp){
#
#   ## "xs" and "ys" must be two vectors of equal length
#   ## with the corresponding y(expression) and x(timepoint)
#   ## values that form the expression plot of interest (one gene).
#   ## "tp" must be the timepoint to impute.
#   ## If the timepoint to be imputed is already present, leave it as is.
#   ## Returns NA if missing the previous or next tp
#
#   if (tp %in% xs){
#
#     idx <- which(xs == tp)
#     imputed <- list(x=xs[idx], y=ys[idx])
#
#   } else {
#
#     before <- which(xs == max(xs[xs < tp]))
#     after <- which(xs == min(xs[xs > tp]))
#
#     if (is.na(ys[before]) | is.na(ys[after])){
#
#       imputed <- list(x=tp, y=NA)
#
#     } else {
#
#       x <- c(xs[c(before, after)])
#       y <- c(ys[c(before, after)])
#
#       imputed <- approx(x, y, xout=tp)
#     }
#   }
#   return(imputed)
# }
# computeArea <- function(eset){
#
#   ## Takes an eset and computes areas.
#   ## pData(eset) must contain a field named "time" with the time-points.
#   ## pData(eset) must have a field named "type" with the grouping variable.
#
#   ## Set needed variables
#   types <- unique(phenoData(eset)$type)
#   type <- phenoData(eset)$type
#   times <- phenoData(eset)$time
#
#   maxminTP <- max(sapply(types,function(x) min(times[type==x])))
#   minmaxTP <- min(sapply(types,function(x) max(times[type==x])))
#   mybreaks <- seq(maxminTP, minmaxTP, length.out=5)
#   tp1 <- mybreaks[2]
#   tp2 <- mybreaks[3]
#   tp3 <- mybreaks[4]
#
#   xsList <- c()
#   for (type in types){
#     xsList <- c(xsList, list(pData(eset)$time[phenoData(eset)$type == type]))
#   }
#
#   ## Main loop
#   all_areas <- c()
#   for (i in 1:dim(eset)[1]){
#
#     gene <- fData(eset)$geneID[i]
#
#     ysList <- c()
#     for (type in types){
#       ysList <- c(ysList, list(exprs(eset)[i, phenoData(eset)$type == type]))
#     }
#
#     ## Estimate points where needed
#     dfs <- list()
#     for (i in 1:length(xsList)){
#
#       x <- unlist(xsList[[i]])
#       y <- unlist(ysList[[i]])
#
#       points <- as.data.frame(cbind(x, y))
#       midpoints <- points[points$x > maxminTP &
#                             points$x < minmaxTP, ]
#
#       first <- imputePoint(x, y, maxminTP)
#       last <- imputePoint(x, y, minmaxTP)
#       p1 <- imputePoint(x, y, tp1)
#       p2 <- imputePoint(x, y, tp2)
#       p3 <- imputePoint(x, y, tp3)
#
#       impPoints <- rbind(first, last, p1, p2, p3)
#       allpoints <- rbind(midpoints, impPoints)
#       allpoints$x <- as.numeric(allpoints$x)
#       allpoints$y <- as.numeric(allpoints$y)
#
#       ordered <- arrange(allpoints, allpoints$x)
#
#       dfs[[i]] <- ordered
#     }
#
#     ## Calculate minY on estimated DFs
#     minY <- min(sapply(dfs, function(df) min(df$y, na.rm = T)))
#
#     rowareas <- c()
#     for (df in dfs){
#
#       df$y <- df$y - minY
#
#       ## Whole polygon from expression data
#       polDF <- rbind(df,
#                      c(minmaxTP, 0),
#                      c(maxminTP, 0),
#                      c(df[1,]))
#
#       ## Create Polygons
#       leftHalf <- rbind(df[which(df$x <= tp2),],
#                         c(tp2, 0),
#                         c(maxminTP, 0),
#                         c(df[1,]))
#
#       rightHalf <- rbind(df[which(df$x >= tp2),],
#                          c(minmaxTP, 0),
#                          c(tp2, 0),
#                          c(df[df$x == tp2,]))
#
#       mid <- rbind(df[which(df$x >= tp1 & df$x <= tp3),],
#                    c(tp3, 0),
#                    c(tp1, 0),
#                    c(df[df$x == tp1,]))
#
#       sides <- rbind(df[which(df$x <= tp1),],
#                      c(tp1, 0),
#                      c(tp3, 0),
#                      df[which(df$x >= tp3),],
#                      c(minmaxTP, 0),
#                      c(maxminTP, 0),
#                      df[1,])
#
#       pols <- list(leftHalf, rightHalf, mid, sides)
#
#       calcArea <- function(x) {ifelse(any(is.na(x)), NA, Polygon(x)@area)}
#       areas <- unlist(lapply(pols, function(x) calcArea(x)))
#
#       rowareas <- c(rowareas, areas)
#
#       ## Plot polygons (for debugging purposes)
#       ##pol <- Polygon(polDF)
#       ##ps = Polygons(list(pol),1)
#       ##sps = SpatialPolygons(list(ps))
#       ##plot(sps)
#     }
#     all_areas <- c(all_areas, list(rowareas))
#   }
#
#   ## Set row and col names for output
#   areaDF <- do.call(rbind, all_areas)
#   titles <- c("Left", "Right", "Middle", "Sides")
#
#   cols <- c()
#   for (i in types){
#     for (t in titles){
#       name <- paste0(i, "_", t)
#       cols <- c(cols, name)
#     }
#   }
#   colnames(areaDF)  <- cols
#   rownames(areaDF) <- rownames(exprs(eset))
#   return(areaDF)
# }

#### Calls: Compute Areas ####

print("Computing Areas...")
# areasDF <- computeArea(xgene)
# old_areasDF <- areasDF


head(areasDF)

areas_df <- as_tibble(areasDF)
areas_df['Gene_id'] <- rownames(areasDF)
areas_df <- areas_df %>%
  select(Gene_id, everything())

#xout <- cbind(fData(xgene), areasDF)
xout <- full_join(fData(xgene), areas_df)
write.csv(xout, paste0(outdir, "/area_geneLevel.csv"), row.names=F)

#### Calculate aAFC (average Area Fold-Change) ####

myMax <- function(x){
  y <- abs(x)
  if (all(is.na(y))){
    return(list(NA, NA))
  } else {
    pos <- which.max(y)
    times <- c("Left", "Right", "Mid", "Sides")
    return(list(maxVal = x[pos],
                maxTime = times[pos]))
  }
}

## Get number of categories.
n <- length(levels(pData(xgene)$type))

ns <- list()
i <- 1
while (i < n+1){
  ns[[i]] <- 1:4+(4*(i-1))
  i <- i+1
}

## Convert de areasDF into a list of DFs separated by types.
areas <- list()
for (i in 1:length(ns)) {
  areas[[i]] <- areasDF[,ns[[i]]]
}


## Calculate area differences

types <- unique(phenoData(xgene)$type)
type <- phenoData(xgene)$type
times <- phenoData(xgene)$time

maxminTP <- max(sapply(types,function(x) min(times[type==x])))
minmaxTP <- min(sapply(types,function(x) max(times[type==x])))
mybreaks <- seq(maxminTP, minmaxTP, length.out=5)

sets <- 1:length(areas)
combs <- combn(sets, 2)
span <- mybreaks[3] - mybreaks[1]
titles <- c("Left", "Right", "Middle", "Sides")


areaDifs <- list()
for (i in 1:dim(combs)[2]){

  one <- combs[1,i]
  two <- combs[2,i]

  dif1 <- as.data.frame((areas[[one]] - areas[[two]])/span)
  names  <- paste0(colnames(areas[[one]]),
                   "_minus_",
                   colnames(areas[[two]]))

  prefix <- paste0(strsplit(colnames(areas[[one]])[1], "_")[[1]][1],
                   "-",
                   strsplit(colnames(areas[[two]])[1], "_")[[1]][1])
  names <- sapply(titles, function(x) paste(prefix, x, sep = "_"))

  colnames(dif1) <- names

  maxval <- apply(dif1, 1, function(x) myMax(x)[[1]])
  maxtime <- apply(dif1, 1, function(x) myMax(x)[[2]])

  mv <- paste0(prefix, "_MaxVal")
  mt <- paste0(prefix, "_MaxTime")

  dif1[mv] <- maxval
  dif1[mt] <- maxtime

  #dif2 <- -dif1

  ## prefix <- paste0(strsplit(colnames(areas[[two]])[1], "_")[[1]][1],
  ##                  "-",
  ##                  strsplit(colnames(areas[[one]])[1], "_")[[1]][1])
  ## names <- sapply(titles, function(x) paste(prefix, x, sep = "_"))

  ## colnames(dif2) <- names

  areaDifs <- c(areaDifs, list(dif1))#, list(dif2))
}

allDifs <- do.call(cbind, areaDifs)
head(allDifs)

#### Areas: Convert Partitions into life stages ####

timetostage <- function(tp){
  while(tp > 48){
    tp = tp-48
  }
  return(tp)
}

getStage <- function(tp) {
  if ((tp >= 0) & (tp < 26)) {stg = "ring"}
  else if ((tp >= 26) & (tp < 38)) {stg = "troph"}
  else if ((tp >= 38) & (tp <= 48)) {stg = "schizont"}
  return(stg)
}

fracToStage <- function(frac, tps){
  if (is.na(frac)){
    return(NA)
  } else if (frac == "Left"){
    tp = tps[2]
    getStage(tp)
  } else if (frac == "Right"){
    tp = tps[4]
    getStage(tp)
  } else if (frac == "Mid"){
    tp = tps[3]
    getStage(tp)
  } else if (frac == "Sides"){
    tp1 = tps[1]+(span/4)
    tp2 = tps[4]+(span/4)
    stg1 <- getStage(tp1)
    stg2 <- getStage(tp2)
    return(paste0(stg1,"-",stg2))
  }
}


ncombs <- dim(combs)[2]

maxValCols <- seq(5, ncombs*6, 6)
maxTimeCols <- seq(6, ncombs*6, 6)

aMAFC <- cbind(allDifs[,c(maxValCols, maxTimeCols)])

timeCols <- (ncombs+1):(ncombs*2)

tps <- sapply(mybreaks, function(x) timetostage(x))

for (i in timeCols){
  aMAFC[,i] <- sapply(aMAFC[,i], function(x) fracToStage(x, tps))
}

write.csv(allDifs, paste0(outdir, "/areaDiferences_geneLevel.csv"))
write.csv(aMAFC, paste0(outdir, "/aMAFC_geneLevel.csv"))

#### Create Max differences table ####

anot_table <- annot %>%
  select(V2, V4, V5) %>%
  rename(Gene_id=V2, Name=V4, Annot=V5) %>%
  dplyr::filter(!is.na(Gene_id)) %>%
  distinct()

head(anot_table)
head(allDifs)

max_df <- allDifs %>%
  mutate(Gene_id = rownames(allDifs)) %>%
  left_join(anot_table, by='Gene_id') %>%
  mutate(MaxMax = pmax(abs(`12B-10G_MaxVal`), abs(`12B-3D7B_MaxVal`), abs(`10G-3D7B_MaxVal`))) %>%
  arrange(desc(abs(MaxMax))) %>%
  select(Gene_id, Name, Annot, contains('MaxVal'), contains('MaxTime'), contains('-'))

max_df_top <- max_df %>% dplyr::filter(abs(`12B-10G_MaxVal`) > 4 | abs(`12B-3D7B_MaxVal`) > 4 | abs(`10G-3D7B_MaxVal`) > 4)

write.csv(max_df, paste0(outdir, "/all_aMAFC.csv"))
write.csv(max_df_top, paste0(outdir, "/top_aMAFC.csv"))

#### PCA Plots ####

print("Plotting PCA..")
noNA <- xgene[complete.cases(exprs(xgene))]
df <- t(exprs(noNA))
df <- as.data.frame(df)

pca <- prcomp(df)
cmp1 <- format(round(summary(pca)$importance[2,1]*100, 2), nsmall = 2)
cmp2 <- format(round(summary(pca)$importance[2,2]*100, 2), nsmall = 2)

df_pca <- as.data.frame(pca$x)
df_pca$Type <- noNA@phenoData@data$type
df_pca$Time <- noNA@phenoData@data$time

p <- ggplot(df_pca, aes(x=PC1,y=PC2, col = Type, group = Type))
p <- p + geom_point(aes(size= Time))
p <- p + geom_path()
p <- p + scale_x_continuous(name=paste0("PC1: ", cmp1, "%"))
p <- p + scale_y_continuous(name=paste0("PC2: ", cmp2, "%"))

p <- p + theme_classic()
p <- p + theme(text = element_text(size=20))
p

ggsave(p, filename = paste0(figPath, "PCA.svg"), device = "svg")

#### Add gametocyte genes info ####
library(readxl)
gene_lists <- read_excel('../All gene lists_160719.xlsx', sheet = 2)

gam_genes <- gene_lists %>%
  select(`Lopez-Barragan`, Lasonder, Gametocites_Young, contains('gam'))

gam_list <- unique(gam_genes %>% pull())
gam_list[gam_list == 'NA'] <- NA
gam_list <- gam_list[!is.na(gam_list)]

finalDF <- max_df %>%
  select(Gene_id, Name, Annot, contains('MaxVal'), contains('MaxTime'))

finalDF['GamGene'] <- finalDF$Gene_id %in% gam_list

write.csv(finalDF, paste0(outdir, "/final_summary_table.csv"), row.names = F)


fData(xgene)

tibble(finalDF) %>%
  left_join(fData(xgene) %>% select(Gene_id, Variant), by = 'Gene_id') %>%
  filter(abs(`12B-10G_MaxVal`) > 1 |
         abs(`12B-3D7B_MaxVal`) > 1 |
         abs(`10G-3D7B_MaxVal`) > 1) %>%
  filter(is.na(Variant))

#### Var genes plot ####

gene_fam <- read_excel('/mnt/Disc4T/Projects/PhD_Project/Binned_Coverage/cvgfamilylist/Supplementary_table_2_CVG_list_161120_ap.xlsx', sheet = 2)
gene_fam <- gene_fam %>%
  rename(Gene_id = `Gene ID`,
         Gene_name = `Gene Name or Symbol`,
         SubFamily = `Family Detail`) %>%
  # mutate(SubFamily = case_when(SubFamily == 'var pseudo,truncated or -like.' ~ 'var-like',
  #                              TRUE ~ SubFamily)) %>%
  select(Gene_id, Gene_name, Family, SubFamily)

red_df <- as.data.frame(exprs(xgene_red))
red_df <- red_df %>%
  mutate(Gene_id = rownames(red_df)) %>%
  left_join(gene_fam, by='Gene_id')

varPlot <- function(strain){

  vals <- red_df %>% select(contains('tp') & contains(strain))
  max_tps <- colnames(vals)[apply(vals,1,which.max)]
  max_tps <- gsub(paste0(strain, '_'), '', max_tps)

  plot_df <- red_df %>%
    mutate(Max = select(., contains(strain)) %>% do.call(pmax, .)) %>%
    mutate(MaxTime = max_tps) %>%
    select(Gene_id, Max, MaxTime, Family, SubFamily) %>%
    dplyr::filter(Family == 'VAR' & SubFamily != 'var pseudo,truncated or -like.')

  write_csv(plot_df, paste0(outdir, 'max_VAR_redsignal_', strain, '.csv'))

  vars_plot <- ggplot(plot_df, aes(x = Gene_id, y = Max)) +
    geom_bar(stat = 'identity') +
    ggtitle(paste0('VAR gene expression ', strain)) +
    ylab(paste0('Normalized Cy5 signal')) +
    theme_classic() +
    theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
    coord_cartesian(ylim = c(0,60000))

  #print(vars_plot)

  ggsave(paste0(figPath, strain, '_var_genes_red.pdf'), vars_plot, device = 'pdf')

  ## Tp 10 Plot
  plot10tp_df <- red_df %>%
    select(Gene_id, contains(strain) & contains('tp10'), Family, SubFamily) %>%
    set_names(c('Gene_id', 'RedSignal', 'Family', 'SubFamily')) %>%
    dplyr::filter(Family == 'VAR' & SubFamily != 'var pseudo,truncated or -like.')

  write_csv(plot10tp_df, paste0(outdir, 'max_VAR_redsignal_tp10_', strain, '.csv'))


  vars10tp_plot <- ggplot(plot10tp_df, aes(x = Gene_id, y = RedSignal)) +
    geom_bar(stat = 'identity') +
    ggtitle(paste0('VAR gene expression ', strain)) +
    ylab(paste0('Normalized Cy5 signal')) +
    theme_classic() +
    theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
    coord_cartesian(ylim = c(0,60000))

  #print(vars10tp_plot)
  ggsave(paste0(figPath, strain, '_var_genes_red_tp10.pdf'), vars10tp_plot, device = 'pdf')
}


strains <- c('12B', '10G', '3D7B')
for (strain in strains) {varPlot(strain)}

#### Single Gene Plot ####
plot_gene <- function(type, gid, out){
  ## Set df and lims depending on what we are plotting.
  if (type == "gene"){
    df = xgene
    path = "/Ratio/Gene_Level/"

  } else if (type == "gene_red"){
    df = xgene_red
    path = "/Red_Signal/Gene_Level/"

  } else if (type == "probe"){
    df = xprobe
    path = "/Ratio/Probe_Level/"

  } else if (type == "probe_red"){
    df = xprobe_red
    path = "/Red_Signal/Probe_Level/"

  }

  ## Set ylims
  ylim = c(min(exprs(df), na.rm = T), max(exprs(df), na.rm = T))


  ## Plot
  ## Set gene for title or gene and probe for probe-level plots.
  gn <- gsub("[/:;.]", "_", gid)
  if (type %in% c("probe", "probe_red")){
    prb <- paste0("_", gsub("[/:;.]", "_" , gid))
  } else {
    prb <- ""
  }
  title <- paste0(gn, prb)

  ## Set y-axis title depending on ratio/red_signal
  ytitle <- ifelse(type %in% c("gene_red", "probe_red"), 'log2(Cy5)', 'log2(Cy3/Cy5)')

  ## Plot
  graf <- melt(df[fData(df)$Gene_id == gid,])
  graf["Type"] <- xgene@phenoData@data$type
  graf["Time"] <- xgene@phenoData@data$time
  p <- ggplot(graf, aes(x = Time, y = value, col = Type, group = Type))
  p <- p + geom_point(aes(color = Type, size = 2))
  p <- p + geom_line(aes(size = 2))
  p <- p + coord_cartesian(ylim = ylim)
  p <- p + theme_classic()
  # p <- p + ggtitle(title)
  # p <- p + ylab(ytitle)
  p <- p + theme(text = element_text(size = 36))
  p <- p + theme(axis.title.x = element_blank())
  p <- p + theme(axis.title.y = element_blank())
  p <- p + theme(legend.position = 'none')
  ggsave(p, file=out, device = "svg")
  print(p)
}

type <- 'gene'
gid <- 'PF3D7_0302200'
outpath <- '/home/lucas/Documents/BioMalPar_2021/Microarrays/Gene_plots/'
outname <- 'PF3D7_0302200_12b10g.svg'

plot_gene(type, gid, paste0(outpath, outname))

xgene_red_tibble <- as.data.frame(exprs(xgene_red)) %>%
  tibble() %>%
  mutate(Gene_id = rownames(exprs(xgene_red))) %>%
  select(Gene_id, everything())

## First get maxcol then percentile (OLD)
## get_red_percent <- function(strain){

##   ## Subset to strain
##   red_strain <- xgene_red_tibble %>%
##     select(Gene_id, contains(strain))

##   ## Calculate Max col
##   maxred <- red_strain %>%
##     dplyr::select(-Gene_id) %>%
##     mutate(Max_Red = do.call(pmax, (.))) %>%
##     select(Max_Red)

##   ## Create "percentile" function from Max col
##   percentile <- ecdf(maxred$Max_Red)
##   red_pcnt <- percentile(maxred$Max_Red)

##   return(red_pcnt)
## }

## perc_12B <- get_red_percent('12B')
## perc_10G <- get_red_percent('10G')

## red_percent <- tibble('Gene_id' = xgene_red_tibble$Gene_id,
##                       '12B' = perc_12B,
##                       '10G' = perc_10G)


## First percentile per col -> max percentile (NEW)

my_percentile <- function(vector){
  ecdf(vector)(vector)*100
}

xgene_red_tibble <- xgene_red_tibble %>%
  mutate(across(.cols = -Gene_id, .fns = my_percentile, .names = "Perc_{.col}"))

get_max_percent <- function(strain){

  ## Subset to strain
  red_strain <- xgene_red_tibble %>%
    select(Gene_id, contains('Perc') & contains(strain))

  ## Calculate Max col
  maxperc <- red_strain %>%
    dplyr::select(-Gene_id) %>%
    mutate(Max_Perc = do.call(pmax, c(., na.rm = T))) %>%
    select(Max_Perc) %>%
    pull()
  return(maxperc)
}

perc_12B <- get_max_percent('12B')
perc_10G <- get_max_percent('10G')
perc_3D7B <- get_max_percent('3D7B')

red_percent <- tibble('Gene_id' = xgene_red_tibble$Gene_id,
                      '12B' = perc_12B,
                      '10G' = perc_10G,
                      '3D7B' = perc_3D7B)

##hist(red_percent$`12B`)

write_csv(red_percent, paste0(outdir, "/red_percentiles.csv"))
names(red_percent)[2:4] <- paste0('Red_Pcnt_', names(red_percent)[2:4])

max_tibble <- as_tibble(max_df)
max_tibble <- max_tibble %>%
  left_join(red_percent, by='Gene_id', suffix = c('', '_Pcnt'))

maxFC_passe_red <- max_tibble %>%
  filter((`12B-10G_MaxVal` >= 2 & Red_Pcnt_12B > 15) |
         (`12B-10G_MaxVal` <= -2 & Red_Pcnt_10G > 15) |
         (`12B-3D7B_MaxVal` >= 2 & Red_Pcnt_12B > 15) |
         (`12B-3D7B_MaxVal` <= -2 & Red_Pcnt_3D7B > 15) |
         (`10G-3D7B_MaxVal` >= 2 & Red_Pcnt_10G > 15) |
         (`10G-3D7B_MaxVal` <= -2 & Red_Pcnt_3D7B > 15))

#### Get MaxTime for each gene ####

myWhichMax <- function(vect){
  if (all(is.na(vect))){
    return(NA)
  } else {
    return(which.max(vect))
  }
}

exp <- as.data.frame(exprs(xgene))

maxcol <- exp %>%
  apply(1, myWhichMax) %>%
  unlist()

times <- pData(xgene) %>%
  select(time) %>%
  pull()

maxtime <- sapply(maxcol, function(x) times[x])
head(maxtime)

max_time <- tibble(Gene_id = names(maxtime), Max_Time = maxtime)
breaks_df <- tibble(Areas_Breaks = mybreaks)

tibble(Gene_id = names(maxtime), Max_Time = maxtime) %>%
  write_csv(paste0(outdir, 'old_arrays_maxtime.csv'))

tibble(Areas_Breaks = mybreaks) %>%
  write_csv(paste0(outdir, 'old_area_breaks.csv'))

## #### Filter by Max-Time ####

## New approach
## Check which areas does maxtimepoint overlapp -> check if aAFC > th at this areas

point_overlap <- function(point, interval){
  point >= interval[1] & point <= interval[2]
}

areas_df <- as_tibble(allDifs) %>%
  mutate(Gene_id = rownames(allDifs)) %>%
  select(Gene_id, everything())

maxtimes_12B_10G <- c()
maxtimes_12B_3D7B <- c()
maxtimes_10G_3D7B <- c()
gids <- c()
th <- 2
for (gid in max_tibble$Gene_id){

  #gid <- 'PF3D7_0324600'
  ## Create time-regions
  breaks <- breaks_df$Areas_Breaks
  left <- c(breaks[1], breaks[3])
  right <- c(breaks[3], breaks[5])
  mid <- c(breaks[2], breaks[4])
  sides_l <- c(breaks[1], breaks[2])
  sides_r <- c(breaks[4], breaks[5])

  ## Get maxtime
  maxtime <- max_time %>%
    filter(Gene_id == gid) %>%
    pull()
  if (is.na(maxtime)){
    maxtimes_12B_10G <- c(maxtimes_12B_10G, NA)
    maxtimes_12B_3D7B <- c(maxtimes_12B_3D7B, NA)
    maxtimes_10G_3D7B <- c(maxtimes_10G_3D7B, NA)
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
    aFCs <- areas_df %>%
      filter(Gene_id == gid) %>%
      select(contains(names(areas[overlaps])))

    fc_12B_10G <- aFCs %>%
      select(contains('12B-10G')) %>%
      replace(is.na(.), 0) %>%
      as.numeric()
    maxtime_FC_12B_10G <- any(abs(fc_12B_10G) > th)

    fc_12B_3D7B <- aFCs %>%
      select(contains('12B-3D7B')) %>%
      replace(is.na(.), 0) %>%
      as.numeric()
    maxtime_FC_12B_3D7B <- any(abs(fc_12B_3D7B) > th)

    fc_10G_3D7B <- aFCs %>%
      select(contains('10G-3D7B')) %>%
      replace(is.na(.), 0) %>%
      as.numeric()
    maxtime_FC_10G_3D7B <- any(abs(fc_10G_3D7B) > th)

    maxtimes_12B_10G <- c(maxtimes_12B_10G, maxtime_FC_12B_10G)
    maxtimes_12B_3D7B <- c(maxtimes_12B_3D7B, maxtime_FC_12B_3D7B)
    maxtimes_10G_3D7B <- c(maxtimes_10G_3D7B, maxtime_FC_10G_3D7B)
    gids <- c(gids, gid)
  }
}
maxtime_aAFC_df <- tibble(
  Gene_id = gids,
  MaxTime_Filter_12B_10G = maxtimes_12B_10G,
  MaxTime_Filter_12B_3D7B = maxtimes_12B_3D7B,
  MaxTime_Filter_10G_3D7B = maxtimes_10G_3D7B
  )

maxtime_aAFC_df %>%
  count(MaxTime_Filter_12B_10G)

max_tibble <- max_tibble %>%
  left_join(maxtime_aAFC_df)


## final_llcm['MaxTime'] <- maxtime
## final_ll0['MaxTime'] <- maxtime
## final_lls['MaxTime'] <- maxtime

## tibble(finalDF)[,1:4]

## trans_df <- tibble(allDifs) %>%
##   mutate(Gene_id = rownames(allDifs)) %>%
##   select(Gene_id, everything())

## trans_df


## ## Define intervals
## left <- c(mybreaks[1], mybreaks[2])
## right <- c(mybreaks[3], mybreaks[5])
## mid <- c(mybreaks[2], mybreaks[4])
## sides1 <- c(mybreaks[1], mybreaks[2])
## sides2 <- c(mybreaks[4], mybreaks[5])

## checkOverap <- function(v1, v2){
##   v1[1] <= v2[2] & v2[1] <= v1[2]
## }

## getInterval <- function(x, width){
##   ## Will break if interval spans > 48h ("circular interval")
##   x_low <- ifelse(x-width >= mybreaks[1], x-width, mybreaks[1])
##   x_high <- ifelse(x+width <= mybreaks[5], x+width, mybreaks[5])
##   max_interval <- c(x_low, x_high)
##   return(max_interval)
## }

## ## Set width arround maxtime
## width <- 10

## ## Check wether each interval overlaps with MaxTime
## trans_df['Left_in_MaxTime'] <- sapply(maxtime, function(x) checkOverap(getInterval(x, width), left))

## final_llcm['Left_in_MaxTime'] <- sapply(maxtime, function(x) checkOverap(getInterval(x, width), left))
## final_llcm['Right_in_MaxTime'] <- sapply(maxtime, function(x) checkOverap(getInterval(x, width), right))
## final_llcm['Mid_in_MaxTime'] <- sapply(maxtime, function(x) checkOverap(getInterval(x, width), mid))
## final_llcm['Sides_in_MaxTime'] <- sapply(maxtime, function(x) checkOverap(getInterval(x, width), sides1) | checkOverap(getInterval(x, width), sides2))

## final_ll0['Left_in_MaxTime'] <- sapply(maxtime, function(x) checkOverap(getInterval(x, width), left))
## final_ll0['Right_in_MaxTime'] <- sapply(maxtime, function(x) checkOverap(getInterval(x, width), right))
## final_ll0['Mid_in_MaxTime'] <- sapply(maxtime, function(x) checkOverap(getInterval(x, width), mid))
## final_ll0['Sides_in_MaxTime'] <- sapply(maxtime, function(x) checkOverap(getInterval(x, width), sides1) | checkOverap(getInterval(x, width), sides2))

## final_lls['Left_in_MaxTime'] <- sapply(maxtime, function(x) checkOverap(getInterval(x, width), left))
## final_lls['Right_in_MaxTime'] <- sapply(maxtime, function(x) checkOverap(getInterval(x, width), right))
## final_lls['Mid_in_MaxTime'] <- sapply(maxtime, function(x) checkOverap(getInterval(x, width), mid))
## final_lls['Sides_in_MaxTime'] <- sapply(maxtime, function(x) checkOverap(getInterval(x, width), sides1) | checkOverap(getInterval(x, width), sides2))

## head(final_llcm)

## ## Check wether max falls in MaxTime (contrast by contrast)
## final_llcm <- final_llcm %>%
##   mutate(Max_in_MaxTime = case_when(MaxFC_Interval == "left" & Left_in_MaxTime ~ TRUE,
##                                     MaxFC_Interval == "right" & Right_in_MaxTime ~ TRUE,
##                                     MaxFC_Interval == "mid" & Mid_in_MaxTime ~ TRUE,
##                                     MaxFC_Interval == "sides" & Sides_in_MaxTime ~ TRUE,
##                                     TRUE ~ FALSE))

## final_ll0 <- final_ll0 %>%
##   mutate(Max_in_MaxTime = case_when(MaxFC_Interval == "left" & Left_in_MaxTime ~ TRUE,
##                                     MaxFC_Interval == "right" & Right_in_MaxTime ~ TRUE,
##                                     MaxFC_Interval == "mid" & Mid_in_MaxTime ~ TRUE,
##                                     MaxFC_Interval == "sides" & Sides_in_MaxTime ~ TRUE,
##                                     TRUE ~ FALSE))

## final_lls <- final_lls %>%
##   mutate(Max_in_MaxTime = case_when(MaxFC_Interval == "left" & Left_in_MaxTime ~ TRUE,
##                                     MaxFC_Interval == "right" & Right_in_MaxTime ~ TRUE,
##                                     MaxFC_Interval == "mid" & Mid_in_MaxTime ~ TRUE,
##                                     MaxFC_Interval == "sides" & Sides_in_MaxTime ~ TRUE,
##                                     TRUE ~ FALSE))

f_path <- '/mnt/Disc4T/Projects/PhD_Project/Paper/Paper_Analysis/Data_Files/Duplication_Deletion_Regions_Mean_Separate_DuplDel/Crossed_with_genes/'

file_list <- c('1.2B_in_sort_q5_noDup_rpkm_normInput_bs10_smth200_bymean_dupl_del_fact_1.75up_0.1dw_minlen500_mergelen_200_filtered_genes.tsv',
               '10G_in_sort_q5_noDup_rpkm_normInput_bs10_smth200_bymean_dupl_del_fact_1.75up_0.1dw_minlen500_mergelen_200_filtered_genes.tsv'
)

dupl_dl_12B <- read_tsv(paste0(f_path, file_list[1]), col_names = F) %>%
  select(X1) %>% pull()

dupl_dl_10G <- read_tsv(paste0(f_path, file_list[2]), col_names = F) %>%
  select(X1) %>% pull()

## Load Latest Annot
info_df <- read_csv(paste0(infiles, 'info_df.csv')) %>%
  dplyr::filter(!is.na(Gene_id))

final_df <- max_tibble %>%
  mutate(Gene_id = ifelse(Gene_id == 'PF3D7_0935400_as', 'PF3D7_0935390', Gene_id)) %>%
  mutate(Not_Plasmodium = !Gene_id %in% info_df$Gene_id)

final_df <- final_df %>%
  select(-Name, -Annot) %>%
  left_join(info_df, by='Gene_id')

colnames(final_df)

all_df <- final_df %>%
  select(
    Gene_id,
    contains('MaxVal'),
    contains('MaxTime'),
    contains('Red'),
    contains('Filter'),
  )

## Set thresholds
red_th <- 15

## 12B vs 10G

final_df %>%
  select(
    Gene_id,
    `12B-10G_MaxVal`,
    `12B-10G_MaxTime`,
    Red_Pcnt_12B,
    Red_Pcnt_10G,
    MaxTime_Filter_12B_10G
  ) %>%
  mutate(PassRed = ifelse(
           `12B-10G_MaxVal` >= 0,
           Red_Pcnt_12B >= red_th,
           Red_Pcnt_10G >= red_th
         )) %>%
  rename(PassMaxtime = MaxTime_Filter_12B_10G) %>%
  mutate(PassDuplDel = !Gene_id %in% dupl_dl_12B & !Gene_id %in% dupl_dl_10G) %>%
  rowwise() %>%
  mutate(PassAll = all(c(PassRed, PassMaxtime, PassDuplDel))) %>%
  write_tsv(paste0(outdir, '12B_10G_final_df.tsv'))

## 12B vs 3D7B

final_df %>%
  select(
    Gene_id,
    `12B-3D7B_MaxVal`,
    `12B-3D7B_MaxTime`,
    Red_Pcnt_12B,
    Red_Pcnt_3D7B,
    MaxTime_Filter_12B_3D7B
  ) %>%
  mutate(PassRed = ifelse(
           `12B-3D7B_MaxVal` >= 0,
           Red_Pcnt_12B >= red_th,
           Red_Pcnt_3D7B >= red_th
         )) %>%
  rename(PassMaxtime = MaxTime_Filter_12B_3D7B) %>%
  mutate(PassDuplDel = !Gene_id %in% dupl_dl_12B) %>%
  rowwise() %>%
  mutate(PassAll = all(c(PassRed, PassMaxtime, PassDuplDel))) %>%
  write_tsv(paste0(outdir, '12B_3D7B_final_df.tsv'))

## 10G vs 3D7B

final_df %>%
  select(
    Gene_id,
    `10G-3D7B_MaxVal`,
    `10G-3D7B_MaxTime`,
    Red_Pcnt_10G,
    Red_Pcnt_3D7B,
    MaxTime_Filter_10G_3D7B
  ) %>%
  mutate(PassRed = ifelse(
           `10G-3D7B_MaxVal` >= 0,
           Red_Pcnt_10G >= red_th,
           Red_Pcnt_3D7B >= red_th
         )) %>%
  rename(PassMaxtime = MaxTime_Filter_10G_3D7B) %>%
  mutate(PassDuplDel = !Gene_id %in% dupl_dl_10G) %>%
  rowwise() %>%
  mutate(PassAll = all(c(PassRed, PassMaxtime, PassDuplDel))) %>%
  write_tsv(paste0(outdir, '10G_3D7B_final_df.tsv'))


###################################################3


maxFC_pass_red <- final_df %>%
  filter((`12B-10G_MaxVal` >= 2 & Red_Pcnt_12B > 15) |
         (`12B-10G_MaxVal` <= -2 & Red_Pcnt_10G > 15) |
         (`12B-3D7B_MaxVal` >= 2 & Red_Pcnt_12B > 15) |
         (`12B-3D7B_MaxVal` <= -2 & Red_Pcnt_3D7B > 15) |
         (`10G-3D7B_MaxVal` >= 2 & Red_Pcnt_10G > 15) |
         (`10G-3D7B_MaxVal` <= -2 & Red_Pcnt_3D7B > 15))

## 12B vs 10G

difs_12B_10G <- final_df %>%
  filter(
  ((`12B-10G_MaxVal` >= 2 & Red_Pcnt_12B > 15) |
   (`12B-10G_MaxVal` <= -2 & Red_Pcnt_10G > 15)) &
  MaxTime_Filter_12B_10G
  ) %>%
  mutate(Dupl_Del = Gene_id %in% dupl_dl_12B | Gene_id %in% dupl_dl_10G) %>%
    dplyr::filter(!Not_Plasmodium & !Is_tRNA) %>%
  select(Gene_id, Name, Annot,
         `12B-10G_MaxVal`,
         Red_Pcnt_12B, Red_Pcnt_10G,
         MaxTime_Filter_12B_10G,
         Variant, Gam_specific, Dupl_Del)

write_csv(difs_12B_10G, paste0(outdir, '12B_vs_10G_log2FC2_red15_maxtime.csv'))


## 12B vs 3D7B

difs_12B_3D7B <- final_df %>%
  filter(
  ((`12B-3D7B_MaxVal` >= 2 & Red_Pcnt_12B > 15) |
   (`12B-3D7B_MaxVal` <= -2 & Red_Pcnt_3D7B > 15)) &
  MaxTime_Filter_12B_3D7B
  ) %>%
  mutate(Dupl_Del = Gene_id %in% dupl_dl_12B) %>%
  select(Gene_id, Name, Annot,
         `12B-3D7B_MaxVal`,
         Red_Pcnt_12B, Red_Pcnt_3D7B,
         MaxTime_Filter_12B_3D7B,
         Variant, Gam_specific, Dupl_Del)

write_csv(difs_12B_3D7B, paste0(outdir, '12B_vs_3D7B_log2FC2_red15_maxtime.csv'))


## 10G vs 3D7B

difs_10G_3D7B <- final_df %>%
  filter(
  ((`10G-3D7B_MaxVal` >= 2 & Red_Pcnt_10G > 15) |
   (`10G-3D7B_MaxVal` <= -2 & Red_Pcnt_3D7B > 15)) &
  MaxTime_Filter_10G_3D7B
  ) %>%
  mutate(Dupl_Del = Gene_id %in% dupl_dl_10G) %>%
  select(Gene_id, Name, Annot,
         `10G-3D7B_MaxVal`,
         Red_Pcnt_10G, Red_Pcnt_3D7B,
         MaxTime_Filter_10G_3D7B,
         Variant, Gam_specific, Dupl_Del)


write_csv(difs_10G_3D7B, paste0(outdir, '10G_vs_3D7B_log2FC2_red15_maxtime.csv'))

write_csv(final_df, paste0(outdir, 'old_arrays_final_df.csv'))

library(eulerr)

A <- difs_12B_10G$Gene_id
B <- difs_12B_3D7B$Gene_id
C <- difs_10G_3D7B$Gene_id

AB <- intersect(A, B)
AC <- intersect(A, C)
BC <- intersect(B, C)

ABC <- intersect(AB, C)

abc <- length(ABC)
ab <- length(AB[!AB %in% ABC])
ac <- length(AC[!AC %in% ABC])
bc <- length(BC[!BC %in% ABC])

a <- length(A) -ab -ac -abc
b <- length(B) -ab -bc -abc
c <- length(C) -ac -bc -abc

fit <- euler(c(A=a, B=b, C=c, "A&B"=ab, "A&C"=ac, "B&C"=bc, "A&B&C" = abc))

scales::viridis_pal()(3)

d <- plot(fit, fills = list(fill = c('#440154FF', "#21908CFF", "#FDE725FF"), alpha = 0.5),
          edges = list(lwd = 0.1),
          quantities = list(quantities = T),
          labels = list(labels=c("1.2B vs 10G", "1.2B vs 3D7B", "10G vs 3D7B")))

ggsave(d, filename = paste0(figPath, "Difs_Venn.pdf"), device = "pdf",
       width = 15, height = 15, units = 'cm')

plot(d)
print(fit)

#### Save environtment ####
##save.image(file = "array_12B10G3D7B_VariantomeData_work_space.RData")
##setwd('/mnt/Disc4T/Projects/PhD_Project/Microarrays/New_Old_separate_approach/Old_Arrays/')
##load('array_12B10G3D7B_VariantomeData_work_space.RData')
##head(finalDF)
