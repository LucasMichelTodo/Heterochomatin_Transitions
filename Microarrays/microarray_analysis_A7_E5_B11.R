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

############## Experiment Setup: MODIFY THIS PART #######################
##*********************************************************************##
##*********************************************************************##

setwd('/mnt/Disc4T/Projects/PhD_Project/Microarrays/New_Old_separate_approach/New_Arrays/')
datadir <- "./RawData/"
outdir <- "./R_results_NewArray/"

annot <- read.csv("./Files/array_anotation.csv", sep="\t", header = F)
## Add new GDV1as annotation
levels(annot$V2) <- c(levels(annot$V2), 'PF3D7_0935390')
annot[annot$V2 == 'PF3D7_0935400_as' & !is.na(annot$V2),]$V2 <- 'PF3D7_0935390'

infiles <- "./Files/"

sample_names <- sub("\\.txt$", "", list.files(datadir, pattern="\\.txt$"))
nsamples <- length(sample_names)

times <- as.integer(factor(sub("^.+_tp", "", sample_names)))
types <- sub("_tp.+", "", sample_names)

array_list <- lapply(list.files(datadir, pattern="\\.txt$", full.names=TRUE), read.table, sep="\t", stringsAsFactors=FALSE, skip=9, header=TRUE)
names(array_list) <- sample_names

run_all_plots <- "no"

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
dir.create(paste0(outdir, "/Plots/Ratio/Gene_Level"))
dir.create(paste0(outdir, "/Plots/Red_Signal"))
dir.create(paste0(outdir, "/Plots/Red_Signal/Gene_Level"))
dir.create(paste0(outdir, "/Plots/Red_Signal/Probe_Level"))
figPath <- paste0(outdir, "/Plots/")

#### Load Array annotation, Gene-list and Variant Genes ####

print("Loading Array Annotation...")
gene_list <- readLines(paste0(infiles, "gene_list.txt"))
gene_list <- c(gene_list, 'PF3D7_0935390')

cvgs <- read.csv2(paste0(infiles, 'taula_CVG_final.csv'), stringsAsFactors = F)
cvgs <- cvgs %>%
  select(Gene_id = Gene.ID, Variant = Final.Customized) %>%
  mutate(Variant = ifelse(Variant == 'YES', TRUE, FALSE)) %>%
  add_row(Gene_id = 'PF3D7_0935390', Variant = TRUE)

annot_df <- read_csv(paste0(infiles, 'info_df.csv')) %>%
  select(Gene_id, Name, Annot) %>%
  dplyr::filter(!is.na(Gene_id)) %>%
  left_join(cvgs, by = 'Gene_id') %>%
  replace_na(list(Variant = FALSE)) %>%
  select(Gene_id, Name, Variant, Annot)

#### Create Probe-DF ####

print("Creating Probe DF...")
probe_df <- array_list[[1]][,c(7,11,14,15)]
probe_df_geo <- array_list[[1]][,c(7,11,14,15)]

getCols <- function(df){
  return(df[,c(11,14,15)])
}

goodCols <- lapply(array_list[2:nsamples], function(x) getCols(x))
df <- do.call("cbind", goodCols)

probe_df <- cbind(probe_df, df)
probe_df_geo <- cbind(probe_df_geo, df)
dim(probe_df)
dim(probe_df_geo)

probe_df["Gene_id"] <- annot$V2
probe_df <- probe_df %>%
  left_join(annot_df, by = 'Gene_id')

probe_df_geo["Gene_id"] <- annot$V2
probe_df_geo <- probe_df_geo %>%
  left_join(annot_df, by = 'Gene_id')

## probe_df["name"] <- annot$V4
## probe_df["Annot"] <- annot$V5

probe_df["Annot"] <- gsub("Plasmodium", "Pl.", probe_df$Annot)
probe_df["Annot"] <- gsub("protein", "prot.", probe_df$Annot)
probe_df["Annot"] <- gsub("membrane", "memb.", probe_df$Annot)
probe_df["Annot"] <- gsub("conserved", "cvd.", probe_df$Annot)
probe_df["Annot"] <- gsub("function", "func.", probe_df$Annot)
probe_df["Annot"] <- gsub("unknown", "ukwn.", probe_df$Annot)
probe_df["Annot"] <- gsub("exported", "xptd.", probe_df$Annot)
probe_df["Annot"] <- gsub("pseudogene", "pseudo", probe_df$Annot)
probe_df["Annot"] <- gsub("putative", "put.", probe_df$Annot)
probe_df["Annot"] <- gsub("%2C", "", probe_df$Annot)

## Remove probes that map to multiple genes.
probe_df <- probe_df[annot$V3 != "drop" | is.na(annot$V3),]

probe_df_geo['Drop_Mulimapper'] <- annot$V3 == "drop"
head(probe_df_geo)


## ## Add Variant Genes information
## varlist <- dplyr::filter(cvgs, Variant) %>% select(Gene_id)
## probe_df["Variant"] <- probe_df$Gene_id %in% varlist

#### Group Columns ####

signalCols <- nsamples*3+1
allcols <- dim(probe_df)[2]

ratioCols <- seq(2, signalCols, 3)
redCols <- seq(4, signalCols, 3)
greCols <- seq(3, signalCols, 3)

infoCols <- c(1, (signalCols+1):allcols)

#### Remove low expression probes ####

print("Removing low expression probes...")
medians <- list()
for (i in 1:nsamples){
  redm <- median(sort(probe_df[probe_df$Gene_id %in% gene_list, redCols][,i])[1:100])
  grenm <- median(sort(probe_df[probe_df$Gene_id %in% gene_list, greCols][,i])[1:100])
  medians[[i]] <- c(grenm, redm)
}

passTest <- list()
for(i in 1:nsamples){
  g <- probe_df[,greCols][i] < 3*medians[[i]][1]
  r <- probe_df[,redCols][i] < 3*medians[[i]][2]
  all <- g & r
  passTest[[i]] <- !all
}

testDF <- as.data.frame(passTest)
pass <- rowSums(testDF) > 0
write.csv(table(!pass), paste0(outdir, "/NA_probes.csv"))
probe_df[!pass,c(ratioCols)] <- NA

## Mark low expression probes for GEO df
medians_geo <- list()
for (i in 1:nsamples){
  redm <- median(sort(probe_df_geo[probe_df_geo$Gene_id %in% gene_list, redCols][,i])[1:100])
  grenm <- median(sort(probe_df_geo[probe_df_geo$Gene_id %in% gene_list, greCols][,i])[1:100])
  medians_geo[[i]] <- c(grenm, redm)
}

passTest_geo <- list()
for(i in 1:nsamples){
  g <- probe_df_geo[,greCols][i] < 3*medians[[i]][1]
  r <- probe_df_geo[,redCols][i] < 3*medians[[i]][2]
  all <- g & r
  passTest_geo[[i]] <- !all
}

testDF_geo <- as.data.frame(passTest_geo)
pass_geo <- rowSums(testDF_geo) > 0
length(pass_geo)
probe_df_geo['Drop_LowExpress'] <- !pass_geo

#### Array Plots ####

print("Plotting Arrays...")
cols <- rev(brewer.pal(11, 'Spectral'))

arrayPlot <- function(df) {
  df_name <- sample_names[i]
  p1 <- qplot(Col, Row, data=df, color=log2(rMedianSignal)<7) + scale_color_manual(values=c("aliceblue", "black")) + ggtitle(df_name)
  ggsave(p1, filename = paste0(figPath, "Array_Plots/sample_", df_name, "_boolean.jpeg"), device = "jpeg")

  p2 <- qplot(Col, Row, data=df, color=log2(rMedianSignal)) + scale_colour_gradientn(colours = cols) + ggtitle(df_name)
  ggsave(p2, filename = paste0(figPath, "Array_Plots/sample_", df_name, ".jpeg"), device = "jpeg")

  p3 <- qplot(Col, Row, data=df, color=is.na(LogRatio)) + scale_color_manual(values=c("aliceblue", "red")) + ggtitle(df_name)
  ggsave(p3, filename = paste0(figPath, "Array_Plots/sample_", df_name, "_NAs.jpeg"), device = "jpeg")
}

for (i in 1:length(array_list)) {arrayPlot(array_list[[i]])}

#### MA Plots ####

print("Plotting MA Plots...")
myMAplot  <- function(mray){
  df_name <- sample_names[i]
  m_vals <- log2(mray$gProcessedSignal) - log2(mray$rProcessedSignal)
  a_vals <- (log2(mray$gProcessedSignal) + log2(mray$rProcessedSignal))/2
  ma_df <- cbind(a_vals, m_vals)
  p <- ggplot(ma_df, aes(x=a_vals, y=m_vals))
  p <- p + geom_point()
  p <- p + geom_smooth(method = "lm", se=F, color= "red")
  p <- p + geom_hline(yintercept=0, color = "blue", size = 1)
  ggsave(p, filename = paste0(figPath, "MA_Plots/sample_", df_name, "_MA.jpeg"), device = "jpeg")
}

for (i in 1:length(array_list)) {myMAplot(array_list[[i]])}

#### Change to Log2 (Log Ratio cols, originally log10) ####

print("Creating Eset...")
probe_df[,ratioCols] <- log2(10**probe_df[,ratioCols])
probe_df_geo[,ratioCols] <- log2(10**probe_df_geo[,ratioCols])

#### Change to Log2 (Raw signal Cols, originally unlogged) ####

probe_df[,redCols]  <- log2(probe_df[,redCols])
probe_df_geo[,redCols]  <- log2(probe_df_geo[,redCols])

#### Create eSet: xprobe ####

exprsx <- as.matrix(probe_df[,ratioCols])
colnames(exprsx) <- sample_names
fdata <- new("AnnotatedDataFrame", probe_df[,infoCols])
teor_time <- times
type <- types
pdata <- data.frame(type=type, teor_time=teor_time); rownames(pdata) <- sample_names
pdata <- new("AnnotatedDataFrame", pdata)
xprobe <- new("ExpressionSet", exprs=exprsx, featureData=fdata, phenoData=pdata)
save(xprobe,file=paste0(outdir, '/probeLevel.RData'))

exprsx <- as.matrix(probe_df[,redCols])
colnames(exprsx) <- sample_names
xprobe_red <- new("ExpressionSet", exprs=exprsx, featureData=fdata, phenoData=pdata)

write.csv(cbind(fdata@data, exprs(xprobe)), paste0(outdir, "/probeLevel_exp.csv"), row.names=F)
write.csv(cbind(fdata@data, exprs(xprobe_red)), paste0(outdir, "/probeLevel_redSignalexp.csv"), row.names=F)

#### Rename and Summarize ####

myRma <- function(x) {
  if (class(x)=='numeric') {
    ans <- x
  } else {
    ans <- medpolish(x,trace.iter=FALSE,na.rm=TRUE)
    ans <- ans$overall + ans$col
  }
  return(ans)
}

renameGenesAndSummarize <- function(genesToRename.sd,exprsx,geneid,summaryMethod=myRma,type) {

  if (type == "ratio"){
    xgene <- by(exprsx[,ratioCols],geneid,myRma)
  } else if (type == "red"){
    xgene <- by(exprsx[,redCols],geneid,myRma)
  }

  xgene <- do.call('rbind',xgene)

  mysd <- function(x) { ans <- ifelse(sum(!is.na(x))==1,0,sd(x,na.rm=TRUE)); return(ans) }
  sdgene <- aggregate(exprsx[, ratioCols],by=list(geneid),FUN=mysd)

  names(sdgene)[1] <- 'geneid'
  xgene <- data.frame(geneid=rownames(xgene),xgene); rownames(xgene) <- NULL

  fdata <- by(exprsx[,(signalCols+1):allcols],geneid,unique)

  genenames <- names(fdata)
  fdata <- do.call('rbind',fdata)

  fdata <- new("AnnotatedDataFrame", data.frame(fdata))
  rownames(fdata) <- as.character(xgene$geneid)

  exprsxgene <- as.matrix(xgene[,-1])
  rownames(exprsxgene) <- as.character(xgene$geneid);
  colnames(exprsxgene) <- sample_names
  eset <- new("ExpressionSet",exprs=exprsxgene, featureData=fdata, phenoData=pdata)
  return(list(eset=eset,sdgene=sdgene,fdata=fdata,geneid=geneid))
}

geneid <- probe_df$Gene_id
geneid <- as.character(geneid)
genesToRename.sd <- NA

tmp <- renameGenesAndSummarize(genesToRename.sd=genesToRename.sd,exprsx=probe_df,geneid=geneid,summaryMethod=myRma, type="ratio")
xgene <- tmp[['eset']]; sdgene <- tmp[['sdgene']]; fdata <- tmp[['fdata']]; geneid <- tmp[['geneid']]

tmp2 <- renameGenesAndSummarize(genesToRename.sd=genesToRename.sd,exprsx=probe_df,geneid=geneid,summaryMethod=myRma, type="red")
xgene_red <- tmp2[['eset']]; sdgene <- tmp2[['sdgene']]; fdata <- tmp2[['fdata']]; geneid <- tmp2[['geneid']]

#### Estimate times ####

print("Estimating times...")
bozdechPath <- paste0(infiles, 'bozdech_Hb3_clean2.csv')
LemieuxFunctionsPath <- paste0(infiles, 'lemieux_et_al_pipeline_functions.r')

getTimeEstimation <- function(x,dataPath,functionsPath,figuresPath,B=100) {
                                        #  x: the expressionSet for which we want to estimate times (our data).
                                        #  dataPath: path to data that will be used to estimate timepoints (from Bozdech et al)
                                        #  functionsPath: path to the script containing the functions from Lemieux's paper.
                                        #  figuresPath: where we want to save the output plots.
  source(functionsPath)
                                        #  z <- read.csv(dataPath, as.is = T,sep='\t')
  z <- read.csv(dataPath, as.is = T)
                                        #  colnames(z)[1] <- 'Name'
                                        #  oldTime <- as.numeric(as.character(pData(x)$time))
  oldTime <- as.numeric(teor_time)
  x <- exprs(x)
  x <- data.frame(Name=as.character(rownames(x)),x,stringsAsFactors=FALSE); rownames(x) <- NULL
  data <- sync_data(x, z)
  x <- data[[1]]
  z <- data[[2]]
  x <- ordinal(x, use.name = T)
  z <- ordinal(z, use.name = T)
                                        #  z.na <- cbind(z[,1:22], rep(NA, nrow(z)), z[,23:27], rep(NA, nrow(z)), z[,28:56])
  z.na <- cbind(z[,1:22], rep(NA, nrow(z)), z[,23:27], rep(NA, nrow(z)), z[,28:46])
  z <- t(apply(z.na, 1, smooth.missing))
  sigma.epsilon <- 789.9056
  z.smooth <- smooth.ref(z, method = "spline", spar = 0.5)
  z.smooth.hourly <- z.smooth[,ll.par$hourly]
                                        #  sigma.eta <- mean(sd(z[,11:ncol(z)] - z.smooth.hourly, na.rm = T), na.rm=T)
  sigma.eta <- mean(sd(z - z.smooth.hourly, na.rm = T), na.rm=T)
  new.sigma <- sqrt(sigma.eta^2 + sigma.epsilon^2)
  ll <- compute.ll(x = x, z = z.smooth, sigma = new.sigma, bootstrap = T, B = B, sample.rate = 0.50)
  myTimes <- mle(ll)
  png(file.path(figuresPath,'/Time_estim/defaultPlots1.png'))
  plot.ll(ll)
  dev.off()
  png(file.path(figuresPath,'/Time_estim/defaultPlots2.png'))
  plot.mle(ll)
  dev.off()
  png(file.path(figuresPath,'/Time_estim/ownPlots1.png'))
  plot(density(myTimes),main='Estimated times density')
  dev.off()
  png(file.path(figuresPath,'/Time_estim/ownPlots2.png'))
  plot(oldTime, as.numeric(myTimes),xlab='Old times',ylab='Estimated times',xlim=c(-5,50),ylim=c(-5,50))
  abline(0,1,col=2,lwd=2)
  abline(v=oldTime,lwd=0.5,lty=3)
  dev.off()
  return(myTimes)
}

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


estimatedTimes <- getTimeEstimation(xgene,bozdechPath,LemieuxFunctionsPath,file.path(figPath),B=100)
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

                                        # Save ExpressionSet at gene level
save(xgene,file=paste0(outdir, '/geneLevel.RData'))
save(xgene_red,file=paste0(outdir, '/geneLevel_redSignal.RData'))

                                        # Boxplot after summarization
pdf(file.path(figPath,'boxplot_afterSummarization.pdf'))
boxplot(exprs(xgene),main='summarization method: median poslish')
dev.off()

#### Save results in CSVs ####

write.csv(xgene@phenoData@data, file = paste0(outdir, "/experiment_data.csv"))
write.csv(cbind(xgene@featureData@data, exprs(xgene)), paste0(outdir, "/geneLevel_exp.csv"), row.names=F)
write.csv(cbind(xgene_red@featureData@data, exprs(xgene_red)), paste0(outdir, "/geneLevel_redSignal_exp.csv"), row.names=F)

#### FC: Probe Level ####

print("Calculating Fold-Changes...")
filter <- function(x,y) {x==y}
combs <- cross2(sample_names, sample_names, .filter = filter)

fc <- list()
for (i in combs) {
  fc[[paste0(i[[1]], "_",  i[[2]])]] <- exprs(xprobe)[,i[[1]]] - exprs(xprobe)[,i[[2]]]
}

probe_fc <- as.data.frame(fc)
probe_fc <- cbind(fData(xprobe), probe_fc)

write.csv(probe_fc, file = paste0(outdir, "/proveLevel_FC.csv"), row.names = F)

#### FC: Gene Level ####

filter <- function(x,y) {x==y}
combs <- cross2(sample_names, sample_names, .filter = filter)

fc <- list()
for (i in combs) {
  fc[[paste0(i[[1]], "_",  i[[2]])]] <- exprs(xgene)[,i[[1]]] - exprs(xgene)[,i[[2]]]
}

gene_fc <- as.data.frame(fc)
gene_fc <- cbind(fData(xgene), gene_fc)

write.csv(gene_fc, file = paste0(outdir, "/geneLevel_FC.csv"), row.names = F)

#### Areas: Functions ####

imputePoint <- function(xs, ys, tp){

  ## "xs" and "ys" must be two vectors of equal length
  ## with the corresponding y(expression) and x(timepoint)
  ## values that form the expression plot of interest (one gene).
  ## "tp" must be the timepoint to impute.
  ## If the timepoint to be imputed is already present, leave it as is.
  ## Returns NA if missing the previous or next tp

  if (tp %in% xs){

    idx <- which(xs == tp)
    imputed <- list(x=xs[idx], y=ys[idx])

  } else {

    before <- which(xs == max(xs[xs < tp]))
    after <- which(xs == min(xs[xs > tp]))

    if (is.na(ys[before]) | is.na(ys[after])){

      imputed <- list(x=tp, y=NA)

    } else {

      x <- c(xs[c(before, after)])
      y <- c(ys[c(before, after)])

      imputed <- approx(x, y, xout=tp)
    }
  }
  return(imputed)
}

computeArea <- function(eset){

  ## Takes an eset and computes areas.
  ## pData(eset) must contain a field named "time" with the time-points.
  ## pData(eset) must have a field named "type" with the grouping variable.

  ## Set needed variables
  types <- unique(phenoData(eset)$type)
  type <- phenoData(eset)$type
  times <- phenoData(eset)$time

  maxminTP <- max(sapply(types,function(x) min(times[type==x])))
  minmaxTP <- min(sapply(types,function(x) max(times[type==x])))
  mybreaks <- seq(maxminTP, minmaxTP, length.out=5)
  tp1 <- mybreaks[2]
  tp2 <- mybreaks[3]
  tp3 <- mybreaks[4]

  xsList <- c()
  for (type in types){
    xsList <- c(xsList, list(pData(eset)$time[phenoData(eset)$type == type]))
  }

  ## Main loop
  all_areas <- c()
  for (i in 1:dim(eset)[1]){

    gene <- fData(eset)$geneID[i]

    ysList <- c()
    for (type in types){
      ysList <- c(ysList, list(exprs(eset)[i, phenoData(eset)$type == type]))
    }

    ## Estimate points where needed
    dfs <- list()
    for (i in 1:length(xsList)){

      x <- unlist(xsList[[i]])
      y <- unlist(ysList[[i]])

      points <- as.data.frame(cbind(x, y))
      midpoints <- points[points$x > maxminTP &
                          points$x < minmaxTP, ]

      first <- imputePoint(x, y, maxminTP)
      last <- imputePoint(x, y, minmaxTP)
      p1 <- imputePoint(x, y, tp1)
      p2 <- imputePoint(x, y, tp2)
      p3 <- imputePoint(x, y, tp3)

      impPoints <- rbind(first, last, p1, p2, p3)
      allpoints <- rbind(midpoints, impPoints)
      allpoints$x <- as.numeric(allpoints$x)
      allpoints$y <- as.numeric(allpoints$y)

      ordered <- arrange(allpoints, allpoints$x)

      dfs[[i]] <- ordered
    }

    ## Calculate minY on estimated DFs
    minY <- min(sapply(dfs, function(df) min(df$y, na.rm = T)))

    rowareas <- c()
    for (df in dfs){

      df$y <- df$y - minY

      ## Whole polygon from expression data
      polDF <- rbind(df,
                     c(minmaxTP, 0),
                     c(maxminTP, 0),
                     c(df[1,]))

      ## Create Polygons
      leftHalf <- rbind(df[which(df$x <= tp2),],
                        c(tp2, 0),
                        c(maxminTP, 0),
                        c(df[1,]))

      rightHalf <- rbind(df[which(df$x >= tp2),],
                         c(minmaxTP, 0),
                         c(tp2, 0),
                         c(df[df$x == tp2,]))

      mid <- rbind(df[which(df$x >= tp1 & df$x <= tp3),],
                   c(tp3, 0),
                   c(tp1, 0),
                   c(df[df$x == tp1,]))

      sides <- rbind(df[which(df$x <= tp1),],
                     c(tp1, 0),
                     c(tp3, 0),
                     df[which(df$x >= tp3),],
                     c(minmaxTP, 0),
                     c(maxminTP, 0),
                     df[1,])

      pols <- list(leftHalf, rightHalf, mid, sides)

      calcArea <- function(x) {ifelse(any(is.na(x)), NA, Polygon(x)@area)}
      areas <- unlist(lapply(pols, function(x) calcArea(x)))

      rowareas <- c(rowareas, areas)

      ## Plot polygons (for debugging purposes)
      ##pol <- Polygon(polDF)
      ##ps = Polygons(list(pol),1)
      ##sps = SpatialPolygons(list(ps))
      ##plot(sps)
    }
    all_areas <- c(all_areas, list(rowareas))
  }

  ## Set row and col names for output
  areaDF <- do.call(rbind, all_areas)
  titles <- c("Left", "Right", "Middle", "Sides")

  cols <- c()
  for (i in types){
    for (t in titles){
      name <- paste0(i, "_", t)
      cols <- c(cols, name)
    }
  }
  colnames(areaDF)  <- cols
  rownames(areaDF) <- rownames(exprs(eset))
  return(areaDF)
}

#### Calls: Compute Areas ####

print("Computing Areas...")
areasDF <- computeArea(xgene)
xout <- cbind(fData(xgene), areasDF)

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

library(dplyr)

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
  mutate(MaxMax = pmax(abs(`A7-B11_MaxVal`), abs(`A7-E5_MaxVal`), abs(`B11-E5_MaxVal`))) %>%
  arrange(desc(abs(MaxMax))) %>%
  select(Gene_id, Name, Annot, contains('MaxVal'), contains('MaxTime'), contains('-'))

max_df_top <- max_df %>% dplyr::filter(abs(`A7-B11_MaxVal`) > 4 | abs(`A7-E5_MaxVal`) > 4 | abs(`B11-E5_MaxVal`) > 4)

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

ggsave(p, filename = paste0(figPath, "PCA.png"), device = "png", dpi = "retina")

p <- p + theme_classic()
p <- p + theme(text = element_text(size=20))
p

ggsave(p, filename = paste0(figPath, "PCA.svg"), device = "svg")

#### Expression Plots ####

print("Plotting Expression Plots...")
expressionPlot <- function(type){

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

  ## Set number of plots
  if (run_all_plots == "yes"){
    nplots = dim(df)[1]
  } else {
    nplots = 20
  }

  ## Main Loop
  for (i in 1:nplots){

    ## Set gene for title or gene and probe for probe-level plots.
    gn <- gsub("[/:;.]", "_", fData(df)$Gene_id[i])
    if (type %in% c("probe", "probe_red")){
      prb <- paste0("_", gsub("[/:;.]", "_" , fData(xprobe)$ProbeName[i]))
    } else {
      prb <- ""
    }
    title <- paste0(gn, prb)

    ## Set y-axis title depending on ratio/red_signal
    ytitle <- ifelse(type %in% c("gene_red", "probe_red"), 'log2(Cy5)', 'log2(Cy3/Cy5)')

    ## Plot
    graf <- melt(df[i,])
    graf["Type"] <- xgene@phenoData@data$type
    graf["Time"] <- xgene@phenoData@data$time
    p <- ggplot(graf, aes(x = Time, y = value, col = Type, group = Type))
    p <- p + geom_point(aes(color = Type, shape = Type)) + geom_line()
    p <- p + coord_cartesian(ylim = ylim)
    p <- p + ggtitle(title)
    p <- p + ylab(ytitle)
    ggsave(p, file=paste0(figPath, path, gn, prb, ".jpeg"),
           device = "jpeg", width = 14, height = 10, units = "cm")

  }

}

expressionPlot("gene")
expressionPlot("gene_red")
expressionPlot("probe")
expressionPlot("probe_red")

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
red_df <- 2**(red_df)
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

strains <- c('A7', 'E5', 'B11')
for (strain in strains) {varPlot(strain)}


vars_plot

xgene_red_tibble <- as.data.frame(exprs(xgene_red)) %>%
  tibble() %>%
  mutate(Gene_id = rownames(exprs(xgene_red))) %>%
  select(Gene_id, everything())

## By max col -> then percentile (old)
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

## perc_A7 <- get_red_percent('A7')
## perc_E5 <- get_red_percent('E5')
## perc_B11 <- get_red_percent('B11')

## red_percent <- tibble('Gene_id' = xgene_red_tibble$Gene_id,
##                       'A7' = perc_A7,
##                       'E5' = perc_E5,
##                       'B11' = perc_B11)


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

perc_A7 <- get_max_percent('A7')
perc_E5 <- get_max_percent('E5')
perc_B11 <- get_max_percent('B11')

red_percent <- tibble('Gene_id' = xgene_red_tibble$Gene_id,
                      'A7' = perc_A7,
                      'E5' = perc_E5,
                      'B11' = perc_B11)


#hist(red_percent$A7)
write_csv(red_percent, paste0(outdir, "/red_percentiles.csv"))

names(red_percent)[2:4] <- paste0('Red_Pcnt_', names(red_percent)[2:4])

max_tibble <- as_tibble(max_df)
max_tibble <- max_tibble %>%
  left_join(red_percent, by='Gene_id', suffix = c('', '_Pcnt'))

maxFC_pass_red <- max_tibble %>%
  dplyr::filter((`A7-B11_MaxVal` >= 2 & Red_Pcnt_A7 > 15) |
         (`A7-B11_MaxVal` <= -2 & Red_Pcnt_B11 > 15) |
         (`A7-E5_MaxVal` >= 2 & Red_Pcnt_A7 > 15) |
         (`A7-E5_MaxVal` <= -2 & Red_Pcnt_E5 > 15) |
         (`B11-E5_MaxVal` >= 2 & Red_Pcnt_B11 > 15) |
         (`B11-E5_MaxVal` <= -2 & Red_Pcnt_E5 > 15))

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

max_time <- tibble(Gene_id = names(maxtime), Max_Time = maxtime)
breaks_df <- tibble(Areas_Breaks = mybreaks)

tibble(Gene_id = names(maxtime), Max_Time = maxtime) %>%
  write_csv(paste0(outdir, 'new_arrays_maxtime.csv'))

tibble(Areas_Breaks = mybreaks) %>%
  write_csv(paste0(outdir, 'new_area_breaks.csv'))

## xgene_tibble <- as.data.frame(exprs(xgene)) %>%
##   tibble() %>%
##   mutate(Gene_id = rownames(exprs(xgene)))

## xgene_tibble %>%
##   select(contains('12B'))

## trans_df['Left_in_MaxTime'] <- sapply(maxtime, function(x) checkOverap(getInterval(x, width), left))

## #### Filter by Max-Time ####

## New approach
## Check which areas does maxtimepoint overlapp -> check if aAFC > th at this areas

point_overlap <- function(point, interval){
  point >= interval[1] & point <= interval[2]
}

areas_df <- as_tibble(allDifs) %>%
  mutate(Gene_id = rownames(allDifs)) %>%
  select(Gene_id, everything())

maxtimes_A7_B11 <- c()
maxtimes_A7_E5 <- c()
maxtimes_B11_E5 <- c()
gids <- c()
th <- 2
for (gid in max_tibble$Gene_id){

  ## Create time-regions
  breaks <- breaks_df$Areas_Breaks
  left <- c(breaks[1], breaks[3])
  right <- c(breaks[3], breaks[5])
  mid <- c(breaks[2], breaks[4])
  sides_l <- c(breaks[1], breaks[2])
  sides_r <- c(breaks[4], breaks[5])

  ## Get maxtime
  maxtime <- max_time %>%
    dplyr::filter(Gene_id == gid) %>%
    pull()
  if (is.na(maxtime)){
    maxtimes_A7_B11 <- c(maxtimes_A7_B11, NA)
    maxtimes_A7_E5 <- c(maxtimes_A7_E5, NA)
    maxtimes_B11_E5 <- c(maxtimes_B11_E5, NA)
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
      dplyr::filter(Gene_id == gid) %>%
      select(contains(names(areas[overlaps])))

    fc_A7_B11 <- aFCs %>%
      select(contains('A7-B11')) %>%
      replace(is.na(.), 0) %>%
      as.numeric()
    maxtime_FC_A7_B11 <- any(abs(fc_A7_B11) > th)

    fc_A7_E5 <- aFCs %>%
      select(contains('A7-E5')) %>%
      replace(is.na(.), 0) %>%
      as.numeric()
    maxtime_FC_A7_E5 <- any(abs(fc_A7_E5) > th)

    fc_B11_E5 <- aFCs %>%
      select(contains('B11-E5')) %>%
      replace(is.na(.), 0) %>%
      as.numeric()
    maxtime_FC_B11_E5 <- any(abs(fc_B11_E5) > th)

    maxtimes_A7_B11 <- c(maxtimes_A7_B11, maxtime_FC_A7_B11)
    maxtimes_A7_E5 <- c(maxtimes_A7_E5, maxtime_FC_A7_E5)
    maxtimes_B11_E5 <- c(maxtimes_B11_E5, maxtime_FC_B11_E5)
    gids <- c(gids, gid)
  }
}
maxtime_aAFC_df <- tibble(
  Gene_id = gids,
  MaxTime_Filter_A7_B11 = maxtimes_A7_B11,
  MaxTime_Filter_A7_E5 = maxtimes_A7_E5,
  MaxTime_Filter_B11_E5 = maxtimes_B11_E5
  )

max_tibble <- max_tibble %>%
  left_join(maxtime_aAFC_df)

f_path <- '/mnt/Disc4T/Projects/PhD_Project/Paper/Paper_Analysis/Data_Files/Duplication_Deletion_Regions_Mean_Separate_DuplDel/Crossed_with_genes/'

file_list <- c(  'A7K9_in_sort_q5_noDup_rpkm_normInput_bs10_smth200_bymean_dupl_del_fact_1.75up_0.1dw_minlen500_mergelen_200_filtered_genes.tsv',
  'E5K9_in_sort_q5_noDup_rpkm_normInput_bs10_smth200_bymean_dupl_del_fact_1.75up_0.1dw_minlen500_mergelen_200_filtered_genes.tsv',
  'B11_in_sort_q5_noDup_rpkm_normInput_bs10_smth200_bymean_dupl_del_fact_1.75up_0.1dw_minlen500_mergelen_200_filtered_genes.tsv'
)

dupl_dl_A7 <- read_tsv(paste0(f_path, file_list[1]), col_names = F) %>%
  select(X1) %>% pull()

dupl_dl_E5 <- read_tsv(paste0(f_path, file_list[2]), col_names = F) %>%
  select(X1) %>% pull()

dupl_dl_B11 <- read_tsv(paste0(f_path, file_list[3]), col_names = F) %>%
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

## Set thresholds
red_th <- 15

## A7 vs B11

final_df %>%
  select(
    Gene_id,
    `A7-B11_MaxVal`,
    `A7-B11_MaxTime`,
    Red_Pcnt_A7,
    Red_Pcnt_B11,
    MaxTime_Filter_A7_B11
  ) %>%
  mutate(PassRed = ifelse(
           `A7-B11_MaxVal` >= 0,
           Red_Pcnt_A7 >= red_th,
           Red_Pcnt_B11 >= red_th
         )) %>%
  rename(PassMaxtime = MaxTime_Filter_A7_B11) %>%
  mutate(PassDuplDel = !Gene_id %in% dupl_dl_A7 & !Gene_id %in% dupl_dl_B11) %>%
  rowwise() %>%
  mutate(PassAll = all(c(PassRed, PassMaxtime, PassDuplDel))) %>%
  write_tsv(paste0(outdir, 'A7_B11_final_df.tsv'))

## A7 vs E5

final_df %>%
  select(
    Gene_id,
    `A7-E5_MaxVal`,
    `A7-E5_MaxTime`,
    Red_Pcnt_A7,
    Red_Pcnt_E5,
    MaxTime_Filter_A7_E5
  ) %>%
  mutate(PassRed = ifelse(
           `A7-E5_MaxVal` >= 0,
           Red_Pcnt_A7 >= red_th,
           Red_Pcnt_E5 >= red_th
         )) %>%
  rename(PassMaxtime = MaxTime_Filter_A7_E5) %>%
  mutate(PassDuplDel = !Gene_id %in% dupl_dl_A7 & !Gene_id %in% dupl_dl_E5) %>%
  rowwise() %>%
  mutate(PassAll = all(c(PassRed, PassMaxtime, PassDuplDel))) %>%
  write_tsv(paste0(outdir, 'A7_E5_final_df.tsv'))

## B11 vs E5

final_df %>%
  select(
    Gene_id,
    `B11-E5_MaxVal`,
    `B11-E5_MaxTime`,
    Red_Pcnt_B11,
    Red_Pcnt_E5,
    MaxTime_Filter_B11_E5
  ) %>%
  mutate(PassRed = ifelse(
           `B11-E5_MaxVal` >= 0,
           Red_Pcnt_B11 >= red_th,
           Red_Pcnt_E5 >= red_th
         )) %>%
  rename(PassMaxtime = MaxTime_Filter_B11_E5) %>%
  mutate(PassDuplDel = !Gene_id %in% dupl_dl_B11 & !Gene_id %in% dupl_dl_E5) %>%
  rowwise() %>%
  mutate(PassAll = all(c(PassRed, PassMaxtime, PassDuplDel))) %>%
  write_tsv(paste0(outdir, 'B11_E5_final_df.tsv'))


###################################################


## A7 vs B11

difs_A7_B11 <- final_df %>%
  dplyr::filter(
  ((`A7-B11_MaxVal` >= 2 & Red_Pcnt_A7 > 15) |
   (`A7-B11_MaxVal` <= -2 & Red_Pcnt_B11 > 15)) &
  MaxTime_Filter_A7_B11
  ) %>%
  mutate(Dupl_Del = Gene_id %in% dupl_dl_A7 | Gene_id %in% dupl_dl_B11) %>%
  dplyr::filter(!Not_Plasmodium & !Is_tRNA) %>%
  select(Gene_id, Name, Annot,
         `A7-B11_MaxVal`,
         Red_Pcnt_A7, Red_Pcnt_B11,
         MaxTime_Filter_A7_B11,
         Variant, Gam_specific, Dupl_Del)

write_csv(difs_A7_B11, paste0(outdir, 'A7_vs_B11_log2FC2_red15_maxtime.csv'))


## A7 vs E5

difs_A7_E5 <- final_df %>%
  dplyr::filter(
  ((`A7-E5_MaxVal` >= 2 & Red_Pcnt_A7 > 15) |
   (`A7-E5_MaxVal` <= -2 & Red_Pcnt_E5 > 15)) &
  MaxTime_Filter_A7_E5
  ) %>%
  mutate(Dupl_Del = Gene_id %in% dupl_dl_A7 | Gene_id %in% dupl_dl_E5) %>%
  dplyr::filter(!Not_Plasmodium & !Is_tRNA) %>%
  select(Gene_id, Name, Annot,
         `A7-E5_MaxVal`,
         Red_Pcnt_A7, Red_Pcnt_E5,
         MaxTime_Filter_A7_E5,
         Variant, Gam_specific, Dupl_Del)

write_csv(difs_A7_E5, paste0(outdir, 'A7_vs_E5_log2FC2_red15_maxtime.csv'))


## B11 vs E5

difs_B11_E5 <- final_df %>%
  dplyr::filter(
  ((`B11-E5_MaxVal` >= 2 & Red_Pcnt_B11 > 15) |
   (`B11-E5_MaxVal` <= -2 & Red_Pcnt_E5 > 15)) &
  MaxTime_Filter_B11_E5
  ) %>%
  dplyr::filter(!Not_Plasmodium & !Is_tRNA) %>%
  mutate(Dupl_Del = Gene_id %in% dupl_dl_B11 | Gene_id %in% dupl_dl_E5) %>%
  select(Gene_id, Name, Annot,
         `B11-E5_MaxVal`,
         Red_Pcnt_B11, Red_Pcnt_E5,
         MaxTime_Filter_B11_E5,
         Variant, Gam_specific, Dupl_Del)

write_csv(difs_B11_E5, paste0(outdir, 'B11_vs_E5_log2FC2_red15_maxtime.csv'))

write_csv(final_df, paste0(outdir, 'new_arrays_final_df.csv'))

library(eulerr)

A <- difs_A7_E5$Gene_id
B <- difs_A7_B11$Gene_id
C <- difs_B11_E5$Gene_id

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
          labels = list(labels=c("A7 vs E5", "A7 vs B11", "B11 vs E5"), fontsize = 7))

ggsave(d, filename = paste0(figPath, "Difs_Venn.pdf"), device = "pdf")

plot(d)
print(fit)

array_list_tibble <- lapply(array_list, as.tibble)
raw_array <- array_list_tibble[[1]] %>%
  select(ProbeName, FeatureNum)

geo_df <- tibble(probe_df_geo)
geo_df['FeatureNum'] <- raw_array$FeatureNum
geo_df <- geo_df %>%
  mutate(Drop_Gene_id = !(Gene_id %in% final_df$Gene_id)) %>%
  mutate(Excluded_from_Analysis = Drop_Mulimapper | Drop_LowExpress | Drop_Gene_id) %>%
  select(
    ProbeName, FeatureNum, Gene_id,
    contains('Drop'), Excluded_from_Analysis,
    contains('Log')
  ) %>%
  write_tsv(paste0(outdir, 'geo_processed_data.tsv'))

#### Save environtment ####
#save.image(file = "array_A7E5B11_work_space.RData")
##setwd('/mnt/Disc4T/Projects/PhD_Project/Microarrays/New_Old_separate_approach/New_Arrays/')
##load('array_A7E5B11_work_space.RData')
head(final_df)
