####################################################################################################################
#
# Script: other_preprocess.R
# Project: An integrated multi-omic single cell atlas of human B cell identity
# Author: David Glass
# Date: 5-29-20
#
# Purpose: Preprocess fcs files from mass cytometry datasets and generates a csv file per dataset
#
# FCS files:
#   https://flowrepository.org/id/FR-FCM-Z2MC
#
# Pseudocode:
# For each dataset:
#   Read in all fcs files
#   Combine data into a single data.table with metadata
#   Remove non-protein channels
#   Asinh transform
#   Normalize if needed
#   Scale individual channels
#   Normalize light chain
#   Assign cluster
#   Assign isotype
#   Write csv
#
# Instructions:
# Install all required packages (see LIBRARIES and/or Versions below)
# Download the dataset linked above
# Put all of the fcs files into a unique directory
# Put all metadata csv files in the same directory:
#   figure_4_metadata.csv
#   figure_5_biosynthesis_metadata.csv
#   figure_5_metabolism_metadata.csv
#   figure_5_signaling_metadata.csv
#   figure_6_metadata.csv
# In the USER INPUTS section below, assign path variable to the path of the fcs directory
# In the MAIN section at the bottom of the script, if you do not wish to generate a csv of processed data,
#   set csv=FALSE in the process...Data functions. By default no data.table is returned, that can be changed
#   by setting return.dt to TRUE
# If you do not wish to process each dataset, delete/comment out irrelevant process...Data function call
# Run script
# The csv file or data.table can be used in figure_X.R to recreate figures 4-6
#
### NOTE: Processing often takes a few minutes to run
#
# Versions:
# R 3.6.3
# RStudio 1.2.5042
# flowCore_1.52.1 
# FlowSOM_1.18.0
# caret_6.0-86
# FNN_1.1.3
# preprocessCore_1.48.0
# dplyr_1.0.0
# data.table_1.12.8
#
#######################################################################################################################



##### USER INPUTS #####

### Path to folder containing screen fcs files and screen_metadata.csv
path <- "~/example2/"



###### LIBRARIES ######

require(flowCore)
require(FlowSOM)
require(caret)
require(FNN)
require(preprocessCore)
require(dplyr)
require(data.table)



##### FUNCTIONS #####

readFile <- function(p=path, files=meta.dat$filename) {
  # takes in a path with fcs files and returns a list of data.tables of the expression matrix
  # Inputs:
  #   p - character vector with directory storing fcs files
  #   files - character vector of filenames
  # Outputs:
  #   frames - a list of data.tables
  print("Reading files")
  frames <- setNames(vector("list", length(files)), files)
  for (i in 1:length(files)) {
    fcs <- read.FCS(paste0(path, files[i]), transformation = FALSE, emptyValue = FALSE)
    frames[[i]] <- data.table(exprs(fcs)) %>%
      setnames(pData(parameters(fcs))$desc)
  }
  return(frames)
}


combineFiles <- function(frames, md=meta.dat, dump=dump.channels, fa=factors, gated=F) {
  # Combines a list of data.tables into a single data.table with factor columns added
  # Inputs:
  #   frames - a list of data.tables
  #   md - data.table of metadata
  #   dump - vector of dump channel names
  #   fa - vector of factors
  #   gated - logical, does the dataset have ungated total B cells and gated subpopulations
  # Outputs:
  #   frame - a single data.table
  print("Combining files")
  all.cols <- unique(unlist(lapply(frames, colnames))) %>%
    .[!. %in% dump] %>%
    c(fa)
  frame <- data.table(matrix(nrow=0, ncol=length(all.cols))) %>%
    setnames(all.cols)
  for (f in md$filename) {
    frame <- frames[[f]] %>%
      cbind(md[filename==f]) %>%
      .[, all.cols, with=F] %>%
      rbind(frame)
  }
  if (gated) {
    frame <- rbind(frame[gate!="Ungated"], frame[gate=="Ungated"]) %>%
      unique(by=all.cols[!all.cols %in% fa])
  }
  return(frame)
}


asinTransform <- function(dt, exceptions=factors) {
  # asinh transforms data
  # Inputs:
  #   dt - data.table
  #   exceptions - character vector of channels not to transform
  # Outputs:
  #   dt - data.table
  print("Asinh transforming")
  to.transform <- setdiff(colnames(dt), exceptions)
  dt[, (to.transform) := asinh(dt[, to.transform, with=F]/5)]
  return(dt)
}


getPeak <- function(x) {
  # returns value of highest density peak in a vector
  # Inputs:
  #   x - a numeric vector
  # Outputs:
  #   Value of the density peak in x
  d <- density(x)
  return(d$x[which.max(d$y)])
  }


normalizeSurfaceSamples <- function(dt, reference=929, others=c(930, 931)) {
  # Normalizes all cells for each subject
  # Inputs:
  #   dt - data.tabs
  #   reference - ID of reference donor
  #   others - ID of other donors
  # Outputs:
  #   dt - data.table
  print("Normalizing files")
  # Manually identified value below dense positive peak.
  peak.mins <- list(CD35=1, CD21=0.2, IgM=2.2, IgD=3, IgA=3.5, IgG=3, CD185=2, HLA_ABC=1, CD31=0.5,
                    CD45RB=3, CD183=2, CD81=0.5, CD29=0.5, CD40=1, CD39=1, CD44=2, CD27=2, CD62L=2, CD99=2,
                    CD38=0.5, CD82=1, CD20=3, CD32=0.5, IgK=2, IgL=3, CD72=1, CD24=0.5, CD19=2.5, CD73=1.5, CD184=2,
                    CD305=1.5, CD22=1, CD45=3)
  
  # not normalizing CD11c, CD5, CD9, CD79b, CD1c, CD95 because peaks unclear
  
  for (marker in names(peak.mins)) {
    marker.ref <- dt[subject==reference & eval(parse(text=marker))>peak.mins[[marker]], marker, with=F] %>%
      as.matrix() %>%
      getPeak()
    for (id in others) {
      subject.ref <- dt[subject==id & eval(parse(text=marker))>peak.mins[[marker]], marker, with=F] %>%
        as.matrix() %>%
        getPeak()
      adjustment <- subject.ref / marker.ref
      dt[subject==id, eval(marker) := dt[subject==id, eval(marker), with=F] / adjustment]
    }
  }
  return(dt)
}


scaleData <- function(dt, exceptions=factors, percentile=0.999) {
  # scales each channel to percentile
  # Inputs:
  #   dt - data.table
  #   exceptions - character vector of channels not to scale
  #   percentile - percentile to scale each channel to
  # Outputs:
  #   dt - data.table
  print("Scaling")
  findRefs <- function(x) quantile(x, probs=c(percentile), na.rm=T)
  channels <- setdiff(colnames(dt), exceptions)
  refs <- apply(dt[, (channels), with=F], 2, findRefs)
  for (channel in channels) {
    dt[, eval(channel) := (dt[, eval(channel), with=F] / refs[channel])]
  }
  return(dt)
}


normalizeLight <- function(dt, threshold=0.5) {
  # Normalizes the light chains and creates light chain channel
  # Inputs:
  #   dt - data.table
  #   threshold - approximate threshold of positivity
  # Outputs:
  #   dt - data.table
  print("Normalizing light chain")
  igk <- dt[IgK>threshold, IgK] %>%
    as.matrix() %>%
    getPeak()
  igl <- dt[IgL>threshold, IgL] %>%
    as.matrix() %>%
    getPeak()
  ref <- igk / igl
  dt[, IgL:=IgL*ref]
  dt[, `:=`(light=pmax(IgK, IgL), IgK=NULL, IgL=NULL)]
  return(dt)
}


somCluster <- function(dtable, channels, ...) {
  # Clusters with a SOM
  # Inputs:
  #   dtable - data.table with all of the data to be clustered
  #   channels - vector of channels to use in clustering
  #   ... other arguments passed on to SOM()
  # Outputs:
  #   cluster assignment column
  RNGkind(sample.kind = "Rounding") # Backwards compatible seed selection
  set.seed(666)
  som.out <- SOM(as.matrix(dtable[, channels, with=F]), ...)
  return(as.factor(som.out$mapping[,1]))
}


metaClustering_hclust <- function(data, nClus){
  # Modified from https://github.com/SofieVG/FlowSOM/blob/master/R/4_metaClustering.R
  # Performs hierachical metaclustering on a matrix
  # Inputs:
  #   data - a matrix of medians
  #   nClus - number of desired clusters
  # Outputs:
  #   Vector with group memberships
  d <- stats::dist(data, method = "minkowski")
  fit <- stats::hclust(d, method="ward.D2")
  stats::cutree(fit, k=nClus)
}


clusterSurface <- function(dt, fa=factors) {
  # Assign metacluster to surface dataset
  # Inputs:
  #   dt - data.table
  #   fa - vector of factor column names
  # Outputs:
  #   all.dat - data.table with meta and cluster assignment
  print("Clustering")
  all.markers <- colnames(dt)[!colnames(dt) %in% fa]
  markers <- all.markers %>% setdiff(c("IgM", "IgD", "IgA", "IgG"))
  dt[, cluster:=somCluster(dtable=dt, channels=all.markers, xdim=13, ydim=13)]

  split.markers <- c("CD45RB", "CD27", "CD305", "CD44", "CD11c")
  split.medians <- dt[, lapply(.SD, median), .SDcols=all.markers, by=.(cluster)]
  split.medians[, meta:=factor(metaClustering_hclust(data=as.matrix(split.medians[, split.markers, with=F]), nClus=2))]
  split.medians[, meta:=recode(meta, "1"="Inexperienced", "2"="Experienced")]
  dt[, meta:=factor(split.medians[order(cluster), meta][dt$cluster])]

  i.dat <- dt[meta=="Inexperienced"] %>%
    .[, cluster:=droplevels(cluster)]
  i.markers <- c("CD38", "CD73", "CD79b")
  i.medians <- i.dat[, lapply(.SD, median), .SDcols=all.markers, by=.(cluster)]
  i.medians[, meta:=factor(metaClustering_hclust(data=as.matrix(i.medians[, i.markers, with=F]), nClus=3))]
  i.medians[, meta:=recode(meta, "1"="Transitional", "2"="73+ Naïve", "3"="73- Naïve")] %>%
    .[, meta:=factor(meta, levels=levels(meta)[c(3, 1, 2)])]
  for (clust in levels(i.dat$cluster)) i.dat[cluster==clust, meta:=i.medians[cluster==clust, meta]]

  e.dat <- dt[meta=="Experienced"] %>%
    .[, cluster:=droplevels(cluster)]
  pb.markers <- c("CD20", "CD268", "CD95", "CD11c")
  e.medians <- e.dat[, lapply(.SD, median), .SDcols=all.markers, by=.(cluster)]
  e.medians[, meta:=factor(metaClustering_hclust(data=as.matrix(e.medians[, pb.markers, with=F]), nClus=4))]
  e.medians[, meta:=recode(meta, "1"="Memory", "4"="Plasma", "2"="95+ Memory", "3"="19++ 11c+ Memory")]
  for (clust in levels(e.dat$cluster)) e.dat[cluster==clust, meta:=e.medians[cluster==clust, meta]]

  a.dat <- e.dat[meta!="Memory"] %>%
    .[, cluster:=droplevels(cluster)]

  m.dat <- e.dat[meta=="Memory"] %>%
    .[, cluster:=droplevels(cluster)]
  m.medians <- m.dat[, lapply(.SD, median), .SDcols=all.markers, by=.(cluster)]
  m.markers <- c("CD45RB")
  m.medians[, meta:=factor(metaClustering_hclust(data=as.matrix(m.medians[, m.markers, with=F]), nClus=2))]
  m.medians[, meta:=recode(meta, "1"="Memory", "2"="RB- Memory")] %>%
    .[, meta:=factor(meta, levels=levels(meta)[c(2, 1)])]
  for (clust in levels(m.dat$cluster)) m.dat[cluster==clust, meta:=m.medians[cluster==clust, meta]]

  no.rb.dat <- m.dat[meta!="Memory"] %>%
    .[, cluster:=droplevels(cluster)]

  rb.dat <-  m.dat[meta=="Memory"] %>%
    .[, cluster:=droplevels(cluster)]
  rb.medians <- rb.dat[, lapply(.SD, median), .SDcols=all.markers, by=.(cluster)]
  rb.markers <- c("CD27", "CD73")
  rb.medians[, meta:=factor(metaClustering_hclust(data=as.matrix(rb.medians[, rb.markers, with=F]), nClus=3))]
  rb.medians[, meta:=recode(meta, "1"="27- Memory", "3"="RB+ 27+ 73- Memory", "2"="RB+ 27+ 73+ Memory")] %>%
    .[, meta:=factor(meta, levels=levels(meta)[c(2, 1, 3)])]
  for (clust in levels(rb.dat$cluster)) rb.dat[cluster==clust, meta:=rb.medians[cluster==clust, meta]]

  all.dat <- rbind(i.dat, rb.dat, no.rb.dat, a.dat)
  all.dat[, meta:=droplevels(meta)]
  all.dat[, meta:=factor(meta, levels(meta)[c(2, 1, 3, 9, 10, 8, 7, 4, 5, 6)])]
  return(all.dat)
}


assignIgSurface <- function(dt) {
  # Assigns isotype class label
  # Inputs:
  #   dt - data.table
  # Outputs:
  #   dt - data.table
  print("Identifying B cell isotype")
  dt[IgA>0.42 & IgA>IgD, isotype:="IgA"]
  dt[IgD>0.4 & IgD>IgA, isotype:="IgD"]
  dt[IgG>0.54 & IgG>IgM, isotype:="IgG"]
  dt[IgM>0.46 & IgM>IgG & IgD<0.4, isotype:="IgM"]
  dt[IgM>0.46 & IgM>IgG & IgD>0.4 & IgD>IgA, isotype:="IgMD"]
  dt[is.na(isotype), isotype:="ND"]
  dt[, isotype:=factor(isotype, levels=c("IgD", "IgMD", "IgM", "IgG", "IgA", "ND"))]
  return(dt)
}


clusterMetabolism <- function(dt) {
  # Assign metacluster to metabolism dataset
  # Inputs:
  #   dt - data.table
  # Outputs:
  #   dt - data.table
  print("Clustering")
  metabolism <- c("MCT1", "AHR", "HK2", "AMPK_p", "ATPA5", "PPARa", "GAPDH", "PFK2", "CytC", "LDHA", "IDH1")
  markers <- c("CD11c", "CD184", "CD19", "CD20", "CD22", "CD23", "CD24", "CD27", "CD305",
                      "CD32", "CD38", "CD39", "CD45", "CD45RB", "CD69", "CD72", "CD73", "CD82", "CD9",
                      "CD95", "IgA", "IgD", "IgG", "IgM", "light")
  clusters <- c("Transitional", "73- Naïve", "73+ Naïve", "27- Memory", "RB+ 27+ 73- Memory",
                "RB+ 27+ 73+ Memory", "RB- Memory", "95+ Memory", "19++ 11c+ Memory", "Plasma")
  all.markers <- c(markers, metabolism)
  
  dt[, cluster:=somCluster(dtable=dt, channels=markers)]
  medians <-dt[, lapply(.SD, median), .SDcols=all.markers, by=.(cluster)]
  medians[CD27>0.75, meta:="Plasma"]
  medians[is.na(meta) & CD11c>0.6, meta:="19++ 11c+ Memory"]
  medians[is.na(meta) & (CD45RB>0.5 | CD27>0.1), meta:="Memory"]
  medians[is.na(meta), meta:="Inexperienced"]
  medians[meta=="Inexperienced" & CD38>0.5, meta:="Transitional"]
  medians[meta=="Inexperienced" & CD73>0.3, meta:="73+ Naïve"]
  medians[meta=="Inexperienced", meta:="73- Naïve"]
  medians[meta=="Memory" & CD95>0.2, meta:="95+ Memory"]
  medians[meta=="Memory" & CD45RB<0.5, meta:="RB- Memory"]
  medians[meta=="Memory" & CD27<0.25, meta:="27- Memory"]
  medians[meta=="Memory" & CD73>0.2, meta:="RB+ 27+ 73+ Memory"]
  medians[meta=="Memory", meta:="RB+ 27+ 73- Memory"]

  medians[, meta:=factor(meta, levels=clusters)]
  for (clust in medians$cluster) dt[cluster==clust, meta:=medians[cluster==clust, meta]]
  return(dt)
}


assignIgMetabolism <- function(dt) {
  # Assigns isotype class label
  # Inputs:
  #   dt - data.table
  # Outputs:
  #   dt - data.table
  print("Identifying B cell isotype")
  dt[IgA>0.4 & IgA>IgD, isotype:="IgA"]
  dt[IgD>0.25 & IgD>IgA, isotype:="IgD"]
  dt[IgM>0.3 & IgM>IgG & IgD<0.25, isotype:="IgM"]
  dt[IgM>0.3 & IgM>IgG & IgD>0.25 & IgD>IgA, isotype:="IgMD"]
  dt[is.na(isotype) & IgG>0.15, isotype:="IgG"]
  dt[is.na(isotype), isotype:="ND"]
  dt[, isotype:=factor(isotype, levels=c("IgD", "IgMD", "IgM", "IgG", "IgA", "ND"))]
  return(dt)
}


clusterBiosynthesis <- function(dt) {
  # Assign metacluster to biosynthesis dataset
  # Inputs:
  #   dt - data.table
  # Outputs:
  #   dt - data.table
  print("Clustering")
  markers <- c("IgM", "IgD", "IgA", "IgG", "CD9", "CD39", "CD45RB", "CD72", "CD45", "CD19", "CD23",
               "CD38", "CD82", "CD32", "CD11c", "CD27", "CD95", "CD24", "CD73", "CD184", "CD305",
               "CD20", "CD22", "light")
  responses <- c("BrdU", "puromycin")
  all.markers <- c(markers, responses)
  clusters <- c("Transitional", "73- Naïve", "73+ Naïve", "27- Memory", "RB+ 27+ 73- Memory",
                "RB+ 27+ 73+ Memory", "RB- Memory", "95+ Memory", "19++ 11c+ Memory", "Plasma")
  
  dt[, cluster:=somCluster(dtable=dt, channels=markers)]
  medians <-dt[, lapply(.SD, median), .SDcols=all.markers, by=.(cluster)]
  medians[CD27>0.75, meta:="Plasma"]
  medians[is.na(meta) & CD11c>0.55, meta:="19++ 11c+ Memory"]
  medians[is.na(meta) & (CD45RB>0.7 | CD27>0.1), meta:="Memory"]
  medians[is.na(meta), meta:="Inexperienced"]
  medians[meta=="Inexperienced" & CD38>0.45, meta:="Transitional"]
  medians[meta=="Inexperienced" & CD73>0.3, meta:="73+ Naïve"]
  medians[meta=="Inexperienced", meta:="73- Naïve"]
  medians[meta=="Memory" & CD95>0.4, meta:="95+ Memory"]
  medians[meta=="Memory" & CD45RB<0.6, meta:="RB- Memory"]
  medians[meta=="Memory" & CD27<0.33, meta:="27- Memory"]
  medians[meta=="Memory" & CD73>0.33, meta:="RB+ 27+ 73+ Memory"]
  medians[meta=="Memory", meta:="RB+ 27+ 73- Memory"]
  medians[, meta:=factor(meta, levels=clusters)]
  for (clust in medians$cluster) dt[cluster==clust, meta:=medians[cluster==clust, meta]]
  return(dt)
}


assignIgBiosynthesis <- function(dt) {
  # Assigns isotype class label
  # Inputs:
  #   dt - data.table
  # Outputs:
  #   dt - data.table
  print("Identifying B cell isotype")
  dt[IgA>0.4 & IgA>IgD, isotype:="IgA"]
  dt[IgD>0.35 & IgD>IgA, isotype:="IgD"]
  dt[IgG>0.37 & IgG>IgM, isotype:="IgG"]
  dt[IgM>0.3 & IgM>IgG & IgD<0.35, isotype:="IgM"]
  dt[IgM>0.3 & IgM>IgG & IgD>0.35 & IgD>IgA, isotype:="IgMD"]
  dt[is.na(isotype), isotype:="ND"]
  dt[, isotype:=factor(isotype, levels=c("IgD", "IgMD", "IgM", "IgG", "IgA", "ND"))]
  return(dt)
}


clusterSignaling <- function(dt) {
  # Assign metacluster to signaling dataset
  # Inputs:
  #   dt - data.table
  # Outputs:
  #   dt - data.table
  print("Clustering")
  clusters <- c("Transitional", "73- Naïve", "73+ Naïve", "27- Memory", "RB+ 27+ 73- Memory",
                "RB+ 27+ 73+ Memory", "RB- Memory", "95+ Memory", "19++ 11c+ Memory", "Plasma")
  
  signaling <- c("p_PLCg2__pY759_", "p_ZAP70_Syk__Y319_Y352_", "p_p38__T180_Y182_", "light")
  ig <- c("IgA", "IgD", "IgG", "IgM", "light")
  non.signaling <- c("CD11c", "CD184", "CD19", "CD20", "CD22", "CD23", "CD24", "CD27", "CD305", "CD32",
                     "CD38", "CD39", "CD45", "CD45RB", "CD72", "CD73", "CD82", "CD9", "CD95")
  markers <- c(non.signaling, ig)
  all.markers <- unique(c(signaling, ig, non.signaling))
  
  dat.0 <- copy(dt[dose==0])
  dat.0[, cluster:=somCluster(dtable=dat.0, channels=markers)]
  medians.0 <-dat.0[, lapply(.SD, median), .SDcols=markers, by=.(cluster)]

  medians.0[CD27>0.75, meta:="Plasma"]
  medians.0[is.na(meta) & CD11c>0.63, meta:="19++ 11c+ Memory"]
  medians.0[is.na(meta) & (CD45RB>0.6 | CD27>0.05), meta:="Memory"]
  medians.0[is.na(meta), meta:="Inexperienced"]
  medians.0[meta=="Inexperienced" & CD24>0.55, meta:="Transitional"]
  medians.0[meta=="Inexperienced" & CD73>0.33, meta:="73+ Naïve"]
  medians.0[meta=="Inexperienced", meta:="73- Naïve"]
  medians.0[meta=="Memory" & CD95>0.2, meta:="95+ Memory"]
  medians.0[meta=="Memory" & CD45RB<0.6, meta:="RB- Memory"]
  medians.0[meta=="Memory" & CD27<0.14, meta:="27- Memory"]
  medians.0[meta=="Memory" & CD73>0.3, meta:="RB+ 27+ 73+ Memory"]
  medians.0[meta=="Memory", meta:="RB+ 27+ 73- Memory"]
  medians.0[, meta:=factor(meta, levels=clusters)]
  for (clust in medians.0$cluster) dat.0[cluster==clust, meta:=medians.0[cluster==clust, meta]]

  # 10-fold classification for each k - best average
  # k 2-20
  # RNGkind(sample.kind = "Rounding") # Backwards compatible seed selection
  # set.seed(666)
  # idx <- createFolds(dat.0$meta, k=10)
  # errors <- matrix(nrow=10, ncol=20)
  # for (i in 1:10) {
  #   train <- dat.0[!idx[[i]]]
  #   test <- dat.0[idx[[i]]]
  #   for (k in 1:20) {
  #     pred <- knn(train=train[, non.signaling, with=F], test=test[, non.signaling, with=F], cl=train$meta, k=k)
  #     errors[i, k] <- mean(pred!=test$meta)
  #   }
  # }
  # k <- apply(errors, 2, mean) %>%
  #   which.min()
  k <- 16L # as selected through commented out cross-validation code

  train <- as.matrix(dat.0[, non.signaling, with=F])
  test <- as.matrix(dt[dose %in% c(0.5, 1), non.signaling, with=F])
  cl <- dat.0$meta
  classes <- FNN::knn(train=train, test=test, cl=cl, k=k)

  dt[dose %in% c(0.5, 1), meta:=classes]
  dat.0[, cluster:=NULL]

  dt <- rbind(dat.0, dt[dose %in% c(0.5, 1)])
  return(dt)
}


batchCorrectTissue <- function(dt, fa=factors, percentile=c(0.999)) {
  # Batch corrects run 1 and 2 from tissue run.
  #   Evaluates 99.9th percentile of each channel from PB sample from run 1.
  #   Does the same for same number of PB cells from run 2.
  #   Identifies a multiplier for each channel to normalize these values (always normalizing downward)
  #   Applies the same multiplier to all the data from that run
  # Inputs:
  #   dt - data.table
  #   fa - character vector of channels not to transform
  #   percentile - percentile for normaliation
  # Outputs:
  #   dt - data.table
  print("Batch correcting")
  n.sample <- nrow(dt[run==1 & tissue=="PB"])
  channels <- setdiff(colnames(dt), fa)
  RNGkind(sample.kind = "Rounding") # Backwards compatible seed selection
  set.seed(666)
  temp <- dt[run==2 & tissue=="PB"] %>% 
    .[, .SD[sample(.N, n.sample)]] %>%
    rbind(dt[run==1 & tissue=="PB"])
  refs <- temp[, lapply(.SD, quantile, probs=percentile), .SDcols=channels, by=run]
  
  for (channel in channels) {
    multiplier <- refs[run==1, eval(parse(text=channel))] / refs[run==2, eval(parse(text=channel))]
    if (multiplier > 1) dt[run==1, eval(channel):=dt[run==1, eval(channel), with=F]  / multiplier]
    else dt[run==2, eval(channel):=dt[run==2, eval(channel), with=F]  * multiplier]
  }
  return(dt)
}


assignIgTissue <- function(dt) {
  # Classifies B cells by isotype
  # Inputs:
  #   dt - data.table
  # Outputs:
  #   dt - data.table
  print("Identifying B cell isotype")
  dt[, isotype:="ND"]
  dt[IgA>0.4 & IgA>IgG & IgA>(IgD+0.2) & IgA>(2*IgM-0.3), isotype:="IgA"]
  dt[isotype=="ND" & IgG>0.4 & IgG>IgA & IgG>(IgM+0.05) & IgG>(2*IgD), isotype:="IgG"]
  dt[isotype=="ND" & IgM>0.35 & IgM>(0.5*IgA+0.15) & IgM>(IgG-0.05), isotype:="IgM"]
  dt[isotype=="ND" & IgD>0.20 & IgD>(IgA-0.2) & IgD>(0.5*IgG), isotype:="IgD"]
  dt[isotype=="IgM" & IgD>0.20 & IgD>(IgA-0.2) & IgD>(0.5*IgG), isotype:="IgMD"]
  dt[, isotype:=factor(isotype, levels=c("IgD", "IgMD", "IgM", "IgG", "IgA", "ND"))]
  return(dt)
}


clusterTissue <- function(dt, fa=factors) {
  # Calls somCluster and then metaclusters into ten populations
  # Inputs:
  #   dt - data.table
  #   fa - factor channel names
  # Outputs:
  #   dt - data.table with added meta.cluster column
  print("Clustering data")
  clusters <- c("Transitional", "73- Naïve", "73+ Naïve", "27- Memory", "RB+ 27+ 73- Memory",
                "RB+ 27+ 73+ Memory", "RB- Memory", "95+ Memory", "19++ 11c+ Memory", "Plasma",
                "Germinal Center", "39+ Tonsilar")
  all.markers <- setdiff(colnames(dt), c(fa, "isotype"))
  
  dt[, cluster:=somCluster(dt, channels=all.markers, xdim=15, ydim=15)]
  medians <- dt[, lapply(.SD, median), .SDcols=all.markers, by=cluster] %>%
    .[order(cluster)]
  medians[CD38>0.5 & CD27>0.1, meta:="Plasma"]
  medians[is.na(meta) & CD39>0.45 & CD73<0.1, meta:="39+ Tonsilar"] 
  medians[is.na(meta) & CD19<0.6 & CD11c>0.4, meta:="Non-B"]  
  medians[is.na(meta) & CD72>0.5 & CD11c>0.4, meta:="19++ 11c+ Memory"]
  medians[is.na(meta) & CD24>0.5 & CD38>0.4, meta:="Transitional"]
  medians[is.na(meta) & CD24<0.25 & CD38>0.4, meta:="Germinal Center"]
  medians[is.na(meta) & CD95>0.39 & CD72<0.4, meta:="95+ Memory"]
  medians[is.na(meta) & (CD45RB>0.45 | CD27>0.03 | IgG+IgA>0.4), meta:="Memory"]
  medians[is.na(meta), meta:="Inexperienced"]
  medians[meta=="Inexperienced" & CD73>0.2, meta:="73+ Naïve"]
  medians[meta=="Inexperienced", meta:="73- Naïve"]
  medians[meta=="Memory" & CD45RB<0.53, meta:="RB- Memory"]
  medians[meta=="Memory" & CD27<0.15, meta:="27- Memory"]
  medians[meta=="Memory" & CD73>0.3, meta:="RB+ 27+ 73+ Memory"]
  medians[meta=="Memory", meta:="RB+ 27+ 73- Memory"]
  for (clust in medians$cluster) dt[cluster==clust, meta:=medians[cluster==clust, meta]]
  dt <- dt[meta!="Non-B"]
  dt[, meta:=factor(meta, levels=clusters)]
  return(dt)
}


processSurfaceData <- function(csv, return.dt, mp=paste0(path, "processed_surface.csv")) {
  # Processes data for surface dataset
  # Inputs:
  #   csv - logical, true if csv should be written
  #   return.dt - logical, true if data.table should be returned
  #   mp - character vector of location to write csv file
  # Outputs:
  #   dt - data.table
  meta.dat <<- fread(paste0(path, "figure_4_metadata.csv")) %>%
    .[, lapply(.SD, as.factor)]
  factors <<- colnames(meta.dat)[-1]
  dump.channels <<- c("Time", "Event_length", "beadDist", "DNA_1", "DNA_2", "Viability", "Lineage")
  
  dt <- readFile() %>%
    combineFiles(gated=T) %>%
    asinTransform() %>%
    normalizeSurfaceSamples() %>%
    scaleData() %>%
    normalizeLight() %>%
    clusterSurface() %>%
    assignIgSurface()
  rm(meta.dat, factors, dump.channels, envir=globalenv())
  if (csv) fwrite(dt, mp)
  if (return.dt) return(dt)
  return()
}


processMetabolismData <- function(csv, return.dt, mp=paste0(path, "processed_metabolism.csv")) {
  # Processes data for surface dataset
  # Inputs:
  #   csv - logical, true if csv should be written
  #   return.dt - logical, true if data.table should be returned
  #   mp - character vector of location to write csv file
  # Outputs:
  #   dt - data.table
  meta.dat <<- fread(paste0(path, "figure_5_metabolism_metadata.csv")) %>%
    .[, lapply(.SD, as.factor)]
  factors <<- colnames(meta.dat)[-1]
  dump.channels <<- c("Time", "Event_length", "BC102", "DNA_3", "BC104", "BC105",
                      "BC106", "BC108", "BC110", "lineage", "barium", "dead")
  dt <- readFile() %>%
    combineFiles() %>%
    asinTransform() %>%
    scaleData() %>%
    normalizeLight() %>%
    clusterMetabolism() %>%
    assignIgMetabolism()
  rm(meta.dat, factors, dump.channels, envir=globalenv())
  if (csv) fwrite(dt, mp)
  if (return.dt) return(dt)
  return()
}


processBiosynthesisData <- function(csv, return.dt, mp=paste0(path, "processed_biosynthesis.csv")) {
  # Processes data for surface dataset
  # Inputs:
  #   csv - logical, true if csv should be written
  #   return.dt - logical, true if data.table should be returned
  #   mp - character vector of location to write csv file
  # Outputs:
  #   dt - data.table
  meta.dat <<- fread(paste0(path, "figure_5_biosynthesis_metadata.csv")) %>%
    .[, lapply(.SD, as.factor)]
  factors <<- colnames(meta.dat)[-1]
  dump.channels <<- c("Time", "Event_length", "Barcode_1", "DNA_3", "Barcode_2", "Barcode_3", "CD79b",
                      "Barcode_4", "Barcode_5", "Barcode_6", "biotin", "barium", "Bead_6", "Bead_7")
  dt <- readFile() %>%
    combineFiles() %>%
    asinTransform() %>%
    scaleData() %>%
    normalizeLight() %>%
    clusterBiosynthesis() %>%
    assignIgBiosynthesis()
  rm(meta.dat, factors, dump.channels, envir=globalenv())
  if (csv) fwrite(dt, mp)
  if (return.dt) return(dt)
  return()
}


processSignalingData <- function(csv, return.dt, mp=paste0(path, "processed_signaling.csv")) {
  # Processes data for surface dataset
  # Inputs:
  #   csv - logical, true if csv should be written
  #   return.dt - logical, true if data.table should be returned
  #   mp - character vector of location to write csv file
  # Outputs:
  #   dt - data.table
  meta.dat <<- fread(paste0(path, "figure_5_signaling_metadata.csv")) %>%
    .[, lapply(.SD, as.factor)]
  factors <<- colnames(meta.dat)[-1]
  dump.channels <<- c("Time", "Event_length", "Barcode_1", "DNA_3", "Barcode_2", "Barcode_3",
                      "Barcode_4", "Barcode_5", "Barcode_6", "biotin", "barium", "Bead_5", "DNA_1",
                      "DNA_2")
  dt <- readFile() %>%
    combineFiles() %>%
    asinTransform() %>%
    scaleData() %>%
    .[, `:=`(light=IgK, IgK=NULL, IgL=NULL)] %>%
    clusterSignaling()
  rm(meta.dat, factors, dump.channels, envir=globalenv())
  if (csv) fwrite(dt, mp)
  if (return.dt) return(dt)
  return()
}


processTissueData <- function(csv, return.dt, mp=paste0(path, "processed_tissue.csv")) {
  # Processes data for surface dataset
  # Inputs:
  #   csv - logical, true if csv should be written
  #   return.dt - logical, true if data.table should be returned
  #   mp - character vector of location to write csv file
  # Outputs:
  #   dt - data.table
  meta.dat <<- fread(paste0(path, "figure_6_metadata.csv")) %>%
    .[, lapply(.SD, as.factor)]
  factors <<- colnames(meta.dat)[-1]
  dump.channels <<- c("Time", "Event_length", "DNA_3", "Barcode_A", "Barcode_B", "barium", "Bead_1", "biotin",
                      "Barcode_C", "Barcode_D", "Barcode_E", "Barcode_F", "DNA_1", "DNA_2", "Viability", "CD25")
  dt <- readFile() %>%
    combineFiles() %>%
    asinTransform() %>%
    batchCorrectTissue() %>%
    scaleData() %>%
    normalizeLight(threshold=0.35) %>%
    assignIgTissue() %>%
    clusterTissue()
  rm(meta.dat, factors, dump.channels, envir=globalenv())
  if (csv) fwrite(dt, mp)
  if (return.dt) return(dt)
  return()
}



##### MAIN #####

processSurfaceData(csv=TRUE, return.dt=FALSE)
processMetabolismData(csv=TRUE, return.dt=FALSE)
processBiosynthesisData(csv=TRUE, return.dt=FALSE)
processSignalingData(csv=TRUE, return.dt=FALSE)
processTissueData(csv=TRUE, return.dt=FALSE)