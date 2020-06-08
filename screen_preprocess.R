####################################################################################################################
#
# Script: screen_preprocess.R
# Project: An integrated multi-omic single cell atlas of human B cell identity
# Author: David Glass
# Date: 5-27-20
#
# Purpose: Preprocess fcs files from the screen and generates a csv file
#
# FCS files:
#   https://flowrepository.org/id/FR-FCM-Z2MA
#
# Pseudocode:
#   Read in all fcs files
#   Combine data into a single data.table with metadata
#   Remove non-protein channels
#   Identify and keep markers expressed on B cells
#   Asinh transform
#   Quantile normalize conserved channels
#   Scale individual channels
#   Write csv
#
# Instructions:
# Install all required packages (see LIBRARIES and/or Versions below)
# Download the dataset linked above
# Put all the fcs files into a unique directory
# Put "screen_metadata.csv" in the same directory
# In the USER INPUTS section below, assign path variable to the path of the fcs directory
# In the MAIN section at the bottom of the script, if you do not wish to generate a csv of processed data,
#   set csv=FALSE in the processData function (a data.table will be returned)
# By default, the csv is written to the fcs file path, you can change that by assigning mp in the processData function
#   to a different directory
# Run script
# The csv file or data.table can be used in screen_figures.R to recreate figures 1-3
#
### NOTE: Processing often takes a few minutes to run
#
# Versions:
# R 3.6.3
# RStudio 1.2.5042
# flowCore_1.52.1 
# preprocessCore_1.48.0
# dplyr_1.0.0
# data.table_1.12.8
#
#######################################################################################################################



##### USER INPUTS #####

### Path to folder containing screen fcs files and screen_metadata.csv
path <- "~/example/"



###### LIBRARIES ######

require(flowCore)
require(preprocessCore)
require(dplyr)
require(data.table)



##### INPUTS #####

meta.dat <- fread(paste0(path, "screen_metadata.csv")) %>%
  .[, lapply(.SD, as.factor)]
factors <- colnames(meta.dat)[-1]
dump.channels <- c("Time", "Event_length", "File Number", "beadDist", "DNA_1", "DNA_2", "Viability") %>% #non-protein
  c("CD100", "CD266", "TRA_1_81", "CD294", "CD209", "TCRa_b", "CD137", "CD336",
    "TIM_4", "integrin_b5", "SUSD2", "CD324", "CD286", "CD69", "FceR1a", "CD127", "integrin_b7") #technical noise
conserved.channels <- c("CD45", "CD19", "CD24", "CD38", "CD27", "IgM", "IgD")



##### FUNCTIONS #####

readFile <- function(p=path) {
  # takes in a path with fcs files and returns a list of data.tables of the expression matrix
  # Inputs:
  #   p - character vector with directory storing fcs files
  # Outputs:
  #   frames - a list of data.tables
  print("Reading files")
  files <- list.files(path=p, pattern=".fcs", full.names=FALSE, recursive=FALSE)
  frames <- setNames(vector("list", length(files)), files)
  for (i in 1:length(files)) {
    fcs <- read.FCS(paste0(path, files[i]), transformation = FALSE, emptyValue = FALSE)
    frames[[i]] <- data.table(exprs(fcs)) %>%
      setnames(pData(parameters(fcs))$desc)
  }
  return(frames)
}


combineFiles <- function(frames, fa=factors, meta=meta.dat, du=dump.channels, con=conserved.channels) {
  # combines a list of data.tables into a single data.table with factor columns added
  # removes duplicates caused by having parent and daughter gates together
  # Inputs:
  #   frames - a list of data.tables
  #   fa - vector of factor column names
  #   meta - data.table with metadata for the dataset being used
  #   du - vector of channel names to remove
  #   con - conserved channel names
  # Outputs:
  #   dt - a single data.table with all data for that dataset
  print("Combining files")
  all.cols <- lapply(frames, colnames) %>%
    unlist() %>%
    unique() %>%
    c(fa)
  dt <- matrix(nrow=0, ncol=length(all.cols)) %>%
    data.table() %>%
    setnames(all.cols)
  
  for (s in levels(meta$subject)) {
    for (p in levels(meta$panel)) {
      files <- meta[subject==s & panel==p, filename] %>%
        as.vector()
      for (f in files) {
        frames[[f]] <- cbind(frames[[f]], meta[filename==f])
      }
      dt <- rbindlist(frames[files]) %>%
        .[order(gate)] %>%
        unique(by=c(con, "Time")) %>%
        rbind(dt, fill=T)
    }
  }
  return(dt[, setdiff(all.cols, du), with=F])
}


keepPositives <- function(dt, fa=factors, percentile=0.999, threshold=40) {
  # removes channels from dt that are below threshold at percentile
  # Inputs:
  #   dt - data.table of values
  #   fa - factor character names
  #   percentile - vector of percentile to evaluate positivity
  #   threshold - minimum number of counts at percentile to be considered positive
  # Outputs:
  #   dt - data.table with only positive protein columns 
  print("IDing positive markers")
  cols <- setdiff(colnames(dt), fa)
  percentiles <- dt[, lapply(.SD, quantile, probs=c(percentile), na.rm=T), .SDcols=cols] %>%
    as.list() %>%
    unlist()
  keep <- names(percentiles)[percentiles >= 40] %>% c(factors)
  return(dt[, keep, with=F])
}


asinTransform <- function(dt, exceptions=factors) {
  # asinh transforms data
  # Inputs:
  #   dt - data.table
  #   exceptions - character vector of channels not to transform
  # Outputs:
  #   dat - data.table
  print("Asinh transforming")
  to.transform <- setdiff(colnames(dt), exceptions)
  dt[, (to.transform) := asinh(dt[, to.transform, with=F]/5)]
  return(dt)
}


quantileNormalize <- function(dt, con=conserved.channels) {
  # Quantile normalizes values for the conserved channels for each patient
  # Inputs:
  #   dt - data.table of expression data
  #   con - vector of conserved column names
  # Outputs:
  #   dt - data.table of expression data
  print("Quantile normalizing conserved channels")
  for (s in unique(dt$subject)) {
    nrows <- max(table(dt[subject==s, panel]))
    for (c in con) {
      m <- matrix(ncol=length(unique(dt$panel)), nrow=nrows)
      for (i in seq(unique(dt$panel))) {
        m[1:nrow(dt[subject==s & panel==unique(dt$panel)[i]]),i] <-
          as.matrix(dt[subject==s & panel==unique(dt$panel)[i], c, with=F])
      }
      m <- normalize.quantiles(m)
      for (i in seq(unique(dt$panel))) {
        dt[subject==s & panel==unique(dt$panel)[i], eval(c):=m[1:nrow(dt[subject==s & panel==unique(dt$panel)[i]]),i]]
      }
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


processData <- function(csv, return.dt, mp=paste0(path, "processed_screen.csv")) {
  # Processes data for a dataset
  # Inputs:
  #   csv - logical, true if csv should be written
  #   return.dt - logical, true if data.table should be returned
  #   mp - character vector of location to write csv file
  # Outputs:
  #   dt - data.table
  dt <- readFile() %>%
    combineFiles() %>%
    keepPositives() %>%
    asinTransform() %>%
    quantileNormalize() %>%
    scaleData()
  if (csv) fwrite(dt, mp)
  if (return.dt) return(dt)
  return()
}



##### MAIN #####

dat <- processData(csv=TRUE, return.dt=FALSE)