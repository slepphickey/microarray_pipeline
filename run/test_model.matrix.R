# This is the third script used to analyze microarray data
# time the script
tic <- as.integer(as.POSIXct( Sys.time() ))

#load the libraries
library(oligo)
library(stringr)
library(affycoretools)
library(tidyverse)
library(purrr)
library(qvalue)


args <- commandArgs(TRUE)
# args[1] is the path to the directory containing the downloaded-files dir
# args[2] name for the output folder and output file prefixes
# args[3] is the path to a sample_table.csv file containing rows for the samples to be included in DE analysis
# args[4] is the path to BrainArrayMapping_withGLP.csv
# args[5] is  a comma separated list of the groups each sample belongs to [optional]
# ie for 6 samples "treated, treated, treated, untreated, untreated, untreated"
# if args[5] is not included this information must be added to arg[3] indepenently as a column named "Condition"

# this is a function that used biomart to map mouse genes to human genes one-to-one homologs
# https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/

#prefix = "test_human"
#dir.name = "/mnt/scratch/hickeys6/test_human"
#meta = read.csv("/mnt/scratch/hickeys6/downloaded-files/hum_diffgpl.csv")
#BAtoGPL = read.csv("/mnt/home/hickeys6/Documents/Pipelines/microarray/data/BrainArrayMapping_withGLP.csv")
print(paste("number of args", length(args)))

prefix = args[2]

dir.name <- paste0(args[1], "/", prefix)
if(!dir.exists(dir.name)) {
  dir.create(dir.name)
}

meta = read.csv(args[3])
BAtoGPL = read.csv(args[4])

if (length(args < 5)) {
  Conditions = unique(meta$Condition)
} else {
  cond = unlist(str_split(args[5], ", "))
  meta$Condition  = cond
}

GSEs = unique(meta$series_id)
GPLs = unique(meta$gpl)

library(limma)

meta$CELfile = basename(meta$CELfile)
f = factor(meta$Condition, levels = Conditions)

print(paste("data comes from", length(GSEs), "GSE(s)"))
print(paste("data comes from", length(GPLs), "platform(s)"))

#if(length(GSEs) > 1 & length(GPLs) > 1) {
#design.mat = model.matrix(~0 + f + series_id + gpl, meta)
#print("removing batch effects caused by GSE and GPL")

#} else if (length(GSEs) > 1 & length(GPLs) == 1) {
#design.mat = model.matrix(~0 + f + series_id, meta)
#print("removing batch effect caused by GSE")

#} else if (length(GSEs) == 1 & length(GPLs) > 1) {

print("removing batch effect caused by GPL")
design.mat = stats::model.matrix(~0 + f + gpl, meta)

# } else {
#design.mat = model.matrix(~0 + f, meta)
#print("no batch effect caused by GSE or GPL")
#}

colnames(design.mat)[1:length(Conditions)] = Conditions
print("design.mat done")
