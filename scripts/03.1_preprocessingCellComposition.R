###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: preprocess the data for cell composition

# tissue: saliva
# method: preprocessing Noob (minfi)

# input: .idat files from each cohort / array type processed separately: IVF data (GSE196432), FLEHS 450k (GSE110128)
#         samples deemed poor quality by RnBeads processing are not included

# output: noob normalised beta values for each sample at each CpG site

############################################################################################################################ preprocessing script for epigenetic age estimation

# load packages
message("loading packages")

suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))

# set data directories

message("loading files")

data.dir <- "idat/"

annotation.file <- "annotationExcluded.csv"

# load the annotation file

annotation <- read.csv(annotation.file) %>% select(-X)

annotation$Basename = as.character(annotation$Basename)

# load the idat files

RGSet = read.metharray.exp(data.dir, targets = annotation)

# process the data using preprocessNoob

message("preprocessing")

MSet = preprocessNoob(RGSet)

message("converting to ratio set")

ratioSet = ratioConvert(MSet)

# extract only the beta values

betas = getBeta(ratioSet)

# write the betas file out

write.csv(betas, file = "betas_Noob.csv")



