###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: applcation of a mixed effects linear model using the variancePartition package

# input: preprocessed beta values generated using 01_preprocessingSWAN.R
#        conducted to compare G3 and K-SICM within the IVF cohort 
#       and once to compare all IVF to all naturally conceived individuals

# targeted analysis: for the targeted analysis, sites were restricted to those contained in the following documents:
#                     imprinted sites: manifestImprinted.csv (source: Ginjala    V. Gene imprinting gateway. Genome Biology 2001.)
#                     birthweight associated sites: manifestBirthweight.csv (source: DOI: 10.1038/s41467-019-09671-3)

# output: (multiple testing corrected) statistics per CpG site (.csv)
#         model design matrix (.csv)

###########################################################################################################################

# load packages

message("loading packages")

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(stringr))

# set up the parallel environment

param = SnowParam(workers = 30, type = "SOCK", progressbar = TRUE)
register(param)

# set data directories

message("loading files")

betas.file <- "betasSWAN.csv"

annotation.file <- "annotationExcluded.csv"

# load the data files

annotation <- read.csv(annotation.file) %>% select(-X)

betas <- fread(betas.file)

# select only the columns relating to the QC passing samples

betas <- betas %>% select(ID, annotation$Sample_ID) %>% 
	column_to_rownames("ID")

# adjust the sample plate variable in the annotation file

annotation$samplePlate <- if_else(annotation$Sample_Plate == "WG5839006-BCD", "a", "b")


# for consistency of naming rownames of annotation file should correspond to colnames of betas file

annotation <- annotation %>% column_to_rownames("Sample_ID")

# create the complex model using the dream package (variancePartition)

# fixed effects = age, gender, culture medium, cell composition
# random effects = sample plate

# specify the model

model <- ~ Culturemedium + Age + Gender + Leukocytes + Epithelial.cells + (1 | samplePlate)

# calculate the weights (this is only possible on sites containing no NA values)

message("calculating weights")

weights <- voomWithDreamWeights(na.omit(betas), model, annotation)

message("model fitting")

fit.complex <- dream(weights, model, annotation)

message("generating output")

design <- fit.complex$design

complex.table <- topTable(fit.complex, coef = "Culturemediumvg3", number = Inf)

# save the resultant objects as .csv files

write.csv(x = design, file = "dmpsDREAMDesign.csv")

write.csv(x = complex.table, file = "dmpsDREAMFit.csv")

