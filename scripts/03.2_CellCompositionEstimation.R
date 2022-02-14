###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: Estimate cellular composition of umbilical cord blood samples

# tissue: saliva
# method: Houseman algorithm
# reference: ewastools - saliva (Middleton et al, 2020 (preprint))

# input: .idat files from each cohort / array type processed separately: IVF data (GSE196432), FLEHS 450k (GSE110128)

# output: cell composition estimnates for the specified cell types (.csv) per sample
#         the estimates are combined with the sample annotation file for downstream processing

###########################################################################################################################

# load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ewastools))
suppressPackageStartupMessages(library(tidyr))
# set the data directories

base.dir <- "output/"

annotation.file <- file.path(base.dir, "annotation.csv")

betas.file <- file.path(base.dir, "betas_noob.csv")

# load the annotation file

annotation <- fread(annotation.file)

betas <- fread(betas.file) %>%
  column_to_rownames("V1")

# calculate the cell composition
cells <- estimateLC(betas, "Saliva") %>% as.data.frame()

rownames(cells) <- colnames(betas)

# write the cell composition file 

write.csv(x = cells, file = file.path(base.dir, "composition.csv"))

# add the cell composition to the annotation file

annotation <- cells %>% rownames_to_column("barcode") %>%
  full_join(annotation, ., by = "barcode")

# write out the annotation file containing the cell composition

write.csv(x = annotation, file = file.path(base.dir, "annotation_cells.csv"))
