###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: apply the iEVORA algorithm to preprocessed beta values

# input: preprocessed beta values generated using 01_preprocessingSWAN.R
#        conducted to compare G3 and K-SICM within the IVF cohort 

###########################################################################################################################

# load packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(matrixTests))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))

# load the data

annotation = fread("annotationExcluded.csv") %>%
  select(-V1)

betas = fread("betas_SWAN.csv") 

betas = betas %>%
  select(ID, all_of(annotation$Sample_ID)) %>%
  column_to_rownames("ID") %>%
  .[ , order(colnames(.))] %>%
  as.matrix()

# ensure that the annotation has the same sample order & create grouping variable

annotation = annotation %>% mutate(groups = if_else(Culturemedium == "vg3", 0, 1)) %>%
  arrange(Sample_ID)

# group 0 = vg3
# group 1 = c123 (K-SICM)

# apply row_iEVORA

result = row_ievora(betas, annotation$groups, cutT = 0.05, cutBfdr = 0.001)

result = result %>% 
  rownames_to_column("ID") %>%
  arrange(rank)

result = result %>% dplyr::rename("obs.G3" = obs.0, "obs.K-SICM" = obs.1, 
                           "mean.G3" = mean.0, "mean.K-SICM" = mean.1, 
                           "var.G3" = var.0, "var.K-SICM" = var.1)

result = result %>% mutate(var_diff = var.G3 - `var.K-SICM`)

# of the significant diffVar sites see how many are more variable with which culture medium

sig = result %>% filter(significant == TRUE)

table(sig$var.G3 > sig$`var.K-SICM`)

# annotate the significant sites with their gene names

manifest = fread("manifest_EPIC_selected_combined.csv") %>%
  filter(ID %in% sig$ID) %>%
  select(ID, custom_Name)

sig = full_join(sig, manifest, by = "ID")

# check the overlap with imprinting genes

imanifest = fread("manifest_EPIC_imprinted.csv") %>%
  select(-V1)

sig = sig %>% mutate(imp = ID %in% imanifest$ID)

# write the significant files as a .csv

write.csv(sig, "suppFile1.csv")



