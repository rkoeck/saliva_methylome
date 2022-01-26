###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: Visualisation of the output of 05_EWAS.R and 06_RegionStatistics.R

# input: preprocessed beta values corresponding to the conducted statistical testing (output from 01preprocessingSWAN.R)
#        statistical testing results generated using 05_EWAS.R or 06_RegionStatistics.R
#        manifest file containing sites of interest (imprinted sites/genes & birth-weight associated sites & genes)

# output: volcano plot(s)

###########################################################################################################################

#load packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(stringr))

#set the data directories

annotation.file = "annotationExcluded.csv"

betas.file = "betas_SWAN.csv"

betas.stats = "dmpsDREAMFit.csv"

# load the data files

annotation = fread(annotation.file) %>% select(-V1)

stats = fread(betas.stats) %>% rename("ID" = V1)

betas = fread(betas.file) 

betas = betas %>% select(ID, Chromosome, Start, End, all_of(annotation$Sample_ID)) %>%
  filter(ID %in% stats$ID)

# calculate the summary statistics for the groups

c123 = annotation %>% filter(Culturemedium == "c123") %>% .$Sample_ID

vg3 = annotation %>% filter(Culturemedium == "vg3") %>% .$Sample_ID

betas$meanC = betas %>% select(all_of(c123)) %>% rowMeans(na.rm = T)
betas$sdC = betas %>% select(all_of(c123)) %>% as.matrix() %>% rowSds()

betas$meanV = betas %>% select(all_of(vg3)) %>% rowMeans(na.rm = T)
betas$sdV = betas %>% select(all_of(vg3)) %>% as.matrix() %>% rowSds()

# mean difference - use vg3 as the reference group

betas$meanDiff = betas$meanV - betas$meanC

# deselct the raw data columns

summary = betas %>% select(-all_of(annotation$Sample_ID))

# join the summary stats and the dream stats together, ensuring the order according to p value is maintained

full = full_join(summary, stats, by = "ID") %>% arrange(P.Value)

# make a column for FDR significance (threshold <0.1)

full$significant = full$adj.P.Val < 0.1

# save the summary statistics for the sites to use for the visualisation of targeted analyses

write.csv(summary, file = "summarySites.csv")


############################## visualisation #################################################

#visualise the results as a volcano plot

plot = full %>% 
  ggplot(aes(x = meanDiff, y = -(log10(P.Value)), colour = significant)) +
  geom_point() +
  #ggtitle("Differentially methylated positions") +
  xlab("mean difference (mean G3 - mean K-SICM)") +
  ylab("-log10 p-value") +
  labs(colour = "FDR significance") +
  scale_color_manual(labels = c("Not sig (>0.1)", "<0.1"), values = c("grey", "red")) +
  theme_classic() +
  theme(legend.position = "bottom") +
  scale_y_continuous(position = "right") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = c( -0.2, -0.1, 0, 0.1, 0.2 ), 
                     limits = c(-0.25, 0.2))

ggsave(filename = "/Users/user/surfdrive/PhD/cultureMedia/cohort1/output/exclusionRerun/dmpVolcano_ESHRE.png", plot = plot,
       width = 10, height = 13, unit = "cm")

####################### load the data for imprinted sites and birthweight sites to get position IDs ##############

complicated_bwt= fread("dmpsComplexBirthwtDreamEWAS.csv") %>%
  rename("ID" = V1)

bwt = complicated_bwt$ID

omplicated_imp = fread("dmpsComplexImprintingDream.csv") %>%
  rename("ID" = V1)

imp = omplicated_imp$ID

# add columns to the original data to signify which positions are bwt-assoc/imp

full[ ,"imp"] = full$ID %in% imp
full[ , "bwt"] = full$ID %in% bwt

# volcano plot highlighting imprinted genes

plot.imprinted = full %>% ggplot(aes(x = meanDiff, y = -(log10(P.Value)))) +
  geom_point(colour = "grey", size = 0.4) +
  geom_point(data = function(x){filter(x, imp)}, colour = "deeppink4", size = 0.4) +
  xlab("mean difference \n (mean G3 - mean K-SICM)") +
  ylab("-log10 p-value") +
  theme_classic(base_size = 12) +
  scale_x_continuous(breaks = c( -0.2, -0.1, 0, 0.1, 0.2 ), 
                     limits = c(-0.25, 0.2))

saveRDS(plot.imprinted, file = file.path("dmpComplexImprinted.rds"))

# volcano plot highlighting birth weight-associated genes

plot.bwt = full %>% ggplot(aes(x = meanDiff, y = -(log10(P.Value)))) +
  geom_point(colour = "grey", size = 0.4) +
  geom_point(data = function(x){filter(x, bwt)}, colour = "chartreuse4", size = 0.4) +
  xlab("mean difference \n (mean G3 - mean K-SICM)") +
  ylab("-log10 p-value") +
  theme_classic(base_size = 12) +
  scale_x_continuous(breaks = c( -0.2, -0.1, 0, 0.1, 0.2 ), 
                     limits = c(-0.25, 0.2))

saveRDS(plot.bwt, file = file.path("dmpComplexBwt.rds"))
