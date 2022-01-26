###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: visualisation of the cell count esimates generated using the 03.2_CellCompositionEstimation.R script

# input: cell composition file generated as output from 03.2_CellCompositionEstimation.R
#        visualisation conducted separately for IVF samples and naturally conceived samples 

# output: violin plot overlaid with box plot & individual data points

###########################################################################################################################

# load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ewastools))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))

# set the data directories

base.dir <- "output/"

annotation <- fread(file.path(base.dir, "/annotationExcluded.csv"))

# visualise the cell compositions

annotation.long <- annotation %>% select(Sample_ID, Culturemedium, Leukocytes, Epithelial.cells) %>%
  gather(., cellType, percentage, Leukocytes:Epithelial.cells)

annotation.long$cellType <- if_else(annotation.long$cellType == "Leukocytes", "Leukocytes", "Epithelial cells" )

annotation.long$Culturemedium = if_else(annotation.long$Culturemedium == "vg3", "G3", "K-SICM")

# plot cell type by culture medium

plot =  annotation.long %>% ggplot(aes(x = Culturemedium, y = percentage)) +
  facet_wrap(~cellType) +
  geom_violin(aes(fill = Culturemedium)) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_point(size = 0.2, position = position_jitter(width = 0.1)) +
  theme_classic() +
  ylab("Proportional composition") +
  scale_fill_manual(values = c("K-SICM" = "darkorchid", "G3" = "darkorange1")) +
  scale_colour_manual(values = c("K-SICM" = "darkorchid", "G3" = "darkorange1" )) +
  xlab("Culture medium") +
  theme(legend.position = "none")
  

saveRDS(plot, file = "composition.rds")

ggsave(filename = file.path(base.dir, "cellComp.png"), plot = comp,
       width = 12, height = 9, units = "cm")

