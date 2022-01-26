###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: sex prediction the sEst package & visualisation of the results

# input: detection p-values (.csv) and raw beta values (.csv) generated in the script 02.1_preprocessingSest.R
#         each cohort is run separately as described in 02.1_preprocessingSest.R

# output: updated annotation file containing pedicted sex columns (.csv)
#         scatter plot of the results

###########################################################################################################################

# sex prediction for cohort 1 using sEst

# load packages

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(sest))
suppressPackageStartupMessages(library(tibble))

# load the data

betas = fread("rawBetas.csv")

pvals = fread("detPvals.csv")

ann = fread("annotation.csv")

# correctly format the betas and pvals for the sest tool

betas = betas %>% column_to_rownames("V1") %>% as.matrix()

pvals = pvals %>% column_to_rownames("V1") %>% as.matrix()

# conduct sex estimation

sex = estimateSex(beta.value = betas, detecP = pvals)

# compare the annotated gender to the predicted gender

ann = sex$test %>% rownames_to_column("Basename") %>% 
  select(Basename, predicted.X, predicted.Y, predicted, X.PC1, Y.PC1) %>%
  full_join(ann, ., by = "Basename")

ann$Sex = if_else(ann$Gender == "Male", "M", "F")

mismatch = ann %>% filter(Sex != predicted)

# remove samples deemed poor quality by greedycut

poor = c("I-088", "I-052", "I-033", "I-041", "I-120")

ann2 = ann %>% filter(!Sample_ID %in% poor)

# adjust the X.PC1 threshold to 0.05

ann2$predicted.X2 = if_else(ann2$X.PC1 < 0.05, "M", "F")
ann2$predicted2 = if_else(ann2$predicted.X2 == "M" & ann2$predicted.Y == "M", "M",
                          if_else(ann2$predicted.X2 == "F" & ann2$predicted.Y == "F", "F", "N"))

ann2 %>% ggplot(aes(x = Gender, y = `Predicted Male Probability`, colour = predicted2)) +
  geom_point(position = position_jitter(width = 0.4), alpha = 0.8) +
  theme_classic()

# visualise the results with the adjusted threshold

plot = ann2 %>% ggplot(aes(x = X.PC1, y = Y.PC1, colour = predicted2, shape = Gender)) +
  geom_point() +
  theme_classic() +
  scale_shape_manual(values = c("F" = 16, "M" = 17), labels = c("Female", "Male")) +
  scale_color_manual(values = c("F" = "plum2", "M" = "skyblue2", "N" = "grey"), labels = c("Female", "Male", "Not specified")) +
  labs(shape = "Recorded sex", colour = "Predicted Sex")

saveRDS(plot, file = "sEST.rds")

# samples with a mismatched sex prediction and recorded sex are removed from downstream analyses
