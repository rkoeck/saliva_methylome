###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: visualisation of the number of outliers per sample & determine if the number of outliers is associated with any 
#                 sample features

# input: output from 08_Outliers.R (summarised total number of outliers per sample)
#        conducted to compare G3 and K-SICM within the IVF cohort 

# output: multi-faceted plot containing scatter plot of hypo to hyper methylation outliers in the centre with adjacent distribution summaries

# source code origin: https://rdrr.io/github/jorvlan/openvis/src/R/R_rainclouds.R

###########################################################################################################################

# load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(Rmisc))
suppressPackageStartupMessages(library(jmuOutlier))

source("rainCloudPlot.R")

# plot formatting

plot.format = theme(text = element_text(family = "Helvetica", size = 7),
                    plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
                    axis.text.x = element_text(size = 7),
                    legend.title = element_text(size = 7, face = "bold"),
                    legend.text = element_text(size = 7),
                    axis.text.y = element_text(size = 7),
                    axis.title.x = element_text(size = 7),
                    axis.title.y = element_text(size = 7),
                    plot.margin = margin())

# load the data

fig.dir = "figures/"

annotation = fread("annotationExcluded.csv") %>% select(-V1)

outliers = fread("sampleThresholdOutliers.csv",
                 header = TRUE) %>% 
  column_to_rownames("V1") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("Sample_ID") %>%
  mutate(Sample_ID = as.character(Sample_ID))

# merge the outlier values with the annotation data

outliers = full_join(annotation, outliers, by = "Sample_ID")

outliers$Culturemedium = if_else(outliers$Culturemedium == "vg3", "G3", "K-SICM")

outliers$Culturemedium = factor(outliers$Culturemedium, levels = c("K-SICM", "G3"))

# visualise the results as a scatter plot with density information summarised at the side 

scatter = outliers %>% ggplot(aes(x = hypo, y = hyper, colour= Culturemedium)) +
  geom_point(size = 0.3) +
  scale_x_log10(breaks = c(20, 100, 1000, 10000), labels = c("20", "100", "1000", "10000"), limits = c(20, 40000), expand = expansion(add = c(0.1, 0.3))) +
  scale_y_log10(breaks = c(20, 100, 1000, 10000), labels = c("20", "100", "1000", "10000"), limits = c(20, 40000), expand = expansion(add = c(0.1, 0.3))) +
  theme_classic() +
  scale_fill_manual(values = c("K-SICM" = "darkorchid", "G3" = "darkorange1")) +
  scale_colour_manual(values = c("K-SICM" = "darkorchid", "G3" = "darkorange1" )) +
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  plot.format +
  theme(legend.title.align = 0.5,
        legend.background = element_rect(color = "black")) +
  ylab("Hypermethylation outliers per sample") +
  xlab("Hypomethylation outliers per sample")

x = outliers %>% ggplot(aes(x = Culturemedium, y = hypo, fill = Culturemedium, colour = Culturemedium)) +
  geom_flat_violin(position = position_nudge(x = -0.2, y = 0), adjust = 0.75, trim = TRUE) +
  geom_boxplot(aes(x = as.numeric(Culturemedium) -0.2, y = hypo), outlier.shape = NA, alpha = 0.5, width = 0.1, colour = "black", lwd = 0.3) +
  theme_classic() + 
  scale_y_log10(breaks = c(20, 100, 1000, 10000), labels = c("20", "100", "1000", "10000"), limits = c(20, 40000), expand = expansion(add = c(0.1, 0.3))) +
  coord_flip() +
  scale_fill_manual(values = c("K-SICM" = "darkorchid", "G3" = "darkorange1")) +
  scale_colour_manual(values = c("K-SICM" = "darkorchid", "G3" = "darkorange1")) +
  scale_shape_manual(values = c(19, 15, 1, 17)) +
  xlab(NULL) +
  ylab("log10(number of outliers per sample)") +
  plot.format +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank())

y = outliers %>% ggplot(aes(x = Culturemedium, y = hyper, fill = Culturemedium, colour = Culturemedium)) +
  geom_flat_violin(position = position_nudge(x = -0.2, y = 0), adjust = 0.75, trim = TRUE) +
  geom_boxplot(aes(x = as.numeric(Culturemedium) -0.2), outlier.shape = NA, alpha = 0.5, width = 0.1, colour = "black", lwd = 0.3) +
  theme_classic() +
  scale_y_log10(breaks = c(20, 100, 1000, 10000), labels = c("20", "100", "1000", "10000"), limits = c(20, 40000), expand = expansion(add = c(0.1, 0.3))) +
  scale_fill_manual(values = c("K-SICM" = "darkorchid", "G3" = "darkorange1")) +
  scale_colour_manual(values = c("K-SICM" = "darkorchid", "G3" = "darkorange1")) +
  scale_shape_manual(values = c(19, 15, 1, 17)) +
  xlab(NULL) +
  ylab("log10(number of outliers per sample)") +
  #theme(legend.position = "none") +
  plot.format +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) 

pdf(file = file.path(fig.dir, "altOutliers.pdf"),
    width = 3.34, height = 3.34)
x + guide_area() + scatter + y + 
  plot_layout(ncol = 2,
              heights = c(1,3),
              widths = c(3, 1),
              guides = "collect")
dev.off()

################################# look at the characteristics of children with a very high number of outliers ##############

upperQ = quantile(outliers$outlier, probs = 0.75)

high = outliers %>% filter(outlier >= upperQ)

demographics = fread("Saliva_Annotation.txt") %>%
  filter(Sample_ID %in% annotation$Sample_ID)

dem.high = demographics %>% filter(Sample_ID %in% high$Sample_ID)

diag = demographics %>% filter(UrologicProbl == 1 | AllergicProblems == 1 | Autism_related_disorder == 1)

diag.outliers = outliers %>% filter(Sample_ID %in% diag$Sample_ID)


############################## look for assoctiations between characteristics and number of outliers #####################

full = full_join(demographics, outliers, by = "Sample_ID")

cor(full$hypo, full$Age.x, method = "pearson")
cor(full$hyper, full$Age.x, method = "pearson")
cor(full$outlier, full$Age.x, method = "pearson")

cor(full$hypo, full$W_9J, method = "pearson")
cor(full$hyper, full$W_9J, method = "pearson")
cor(full$outlier, full$W_9J, method = "pearson")

cor(full$hypo, full$Geboortegewicht, method = "pearson")
cor(full$hyper, full$Geboortegewicht, method = "pearson")
cor(full$outlier, full$Geboortegewicht, method = "pearson")

cor(full$hypo, full$Leukocytes, method = "pearson")
cor(full$hyper, full$Leukocytes, method = "pearson")
cor(full$outlier, full$Leukocytes, method = "pearson")

cor(full$hypo, full$Epithelial.cells, method = "pearson")
cor(full$hyper, full$Epithelial.cells, method = "pearson")
cor(full$outlier, full$Epithelial.cells, method = "pearson")

significances = data.frame()

for(v in c("Age.x", "W_9J", "Geboortegewicht", "Leukocytes", "Epithelial.cells")){
  
  for(o in c("hypo", "hyper", "outlier")){
    
    variable = v
    
    outlier = o
    
    data = full[ , c(..variable, ..outlier)] %>% as.matrix()
    
    test = perm.cor.test(x = data[,1], y = data[,2], alternative = "two.sided", method = "pearson", num.sim = 10000)
    
    significances[v, o] = test$p.value
    
  }
  
}

# for sample plate use wilcoxon signed rank test

x = full %>% filter(Sample_Plate.x == "WG5839006-BCD") %>% select(outlier)

x = x$outlier

y = full %>% filter(Sample_Plate.x != "WG5839006-BCD") %>% select(outlier) 

y = y$outlier 

test = wilcox.test(x = x, y = y, alternative = "two.sided")

test$p.value


x = full %>% filter(Sample_Plate.x == "WG5839006-BCD") %>% select(hypo)

x = x$hypo

y = full %>% filter(Sample_Plate.x != "WG5839006-BCD") %>% select(hypo) 

y = y$hypo

test = wilcox.test(x = x, y = y, alternative = "two.sided")

test$p.value


x = full %>% filter(Sample_Plate.x == "WG5839006-BCD") %>% select(hyper)

x = x$hyper

y = full %>% filter(Sample_Plate.x != "WG5839006-BCD") %>% select(hyper) 

y = y$hyper

test = wilcox.test(x = x, y = y, alternative = "two.sided")

test$p.value
