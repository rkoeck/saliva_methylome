###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: PCA to explore global methylation & determine associations with covariates

# input: swan normalised beta values (output from 01_preprocessingSWAN.R), poor quality samples have been removed
#        conducted on each cohort separately to identify batch effects 
#        conducted on the full dataset (IVF & naturally conceived data)

# output: PCA plot
#         statistics for associations between sample characteristics & PCs
#         heatmap showing associations between sample characteristics & PCs

###########################################################################################################################

#load packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(jmuOutlier))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(tidyr))

# set data directories
annotation.file = "annotationExcluded.csv"

betas.file = "betas_SWAN.csv"

# load the data files

annotation <- fread(annotation.file) %>% 
  select(-V1)

betas <- fread(betas.file) 

betas = betas %>% 
  select(ID, all_of(annotation$Sample_ID)) %>% column_to_rownames("ID")

betas.pca <- betas %>% na.omit() %>%
  as.matrix() %>%
  t()

# conduct a principal component analysis
pca <- prcomp(betas.pca, center = T, scale. = F)

# to see % of variance explained by each PC

summary(pca)

# plot the PCA

coords <- as.data.frame(pca$x) %>% rownames_to_column("Sample_ID") %>%
  full_join(annotation, .[ , 1:15], by = "Sample_ID")

coords$Culturemedium = if_else(coords$Culturemedium == "vg3", "G3", "K-SICM")

pca.plt <- coords %>% ggplot(aes(x = PC1, y = PC2, colour = Culturemedium)) + 
  geom_point(size = 0.6) +
  scale_colour_manual(values = c("K-SICM" = "darkorchid", "G3" = "darkorange1"), name = "Culture medium") + 
  theme_classic() +
  xlab("PC1 (32%)") +
  ylab("PC2 (18%)")
  
saveRDS(pca.plt, file = file.path("globalPCA.rds"))

ggsave(filename = "globalPCA_ESHRE.png", plot = pca.plt, 
       width = 15, height = 10, units = "cm",
       dpi = 320)


# determine which (if any) PCs are associated with sample features

plate = list()

for(n in 1:8){
  
  pc = paste0("PC", n)
  
  x = coords %>% filter(Sample_Plate == "WG5839006-BCD") %>% select(all_of(pc))
  
  x = x[,1]
  
  y = coords %>% filter(Sample_Plate != "WG5839006-BCD") %>% select(all_of(pc)) %>% as.vector()
  
  y = y[,1]
  
  test = wilcox.test(x = x, y = y, alternative = "two.sided")
  
  plate[[pc]] = test$p.value
  
}

gender = list()

for(n in 1:8){
  
  pc = paste0("PC", n)
  
  x = coords %>% filter(Gender == "M") %>% select(all_of(pc))
  
  x = x[,1]
  
  y = coords %>% filter(Gender == "F") %>% select(all_of(pc)) %>% as.vector()
  
  y = y[,1]
  
  test = wilcox.test(x = x, y = y, alternative = "two.sided")
  
  gender[[pc]] = test$p.value
  
}


medium = list()

for(n in 1:8){
  
  pc = paste0("PC", n)
  
  x = coords %>% filter(Culturemedium == "vg3") %>% select(all_of(pc))
  
  x = x[,1]
  
  y = coords %>% filter(Culturemedium == "c123") %>% select(all_of(pc)) %>% as.vector()
  
  y = y[,1]
  
  test = wilcox.test(x = x, y = y, alternative = "two.sided")
  
  medium[[pc]] = test$p.value
  
}

age = list()
cor.age = list()

for(n in 1:8){
  
  pc = paste0("PC", n)
  
  data = coords[ , c("Age", pc)]
  
  test = perm.cor.test(x = data[,1], y = data[,2], alternative = "two.sided", method = "pearson", num.sim = 10000)
  
  age[[pc]] = test$p.value
  
  cor.age[[pc]] = cor(data[,1], data[,2], method = "pearson")
  
}


leukocytes = list()

cor.leuk = list()

for(n in 1:8){
  
  pc = paste0("PC", n)
  
  data = coords[ , c("Leukocytes", pc)]
  
  test = perm.cor.test(x = data[,1], y = data[,2], alternative = "two.sided", method = "pearson", num.sim = 10000)
  
  leukocytes[[pc]] = test$p.value
  
  cor.leuk[[pc]] = cor(data[,1], data[,2], method = "pearson")
  
}

epis = list()

cor.epis = list()

for(n in 1:8){
  
  pc = paste0("PC", n)
  
  data = coords[ , c("Epithelial.cells", pc)]
  
  test = perm.cor.test(x = data[,1], y = data[,2], alternative = "two.sided", method = "pearson", num.sim = 10000)
  
  epis[[pc]] = test$p.value
  
  cor.epis[[pc]] = cor(data[,1], data[,2], method = "pearson")
  
}

sentID = list()

for(n in 1:8){
  
  pc = paste0("PC", n)
  
  data = coords[ , c("Sentrix_ID", pc)]
  
  test = kruskal.test(data[,1]~data[,2])
  
  sentID[[pc]] = test$p.value
  
}

sentPos = list()

for(n in 1:8){
  
  pc = paste0("PC", n)
  
  data = coords[ , c("Sentrix_Position", pc)]
  
  test = kruskal.test(data[,1]~data[,2])
  
  sentPos[[pc]] = test$p.value
  
}

associations = data.frame(row.names = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")) %>%
  rownames_to_column("Principal_component") %>%
  mutate(Sample_plate = unlist(plate), 
         Sentrix_ID = unlist(sentID),
         Sentrix_Position = unlist(sentPos),
         Sex = unlist(gender),
         Age = unlist(age),
         Culture_Medium = unlist(medium),
         Leukocytes = unlist(leukocytes),
         Epithelial_Cells = unlist(epis))


# plot the associations as a heatmap 

t = to.plot %>% as.data.frame() %>% rownames_to_column("V1")

t2 = gather(t, key = "PC", value = "log.pvals", PC1:PC8)

t2$V1 = factor(t2$V1, levels = rev(rownames(to.plot)), ordered = T)

t2 = t2 %>% mutate(p.vals = 10^-log.pvals, text = p.vals <= 0.05, pval.text = formatC(p.vals, format = "e", digits = 1))

plot = t2 %>% ggplot(aes(x = PC, y = V1, fill = log.pvals)) + 
  geom_tile(colour = "white") + 
  scale_fill_gradientn(colours = c("white", "lightblue", "dodgerblue2"), 
                       values = scales::rescale(c(0, -log10(0.05), 8)), limits = c(0,5), oob = scales::squish) +
  geom_text(data = function(x){filter(x, text)}, aes(label = pval.text), size = 2, angle = 20) +
  theme_classic() +
  labs(fill = "-log10 \n p-value") +
  scale_y_discrete(labels = c("Sample_plate" = "Sample plate", "Sentrix_ID" = "Sentrix ID", "Sentrix_Position" = "Sentrix position",
                              "Culture_Medium" = "Culture medium", "Epithelial_Cells" = "Epithelial cells")) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(), 
        text = element_text(family = "Helvetica", size = 7),
        plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 7),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 7),
        axis.text.y = element_text(size = 7))

saveRDS(plot, file = file.path("associationsHeatmap.rds"))