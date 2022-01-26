###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: script to compose (multi-panel) figures, make them uniform and save them in an appropriate format for submission

# source code origin: https://rdrr.io/github/jorvlan/openvis/src/R/R_rainclouds.R

###########################################################################################################################

# load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(cowplot))

source("rainCloudPlot.R")

# set directory

fig.dir = "figures/"

# function to generalise formatting (fonts etc)

plot.format = theme(text = element_text(family = "Helvetica", size = 7),
                    plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
                    axis.text.x = element_text(size = 7),
                    legend.title = element_text(size = 7, face = "bold"),
                    legend.text = element_text(size = 7),
                    axis.text.y = element_text(size = 7),
                    axis.title.x = element_text(size = 7),
                    axis.title.y = element_text(size = 7))
  
  
  
####################### figure 1: global methylation #############################

p1 = readRDS(file.path(fig.dir, "globalPCA.rds")) + plot.format +
  theme(legend.position = "bottom",
        legend.margin = margin(t = -10))

pdf(width = 3.34, height = 3, file = file.path(fig.dir, "Figure1PCA.pdf"))
p1
dev.off()

jpeg(file.path(fig.dir, "Figure1.jpg"),
     width = 3.34, height = 3, units = "in",
     res = 900)
p1
dev.off()

################### figure 2: technical data processing ###################
p2 = readRDS(file.path(fig.dir, "sEST.rds")) + plot.format +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.key.height = unit(0.01, "cm"))

p3 = readRDS(file.path(fig.dir, "composition.rds")) + plot.format +
  theme(legend.position = "none")

p4 = readRDS(file.path(fig.dir, "associationsHeatmap.rds")) 

top = plot_grid(p2, p3, ncol = 2, rel_widths = c(1,1),
                labels = c("A", "C"),
                label_size = 8,
                label_fontface = "bold",
                scale = 0.95)

pdf(width = 6.69, height =6, file = file.path(fig.dir, "Figure2.pdf"))
plot_grid(top, p4, nrow = 2,
          rel_heights = c(1, 0.9),
          labels = c("", "B"),
          label_size = 8,
          scale = 0.95)
dev.off()


##################### figure 3: DNA methylation at individual CpGs #############

p5 = readRDS(file.path(fig.dir, "dmpComplexImprinted.rds")) + plot.format

p6 = readRDS(file.path(fig.dir, "dmpComplexBwt.rds")) + plot.format

jpeg(file.path(fig.dir, "Figure3ab.jpg"),
     width = 4.5, height = 3, units = "in",
     res = 900)
p5 + p6 + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 8, face = "bold"))
dev.off()


pdf(width = 4.5, height = 3, file = file.path(fig.dir, "Figure3ab.pdf"))
p5 + p6 + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 8, face = "bold"))
dev.off()

##################### figure 3: regional DNA methylation ######################

p7 = readRDS(file.path(fig.dir, "genesComplexImprinted.rds")) + ggtitle("Genes") + plot.format + theme(legend.position = "none")

p8 = readRDS(file.path(fig.dir, "promoterComplex.rds")) + plot.format + theme(legend.position = "none") +
  xlab("mean difference \n (mean G3 - mean K-SICM)")

p9 = readRDS(file.path(fig.dir, "islandComplicated.rds")) + plot.format + theme(legend.position = "none")

jpeg(file.path(fig.dir, "Figure3cde.jpg"),
     width = 6, height = 3, units = "in",
     res = 900)
p7 + p8 + p9 + plot_annotation(tag_levels = list(c("C", "D", "E"))) & theme(plot.tag = element_text(size = 8, face = "bold"))
dev.off()

pdf(width = 10, height = 5, file = file.path(fig.dir, "Figure3cde.pdf"))
p7 + p8 + p9 + plot_annotation(tag_levels = list(c("C", "D", "E"))) & theme(plot.tag = element_text(size = 8, face = "bold"))
dev.off()

################## figure 4: outliers per sample ##########################

# figure finalised and saved as .pdf in the script from which it originates

################# Supplementary figure NC data processing ################

partA = readRDS("sESTNC.rds") +
  plot.format +
  ylab("Y.PC1") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2, byrow = T),
         shape = guide_legend(nrow = 2, byrow = T))

partB = readRDS("compositionNC.rds") +
  plot.format +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 3, byrow = T))

partC = readRDS("associationsHeatmapNC.rds") +
  plot.format +
  xlab("") +
  ylab("")

top = plot_grid(partA, partB, ncol = 2,
                rel_widths = c(1,1),
                labels = c("A", "B"),
                label_size = 8,
                label_fontface = "bold",
                scale = 0.95)

full = plot_grid(top, partC, nrow = 2,
                 rel_heights = c(1.3,1),
                 labels = c("", "C"),
                 label_size = 8,
                 label_fontface = "bold",
                 scale = c(1, 0.95))

pdf(width = 6.69, height = 6, 
    file = "ncProcessing.pdf")
full
dev.off() 

jpeg("ncProcessing.jpg",
     width = 6.69, height = 6, units = "in",
     res = 900)
full
dev.off()


################ Supplementary figure NC & IVF comparison ##################

plotA = readRDS("globalPCANCIVF.rds") +
  plot.format 

plotB = readRDS("associationsHeatmapAllNCIVF.rds") +
  plot.format +
  ylab("") + xlab("")

plotC = readRDS("densityNCIVF.rds")

plotC = plotC + plot.format + theme(legend.position = "none")

plotD = readRDS("volcanoDmpsNCIVF.rds") +
  plot.format

left = plot_grid(plotA, plotC, nrow = 2,
                rel_heights = c(2, 1),
                labels = c("A", "C"),
                label_size = 8,
                label_fontface = "bold",
                scale = 0.95)

rightbottom = plot_grid(NULL, plotD, NULL, ncol = 3,
                        rel_widths = c(0.5, 1, 0.5),
                        labels = c("", "D", ""),
                        label_size = 8,
                        label_fontface = "bold")

right = plot_grid(plotB, rightbottom, nrow = 2,
                rel_widths = c(1, 1),
                labels = c("B", ""),
                label_size = 8,
                label_fontface = "bold",
                scale = 0.95)

all = plot_grid(left, right, ncol = 2,
                rel_widths = c(1, 2),
                labels = "")

pdf(width = 6.69, height =8, 
    file = "ncIVFComp2.pdf")
all
dev.off() 

jpeg("ncIVFComp2.jpg",
     width = 6.69, height = 8, units = "in",
     res = 900)
all
dev.off()
