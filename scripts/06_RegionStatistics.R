###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: statistics testing of CpG sites aggregated into regions (genes, promoters and CpG islands)

# input: preprocessed beta values generated using 01_preprocessingSWAN.R to calculate average methylation values per gene
#        aggregated beta values per promoter and CpG island extracted frmo the output of 01_preprocessingSWAN.R
#        conducted to compare G3 and K-SICM within the IVF cohort (including and excluding pregnancy complications)
#        and once to compare all IVF to all naturally conceived individuals for each region type (genes, promoters, CpG islands)

# targeted analysis: for the targeted analysis, genes were restricted to those contained in the following documents:
#                     imprinted sites: manifestImprinted.csv (source: Ginjala    V. Gene imprinting gateway. Genome Biology 2001.)

# output: (multiple testing corrected) statistics per region - gene/promoter/CGI (.csv)

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

param = SnowParam(workers = 15, type = "SOCK", progressbar = TRUE)
register(param)

# set data directories

message("loading files")

annotation.file = "annotationExcluded.csv"

betas.file =  "betasSWAN.csv"

manifest.file <- "manifest_EPIC_selected_combined.csv"

# load the data files

annotation <- read.csv(annotation.file) %>% select(-X)

betas <- fread(betas.file)

manifest <- fread(manifest.file) %>% select(ID, CHR, custom_Name, custom_Accession, custom_Group)

# select only the columns relating to the QC passing samples

betas <- betas %>% select(ID, annotation$Sample_ID) %>%
        column_to_rownames("ID")

# adjust the sample plate variable in the annotation file

annotation$samplePlate <- if_else(annotation$Sample_Plate == "WG5839006-BCD", "a", "b")

# for consistency of naming rownames of annotation file should correspond to colnames of betas file

annotation <- annotation %>% column_to_rownames("Sample_ID")

# filter the manifest file to contain only probes passing QC

message("grouping by gene name")

manifest <- manifest %>% filter(ID %in% rownames(betas))

sites_per_gene <- manifest %>%
        filter(nchar(custom_Name) > 0) %>%
        group_by(custom_Name, CHR) %>%
        summarise(sites = n())

chroms_per_gene <- sites_per_gene %>%
        group_by(custom_Name) %>%
        summarise(chrs = n())

#remove genes that aren't annotated to a unique chromosome
chroms_per_gene <- chroms_per_gene %>% filter(chrs == 1)

# remove genes that aren't represented by 3 or more individual CpG sites
sites_per_gene <- sites_per_gene %>% filter(custom_Name %in% chroms_per_gene$custom_Name) %>% filter(sites >=3)

#annotate the betas file to contain the gene names

betas.annotated <- betas %>%
        rownames_to_column("ID") %>%
        full_join(manifest, . , by = "ID")

# remove CpG sites that are not relevant for the genes analysis and then aggregate the methylation values by gene

betas.genes <- betas.annotated %>% filter(custom_Name %in% sites_per_gene$custom_Name)

betas.genes <- betas.genes %>%
        select(-ID, -CHR, - custom_Accession, -custom_Group) %>%
        group_by(custom_Name) %>%
        summarise_all(., list(~mean(.,na.rm = T)))

# reformat the betas.genes object to be suitable for statistical modelling

betas.genes <- betas.genes %>% as.data.frame() %>% column_to_rownames("custom_Name")

################ FIT THE STATISTICAL  MODEL ###################################

# create the complex model using the dream package (variancePartition)

# fixed effects = age, gender, culture medium, cell composition
# random effects = sample plate

# specify the model

model <- ~ Culturemedium + Age + Gender + Leukocytes + Epithelial.cells + (1 | samplePlate)

# calculate the weights (this is only possible on sites containing no NA values)

message("calculating weights")

weights <- voomWithDreamWeights(na.omit(betas.genes), model, annotation)

message("model fitting")

fit.complex <- dream(weights, model, annotation)

message("generating output")

complex.table <- topTable(fit.complex, coef = "Culturemediumvg3", number = Inf)

# save the resultant objects as .csv files

write.csv(x = complex.table, file = "dmrsGenesCustom.csv")