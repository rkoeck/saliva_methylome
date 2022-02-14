###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: preprocessing of raw .idat files using RnBeads

# input: data from each cohort / array type processed separately: IVF data (GSE196432), FLEHS 450k (GSE110128)

# output:  methylation beta value per site per sample (format: .csv), 
#         aggregated methylation values for promoters and CpG islands per sample (format: .csv)

###########################################################################################################################


rm(list=ls(all=T))

#attach libraries
suppressPackageStartupMessages(library(RnBeads))

#set up parallel processing

num.cores <- 30

parallel.setup(num.cores)

#set file directories

data.dir <- "Projects/CultMed"

idat.dir <- file.path(data.dir, "idat/")

sample.annotation <- file.path("annotation.txt")

analysis.dir <- "cohort1"

report.dir <- file.path(analysis.dir, "preprocessing")

#name the Cohort

Cohort = "CultureMedia_Cohort1"

#set the options for RnBeads

rnb.options(analysis.name = "CultureMedia_Cohort1_preprocessing", #name the analysis
            logging = TRUE, #creates a log in the automatic run of the pipeline
            assembly = "hg19", #assembly genome
            analyze.sites = TRUE, #analyse per site/probe - always done for preprocessing steps
            identifiers.column = "Sample_ID", # column name in table of phenotype information to use as sample identifiers otherwise rownames
            gz.large.files = TRUE, #large files should be compressed
            import = TRUE, #carry out import module, only false if using previous RnBSet
            import.sex.prediction = TRUE, #does sex prediction when data is imported
            qc = TRUE, # QC module carried out
            preprocessing = TRUE, #preprocessing steps are completed
            normalization = TRUE, #normalisation is never carried out on sequencing data
            normalization.method = "illumina", #select normalisation method
            normalization.background.method = "none", #no background normalisation is carried out
            filtering.context.removal = c("CC", "CAG", "CAH", "CTG", "CTG"), #only retain CpG probes
            filtering.snp = "any", #remove all probes that overlap with SNPs
            filtering.greedycut = TRUE, #run greedycut filtering to remove low quality probes and samples
            filtering.sex.chromosomes.removal = TRUE, #all probes on sex chromosomes are removed
            filtering.missing.value.quantile = 0.05, # proportion of samples that must have value for probe for it to be included
            imputation.method = "none", #no imputation to be carried out,
            inference = FALSE, #no covariate inference to be done
            exploratory = TRUE, #carry out some steps of the exploratory module
	    exploratory.columns = c("Sample_Plate", "Age", "Gender", "Culturemedium",
				    "Smoking_mother", "Smoking_father"),
	    exploratory.intersample = TRUE,
	    exploratory.deviation.plots = NULL,	
            differential = FALSE, #differential module not to be carried out,
            export.to.bed = FALSE, #export data to bed file
            export.to.trackhub = NULL, #disable export to trackhub (disc full??)
            export.to.csv = TRUE #methylation values are exported to csv files
)

options(fftempdir= analysis.dir)
options(ffcaching="ffeachflush")

rnb.run.analysis(dir.reports=report.dir, sample.sheet=sample.annotation,data.dir=idat.dir, data.type="infinium.idat.dir")

# samples deemed poor quality by this process are removed from subsequent analyses
