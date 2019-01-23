########################################################
# This script is for performing an ANCOM (Analysis of
# Composition of Microbiomes (ANCOM) to detect 
# differentially abundant taxa in microbial surveys.
# Created by Rebecca Maher
# Created on 10/26/18
# Edited on 01/23/19
########################################################

# clear workspace-----------------------------
rm(list=ls())

# load libraries
library(exactRankTests)
library(nlme)
library(stats)

# load data
# OTU data or taxa data: This should be a data frame with each
# sample in rows and OTUs (or taxa) in columns. The first 
# column should be the sample identifier with column name
# "Sample.ID"

OTUdat <- read.csv(file = "~/data/otu-table-no-mitochondria-no-chloroplast-min2-names-wtax_L5.csv") # by family

# Metadata: Dataframe with the first columns being the sample 
# identifier with column name "Sample.ID"

Vardat <- read.csv(file = "~/data/map.csv")
colnames(Vardat)[1] <- "Sample.ID"

# Subset map to only includ ACR samples from March and May with no NAs
Vardat <- Vardat[which(Vardat$species == "ACR" & Vardat$time == c("T1", "T2") & Vardat$Sample.ID != "March2016.acr35"),]


# ANCOM test
# First must run the function ANCOM.main from the ANCOM_updated_code.R files

longitudinal_comp_test = ANCOM.main(OTUdat = otudf, 
                                    Vardat = Vardat, 
                                    adjusted = F,
                                    repeated = F,
                                    main.var = "time",
                                    adj.formula = "treatment",
                                    repeat.var = NULL,
                                    longitudinal = FALSE,
                                    random.formula = NULL,
                                    multcorr = 2,
                                    sig=0.05,
                                    prev.cut = 0.90)

longitudinal_comp_test$W.taxa

# Working with count otu table
otudf <- read.csv(file = "/Users/Becca/Box Sync/RAPID/ANCOM/otu_table.csv")
