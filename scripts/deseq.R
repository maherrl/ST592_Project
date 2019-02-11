########################################################
###           DESeq2 for RAPID                       ###
########################################################

rm(list=ls())


library(phyloseq)
library("DESeq2"); packageVersion("DESeq2")
library(dplyr)
library(data.table)
library("BiocParallel")
library('purrr')
library('tibble')

# Functions
# Function fast_melt for filtering later on
fast_melt = function(physeq,
                     includeSampleVars = character(),
                     omitZero = FALSE){
  require("phyloseq")
  require("data.table")
  # supports "naked" otu_table as `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "TaxaID")
  # Enforce character TaxaID key
  otudt[, TaxaIDchar := as.character(TaxaID)]
  otudt[, TaxaID := NULL]
  setnames(otudt, "TaxaIDchar", "TaxaID")
  # Melt count table
  mdt = melt.data.table(otudt, 
                        id.vars = "TaxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  if(omitZero){
    # Omit zeroes and negative numbers
    mdt <- mdt[count > 0]
  }
  # Omit NAs
  mdt <- mdt[!is.na(count)]
  # Calculate relative abundance
  mdt[, RelativeAbundance := count / sum(count), by = SampleID]
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "TaxaID")
    # Enforce character TaxaID key
    taxdt[, TaxaIDchar := as.character(TaxaID)]
    taxdt[, TaxaID := NULL]
    setnames(taxdt, "TaxaIDchar", "TaxaID")
    # Join with tax table
    setkey(taxdt, "TaxaID")
    setkey(mdt, "TaxaID")
    mdt <- taxdt[mdt]
  }
  # includeSampleVars = c("DaysSinceExperimentStart", "SampleType")
  # includeSampleVars = character()
  # includeSampleVars = c()
  # includeSampleVars = c("aksjdflkas") 
  wh.svars = which(sample_variables(physeq) %in% includeSampleVars)
  if( length(wh.svars) > 0 ){
    # Only attempt to include sample variables if there is at least one present in object
    sdf = as(sample_data(physeq), "data.frame")[, wh.svars, drop = FALSE]
    sdt = data.table(sdf, keep.rownames = TRUE)
    setnames(sdt, "rn", "SampleID")
    # Join with long table
    setkey(sdt, "SampleID")
    setkey(mdt, "SampleID")
    mdt <- sdt[mdt]
  }
  setkey(mdt, "TaxaID")
  return(mdt)
}


# First load unrarefied data (qd.Rdata)
# Subset to only ACR and POC
qd <- subset_samples(qd, species !="POR")
qd <- subset_samples(qd, species !="POC")

# Filter by prevalence and total counts
qdt = fast_melt(qd)
prevdt = qdt[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = TaxaID]
keepTaxa = prevdt[(Prevalence >=3 & TotalCounts >10), TaxaID]
qd = prune_taxa(keepTaxa,qd)


des <- phyloseq_to_deseq2(qd, ~1)
design(des) <- ~treatment

cts = counts(des)
geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(sum(log(row[row != 0]))))
des = estimateSizeFactors(des, geoMeans=geoMeans)
des_acr_trt <- DESeq(des,betaPrior=F)

res_des_acr_trt <- results(des_acr_trt)
summary(res_des_acr_trt) # summary only presents results for one of the comparisons
resultsNames(res_des_acr_trt) # Empty character, not sure why, but the comparisons still exist and can be called by name

head(coef(des_acr_int))

## The following code is for extracting comparisons from the results object
## This is for all the results (time*treatment) but can be adjusted for just time or just treatment

alpha = 0.05
res1 <- results(des_acr_int, name = "time_T1_vs_T0", cooksCutoff = FALSE)
res1 = res1[which(res1$padj < alpha),]
res1 = cbind(as(res1, "data.frame"))
res1$comparison <- "time_T1_vs_T0"
res2 <- results(des_acr_int, name = "time_T2_vs_T0", cooksCutoff = FALSE)
res2 = res2[which(res2$padj < alpha),]
res2 = cbind(as(res2, "data.frame"))
res2$comparison <- "time_T2_vs_T0"
res3 <- results(des_acr_int, name = "time_T3_vs_T0", cooksCutoff = FALSE)
res3 = res3[which(res3$padj < alpha),]
res3 = cbind(as(res3, "data.frame"))
res3$comparison <- "time_T3_vs_T0"
res4 <- results(des_acr_int, name = "time_T4_vs_T0", cooksCutoff = FALSE)
res4 = res4[which(res4$padj < alpha),]
res4 = cbind(as(res4, "data.frame"))
res4$comparison <- "time_T4_vs_T0"
res5 <- results(des_acr_int, name = "treatment_N_vs_C", cooksCutoff = FALSE)
res5 = res5[which(res5$padj < alpha),]
res5 = cbind(as(res5, "data.frame"))
res5$comparison <- "treatment_N_vs_C"
res6 <- results(des_acr_int, name = "treatment_U_vs_C", cooksCutoff = FALSE)
res6 = res6[which(res6$padj < alpha),]
res6 = cbind(as(res6, "data.frame"))
res6$comparison <- "treatment_U_vs_C"
res7 <- results(des_acr_int, name = "timeT1.treatmentN", cooksCutoff = FALSE)
res7 = res7[which(res7$padj < alpha),]
res7 = cbind(as(res7, "data.frame"))
res7$comparison <- "timeT1.treatmentN"
res8 <- results(des_acr_int, name = "timeT2.treatmentN", cooksCutoff = FALSE)
res8 = res8[which(res8$padj < alpha),]
res8 = cbind(as(res8, "data.frame"))
res8$comparison <- "timeT2.treatmentN"
res9 <- results(des_acr_int, name = "timeT3.treatmentN", cooksCutoff = FALSE)
res9 = res9[which(res9$padj < alpha),]
res9 = cbind(as(res9, "data.frame"))
res9$comparison <- "timeT3.treatmentN"
res10 <- results(des_acr_int, name = "timeT4.treatmentN", cooksCutoff = FALSE)
res10 = res10[which(res10$padj < alpha),]
res10 = cbind(as(res10, "data.frame"))
res10$comparison <- "timeT4.treatmentN"
res11 <- results(des_acr_int, name = "timeT1.treatmentU", cooksCutoff = FALSE)
res11 = res11[which(res11$padj < alpha),]
res11 = cbind(as(res11, "data.frame"))
res11$comparison <- "timeT1.treatmentU"
res12 <- results(des_acr_int, name = "timeT2.treatmentU", cooksCutoff = FALSE)
res12 = res12[which(res12$padj < alpha),]
res12 = cbind(as(res12, "data.frame"))
res12$comparison <- "timeT2.treatmentU"
res13 <- results(des_acr_int, name = "timeT3.treatmentU", cooksCutoff = FALSE)
res13 = res13[which(res13$padj < alpha),]
res13 = cbind(as(res13, "data.frame"))
res13$comparison <- "timeT3.treatmentU"
res14 <- results(des_acr_int, name = "timeT4.treatmentU", cooksCutoff = FALSE)
res14 = res14[which(res14$padj < alpha),]
res14 = cbind(as(res14, "data.frame"))
res14$comparison <- "timeT4.treatmentU"

# Combining results into one dataframe and saving
res <- list(res1,res2, res3, res4, res5, res6, res7, res8, res9) %>% map_df(rownames_to_column, 'OTUID')
write.csv(res, file = "/Users/Becca/Box Sync/Winter 2019/Statistical Genomics/ST592_Project/deseq_res_acr_int_try1.csv")
write.csv(res_des_acr_int, file = "/Users/Becca/Box Sync/Winter 2019/Statistical Genomics/ST592_Project/deseq_res_acr_int_all.csv")
write.csv(res_des_acr_time, file = "/Users/Becca/Box Sync/Winter 2019/Statistical Genomics/ST592_Project/deseq_res_acr_time_all.csv")
write.csv(res_des_acr_trt, file = "/Users/Becca/Box Sync/Winter 2019/Statistical Genomics/ST592_Project/deseq_res_acr_trt_all.csv")

#######################################################################################
# Now trying with family table importing as a .csv
#otu <- as.matrix(read.csv(file = "/Users/Becca/Box Sync/RAPID/family-table/otu-table-no-mitochondria-no-chloroplast-min2-names-wtax_L5_t.csv", row.names = "otu.ID"))
#class(otu)

#meta <- read.csv(file = "/Users/Becca/Box Sync/Winter 2019/Statistical Genomics/ST592_Project/data/map-ACRonly.csv", row.names = "sample.id")
#dds <- DESeqDataSetFromMatrix(countData = otu, colData = meta, design = ~1)
# This did not work because the family table from QIIME is proportion data, and DESeq only takes counts
#######################################################################################

## Extracting qd tables for use in ANCOM so the two analyses are comparable and start with the same OTU table.
OTUdat <- as.data.frame(t(otu_table(qd)))
Vardat <- as.data.frame(sample_data(qd))

# Export them here to clean up the tables and add "Sample.ID" as the rownames for ANCOM input.
write.csv(OTUdat, file = "/Users/Becca/Box Sync/Winter 2019/Statistical Genomics/ST592_Project/data/OTUdat_acr_try1.csv")
write.csv(Vardat, file = "/Users/Becca/Box Sync/Winter 2019/Statistical Genomics/ST592_Project/data/Vardat_acr_try1.csv")

