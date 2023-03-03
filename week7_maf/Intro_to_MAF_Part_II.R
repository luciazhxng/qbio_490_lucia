if(!require(survival)) {
  install.packages("survival")
}
library(survival)
if(!require(survminer)) {
  install.packages("survminer")
}
library(survminer)
if(!require(ggplot2)) {
  install.packages("ggplot2")
}
library(ggplot2)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
if (!require("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks")
if (!require("maftools", quietly = TRUE))
  BiocManager::install("maftools")
library(BiocManager)
library(TCGAbiolinks)
library(maftools)
setwd("C:/Users/zhang/OneDrive/college/qbio490/qbio_490_lucia/analysis_data") # set working directory

# read in the clinical data
clinical <- read.csv("C:/Users/zhang/OneDrive/college/qbio490/qbio_490_lucia/analysis_data/brca_clinical_data.csv")
# read in the MAF data
maf_query <- GDCquery(project = 'TCGA-BRCA', data.category = "Simple Nucleotide Variation", access = "open", data.type = "Masked Somatic Mutation", workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
#GDCdownload(maf_query)
maf <- GDCprepare(maf_query)
maf_object <- read.maf(maf = maf, clinicalData = clinical, isTCGA = TRUE)

# change the name of the clinical data to match the MAF
colnames(clinical)[ colnames(clinical) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"
write.csv(clinical, "C:/Users/zhang/OneDrive/college/qbio490/qbio_490_lucia/analysis_data/brca_clinical_data.csv", row.names = FALSE)
# remove rows where the lymph node is NA
maf_object@clinical.data <- maf_object@clinical.data[!is.na(maf_object@clinical.data$lymph_node_examined_count),]
# create a categorical version of lymph node numbers
maf_object@clinical.data$lymph_node_categorical <- ifelse(maf_object@clinical.data$lymph_node_examined_count >= 15, "Greater than 15", "Less than 15") 
maf_object@clinical.data$lymph_node_categorical

less_mask <- ifelse(maf_object@clinical.data$lymph_node_categorical == "Less than 15", T, F)
less_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[less_mask]

# create a MAF subset of just the people with less lymph nodes
less_maf <- subsetMaf(maf = maf_object, tsb = less_barcodes)
# create a MAF subset of just the people with more lymph nodes
more_mask <- ifelse(maf_object@clinical.data$lymph_node_categorical == "Greater than 15", T, F)
more_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[more_mask]
more_maf <- subsetMaf(maf = maf_object, tsb = more_barcodes)
# find the top 10 genes
less_maf.genes = getGeneSummary(less_maf)[1:10]
more_maf.genes = getGeneSummary(more_maf)[1:10]
mdt = merge(less_maf.genes[,.(Hugo_Symbol, MutatedSamples)], more_maf.genes[,.(Hugo_Symbol, MutatedSamples)], by = 'Hugo_Symbol', all = TRUE)
mdt$MutatedSamples.x[is.na(mdt$MutatedSamples.x)] = 0
mdt$MutatedSamples.y[is.na(mdt$MutatedSamples.y)] = 0
mdt$max = apply(mdt[,.(MutatedSamples.x, MutatedSamples.y)], 1, max)
mdt = mdt[order(max, decreasing = TRUE)]
mdt
top_10 <- (c(mdt$Hugo_Symbol))
top_10 <- c("TP53",  "PIK3CA" ,"TTN",    "CDH1", "GATA3",  "MUC16" , "KMT2C",  "MAP3K1", "HMCN1",  "FLG")
jpeg("C:/Users/zhang/OneDrive/college/qbio490/qbio_490_lucia/week7_maf/part_ii_cooncoplot.png")
# create a co-oncoplot between the two mafs
coOncoplot(m1 = less_maf, m2 = more_maf, m1Name = "Patients with Less than 15 Lymph Nodes", m2Name = "Patients with Greater than 15 Lymph Nodes", borderCol = NA, genes = top_10) 
# none of the genes had a big discrepency. TTN has a 4% discrepency. TTN is the largest known protein that plays a role in the elasticity and stility of muscle fibers. Mutations can cause muscular disorders. I'm unsure why there is a discrepency. 
maf_object@data$MUC16_status <- ifelse(maf_object@data$Hugo_Symbol == "MUC16", "Mutated", "Normal")
# merge a MAF with MUC16 and the lymph nodes.
merge <- merge(maf_object@data[, c("MUC16_status", "Tumor_Sample_Barcode")], maf_object@clinical.data[,c( "lymph_node_categorical", "Tumor_Sample_Barcode")], by="Tumor_Sample_Barcode")
# create a contingency table
contig <- table(merge$MUC16_status, merge$lymph_node_categorical)
# plot it
mosaicplot(contig)
# run a fischer test
fisher_test <- fisher.test(contig)
fisher_test
# the p value is 0.03763, which is a significant value. the odds ratio measures the association b/t the two variable, and it is at 0.6368. 
lollipopPlot2(m1=less_maf, m2=more_maf, m1_name="Less than 15", m2_name="Greater than 15", gene = "MUC16")
# create a survival status column depending on whether or not a patient is alive
maf_object@clinical.data$Overall_Survival_Status <- ifelse(maf_object@clinical.data$vital_status == "Alive", TRUE, FALSE)
# create a K-M plot
mafSurvival(maf = maf_object, genes = "MUC16", time = "days_to_last_followup", Status = "Overall_Survival_Status", isTCGA = TRUE)
# with a p-value of 0.165, there does not seem to be a differene. I chose MUC16 this time, and which is a glycoprotein. I'm unsure why there is no difference between the two. 

