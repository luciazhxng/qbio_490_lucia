# Import necessary packages
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
# BiocManager::install(version = "3.16")
if (!require("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks")
if (!require("maftools", quietly = TRUE))
  BiocManager::install("maftools")
library(BiocManager)
library(TCGAbiolinks)
library(maftools)
if (!require("EnhancedVolcano")){
  install.packages("EnhancedVolcano") 
}
library(EnhancedVolcano)
if (!require("DESeq2")){
  install.packages("DESeq2") 
}
library(DESeq2)

# Set directory to where data is stored on local computer
# Query and load data

setwd("C:/Users/zhang/OneDrive/college/qbio490/qbio_490_lucia/analysis_data")
clin_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Clinical", 
  file.type = "xml"
)
# GDCdownload(clin_query)
clinic <- GDCprepare_clinic(clin_query, clinical.info = "patient")
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

maf_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Simple Nucleotide Variation", 
  access = "open", # we only have access to somatic mutations which are open access
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

#GDCdownload(maf_query)

maf <- GDCprepare(maf_query) # as long as it runs, ignore any errors

maf_object <- read.maf(maf = maf, clinicalData = clinic, isTCGA = TRUE)


rna_query <- GDCquery(project ="TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")

# GDCdownload(rna_query)

rna_se<- GDCprepare(rna_query)

# Change working directory to new outputs folder
setwd("C:/Users/zhang/OneDrive/college/qbio490/qbio_490_lucia/midterm_project_zhang")
dir.create("outputs")
setwd("outputs")

# Look at pre versus post menopause patients, check the variables
sum(is.na(clinic$menopause_status))
levels(clinic$menopause_status)
clinic$menopause_status
clinic$menopause_status <- as.factor(clinic$menopause_status)
clinic_menopause <- clinic

# Remove those without a status or are indeterminate
indeterminate_mask <- ifelse(clinic_menopause$menopause_status == "" | clinic_menopause$menopause_status == "Indeterminate (neither Pre or Postmenopausal)", F, T)
clinic_menopause<- clinic_menopause[indeterminate_mask, ]

# Count those who are over six months since LMP as "Post" menopause
clinic_menopause$menopause_status <- ifelse(clinic_menopause$menopause_status == "Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)", "Pre", "Post")

# Convert menopause toa factor
clinic_menopause$menopause_status <- as.factor(clinic_menopause$menopause_status)
levels(clinic_menopause$menopause_status)

# Create box plot comparing ages of pre vs. post menopause patients
jpeg("menopause_status.jpg", quality = 100, width=1000, height=700, res=120)
box_plot <- ggplot(clinic_menopause, aes(x=menopause_status, y=age_at_initial_pathologic_diagnosis))+ geom_boxplot() + scale_x_discrete(guide = guide_axis(n.dodge=2)) +  labs(x="Menopause Status", y="Age at Initial Pathologic Diagnosis")
box_plot
dev.off()

# Create death event column using survival time
clinic_menopause$survival_time <- ifelse(is.na(clinic_menopause$days_to_death), clinic_menopause$survival_time <- clinic_menopause$days_to_last_followup, clinic_menopause$days_to_death) 
clinic_menopause <- clinic_menopause[ifelse(clinic_menopause$survival_time == "-Inf", FALSE, TRUE), ] 

# Create Kaplan-Meier plot
survival_object <- Surv(clinic_menopause$survival_time, event = ifelse(clinic_menopause$vital_status == "Dead", TRUE, FALSE))  
fit_object <- survfit(survival_object~menopause_status, data = clinic_menopause)
survplot <- ggsurvplot(fit_object, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")
jpeg("KM_plot_histo.jpg", quality = 100, width=1000, height=400, res=120)
KM_plot_histo <- survplot$plot + theme_bw() + theme(axis.title = element_text(size=20), axis.text = element_text(size=16), legend.title = element_text(size=14), legend.text = element_text(size=12))
KM_plot_histo
dev.off()

# Create new MAF, clean the clinical data the same way we cleaned previous clinical data
maf_menopause <- maf_object
maf_menopause@clinical.data <- maf_menopause@clinical.data[ifelse(maf_menopause@clinical.data$menopause_status == "" | maf_menopause@clinical.data$menopause_status == "Indeterminate (neither Pre or Postmenopausal)", F, T), ]
maf_menopause@clinical.data$menopause_status <- ifelse(maf_menopause@clinical.data$menopause_status == "Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)", "Pre", "Post")

# Subset MAF into a MAF with pre-menopause and a MAF with post-menopause
mask <- ifelse(maf_menopause@clinical.data$menopause_status == "Pre", T, F)
pre_barcodes <- maf_menopause@clinical.data$Tumor_Sample_Barcode[mask]
pre_maf <- subsetMaf(maf = maf_menopause,
                       tsb = pre_barcodes)
maf_menopause@clinical.data$menopause_status

mask <- ifelse(maf_menopause@clinical.data$menopause_status == "Post", T, F)
post_barcodes <-maf_menopause@clinical.data$Tumor_Sample_Barcode[mask]

post_maf <- subsetMaf(maf = maf_menopause,
                     tsb = post_barcodes)

# Create co-oncoplot for the two subsetted MAFs
jpeg("co-oncoplot.jpg", quality = 100, width=1300, height=700, res=120)

cooncoplot <- coOncoplot(m1 = pre_maf, 
           m2 = post_maf, 
           m1Name = "Pre-menopause", 
           m2Name = "Post-menopause", 
           borderCol = NA)
cooncoplot
dev.off()

# Create clinica variable from rna_se, clean it up by removing unnecessary columns
rna_clinical <- rna_se@colData
rna_clinical <- as.data.frame(rna_clinical)
treatments_mask <- !colnames(rna_clinical) %in% c("treatments", "primary_site", "disease_type")
rna_clinical <- rna_clinical[ , treatments_mask]
row.names(rna_clinical) <- rna_clinical[, "barcode"]

# Convert controlled variables to factors
rna_clinical$gender <- as.factor(rna_clinical$gender)
rna_clinical$ajcc_pathologic_stage <- as.factor(rna_clinical$ajcc_pathologic_stage)

# Create rna_genes and rna_counts variables
rna_genes <-rna_se@rowRanges@elementMetadata
rna_genes <- as.data.frame(rna_genes)
row.names(rna_genes) <- rna_genes$gene_id
rna_counts <- rna_se@assays@data$unstranded
rna_counts <- as.data.frame(rna_counts)

# Remove NA values
na_mask <-  ifelse(is.na(rna_clinical$ajcc_pathologic_stage) == TRUE | is.na(rna_clinical$gender == TRUE), F, T)
rna_clinical <-  rna_clinical[na_mask, ]
rna_counts <- rna_counts[, na_mask ]

# Remove RNA counts that have a sum of less than 10 counts
row_sums <- rowSums(rna_counts)
low_counts_mask <- ifelse(row_sums < 10, F, T)
rna_counts <- rna_counts[low_counts_mask, ]
rna_genes <- rna_genes [low_counts_mask, ]

# Add menopause_status column to rna_clinical using previous clinic variable
clinic$patient <- clinic$Tumor_Sample_Barcode
rna_clinical$menopause_status <- clinic$menopause_status[match(rna_clinical$patient, clinic$patient)]
rna_clinical$menopause_status <- as.factor(rna_clinical$menopause_status)

# Clean rna_clinical$menopause_status the same way we cleaned clinic
indeterminate_mask <- ifelse(rna_clinical$menopause_status == "" | rna_clinical$menopause_status == "Indeterminate (neither Pre or Postmenopausal)", F, T)
rna_clinical<- rna_clinical[indeterminate_mask, ]
rna_counts <- rna_counts[, indeterminate_mask]
rna_clinical$menopause_status <- ifelse(rna_clinical$menopause_status == "Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)", "Pre", "Post")
rna_clinical$menopause_status <- as.factor(rna_clinical$menopause_status)

# Make the column and row names more informative
colnames(rna_counts) <- rna_clinical$barcode
rownames(rna_counts) <- rna_genes$gene_id

# Run differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = rna_counts,
                              colData = rna_clinical,
                              design = ~ajcc_pathologic_stage + menopause_status)
dds_obj <- DESeq(dds)


results <- results(dds_obj, format = "DataFrame", contrast = c("menopause_status", "Pre", "Post"))

# Create volcano plot using results from DESeq
# Match gene ids to the gene names in the results
results@listData$gene_names <- rna_genes$gene_name[match(results@rownames, rna_genes$gene_id)]

# Create a new data frame with only the necessary data for the volcano plot
volcano <- data.frame(results@rownames, results@listData$log2FoldChange, results@listData$pvalue, results@listData$padj, -log10(results@listData$padj), results@listData$gene_names)
sum(duplicated(volcano$ensembl))
colnames(volcano) <- c("ensembl", "log2_fold_change", "p_value", "p_adjusted", "-log10_p_adjusted", "gene_names")
row.names(volcano) <- volcano$ensembl
jpeg("volcano_plot.jpeg", quality = 100, width=1000, height=700, res=120)

volcano_plot <- EnhancedVolcano(volcano, title = "Gene Expression in Pre vs. Post Menopause Patients", lab = volcano$gene_names, x= "log2_fold_change", y="p_adjusted", labSize = 2, boxedLabels = TRUE, drawConnectors = TRUE, max.overlaps = 35, shape = 16, colAlpha = 0.5)

volcano_plot
dev.off()
