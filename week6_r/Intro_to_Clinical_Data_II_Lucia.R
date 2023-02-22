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
clinic <- read.csv("C:/Users/zhang/OneDrive/college/qbio490/qbio_490_lucia/analysis_data/brca_clinical_data.csv")
clin_query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", file.type = "xml")
clinical_drug <- GDCprepare_clinic(query = clin_query,clinical.info = "drug") # get clinical_drug data
clinical_rad <- GDCprepare_clinic(query = clin_query, clinical.info = "radiation") # get clinical_rad data
sum(is.na(clinic$histological_type)) # find number of NAs for histological type 
clinic$histological_type # view histological_type column
# chose histological_type
# a discrete variable

sum(is.na(clinical_drug$prescribed_dose)) 
clinical_drug$prescribed_dose
# total_dose is continuous 
# relating variable to each other: certain histological types of cancer will have higher prescribed doses
# relate to breast cancer: certain types of cancer mean a lower survival rate
# relate to breast cancer: higher dose means lower survival rate
clinic_drug_merge <- merge(clinic, clinical_drug, by="bcr_patient_barcode") # merge clinical_drug and clinic so we can plot them
# clinic_drug_merge$histological_type
clinic_drug_merge <- clinic_drug_merge[ifelse(clinic_drug_merge$histological_type == "Other, specify" | clinic_drug_merge$histological_type == "Mixed Histology (please specify)", FALSE, TRUE), ] # remove these two types of histological type because they aren't specific
# str(clinic_drug_merge$prescribed_dose)
clinic_drug_merge$prescribed_dose <- as.numeric(clinic_drug_merge$prescribed_dose) #set the prescribed dose to a numeric (it's currently a factor) so that you can graph as if it were continuous 
jpeg("C:/Users/zhang/OneDrive/college/qbio490/qbio_490_lucia/week6_r/box_plot.jpg")
box_plot <- ggplot(clinic_drug_merge, aes(x=histological_type, y=prescribed_dose))+ geom_boxplot() + scale_x_discrete(guide = guide_axis(n.dodge=2)) # create boxplot
box_plot
dev.off()

# (clinic_drug_merge$prescribed_dose)
clinic_drug_merge$dose_categorical <- ifelse(clinic_drug_merge$prescribed_dose >= 100, "Greater than 100", "Less than 100") #create new column that splits prescribed dose into two categories

# clinic_drug_merge$dose_categorical
clinic_drug_merge$survival_time <- ifelse(is.na(clinic_drug_merge$days_to_death), clinic_drug_merge$survival_time <- clinic_drug_merge$days_to_last_followup, clinic_drug_merge$days_to_death) # if they died, set the survival time to reflect that; otherwise, setset days to last followup because they're still alive
clinic_drug_merge <- clinic_drug_merge[ifelse(clinic_drug_merge$survival_time == "-Inf", FALSE, TRUE), ] # remove -Inf values
survival_object <- Surv(clinic_drug_merge$survival_time, event = ifelse(clinic_drug_merge$vital_status == "Dead", TRUE, FALSE)) # create survival object, set the death event based on vital status
fit_object <- survfit(survival_object~dose_categorical, data = clinic_drug_merge) # set variable to the categorical dose column
survplot <- ggsurvplot(fit_object, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")
jpeg("C:/Users/zhang/OneDrive/college/qbio490/qbio_490_lucia/week6_r/KM_plot_dose.jpg")
KM_plot_dose <- survplot$plot + theme_bw() + theme(axis.title = element_text(size=20), axis.text = element_text(size=16), legend.title = element_text(size=14), legend.text = element_text(size=12)) # create the plot
KM_plot_dose
dev.off()


# repeat for histological type, except you don't need to change it to a categorical variable because it is categorical by default
fit_object <- survfit(survival_object~histological_type, data = clinic_drug_merge)
survplot <- ggsurvplot(fit_object, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")
jpeg("C:/Users/zhang/OneDrive/college/qbio490/qbio_490_lucia/week6_r/KM_plot_histo.jpg")
KM_plot_histo <- survplot$plot + theme_bw() + theme(axis.title = element_text(size=20), axis.text = element_text(size=16), legend.title = element_text(size=14), legend.text = element_text(size=12))
KM_plot_histo
dev.off()
