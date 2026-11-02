#################
## Make dataset
#################
##
## Descriptuion:  This R script generates the CSV file dataset_cancer_cell.csv, used for data analysis.
##                The source data are available at https://ega-archive.org/studies/EGAS00001007087. 
##                Information about how to request Lifelines data and the conditions of use can be 
##                found at https://www.lifelines-biobank.com/researchers/working-with-us.

# Check and install packages if needed
packages <- c("dplyr", "stringr")

for (pkg in packages) {
  if (!(pkg %in% rownames(installed.packages()))) {
    install.packages(pkg)
  }
}

# Load libraries
library(dplyr)
library(stringr)

## current directory
currDir <- getwd()
setwd(currDir)

## inout data
variants <- read.csv2("./variant_dataset.csv")
pheno <- read.csv2("./phenotype_dataset.csv")

data <- left_join(variants, pheno)

## selecting columns of interest
data <- select(data, indiv_ID, gene, genomic_location, exon, variant_class, cDNA_change, protein_change, age, sex, VAF_baseline, variant_reads_baseline, VAF_FU, variant_reads_FU, VAF_thirdvisit, variant_reads_thirdvisit, date_diff_baseline_FU, included_thirdvisit, date_diff_FU_thirdvisit)

## making total reads variable
data$total_reads_baseline <- round(data$variant_reads_baseline / as.numeric(data$VAF_baseline) *100, digits = 0)
data$total_reads_FU <- round(data$variant_reads_FU / as.numeric(data$VAF_FU) *100, digits = 0)
data$total_reads_thirdvisit <- round(data$variant_reads_thirdvisit / as.numeric(data$VAF_thirdvisit) *100, digits = 0)

## remove variant read variable
data <- select(data, -variant_reads_baseline, -variant_reads_FU, -variant_reads_thirdvisit)

## order colnames
data <- select(data, indiv_ID, gene, genomic_location, exon, variant_class, cDNA_change, protein_change, VAF_baseline, VAF_FU, VAF_thirdvisit, total_reads_baseline, total_reads_FU, total_reads_thirdvisit, date_diff_baseline_FU, included_thirdvisit, date_diff_FU_thirdvisit, age, sex)

## export data
write.csv(data,"./dataset_cancer_cell.csv", row.names = FALSE)


