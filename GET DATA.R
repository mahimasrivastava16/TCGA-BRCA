library(TCGAbiolinks)
library(dplyr)
library(DT)
library(sesameData)
library(tidyverse)
library(SummarizedExperiment)
library(dbplyr)
library(sesame)
library(ggVennDiagram)
library(parallel)
library(doParallel)

#===============================
#----------FUNCTIONS------------
#===============================
# Impute missing values with row mean
impute_row_mean <- function(row) {
  row_mean <- mean(row, na.rm = TRUE)
  row[is.na(row)] <- row_mean
  return(row)
}


missing_threshold <- 0.2

#===============================
#----------QUERY DATA-----------
#===============================

clin_BRCA <- GDCquery_clinic("TCGA-BRCA", "Clinical")

clin_BRCA$deceased <- ifelse(clin_BRCA$vital_status == "Alive", 0, 1)
# create an "overall survival" variable that is equal to days_to_death for dead patients
# and to days_to_last_follow_up for patients who are still alive

clin_BRCA$overall_survival <- ifelse(clin_BRCA$vital_status == "Alive",
                                     clin_BRCA$days_to_last_follow_up,
                                     clin_BRCA$days_to_death)

clin_BRCA_data <- clin_BRCA[c('submitter_id', 'deceased', 'overall_survival')]

dim(clin_BRCA_data)

write.csv(clin_BRCA,'BRCA/BRCA_complete_clinical.csv')

# build a query to retrieve gene expression data ------------
query_BRCA_RNA_TP <- GDCquery(project = 'TCGA-BRCA',
                              data.category = 'Transcriptome Profiling',
                              experimental.strategy = 'RNA-Seq',
                              workflow.type = 'STAR - Counts',
                              
                              access  = 'open')

Output_query_BRCA_RNA_TP <- getResults(query_BRCA_RNA_TP)
BRCA_GE_sample <- Output_query_BRCA_RNA_TP[c('cases.submitter_id')]

# build a query to retrieve Copy Number Variation ------------

query_BRCA_CNV <- GDCquery(project = 'TCGA-BRCA',
                           data.category = 'Copy Number Variation',
                           
                           data.type = 'Gene Level Copy Number',
                           access = 'open')
Output_query_BRCA_CNV <- getResults(query_BRCA_CNV)
BRCA_CNV_sample <- Output_query_BRCA_CNV[c('cases.submitter_id')]


# build a query to retrieve DNA Methylation data ------------

query_BRCA_Meth <- GDCquery(project = 'TCGA-BRCA',
                            data.category = 'DNA Methylation',
                            platform = 'Illumina Human Methylation 450',
                            
                            data.type = 'Methylation Beta Value',
                            access = 'open')


Output_query_BRCA_Meth <- getResults(query_BRCA_Meth)
BRCA_Meth_sample <- Output_query_BRCA_Meth[c('cases.submitter_id')]


# build a query to retrieve miRNA expression data ------------

query_BRCA_ME <- GDCquery(project = 'TCGA-BRCA',
                          data.category = 'Transcriptome Profiling',
                          experimental.strategy = 'miRNA-Seq',
                          workflow.type = 'BCGSC miRNA Profiling',
                          data.type = 'miRNA Expression Quantification',
                          
                          access = 'open')

Output_query_BRCA_ME <- getResults(query_BRCA_ME)
BRCA_ME_sample <- Output_query_BRCA_ME[c('cases.submitter_id')]

# Get COMMON SAMPLES ACROSS ALL OMICS DATA
common_samples <- Reduce(intersect, list(BRCA_GE_sample[[1]], BRCA_CNV_sample[[1]], BRCA_Meth_sample[[1]], BRCA_ME_sample[[1]]))


# Pre-process gene expression data ------------------------------------


GDCdownload(query_BRCA_RNA_TP)
tcga_BRCA_GE <- GDCprepare(query_BRCA_RNA_TP)
dim(tcga_BRCA_GE) # Gets HOW MANY FEATURES THEN SAMPLES
colnames(colData(tcga_BRCA_GE)) # GETS COLUMN NAMES -> BOTH CLINCIAL AND EXPRESSION DATA IS PRESENT
BRCA_matrix_GE <- assay(tcga_BRCA_GE, 'fpkm_unstrand')

# Pre-process -> remove missing values + get high gene variance

missing_percentage_features <- rowMeans(is.na(BRCA_matrix_GE))
selected_features <- which(missing_percentage_features <= missing_threshold)
BRCA_matrix_GE_filtered_features <- BRCA_matrix_GE[selected_features, ]
missing_percentage_samples <- colMeans(is.na(BRCA_matrix_GE_filtered_features))
selected_samples <- which(missing_percentage_samples <= missing_threshold)
BRCA_matrix_GE_filtered <- BRCA_matrix_GE_filtered_features[, selected_samples]


#  By focusing on genes with higher variance, you prioritize those that exhibit more dynamic expression patterns across samples. 
#  This can help identify genes that are likely to be biologically relevant or associated with specific conditions.
BRCA_matrix_GE_filtered <- t(apply(BRCA_matrix_GE_filtered, 1, impute_row_mean))
gene_variances <- apply(BRCA_matrix_GE_filtered, 1, var)
variance_threshold <- quantile(gene_variances, probs = 0.95)
high_var_genes_indices <- which(gene_variances >= variance_threshold)
high_var_genes <- BRCA_matrix_GE_filtered[gene_variances >= variance_threshold, ]


dim(high_var_genes) # GETS HOW MANY GENES REMAIN AFTER FILTERING


BRCA_gene_metadata <- as.data.frame(rowData(tcga_BRCA_GE)) # To get gene name
BRCA_gene_data <- BRCA_gene_metadata[c('gene_id', 'gene_name')]


BRCA_GE <- high_var_genes %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., BRCA_gene_data, by = "gene_id")

# Pre-process -> Take only samples that are shared by all omics data
BRCA_GE$case_id <- gsub('-01.*', '', BRCA_GE$case_id)

BRCA_GE_data<- BRCA_GE %>% 
  filter(case_id %in% common_samples)

# Add clinical information to BRCA_GE_data
BRCA_GE_data <- merge(BRCA_GE_data, clin_BRCA_data, by.x = 'case_id', by.y = 'submitter_id')

# Log (FPKM + 1) the counts
BRCA_GE_data$counts <- log(BRCA_GE_data$counts + 1)


BRCA_RNA <- BRCA_GE_data[,-1]
BRCA_RNA <- BRCA_GE_data[,-2]
BRCA_RNA <- na.omit(BRCA_RNA)

BRCA_RNA <- pivot_wider(BRCA_RNA, names_from = gene_name, values_from = counts, values_fn = mean)

dir_name <- 'BRCA'
if (!dir.exists(dir_name)) {
  dir.create(dir_name)
}

write.csv(BRCA_RNA,'BRCA/BRCA_GE_data.csv')





################# Pre-process Copy Number Variation data ------------------------------------
# Extract the CNV results
Output_query_BRCA_CNV <- getResults(query_BRCA_CNV)
BRCA_CNV_sample <- Output_query_BRCA_CNV[c('cases.submitter_id')]

# Check for duplicate samples
duplicated_samples <- BRCA_CNV_sample[duplicated(BRCA_CNV_sample$cases.submitter_id), ]
cat("Duplicated samples: ", duplicated_samples, "\n")

# Remove duplicate samples from the query results
BRCA_CNV_sample_unique <- BRCA_CNV_sample[!duplicated(BRCA_CNV_sample$cases.submitter_id), ]

# Update the query object with only unique samples 
query_BRCA_CNV$results[[1]] <- Output_query_BRCA_CNV[!duplicated(Output_query_BRCA_CNV$cases.submitter_id), ]

# Now download and prepare the CNV data with unique samples
GDCdownload(query_BRCA_CNV)
tcga_BRCA_CNV <- GDCprepare(query_BRCA_CNV)

# Check dimensions of the prepared data
dim(tcga_BRCA_CNV)

# Filter out the features that had more than 20% missing values across all samples 
# and filtered the samples that had more than 20% missing values across all features

BRCA_matrix_CNV <- assay(tcga_BRCA_CNV, 'copy_number')

# Step 1: Filter features based on missing data threshold
missing_percentage_features <- rowMeans(is.na(BRCA_matrix_CNV))
selected_features <- which(missing_percentage_features <= missing_threshold)
BRCA_matrix_CNV_filtered_features <- BRCA_matrix_CNV[selected_features, ]

# Step 2: Filter samples based on missing data threshold
missing_percentage_samples <- colMeans(is.na(BRCA_matrix_CNV_filtered_features))
selected_samples <- which(missing_percentage_samples <= missing_threshold)
BRCA_matrix_CNV_filtered <- BRCA_matrix_CNV_filtered_features[, selected_samples]

# Step 3: Impute missing values with row mean
BRCA_matrix_CNV_filtered <- t(apply(BRCA_matrix_CNV_filtered, 1, impute_row_mean))

# Step 4: Recalculate gene variances after filtering
gene_variances <- apply(BRCA_matrix_CNV_filtered, 1, var)

# Step 5: Filter high variance genes
variance_threshold <- quantile(gene_variances, probs = 0.95)
high_var_genes_indices <- which(gene_variances >= variance_threshold)

# Step 6: Subset the matrix for high variance genes
BRCA_matrix_CNV_high_var <- BRCA_matrix_CNV_filtered[high_var_genes_indices, ]

# Check the dimensions of the final filtered matrix
dim(BRCA_matrix_CNV_high_var)


write.table(BRCA_matrix_CNV_high_var, 'BRCA/BRCA_CNV_before_clustering_data.csv', col.names = NA, sep = ",")

