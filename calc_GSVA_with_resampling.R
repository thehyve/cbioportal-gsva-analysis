#!/usr/bin/Rscript

###############################################################
## Script to perform GSVA analysis for studies in datahub    ##
## Input: expression file (rows: genes, columns: samples),   ##
##		    meta expression file, geneset file (in gmt),       ##
##		    number of cores to use for analysis, number of     ##
## 		    resamplings to use for analysis, prefix for output ##
## 		    files (with desired path)                          ##
## Output: GSVA score and p-value file (genesets in rows,    ##
##		     samples in columns), meta score and p-value file, ##
##         case list file and the Rdata from the analysis	   ##
## Author: Dionne Zaal - The Hyve									           ##
###############################################################

cat("\n\n---> Load R libraries to perform the analysis: parallel, snow, qusage and GSVA\n\n")
library("parallel")  # To be able to do parallel analysis
library("snow")  # To perform parallel analysis
library("qusage") # To read in geneset file (.gmt)
library("GSVA") # To perform GSVA analysis

# Get arguments from command: 
# [1] name (and path) of expression file
# [2] name (and path) of geneset file
# [3] name of geneset version in cBioPortal
# [4] name (and path) of meta expression file
# [5] prefix for name outputfile
# [6] number of cores to use for the analysis
# [7] number of resamplings to do for the analysis
c_args <- commandArgs(TRUE)
expr_file <- c_args[1]
geneset_file <- c_args[2]
geneset_def_version <- c_args[3]
meta_expression_file <- c_args[4]
prefix_out <- c_args[5]
n_cores <- as.numeric(c_args[6])
n_resampling <- as.numeric(c_args[7])

# Change this text in case you are using different genesets
profile_descrip_meta_gsva_scores <- paste0("GSVA scores for a selection of test gene sets from MSigDB v6.1")

# Load and normalize expression file
cat(paste0("\n\n---> Load expression file ", expr_file, " and geneset file ",  geneset_file, "\n\n"))

# Load genesets, geneset file should be in gmt format
genesets <- read.gmt(geneset_file)  # Will put genesets in named list

# Load expression file
expr <- read.delim(expr_file, sep="\t", header=T, row.names=2, quote="")

# Check if we need to remove the hugo or entrez gene columns
if (length(grep("[A-Z]", genesets[1], perl=TRUE, value=TRUE)) == 0) {
  expr <- expr[,-1] # Use in case of entrez gene ids
} else {
  # in case of hugo symbols:
  expr <- expr[!duplicated(expr[,1]),]
  expr <- expr[!is.na(expr[,1]),]
  row.names(expr) <- expr[,1]
  expr <- expr[,-1]
}

# Normalize expression values and only keep genes with at least a minimum expression of 1
expr_norm <- log2(expr + 1)
mexp <- rowMeans(expr_norm)
expr_norm_high <- expr_norm[mexp > 1, ]

# Calculate original gene set scores
cat("\n\n---> Calculate GSVA scores original genesets with ", n_cores," number of cores\n\n")
full_set_gsva_result <- gsva(as.matrix(expr_norm_high), genesets, method="gsva", parallel.sz=n_cores, parallel.type="SOCK")

# Write scores to file
cat(paste0("\n\n---> Output file original GSVA written to ", prefix_out, "gsva_scores.txt \n\n"))
colnames(full_set_gsva_result$es.obs) <- gsub("\\.", "-", colnames(full_set_gsva_result$es.obs))
scores = cbind(geneset_id = rownames(full_set_gsva_result$es.obs), full_set_gsva_result$es.obs)
write.table(scores, paste0(prefix_out, "data_gsva_scores.txt"), quote = F, sep = "\t", col.names = T, row.names = F)

# If the user indicated that resamplings should be done, do resamplings
if (n_resampling > 0){
  cat("\n\n---> Randomize gene sets and calculate resampling scores")
  # Create empty list for new genesets
  new_genesets <- vector("list", length(genesets))
  names(new_genesets) <- names(genesets)
  
  # Create list to add all resampled scores separatly
  list_scores <- vector("list", n_resampling)
  names(list_scores) <- paste0("resampling_", 1:n_resampling)
  all_genes_geneset <- unique(as.vector(unlist(genesets)))
  
  # source function to resample gene sets
  source("func_sampling_same_dist.R")
  
  for (i in 1:n_resampling){
    cat(paste0("\n\n---> Resampling Number: ", i, " \n"))
    # Resample gene sets from the same distribution with function in separate R file
    new_genesets <- sample_geneset_from_dist(genesets)
    
    # Calculate GSVA scores for resampled gene sets
    gene_resamp_gsva <- gsva(as.matrix(expr_norm_high), new_genesets, method="gsva", parallel.sz=n_cores, parallel.type="SOCK")
    results <- data.frame(matrix(NA, nrow=nrow(full_set_gsva_result$es.obs), ncol=ncol(full_set_gsva_result$es.obs)))
    
    # Save 'resampled' scores to list with dataframes (size of dataframes are the same and therefore can complete be added to list)
    rownames(results) = rownames(full_set_gsva_result$es.obs)
    colnames(results) = colnames(full_set_gsva_result$es.obs)
    results[rownames(gene_resamp_gsva$es.obs),] <- gene_resamp_gsva$es.obs
    list_scores[[i]] <- results
  }

  # Create empty dataframe to add p-values
  pvalues_calc <- data.frame(matrix(NA, nrow=nrow(full_set_gsva_result$es.obs), ncol=ncol(full_set_gsva_result$es.obs)))
  rownames(pvalues_calc) = rownames(full_set_gsva_result$es.obs)
  colnames(pvalues_calc) = colnames(full_set_gsva_result$es.obs)
  
  # Calculating p-values with resampled values
  cat("\n\n---> Calculate pvalues")
  for (sample in 1:ncol(pvalues_calc)){
    for (geneset in 1:nrow(pvalues_calc)){
      # Get all scores from each sample and gene sets for the amount of resamplings
      resampled_scores <- sapply(list_scores, '[', geneset, sample) 
      
      # Check how many of the absolute values from the resampled gene sets are better
      # than the calculated score for the original gene set
      # (-0.4 is still stronger than 0.2 but with previous method was not taken into account)
      pvalues_calc[geneset,sample] <- sum(abs(resampled_scores) >= abs(full_set_gsva_result$es.obs[geneset, sample])) / n_resampling

      # PREVIOUS METHOD USED, CHANGED TO CHECKING THE ABSOLUTE VALUES
      ## Positive and negative scores needed separatly for calculation
      #pos_scores <- bootstrap_scores[which(bootstrap_scores >= 0)]
      #neg_scores <- bootstrap_scores[which(bootstrap_scores < 0)]
      #if (full_set_gsva_result$es.obs[geneset, sample] >= 0) {
      #  pvalues_calc[geneset,sample] <- sum(pos_scores >= full_set_gsva_result$es.obs[geneset, sample]) / (length(pos_scores) +1)
      #} else {
      #  pvalues_calc[geneset,sample] <- sum(neg_scores <= full_set_gsva_result$es.obs[geneset, sample]) / (length(neg_scores) +1)
      #}

      ## In case no scores are above (or below in case of negative) the original, set the pvalue to the lowest possible
      #if (pvalues_calc[geneset,sample] == 0){
      #  pvalues_calc[geneset,sample] = 1/n_bootstrap
      #}
    }
  }
  
  # Write p-value results to file
  cat(paste0("\n\n---> Output file resampling written to ", prefix_out, "data_gsva_pvalues.txt\n\n"))
  colnames(pvalues_calc) <- gsub("\\.", "-", colnames(pvalues_calc))
  pvalues_with_genesetid = cbind(geneset_id = rownames(pvalues_calc), pvalues_calc)
  write.table(pvalues_with_genesetid, paste0(prefix_out, "data_gsva_pvalues.txt"), quote = F, sep = "\t", col.names = T, row.names = F)
  
}


### Write meta files
cat(paste0("\n\n---> Create meta datafiles for gsva scores and pvalues, also create new case list"))
# Get versions R and GSVA
r_version <- as.character(getRversion())
gsva_version <- as.character(packageVersion("GSVA"))

# Read in expression meta file
meta_expression <- read.table(meta_expression_file, sep=":", header=FALSE)
study_id <- trimws(as.character(meta_expression[which(meta_expression[,1] == "cancer_study_identifier"), 2]))
source_stable_id <- trimws(as.character(meta_expression[which(meta_expression[,1] == "stable_id"), 2]))

meta_scores <- paste0("cancer_study_identifier: ", study_id, "
genetic_alteration_type: GENESET_SCORE
datatype: GSVA-SCORE
stable_id: gsva_scores
source_stable_id: ", source_stable_id, "
profile_name: GSVA scores
profile_description: ", profile_descrip_meta_gsva_scores, " calculated with GSVA version ", gsva_version,", R version ", r_version, ". See https://github.com/thehyve/cbioportal-gsva-analysis for documentation and R code.
data_filename: data_gsva_scores.txt
show_profile_in_analysis_tab: true
geneset_def_version: ", geneset_def_version)
write(meta_scores, paste0(prefix_out, "meta_gsva_scores.txt"))

# In case resampling is done make meta file which indicates the amount of
# resaplings that were done, otherwise create a dummy pvalue file and
# indicate this in the meta p-value file
if (n_resampling > 0){
  meta_pvalues <- paste0("cancer_study_identifier: ", study_id, "
genetic_alteration_type: GENESET_SCORE
datatype: P-VALUE
stable_id: gsva_pvalues
source_stable_id: gsva_scores
profile_name: GSVA p-values
profile_description: P-values calculated with permutation test (n=", n_resampling, ") based on GSVA scores of random gene sets. See https://github.com/thehyve/cbioportal-gsva-analysis for documentation and R code.
data_filename: data_gsva_pvalues.txt
geneset_def_version: ", geneset_def_version)
  write(meta_pvalues, paste0(prefix_out, "meta_gsva_pvalues.txt"))
} else {
  # create dummy p-values instead of empty file
  dummy_pvalues <- (full_set_gsva_result$es.obs * 0) + 0.01
  gsva_pvalues <- data.frame("geneset_id" = rownames(full_set_gsva_result$es.obs), dummy_pvalues, check.names = F)
  colnames(gsva_pvalues) <- gsub("\\.", "-", colnames(gsva_pvalues))
  write.table(gsva_pvalues, paste0(prefix_out, "data_gsva_pvalues.txt"), quote = F, sep = "\t", col.names = T, row.names = F)
  
  # create also meta p-value file
  meta_pvalues <- paste0("cancer_study_identifier: ", study_id, "
genetic_alteration_type: GENESET_SCORE
datatype: P-VALUE
stable_id: gsva_pvalues
source_stable_id: gsva_scores
profile_name: GSVA p-values
profile_description: Dummy P-values, no resampling done. See https://github.com/thehyve/cbioportal-gsva-analysis for documentation and R code.
data_filename: data_gsva_pvalues.txt
geneset_def_version: ", geneset_def_version)
  write(meta_pvalues, paste0(prefix_out, "meta_gsva_pvalues.txt"))
}

# Create case list file and if necessary also case lists dir
dir.create(file.path(dirname(paste0(prefix_out, "meta_gsva_pvalues.txt")), "case_lists"))

case_list <- paste0(c("cancer_study_identifier: ", study_id, "
stable_id: ", study_id, "_gsva_scores
case_list_name: Samples with GSVA data
case_list_description: All samples with GSVA data
case_list_category: all_cases_with_gsva_data
case_list_ids:    ", paste0(colnames(full_set_gsva_result$es.obs), collapse = "\t")), collapse = "")
write(case_list, paste0(dirname(paste0(prefix_out, "meta_gsva_pvalues.txt")), "/case_lists/cases_GSVA.txt"))

cat(paste0("\n\n---> Meta files written to ", prefix_out, "meta_gsva_scores.txt, ", prefix_out, "meta_gsva_pvalues.txt and ", prefix_out, "case_lists/cases_GSVA.txt\n\n"))

#save.image(file=paste0(prefix_out, "GSVA_calc.RData"))

#cat(paste0("\n\n---> R environment variables written to ", prefix_out, "GSVA_calc.RData\n\n"))

