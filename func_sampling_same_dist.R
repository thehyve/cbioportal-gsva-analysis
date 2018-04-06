#!/usr/bin/Rscript

## Script with function to carry out resampling of genes in a geneset
## Author: Dionne Zaal - The Hyve

# Function carries out the resampling of genes for geneset taking into 
# account the distribution of genes in the original genesets
sample_geneset_from_dist <- function(genesets){
  # Make empty variables for new genesets
  resampled_genesets <- vector("list", length(genesets))
  names(resampled_genesets) <- names(genesets)
  
  # all genes from genesets are needed
  all_genes <- as.vector(unlist(genesets))
    
    # Create new content for gene sets, drawn from the same distribution of original
    # also with the same size as the original gene set
    for (geneset in 1:length(genesets)){
      # Sample genes from list of original genes
      resampled_genesets[geneset][[1]] <- sample(all_genes, length(genesets[geneset][[1]]), replace=FALSE)
      # Check if genes are double in geneset, until there are the geneset is not accepted
      dupl_gene_check = 0
      while(dupl_gene_check != 1){
        # If all genes in the geneset are unique, remove genes used from list,
        # and continue with next geneset
        if (length(resampled_genesets[geneset][[1]]) == length(unique(resampled_genesets[geneset][[1]]))){
          # Remove genes from full list
          all_genes <- all_genes[-match(resampled_genesets[geneset][[1]], all_genes)]
          dupl_gene_check = 1
        } else {
          ## If there are duplicate genes, determine the amount and try to resample for those
          # First only grab the unique genes for the geneset
          resampled_genesets[geneset][[1]] <- unique(resampled_genesets[geneset][[1]])
          diff = length(genesets[geneset][[1]]) - length(resampled_genesets[geneset][[1]])
          # Continue when new resampling has only unique genes and are not in the current geneset
          uniq_genes = 0
          uniq_geneset = 0
          attempt = 0
          while(uniq_genes != 1 && uniq_geneset != 1){
            add_genes <- sample(all_genes, diff, replace=FALSE)
            if (length(unique(add_genes)) == diff){
              uniq_genes = 1
              # Now check if genes do not occur yet in gene set
              check_genes <- match(add_genes, resampled_genesets[geneset][[1]])
              check_na <- check_genes[!is.na(check_genes)]
              if (length(check_na) == 0){
                # This means new genes were found, add to geneset and remove genes from list
                resampled_genesets[geneset][[1]] <- c(resampled_genesets[geneset][[1]], add_genes)
                all_genes <- all_genes[-match(resampled_genesets[geneset][[1]], all_genes)]
                uniq_geneset = 1
                dupl_gene_check = 1
                break
              } else {
                uniq_genes = 0
              }
            } else {
              uniq_genes = 0
            }
            
            attempt = attempt + 1
            # After 10 attempts to find unique geneset just save the current geneset with only
            # unique genes. It can happen that no unique geneset can be found with large genesets
            # but also with the last genesets that are being assessed, if only duplicate genes are 
            # still in the list it is difficult to find a unique gene
            if (attempt == 10){
              resampled_genesets[geneset][[1]] <- unique(c(resampled_genesets[geneset][[1]], add_genes))
              all_genes <- all_genes[-match(resampled_genesets[geneset][[1]], all_genes)]
              uniq_genes = 1
              uniq_geneset = 1
              dupl_gene_check = 1
            }
          }
        }
      }
    }
  return(resampled_genesets)
  }
