## Resampling gene sets
Gene sets are resampled with the same size as the initial gene sets. Only genes from the original gene sets are used to fill the resampled gene sets, following the same distribution as genes in the original gene sets. When genes are used in a resampled gene set, they will be removed from the full list of genes and will not be used again (resampling without replacement). Within the same gene sets no duplicate genes are allowed, but the same genes can be used in different gene sets.

## P-value calculation
The p-value for the GSVA score indicates the chance of finding the same or a higher GSVA score (closer to -1 or 1) with a geneset of the same size in the same set of samples. A GSVA score closer to -1 or 1 indicates more concordantly activated genes (in one direction, either over- or under-expressed) than the expression of these genes in the other samples ([Hänzelmann _et al_. (2013)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-7)). This permutation test is done by measuring how often the same or higher score results is calculated for a random gene set of the same size as the original. 

![pvalue calculation](https://raw.githubusercontent.com/thehyve/cbioportal-gsva-analysis/master/GSVA-pvalue-method.jpg)
