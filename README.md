# integration_analysis_scripts
Scripts for multi-omics integration

## Unsupervised analysis: *integration_unsupervised.R*

This script performs unsupervised analyses (clustering) from transformed expression data (e.g., log fpkm) and methylation beta values

## Prerequisites
This R script requires the following packages:
- iClusterPlus
- gplots
- lattice

### Usage
```bash
Rscript integration_unsupervised.R [options]
```

| **PARAMETER** | **DEFAULT** | **DESCRIPTION** |
|-----------|--------------:|-------------| 
*-d* | NULL | File with somatic mutation data |
*-C* | NULL | File with copy number variation data |
*-r* | NULL | File with expression data |
*-m* | NULL | File with methylation data (beta values) |
*-k* | 2 | Minimum number of clusters |
*-K* | 6 | Maximum number of clusters |
*-c* | 2 | Number of cores |
*-o* | out | output prefix |
*-h*    |  | Show help message and exit|

For example, one can type
```bash
Rscript integration_unsupervised.R -r expression_matrix.txt -o output/
```

### Details
The script involves 3 steps
- **Data transformation** of methylation beta values, using the logit function
- **Clustering** across a range of LASSO lambda penalties and for each number of clusters *K* using iClusterPlus
- **Selection** of the best lambda value (BIC) for each *K*, and plot of the R^2 as a function of *K* to help the choice of *K*
- **Selection of the top features** differentiating the clusters

### Output
- A figure with R^2 as a function of *K*, and cluster memberships of each sample as a function of *K*

In addition, for each value of *K*:
- an .RData file with clustering results
- a heatmap with the top features for each dataset
- a .txt file with the name of the top features for each dataset

## Regression analysis for unsupervised analysis: *PCA_regression.R*

This script provides functions to perform regression analysis between variables (e.g., batch variables or clinical variables) and latent factors as obtained by PCA or group factor analysis.


