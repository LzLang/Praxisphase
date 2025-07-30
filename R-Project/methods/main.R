# Check and install necessary packages
if (!require('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
  BiocManager::install('limma')
  BiocManager::install("edgeR")
  BiocManager::install("clusterProfiler")
  BiocManager::install("pathview")
  BiocManager::install("enrichplot")

if(!require('tidyverse', quietly = TRUE)) {
  install.packages("tidyverse", type="source")

}
  
#library(DESeq2)
library(tidyverse)
library(edgeR)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library("org.Hs.eg.db", character.only = TRUE)

load_data <- function(outlier = TRUE) {
  if (outlier) {
    # load data
    counts <- readTargets(gzfile(paste(getwd(),'/data/merged_data_R_with_outlier.csv.gz', sep = '')), sep = ',', row.names = 'gene_id') %>% 
      .[, -1] %>% .[, order(names(.))]
    metadata <- read_csv(paste(getwd(), '/data/merged_metadata_with_outlier.csv', sep='')) %>%
      .[order(.$'Run'),] %>%
      column_to_rownames(., var="Run")
    ml_features <- read_csv(paste(getwd(), '/data/ml_features_with_outlier.csv', sep=''))
  }
  else {
    # load data
    counts <- readTargets(gzfile(paste(getwd(),"/data/merged_data_R_without_outlier.csv.gz", sep = "")), sep = ',', row.names = 'gene_id') %>% 
      .[, -1] %>% .[, order(names(.))]
    metadata <- read_csv(paste(getwd(), '/data/merged_metadata_without_outlier.csv', sep='')) %>%
      .[order(.$'Run'),] %>%
      column_to_rownames(., var="Run")
    ml_features <- read_csv(paste(getwd(), '/data/ml_features_without_outlier.csv', sep=''))
  }
  
  assign('counts', counts, envir = .GlobalEnv)  
  assign('metadata', metadata, envir = .GlobalEnv)  
  assign('ml_features', ml_features, envir = .GlobalEnv)
}

differential_expression <- function() {
  # Create DGEList object, filter and calculate normalization factors
  # Important: calcNormFactors doesn't normnalize the data
  # -> just calculates normaliztion factors for use downstream
  expression_df <- DGEList(counts) %>% calcNormFactors(.)
  # Set up design matrix
  # 0: Normal
  # 1: Tumor
  design_matrix <- model.matrix(~., data=metadata)
  #keep <- filterByExpr(expression_df, design_matrix) 
  #expression_df <- expression_df[keep,,keep.lib.size=FALSE] %>% calcNormFactors(.)
  
  # Perform differential expression
  ## Use LogCPM -> Use prior counts to damp down the variacnes of logarithms of low counts
  logCPM <- cpm(expression_df, log=TRUE, prior.count=3)
  ## Use lmFit to test each gene for differential expression between the two groups using a linear model
  fit <- lmFit(logCPM, design_matrix)
  ## Apply empirical Bayes smoothing with the eBayes() function
  fit <- eBayes(fit, trend=TRUE)
  ## Apply multiple testing correction and obtain stats
  ### We are using Benjamini-Hochberg methond
  ### By default, results are orderd by largest B (the log odds value)
  ### -> Most differentially expressed genes should be toward the top
  stats_df <- topTable(fit, coef = ncol(design_matrix), n = Inf)
  print(head(stats_df))
  
  ##########
  stats_df$gene_id <- rownames(stats_df)
  stats_df <- stats_df[,c("gene_id", names(stats_df)[1:6])]
  write.csv(stats_df, file = "output.csv")
  ##########
  #### logFC: log2 fold change
  #### AveExpr: Average expression across all samples, in log2 CPM
  #### t: logFC divided by its standard error
  #### P.Value: Raw p-value (based on t) from test that logFC differs from 0
  #### andj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-Value
  #### B: log-odds that gene is DE (arguably less unsefull than the other columns)
  
  # Mean-variance trend
  ## Counts are transformes to log2 counts per million reads (CPM), where "per million reads" is defined based on the normalization factors we calculated earlier
  ## A linear model is fitted to the log2 CPM for each gene, and the residuals are calculated
  ## A smoothed curve is fitted to the sqrt(residual standard deviation) by average expression (see red line in plot above)
  ## The smoothed curve is used to obtain weights for each gene and sample that are passed into limma along with the log2 CPMs.
  ## https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
  
  png(file="voom.png", width=1920, height = 1080)
  y <- voom(expression_df, design_matrix, plot = T)
  dev.off()
  
  # How many DE genes are there?
  #length(which(stats_df$adj.P.Val<0.05))
  # Mit outlier: 3568
  # Ohne: 5829
  return(stats_df)
}

plot_differential_expression <- function(){
  # Check result by plotting gene
  ## In our case ENSG00000142748.13 has the greatest B value
  top_gene_df <- counts %>%
    # Extract this gene from 'expression_df'
    #dplyr::filter(rownames(.)=='ENSG00000142748.13') %>%
    # Transpose so the gene is a column
    t() %>%
    # Transpose made this a matrix, let's make it back into a data frame
    data.frame() %>%
    rownames_to_column('Run') %>%
    dplyr::inner_join(dplyr::select(
      rownames_to_column(metadata, "Run"),
      Run,
      Tumor_Status
    ))
  head(top_gene_df)
  ggplot(top_gene_df, aes(x = Tumor_Status, y=ENSG00000142748.13,color = Tumor_Status,)) +
    geom_jitter(width = 0.2, height = 0) + # We'll make this a jitter plot
    theme_classic() # This makes some aesthetic changes
}

enrichment_analysis <- function(number = 1000, ml_feature = NULL) {
  # Enrichment Analysis
  #Prepare Input
  if (!is.null(ml_feature)){
    if (ml_feature=='logistic_regression') {
      gene_list <- ml_features$logreg_log2fc[1:number]
      names(gene_list) <- ml_features$logistic_regression[1:number]
    }
    else if (ml_feature=='random_forest_classification') {
      gene_list <- ml_features$rf_log2fc[1:number]
      names(gene_list) <- ml_features$random_forest_classification[1:number]
    }
    else stop("Wrong ml_feature selected! You can choose between logistic_regression and random_forest_classification")
  } else {
    ## We want the log2 fold change
    gene_list <- stats_df$logFC
    ## name the vector
    names(gene_list) <- stats_df$gene_id
    ## sort the list in decreasing order (required for clusterProfiler)
    gene_list <- sort(gene_list, decreasing=TRUE)[1:number]
  }
  
  gene_list = sort(gene_list, decreasing = TRUE)
  
  
  
  # Gene Set Enrichment
  gse <- gseGO(geneList = gene_list,
              ont = "ALL",
              keyType = "ENSEMBL",
              minGSSize = 3,
              maxGSSize = 800,
              pvalueCutoff = 0.05,
              verbose = TRUE,
              OrgDb = "org.Hs.eg.db",
              pAdjustMethod = "none")
  if(!is.null(number)) {gseaplot(gse, by="all", paste("Top", number), geneSetID=1)}
  else {gseaplot(gse, by="all", gse$Description[1], geneSetID=1)}
  return(gse)
}

load_data()
stats_df <- differential_expression()
# 
#enrichment_analysis(100)
#enrichment_analysis(100, 'logistic_regression')
#enrichment_analysis(100, 'random_forest_classification')
#gse <- enrichment_analysis(1000)
gse <- enrichment_analysis(1000, 'logistic_regression')
gse <- enrichment_analysis(1000, 'random_forest_classification')
gseaplot(gse, by="all", paste("Top", 1000), geneSetID=1)

#require(DOSE)
#dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

#edox2 <- pairwise_termsim(gse)
#p1 <- treeplot(edox2)
#p2 <- treeplot(edox2, hclust_method="average")
#aplot::plot_list(p1, p2, tag_levels = 'A')

#clusterProfiler::emapplot(edox2, showCategory = 10)


#ridgeplot(gse) + labs(x = "enrichment distribution")
