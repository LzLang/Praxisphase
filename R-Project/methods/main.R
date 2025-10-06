# Check and install necessary packages
if (!require('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
  BiocManager::install('limma')
  BiocManager::install("edgeR")
  BiocManager::install("clusterProfiler")
  BiocManager::install("pathview")
  BiocManager::install("enrichplot")
  BiocManager::install("sva")
  BiocManager::install("GOSim")
  BiocManager::install("biomaRt")
  BiocManager::install("fgsea")

if(!require('tidyverse', quietly = TRUE)) install.packages("tidyverse", type="source")

if(!require("VennDiagram", quietly = TRUE)) install.packages("VennDiagram", type="source")
  
if(!require("msigdbr", quietly = TRUE)) install.packages("msigdbr", type="source")

install.packages("UpSetR")

library(tidyverse)
library(edgeR)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(sva)
library(VennDiagram)
library(GO.db)
library(biomaRt)
library(fgsea)
library(msigdbr)
library(UpSetR)
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
    batch <- read_csv(gzfile(paste(getwd(),'/data/merged_data_R.csv.gz', sep = ''))) %>%
      column_to_rownames(., var="gene_id")
  }
  else {
    # load data
    counts <- readTargets(gzfile(paste(getwd(),"/data/merged_data_R_without_outlier.csv.gz", sep = "")), sep = ',', row.names = 'gene_id') %>% 
      .[, -1] %>% .[, order(names(.))]
    metadata <- read_csv(paste(getwd(), '/data/merged_metadata_without_outlier.csv', sep='')) %>%
      .[order(.$'Run'),] %>%
      column_to_rownames(., var="Run")
    ml_features <- read_csv(paste(getwd(), '/data/ml_features_without_outlier.csv', sep=''))
    batch <- read_csv(gzfile(paste(getwd(),'/data/merged_data_R.csv.gz', sep = ''))) %>%
      column_to_rownames(., var="gene_id")
  }
  
  assign('counts', as.matrix(counts), envir = .GlobalEnv)  
  assign('metadata', metadata, envir = .GlobalEnv)  
  assign('ml_features', ml_features, envir = .GlobalEnv)
  assign('batch', batch, envir = .GlobalEnv)
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
  
  return(stats_df)
}

over_representation <- function(number=1000, database = 'uniprot', ml_feature=NULL, dataframe=stats_df) {
  #Prepare Input
  if (!is.null(ml_feature)){
    if (ml_feature=='logistic_regression') {
      gene_list <- ml_features$logreg_log2fc
      names(gene_list) <- ml_features$logistic_regression
      gene_list <- gene_list[1:number]
    }
    else if (ml_feature=='random_forest_classification') {
      gene_list <- ml_features$rf_log2fc
      names(gene_list) <- ml_features$random_forest_classification
      gene_list <- gene_list[1:number]
    }
    else stop("Wrong ml_feature selected! You can choose between logistic_regression and random_forest_classification")
  } else {
    ## We want the log2 fold change
    gene_list <- dataframe$logFC
    ## name the vector
    names(gene_list) <- dataframe$gene_id
    ## sort the list in decreasing order (required for clusterProfiler)
  }
  
  gene_list <- names(sort(gene_list, decreasing = TRUE)[1:number])
  
  if(database=='uniprot' || database=='entrez'){
    mart <- useMart('ENSEMBL_MART_ENSEMBL')
    mart <- useDataset('hsapiens_gene_ensembl', mart)
    
    biomart_id <- getBM(
      mart = mart,
      attributes = c('entrezgene_id', 'uniprot_gn_id'),
      filter = 'ensembl_gene_id',
      values = gene_list,
      uniqueRows = TRUE)
    
    if(database=='uniprot') {
      gene_list <- biomart_id$uniprot_gn_id
      keytype <- 'uniprot'
    } else if(database=="entrez"){
      gene_list <- biomart_id$entrezgene_id
      keytype <- 'ncbi-geneid'
    }
    return(enrichKEGG(
      gene = gene_list,
      organism = "hsa",
      keyType = keytype,
      pvalueCutoff = 0.05,
      pAdjustMethod = 'BH',
      universe = names(gene_list),
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    ))
  } else if(database=='GO'){
    return(enrichGO(
      gene = gene_list,
      OrgDb = 'org.Hs.eg.db',
      keyType = 'ENSEMBL',
      ont = 'ALL',
      pvalueCutoff = 0.05,
      pAdjustMethod = 'BH',
      universe = names(gene_list),
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 0.2
    ))
  } else {
    print("Wrong Database!")
    return(NULL)
  }
}

fora_analysis <- function(number=1000, ml_feature=NULL, dataframe = stats_df, name_vector=NULL) {
  #Prepare Input
  if (!is.null(ml_feature)){
    if (ml_feature=='logistic_regression') {
      gene_list <- ml_features$logreg_log2fc
      names(gene_list) <- ml_features$logistic_regression
      gene_list <- gene_list[1:number]
      universe <- ml_features$logistic_regression
    }
    else if (ml_feature=='random_forest_classification') {
      gene_list <- ml_features$rf_log2fc
      names(gene_list) <- ml_features$random_forest_classification
      gene_list <- gene_list[1:number]
      universe <- ml_features$random_forest_classification
    }
    else stop("Wrong ml_feature selected! You can choose between logistic_regression and random_forest_classification")
  } else if(is.null(name_vector)){
    gene_list <- dataframe$logFC
    ## name the vector
    names(gene_list) <- dataframe$gene_id
    universe <- dataframe$gene_id
  } else {
    gene_list <- name_vector
    universe <- dataframe$gene_id
  }
  
  if(is.null(name_vector)) gene_list <- names(sort(gene_list, decreasing = TRUE)[1:number])
  
  gene_sets <- gene_set_df %>% split(x = .$ensembl_gene, f = .$gs_name)
  
  return(fora(genes = gene_list, universe = universe, pathways = gene_sets))
  
 
}

do_enrichment <- function(number=1000, ml_feature=NULL, dataframe = stats_df, name_vector=NULL, name=NULL) {
  if(is.null(name)) throw('Filename vergessen')
  else print(name)
  #Prepare Input
  if (!is.null(ml_feature)){
    if (ml_feature=='logistic_regression') {
      gene_list <- ml_features$logreg_log2fc
      names(gene_list) <- ml_features$logistic_regression
      gene_list <- gene_list[1:number]
      universe <- ml_features$logistic_regression
    }
    else if (ml_feature=='random_forest_classification') {
      gene_list <- ml_features$rf_log2fc
      names(gene_list) <- ml_features$random_forest_classification
      gene_list <- gene_list[1:number]
      universe <- ml_features$random_forest_classification
    }
    else stop("Wrong ml_feature selected! You can choose between logistic_regression and random_forest_classification")
  } else if(is.null(name_vector)){
    gene_list <- dataframe$logFC
    ## name the vector
    names(gene_list) <- dataframe$gene_id
    universe <- dataframe$gene_id
  } else {
    gene_list <- name_vector
    universe <- dataframe$gene_id
  }
  
  if(is.null(name_vector)) gene_list <- names(sort(gene_list, decreasing = TRUE)[1:number])
  
  mart <- useMart('ENSEMBL_MART_ENSEMBL')
  mart <- useDataset('hsapiens_gene_ensembl', mart)
    
  biomart_id <- getBM(
    mart = mart,
    attributes = 'entrezgene_id',
    filter = 'ensembl_gene_id',
    values = gene_list,
    uniqueRows = TRUE)
  universe_id <- getBM(
    mart = mart,
    attributes = 'entrezgene_id',
    filter = 'ensembl_gene_id',
    values = universe,
    uniqueRows = TRUE)
    
  gene_list <- as.character(biomart_id$entrezgene_id)
  universe <- as.character(universe_id$entrezgene_id)
  
  edo <- enrichDO(
    gene = gene_list,
    ont = 'HDO',
    organism = 'hsa',
    pvalueCutoff = 0.05,
    pAdjustMethod = 'BH',
    universe = universe,
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.2,
    readable = FALSE
  )
  
  plot <- barplot(edo, showCategory=20, x='Count', color='p.adjust', title=paste('DO Enrichment Analysis - ', name, sep=''), font.size=8)
  ggsave(filename = paste(getwd(),'/output/DO/', name, '.png', sep = ''))
  
  return(edo)
}

batch_correction <- function(){
  ## Setting up the data
  expression_df <- DGEList(counts) %>% calcNormFactors(.)
  logCPM <- cpm(expression_df, log=TRUE, prior.count=3)
  mod = model.matrix(~as.factor(Tumor_Status), data=metadata)
  mod0 = model.matrix(~1,data=metadata)
  
  ## Applying the sva function to estimate batch and other artifacts
  nsv = num.sv(logCPM, mod, method="leek")
  svobj = sva(logCPM, mod, mod0, n.sv=nsv)
  
  ## Adjusting for surrogate variables using the f.pvalue function
  pValues = f.pvalue(logCPM, mod, mod0)
  qValues = p.adjust(pValues, method = "BH")
  modSV = cbind(mod, svobj$sv)
  mod0SV = cbind(mod0, svobj$sv)
  pValuesSV = f.pvalue(as.matrix(logCPM), modSV, mod0SV)
  qValuesSV = p.adjust(pValuesSV, method = "BH")
  
  ## Adjusting for surrogate variables using the limma package
  fit = lmFit(logCPM, modSV)
  eb <- eBayes(fit)
  sv_df <- topTable(eb, coef=2, adjust.method = "BH", number = Inf)
  sv_df$gene_id <- rownames(sv_df)
  sv_df <- sv_df[,c("gene_id", names(sv_df)[1:6])]
  
  return(sv_df)
}

fora_common_pathways <- function(batch_cor, pathway_count = 100, gene_count = 1000, name='') {
  foraDGE <- fora_analysis(gene_count)
  foraLogreg <- fora_analysis(gene_count, ml_feature = 'logistic_regression')
  foraForest <- fora_analysis(gene_count, ml_feature = 'random_forest_classification')
  foraBatch <- fora_analysis(gene_count, dataframe = batch_cor)
  
  foraDGEPathways <- head(foraDGE$pathway[order(foraDGE$padj, decreasing = FALSE)], pathway_count)
  foraBatchPathways <- head(foraBatch$pathway[order(foraBatch$padj, decreasing = FALSE)], pathway_count)
  foraLogregPathways <- head(foraLogreg$pathway[order(foraLogreg$padj, decreasing = FALSE)], pathway_count)
  foraForestPathways <- head(foraForest$pathway[order(foraForest$padj, decreasing = FALSE)], pathway_count)
  
  commonPathways <- Reduce(intersect, list(
    foraDGEPathways,
    foraLogregPathways,
    foraForestPathways,
    foraBatchPathways
  ))  %>% c(., rep(NA, pathway_count-length(.)))
  
  listInput <- list(DGE = foraDGEPathways, LogReg = foraLogregPathways, RandForest = foraForestPathways, BatchCor = foraBatchPathways)
  
  venn.diagram(
    x = listInput,
    category.names = c("DGE", "LogReg", "RandForest", "BatchCor"),
    filename = paste(getwd(), '/output/Pathways_venn_diagram_',name, '.png', sep=''),
    main = paste(name, sep=''),#+"Pathway-Level Overlap - ", 
    output = FALSE,
    disable.logging = TRUE,
    
    main.cex = 2.1,
    main.fontface = "bold",
    
    # Numbers
    cex = 2,
    fontface = "bold",
    fontfamily = "sans",
    
    
    # Set names
    cat.cex = 1.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.fontfamily = "sans"
  )
  
  plot <- upset(fromList(listInput), sets = c("DGE", "LogReg", "RandForest", "BatchCor"), order.by = "freq")
  pdf(file = paste(getwd(),'/output/upset/', name, '.pdf', sep = ''), onefile = FALSE)
  print(plot)
  dev.off() 
  
  
  output <- list(
    DGEPathways = foraDGEPathways,
    Batch = foraBatchPathways,
    LogregPathways = foraLogregPathways,
    ForestPathways = foraForestPathways,
    CommonPathways = commonPathways
  )
  
  return(output)
}

common_pathways <- function(batch_cor, pathway_count = 100, gene_count = 1000, database=NULL, name=''){
  ####################### Export the patchways
  if(database=='uniprot') db<-'uniprot'
  else if(database=="entrez") db<-'entrez'
  else return(NULL)
  gseDGE <- over_representation(gene_count, db)
  gseLogreg <- over_representation(gene_count, db, 'logistic_regression')
  gseForest <- over_representation(gene_count, db, 'random_forest_classification')
  gseBatch <- over_representation(gene_count, db, dataframe = batch_cor)
  
  ## Pathways
  gseDGEPathways <- head(gseDGE@result[order(gseDGE@result$p.adjust, decreasing = FALSE), c("Description", "ID", "category", "p.adjust")], pathway_count)
  gseLogregPathways <- head(gseLogreg@result[order(gseLogreg@result$p.adjust, decreasing = FALSE), c("Description", "ID", "category", "p.adjust")], pathway_count)
  gseForestPathways <- head(gseForest@result[order(gseForest@result$p.adjust, decreasing = FALSE), c("Description", "ID", "category", "p.adjust")], pathway_count)
  gseBatchPathways <- head(gseBatch@result[order(gseBatch@result$p.adjust, decreasing = FALSE), c("Description", "ID", "category", "p.adjust")], pathway_count)
  
  
  commonPathways <- Reduce(intersect, list(
    paste(gseDGEPathways$Description, gseDGEPathways$ID, gseDGEPathways$category, sep="||"),
    paste(gseLogregPathways$Description, gseLogregPathways$ID, gseLogregPathways$category, sep="||"),
    paste(gseForestPathways$Description, gseForestPathways$ID, gseForestPathways$category, sep="||"),
    paste(gseBatchPathways$Description, gseBatchPathways$ID, gseBatchPathways$category, sep="||")
  )) %>% c(., rep("NA||NA||NA", pathway_count-length(.)))
  
  listInput = list(DGE = gseDGEPathways$ID, LogReg = gseLogregPathways$ID, RandForest = gseForestPathways$ID, BatchCor = gseBatchPathways$ID)
  
  venn.diagram(
    x = listInput,
    category.names = c("DGE", "LogReg", "RandForest", "BatchCor"),
    filename = paste(getwd(), '/output/Pathways_venn_diagram_', name, '.png', sep=''),
    main = paste(name, sep=''), #"Pathway-Level Overlap - ", 
    output = FALSE,
    disable.logging = TRUE,
    
    main.cex = 2.1,
    main.fontface = "bold",
    
    # Numbers
    cex = 2,
    fontface = "bold",
    fontfamily = "sans",
    
    
    # Set names
    cat.cex = 1.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.fontfamily = "sans"
  )
  
  plot <- upset(fromList(listInput), order.by = "freq")
  ggsave(filename = paste(getwd(),'/output/upset/', name, '.png', sep = ''))
  
  output <- list(
    DGEPathways = gseDGEPathways,
    BatchPathways = gseBatchPathways,
    ForestPathways = gseForestPathways,
    LogregPathways = gseLogregPathways,
    CommonPathways = commonPathways
  )
  
  return(output)
}

common_features<- function(batch_cor, number=1000, name= '') {
  gene_list_features <- stats_df$logFC
  names(gene_list_features) <- stats_df$gene_id
  gene_list_features <- names(sort(gene_list_features, decreasing = TRUE)[1:number])
  
  BatchCor <- batch_cor$logFC
  names(BatchCor) <- batch_cor$gene_id
  BatchCor <- names(sort(BatchCor, decreasing = TRUE)[1:number])
  logreg <- ml_features$logreg_log2fc[1:number]
  names(logreg) <- ml_features$logistic_regression[1:number]
  logreg <- names(sort(logreg, decreasing = TRUE))
  randforst <- ml_features$rf_log2fc[1:number]
  names(randforst) <- ml_features$random_forest_classification[1:number]
  randforst <- names(sort(randforst, decreasing = TRUE))
  
  commonFeatures <- Reduce(intersect, list(gene_list_features, logreg, randforst, BatchCor))
  
  listInput <- list(DGE = gene_list_features, LogReg = logreg, RandForest = randforst, BatchCor = BatchCor)
  venn.diagram(
    x = listInput,
    category.names = c("DGE", "LogReg", "RandForest", "BatchCor"),
    filename = paste(getwd(), '/output/Features_venn_diagram_',name, '.png', sep=''),
    main = paste(name, sep=''), #"Gene-Level Overlap - ", 
    output = FALSE,
    disable.logging = TRUE,
    
    main.cex = 2.1,
    main.fontface = "bold",
    
    # Numbers
    cex = 2,
    fontface = "bold",
    fontfamily = "sans",
    
    
    # Set names
    cat.cex = 1.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.fontfamily = "sans"
  )
  
  plot <- upset(fromList(listInput), order.by = "freq")
  ggsave(filename = paste(getwd(),'/output/upset/', name, '.png', sep = ''))
  
  return(commonFeatures)
}

common_Intersection <- function(commonWith, commonWithout, method="kegg"){
  max_len <- max(
    length(commonWith),
    length(commonWithout)
  )
  length(commonWith) <- max_len
  length(commonWithout) <- max_len
  
  commonIntersection = Reduce(intersect, list(commonWith, commonWithout))
  if(length(commonIntersection)==0){
    if(method=="kegg") commonIntersection <- rep("NA||NA||NA", max_len-length(commonIntersection))
    else commonIntersection <- rep(NA, max_len-length(commonIntersection))
  } 
  return(commonIntersection)
}

get_ontologie <- function(dataframe) {
  dataframe <- dataframe %>% 
    rowwise() %>% 
    mutate(
      Ancestor = paste(
        switch(
          ONTOLOGY, 
          BP = as.character(GOBPANCESTOR[[ID]]), 
          MF = as.character(GOMFANCESTOR[[ID]]), 
          CC = as.character(GOCCANCESTOR[[ID]]), 
          NA
        ), collapse = ", "),
    Parent = paste(
      switch(
        ONTOLOGY, 
        BP = as.character(GOBPPARENTS[[ID]]), 
        MF = as.character(GOMFPARENTS[[ID]]), 
        CC = as.character(GOCCPARENTS[[ID]]), 
        NA
      ), collapse = ", "),
    Children = paste(
      switch(
        ONTOLOGY, 
        BP = as.character(GOBPCHILDREN[[ID]]), 
        MF = as.character(GOMFCHILDREN[[ID]]), 
        CC = as.character(GOCCCHILDREN[[ID]]), 
        NA
      ), collapse = ", ")
    )
  
  return(dataframe)
}

get_ontology_counts <- function(ontology_list){
  vector_counts <- unlist(strsplit(ontology_list, ", "))
  vector_counts <- vector_counts[!is.na(vector_counts) & vector_counts != "" & vector_counts != "all" & vector_counts != "NA"]
  ontology_counts <- table(vector_counts)
  if(length(ontology_counts)==0) return(data.frame(GO_ID=NA, Count=NA))
  else {
    counts_df <- as.data.frame(ontology_counts)
    colnames(counts_df) <- c("GO_ID", "Count")
    counts_df <- counts_df[order(-counts_df$Count), ]
    
    return(counts_df)
  }
}

export_common <- function(commonWith, commonWithout, commonIntersection, BooleanPathways, name="Common", method="kegg") {
  uncommonCommon <- setdiff(union(commonWith, commonWithout), commonIntersection)
  uncommonCommonWith <- setdiff(commonWith, commonIntersection)
  uncommonCommonWithout <- setdiff(commonWithout, commonIntersection)
  
  max_len <- max(
    length(commonWith),
    length(commonWithout),
    length(commonIntersection),
    length(uncommonCommon),
    length(uncommonCommonWith),
    length(uncommonCommonWithout)
  )
  if(length(uncommonCommon)==0 || length(uncommonCommonWith)==0 || length(uncommonCommonWithout)==0){
    if(method=="kegg"){
      uncommonCommon <- rep("NA||NA||NA", max_len-length(uncommonCommon))
      uncommonCommonWith <- rep("NA||NA||NA", max_len-length(uncommonCommonWith))
      uncommonCommonWithout <- rep("NA||NA||NA", max_len-length(uncommonCommonWithout))
    } else {
      uncommonCommon <- rep(NA, max_len-length(uncommonCommon))
      uncommonCommonWith <- rep(NA, max_len-length(uncommonCommonWith))
      uncommonCommonWithout <- rep(NA, max_len-length(uncommonCommonWithout))
    }
  }
  length(commonWith) <- max_len
  length(commonWithout) <- max_len
  length(commonIntersection) <- max_len
  length(uncommonCommon) <- max_len
  length(uncommonCommonWith) <- max_len
  length(uncommonCommonWithout) <- max_len
  
  if(BooleanPathways && method=="kegg") {
    commonWith <- data.frame(do.call(rbind, strsplit(commonWith, "\\|\\|")))
    colnames(commonWith) <- c("Description", "ID", "ONTOLOGY")
    commonWithout <- data.frame(do.call(rbind, strsplit(commonWithout, "\\|\\|")))
    colnames(commonWithout) <- c("Description", "ID", "ONTOLOGY")
    commonIntersection <- data.frame(do.call(rbind, strsplit(commonIntersection, "\\|\\|")))
    colnames(commonIntersection) <- c("Description", "ID", "ONTOLOGY")
    uncommonCommon <- data.frame(do.call(rbind, strsplit(uncommonCommon, "\\|\\|")))
    colnames(uncommonCommon) <- c("Description", "ID", "ONTOLOGY")
    uncommonCommonWith <- data.frame(do.call(rbind, strsplit(uncommonCommonWith, "\\|\\|")))
    colnames(uncommonCommonWith) <- c("Description", "ID", "ONTOLOGY")
    uncommonCommonWithout <- data.frame(do.call(rbind, strsplit(uncommonCommonWithout, "\\|\\|")))
    colnames(uncommonCommonWithout) <- c("Description", "ID", "ONTOLOGY")
    
    commonWith <- get_ontologie(commonWith) %>% mutate(across(everything(), ~na_if(., "NA")))
    commonWithout <- get_ontologie(commonWithout) %>% mutate(across(everything(), ~na_if(., "NA")))
    commonIntersection <- get_ontologie(commonIntersection) %>% mutate(across(everything(), ~na_if(., "NA")))
    uncommonCommon <- get_ontologie(uncommonCommon) %>% mutate(across(everything(), ~na_if(., "NA")))
    uncommonCommonWith <- get_ontologie(uncommonCommonWith) %>% mutate(across(everything(), ~na_if(., "NA")))
    uncommonCommonWithout <- get_ontologie(uncommonCommonWithout) %>% mutate(across(everything(), ~na_if(., "NA")))
  }
  
  output <- list(
    CommonWith = commonWith,
    CommonWithout = commonWithout,
    CommonIntersection = commonIntersection,
    Uncommon = uncommonCommon,
    UncommonWith = uncommonCommonWith,
    UncommonWithout = uncommonCommonWithout
  )
  write.csv(output, file=paste(getwd(),'/output/', name, '.csv', sep = ''), row.names = F)
  
  if(BooleanPathways) return(output)
}

export_ontology_counts <- function(dataframe, name="Counts"){
  cols <- list(
    CommonWithAncestor=dataframe$CommonWith$Ancestor, CommonWithParent=dataframe$CommonWith$Parent, CommonWithChildren=dataframe$CommonWith$Children,
    CommonWithoutAncestor=dataframe$CommonWithout$Ancestor, CommonWithoutParent=dataframe$CommonWithout$Parent, CommonWithoutChildren=dataframe$CommonWithout$Children,
    CommonIntersectionAncestor=dataframe$CommonIntersection$Ancestor, CommonIntersectionParent=dataframe$CommonWith$Parent, CommonIntersectionChildren=dataframe$CommonWith$Children,
    UncommonAncestor=dataframe$Uncommon$Ancestor, UncommonParent=dataframe$Uncommon$Parent, UncommonChildren=dataframe$Uncommon$Children,
    UncommonWithAncestor=dataframe$UncommonWith$Ancestor, UncommonWithParent=dataframe$UncommonWith$Parent, UncommonWithChildren=dataframe$UncommonWith$Children,
    UncommonWithoutAncestor=dataframe$UncommonWithout$Ancestor, UncommonWithoutParent=dataframe$UncommonWithout$Parent, UncommonWithoutChildren=dataframe$UncommonWithout$Children
  )
  
  # apply over all rows
  results <- lapply(cols, get_ontology_counts)
  
  # set name
  #names(results) <- gsub("\\.", "", cols)
  
  max_len <- max(sapply(results, nrow))
  
  results <- lapply(results, function(dataframe, max_len){
    n_missing <- max_len - nrow(dataframe)
    if(n_missing >0) {
      dataframe <- rbind(dataframe, data.frame(
        GO_ID = rep(NA, n_missing),
        Count = rep(NA, n_missing)
      ))
    }
    return(dataframe)
  }, max_len=max_len)
  
  write.csv(results, file=paste(getwd(),'/output/', name, '.csv', sep = ''), row.names = F)
  return(results)
}


####################


gene_number <- 1000
pathway_number <- 100
top_pathway <- 50
gene_set_df <- msigdbr(species = 'Homo sapiens')

### Preparation
## With Outlier
load_data()
stats_df <- differential_expression()
# Gene List
gene_list_with <- stats_df$logFC
names(gene_list_with) <- stats_df$gene_id
gene_list_with <- names(sort(gene_list_with, decreasing = TRUE)[1:gene_number])
# Batch Correction
batch_cor_with <- batch_correction()
BatchCorWith <- batch_cor_with$logFC
names(BatchCorWith) <- batch_cor_with$gene_id
BatchCorWith <- names(sort(BatchCorWith, decreasing = TRUE)[1:gene_number])
# Logistic Regression
logreg_with <- ml_features$logreg_log2fc[1:gene_number]
names(logreg_with) <- ml_features$logistic_regression[1:gene_number]
logreg_with <- names(sort(logreg_with, decreasing = TRUE))
# Random Forest
randforst_with <- ml_features$rf_log2fc[1:gene_number]
names(randforst_with) <- ml_features$random_forest_classification[1:gene_number]
randforst_with <- names(sort(randforst_with, decreasing = TRUE))
# Common Pathways
commonPathwaysWith <- fora_common_pathways(batch_cor_with, pathway_number, gene_number, 'Pathway-Level - With Outlier')
# Common Features
commonFeaturesWith <- common_features(batch_cor_with, name='Feature-Level - With Outlier')

## Without Outlier
load_data(FALSE)
stats_df <- differential_expression()
# Gene List
gene_list_without <- stats_df$logFC
names(gene_list_without) <- stats_df$gene_id
gene_list_without <- names(sort(gene_list_without, decreasing = TRUE)[1:gene_number])
# Batch Correction
batch_cor_without <- batch_correction()
BatchCorWithout <- batch_cor_without$logFC
names(BatchCorWithout) <- batch_cor_without$gene_id
BatchCorWithout <- names(sort(BatchCorWithout, decreasing = TRUE)[1:gene_number])
# Logistic Regression
logreg_without <- ml_features$logreg_log2fc[1:gene_number]
names(logreg_without) <- ml_features$logistic_regression[1:gene_number]
logreg_without <- names(sort(logreg_without, decreasing = TRUE))
# Random Forest
randforst_without <- ml_features$rf_log2fc[1:gene_number]
names(randforst_without) <- ml_features$random_forest_classification[1:gene_number]
randforst_without <- names(sort(randforst_without, decreasing = TRUE))
# Common Pathways
commonPathwaysWithout <- fora_common_pathways(batch_cor_without, pathway_number, gene_number, 'Pathway-Level - Without Outlier')
# Common Features
commonFeaturesWithout <- common_features(batch_cor_without, name='Feature-Level - Without Outlier')

##########

## Pathway and Feature Intersection
commonPathwayIntersection <- common_Intersection(commonPathwaysWith$CommonPathways, commonPathwaysWithout$CommonPathways, method = "fora")
commonFeatureIntersection <- common_Intersection(commonFeaturesWith, commonFeaturesWithout)

## Export
outputPathways = export_common(commonPathwaysWith$CommonPathways, commonPathwaysWithout$CommonPathways, commonPathwayIntersection, BooleanPathways = TRUE, name='Common_Pathways', method = "fora")
export_common(commonFeaturesWith, commonFeaturesWithout, commonFeatureIntersection, BooleanPathways = FALSE, name='Common_Features')
Genes <- data.frame(
  DGEWith = gene_list_with,
  BatchWith = BatchCorWith,
  LogRegWith = logreg_with,
  RandForestWith = randforst_with,
  DGEWithout = gene_list_without,
  BatchWithout = BatchCorWithout,
  LogRegWithout = logreg_without,
  RandForestWithout = randforst_without
)
write.csv(Genes, file=paste(getwd(),'/output/', 'Genes', '.csv', sep = ''), row.names = F)

####################

### Downstream Analysis
# Prepare and Count Features
all_features_with <- c(gene_list_with, BatchCorWith, logreg_with, randforst_with)
all_features_without <- c(gene_list_without, BatchCorWithout, logreg_without, randforst_without)
all_features <- union(all_features_with, all_features_without)
counts_with <- table(all_features_with)
counts_without <- table(all_features_without)

get_unique <- function(lst, counts) {
  return(lst[counts[lst] == 1])
}

# Unique Features per method
DGEWith <- get_unique(gene_list_with, counts_with)
DGEWithout <- get_unique(gene_list_without, counts_without)
BatchWith <- get_unique(BatchCorWith, counts_with)
BatchWithout <- get_unique(BatchCorWithout, counts_without)
LogRegWith <- get_unique(logreg_with, counts_with)
LogRegWithout <- get_unique(logreg_without, counts_without)
RandWith <- get_unique(randforst_with, counts_with)
RandWithout <- get_unique(randforst_without, counts_without)

# Overlaps and Differences
common_with <- Reduce(intersect, list(gene_list_with, BatchCorWith, logreg_with, randforst_with))
common_without <- Reduce(intersect, list(gene_list_without, BatchCorWithout, logreg_without, randforst_without))
complete_overlap <- intersect(common_with, common_without)
complete_difference <- setdiff(union(common_with, common_without), complete_overlap)
DGE_overlap <- intersect(DGEWith, DGEWithout)
DGE_difference <- setdiff(union(DGEWith, DGEWithout), DGE_overlap)
BatchCor_overlap <- intersect(BatchWith, BatchWithout)
BatchCor_difference <- setdiff(union(BatchWith, BatchWithout), BatchCor_overlap)
LogReg_overlap <- intersect(LogRegWith, LogRegWithout)
LogReg_difference <- setdiff(union(LogRegWith, LogRegWithout), LogReg_overlap)
RandForest_overlap <- intersect(RandWith, RandWithout)
RandForest_difference <- setdiff(union(RandWith, RandWithout), RandForest_overlap)

##########

## ORA - Overlap/Differemce between the methods
#ora_complete_overlap <- fora_analysis(length(complete_overlap), name_vector = complete_overlap)
#pathways_complete_overlap <- head(ora_complete_overlap$pathway[order(ora_complete_overlap$padj, decreasing = FALSE)], 50)
#ora_complete_difference <- fora_analysis(length(complete_difference), name_vector = complete_difference)
#pathways_complete_difference <- head(ora_complete_difference$pathway[order(ora_complete_difference$padj, decreasing = FALSE)], 50)
#ora_DGE_difference <- fora_analysis(length(DGE_difference), name_vector = DGE_difference)
#pathways_DGE_difference <- head(ora_DGE_difference$pathway[order(ora_DGE_difference$padj, decreasing = FALSE)], 50)
#ora_BatchCor_difference <- fora_analysis(length(BatchCor_difference), name_vector = BatchCor_difference)
#pathways_BatchCor_difference <- head(ora_BatchCor_difference$pathway[order(ora_BatchCor_difference$padj, decreasing = FALSE)], 50)
#ora_LogReg_overlap <- fora_analysis(length(LogReg_overlap), name_vector = LogReg_overlap)
#pathways_LogReg_overlap <- head(ora_LogReg_overlap$pathway[order(ora_LogReg_overlap$padj, decreasing = FALSE)], 50)
#ora_LogReg_difference <- fora_analysis(length(LogReg_difference), name_vector = LogReg_difference)
#pathways_LogReg_difference <- head(ora_LogReg_difference$pathway[order(ora_LogReg_difference$padj, decreasing = FALSE)], 50)
#ora_RandForest_overlap <- fora_analysis(length(RandForest_overlap), name_vector = RandForest_overlap)
#pathways_RandForest_overlap <- head(ora_RandForest_overlap$pathway[order(ora_RandForest_overlap$padj, decreasing = FALSE)], 50)
#ora_RandForest_difference <- fora_analysis(length(RandForest_difference), name_vector = RandForest_difference)
#pathways_RandForest_difference <- head(ora_RandForest_difference$pathway[order(ora_RandForest_difference$padj, decreasing = FALSE)], 50)

#pathway_output <- data.frame(
#  complete_overlap = pathways_complete_overlap,
#  complete_difference = pathways_complete_difference,
#  DGE_difference = pathways_DGE_difference,
#  BatchCor_difference = pathways_BatchCor_difference,
#  LogReg_overlap = pathways_LogReg_overlap,
#  LogReg_difference = pathways_LogReg_difference,
#  RandForest_overlap = pathways_RandForest_overlap,
#  RandForest_difference = pathways_RandForest_difference
#)
#write.csv(pathway_output, file=paste(getwd(),'/output/', 'pathways', '.csv', sep = ''), row.names = F)

##########

## DO Analysis - Per Method
# With Outlier
#DO_DGE_with <- do_enrichment(1000, name = 'DGE_with')
#DO_BatchCor_with <- do_enrichment(length(BatchCorWith), name_vector = BatchCorWith, name = 'BatchCor_with')
#DO_LogReg_with <- do_enrichment(1000, ml_feature = 'logistic_regression', name = 'LogReg_with')
#DO_RandForest_with <- do_enrichment(1000, ml_feature = 'random_forest_classification', name = 'RandForest_with')
#commonFeaturesWith <- common_features(batch_cor_with, name='With Outlier')

# Without Outlier
#DO_DGE_without <- do_enrichment(1000, name = 'DGE_without')
#DO_BatchCor_without <- do_enrichment(length(BatchCorWith), name_vector = BatchCorWithout, name = 'BatchCor_without')
#DO_LogReg_without <- do_enrichment(1000, ml_feature = 'logistic_regression', name = 'LogReg_without')
#DO_RandForest_without <- do_enrichment(1000, ml_feature = 'random_forest_classification', name = 'RandForest_without')
#commonFeaturesWithout <- common_features(batch_cor_without, name='Without Outlier')

####

## DO Analysis - Overlaps and Differences
#DO_common_with <- do_enrichment(length(common_with), name_vector = common_with, name = 'common_with')
#DO_common_without <- do_enrichment(length(common_without), name_vector = common_without, name = 'common_without')
#DO_complete_overlap <- do_enrichment(length(complete_overlap), name_vector = complete_overlap, name = 'complete_overlap')
#DO_complete_difference <- do_enrichment(length(complete_difference), name_vector = complete_difference, name = 'complete_difference')
#DO_DGE_difference <- do_enrichment(length(DGE_difference), name_vector = DGE_difference, name = 'DGE_difference')
#DO_BatchCor_difference <- do_enrichment(length(BatchCor_difference), name_vector = BatchCor_difference, name = 'BatchCor_difference')
#DO_LogReg_overlap <- do_enrichment(length(LogReg_overlap), name_vector = LogReg_overlap, name = 'LogReg_overlap')
#DO_LogReg_difference <- do_enrichment(length(LogReg_difference), name_vector = LogReg_difference, name = 'LogReg_difference')
#DO_RandForest_overlap <- do_enrichment(length(RandForest_overlap), name_vector = RandForest_overlap, name = 'RandForest_overlap')
#DO_RandForest_difference <- do_enrichment(length(RandForest_difference), name_vector = RandForest_difference, name = 'RandForest_difference')

###

## DO Analysis - Results
#get_do_results <- function(enrichResult, number=10) {
#  return(head(enrichResult@result$Description[order(enrichResult@result$p.adjust)], number))
#}

#DO_DGE_with_result <- get_do_results(DO_DGE_with)
#DO_DGE_without_result <- get_do_results(DO_DGE_without)
#DO_BatchCor_with_result <- get_do_results(DO_BatchCor_with)
#DO_BatchCor_without_result <- get_do_results(DO_BatchCor_without)
#DO_LogReg_with_result <- get_do_results(DO_LogReg_with)
#DO_LogReg_without_result <- get_do_results(DO_LogReg_without)
#DO_RandForest_with_result <- get_do_results(DO_RandForest_with)
#DO_RandForest_without_result <- get_do_results(DO_RandForest_without)

#DO_common_with_result <- get_do_results(DO_common_with)
#DO_common_without_result <- get_do_results(DO_common_without)
#DO_complete_overlap_result <- get_do_results(DO_complete_overlap)
#DO_complete_difference_result <- get_do_results(DO_complete_difference)
#DO_DGE_difference_result <- get_do_results(DO_DGE_difference)
#DO_BatchCor_difference_result <- get_do_results(DO_BatchCor_difference)
#DO_LogReg_overlap_result <- get_do_results(DO_LogReg_overlap)
#DO_LogReg_difference_result <- get_do_results(DO_LogReg_difference)
#DO_RandForest_overlap_result <- get_do_results(DO_RandForest_overlap)
#DO_RandForest_difference_result <- get_do_results(DO_RandForest_difference)

#DO_output = data.frame(
#  DGE_with = DO_DGE_with_result,
#  DGE_without = DO_DGE_without_result,
#  BatchCor_with = DO_BatchCor_with_result,
#  BatchCor_without = DO_BatchCor_without_result,
#  LogReg_with = DO_LogReg_with_result,
#  LogReg_without = DO_LogReg_without_result,
#  RandForest_with = DO_RandForest_with_result,
#  RandForest_without = DO_RandForest_without_result,
#  common_with = DO_common_with_result,
#  common_without = DO_common_without_result,
#  complete_overlap = DO_complete_overlap_result,
#  complete_difference = DO_complete_difference_result,
#  DGE_difference = DO_DGE_difference_result,
#  BatchCor_difference = DO_BatchCor_difference_result,
#  LogReg_overlap = DO_LogReg_overlap_result,
#  LogReg_difference = DO_LogReg_difference_result,
#  RandForest_overlap = DO_RandForest_overlap_result,
#  RandForest_difference = DO_RandForest_difference_result
#)
#write.csv(DO_output, file=paste(getwd(),'/output/', 'DO', '.csv', sep = ''), row.names = F)

#################### Not further used

# keggg/uniprot/GO
# With Outlier
#load_data()
#stats_df <- differential_expression()
#batch_cor_with <- batch_correction()
#commonPathwaysWith <- common_pathways(batch_cor_with, 100, 1000, 'uniprot', 'with')

# Without Outlier
#load_data(FALSE)
#stats_df <- differential_expression()
#batch_cor_without <- batch_correction()
#commonPathwaysWithout <- common_pathways(batch_cor_without, 100, 1000, 'uniprot', 'without')