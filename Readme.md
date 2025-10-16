# Practical Phase Laszlo Lang
Subject: Robustness of machine learning against Batch effects in RNA-seq data <br />
Supervised by Dr. Maximilien Sprang - The Mayer Lab Mainz <br />
TH-Bingen - Bachelor's Degree Program "Angewandte Bioinformatik" - 2025

---

RNA sequencing (RNA-seq) enables the systematic study of gene expression and provides insights into biological processes at the molecular level. However, RNA-seq data are often affected by batch effects.
These arise from systematic technical variation, for example, due to sequencing runs, sample preparation, or laboratory conditions. Such variation can obscure true biological signals and complicate downstream
analyses. Machine learning (ML) methods are powerful tools for detecting complex patterns in high-dimensional data. Yet, in the context of RNA-seq, it remains unclear whether ML models capture genuine
biological variation or adapt to technical artifacts. <br />
In this project, we investigate the robustness of selected ML methods against batch effects in RNA-seq data. Classical differential expression analysis is compared with ML-based approaches, and the resulting
gene lists are interpreted in biological context using over-representation analysis (ORA) and correlation analysis. This allows us to highlight differences between methods and assess which biological pathways
remain consistently detectable despite batch effects.

---

## Chapters:
- [Python](#python)
  - [Packages](#packages)
  - [Functions to prepare the DataFrame](#functions-to-prepare-the-dataframe)
    - [prepare_dataframe](#prepare_dataframe)
    - [read_excel_file](#read_excel_file)
    - [merge_dfs](#merge_dfs)
    - [label_group](#label_group)
    - [calc_log2fc](#calc_log2fc)
  - [Functions for basic analysis and visualization](#functions-for-basic-analysis-and-visualization)
    - [heatmap](#heatmap)
    - [pca](#pca)
    - [distribution_plot](#distribution_plot)
    - [mean_expression](#mean_expression)
  - [Functions for machine learning methods](#functions-for-machine-learning-methods)
    - [logistic_regression](#logistic_regression)
    - [decision_tree](#decision_tree)
    - [random_forest_classification](#random_forest_classification)
  - [Functions for plotting results](#functions-for-plotting-results)
    - [plot_logistic_regression](#plot_logistic_regression)
    - [plot_decision_tree](#plot_decision_tree)
    - [plot_feature_importance](#plot_feature_importance)
    - [plot_perm_importance](#plot_perm_importance)
- [R](#r)
  - [Packages](#packages-1)
  - [Functions to prepare the Data](#functions_to_prepare_the_data)
    - [load_data](#load_data)
    - [differential_expression](#differential_expression)
    - [batch_correction](#batch_correction)
  - [Functions for analysis and visualization](#functions_for_analysis_and_visualization)
    - [over_representation](#over_representation)
    - [fora_analysis](#fora_analysis)
    - [do_enrichment](#do_enrichment)
  - [Functions for common analysis](#functions_for_common_analysis)
    - [fora_common_pathways](#fora_common_pathways)
    - [common_pathways](#common_pathways)
    - [common_features](#common_features)
    - [common_Intersection](#common_Intersection)
  - [Functions for export](#functions_for_export)
    - [get_ontologie](#get_ontologie)
    - [get_ontology_counts](#get_ontology_counts)
    - [export_common](#export_common)
    - [export_ontology_counts](#export_ontology_counts)

---

## Python
### Packages
- `os`
- `time`
- `pandas`
- `seaborn`
- `matplotlib.pyplot`
- `numpy`
- `plotly.express`
- `functools` &#8594; `reduce`
- `sklearn.decomposition` &#8594; `PCA`
- `sklearn.preprocessing` &#8594; `StandardScaler`
- `sklearn.linear_model` &#8594; `LogisticRegression`
- `sklearn.ensemble` &#8594; `RandomForestClassifier`
- `sklearn.model_selection` &#8594; `train_test_split`
- `sklearn.metrics` &#8594; `confusion_matrix`, `accuracy_score`, `classification_report`, `roc_curve`, `roc_auc_score`, `f1_score`
- `sklearn.inspection` &#8594; `permutation_importance`
- `sklearn.tree` &#8594; `DecisionTreeClassifier`, `plot_tree`
- `typing` &#8594; `Literal`
- `scipy.stats` &#8594; `zscore`
- `collections` &#8594; `Counter`

---

### Functions to prepare the DataFrame
#### prepare_dataframe
Input: The dataframe to be transformed <br />
Output: A dataframe containing log2 tarnsformed expression data <br />
Code: [Click](https://github.com/LzLang/Praxisphase/blob/42be1f484280a211f4093c8382eb2e0b888e08c8/methods/main.py#L22)


#### read_excel_file
Reads Excel files from each subdirectory in the specified directory. <br />
Input: The path to the directory containing folders with Excel files. <br />
Output: A list of DataFrames read from the Excel files or none, if f an error occurs while reading the files. <br />
Code: [Click](https://github.com/LzLang/Praxisphase/blob/42be1f484280a211f4093c8382eb2e0b888e08c8/methods/main.py#L80)


#### merge_dfs
Merges a list of DataFrames into a single DataFrame. <br />
Input: A list of DataFrames to be merged. <br />
Output: A single DataFrame containing all merged data or none, if no DataFrames are provided or if an error occurs. <br />
Code: [Click](https://github.com/LzLang/Praxisphase/blob/42be1f484280a211f4093c8382eb2e0b888e08c8/methods/main.py#L117)


#### label_group
Labels the group based on keywords indicating tumor or normal status. <br />
Input: The group to label. <br />
Output: The label for the group ("Normal", "Tumor", or "Other"). <br />
Code: [Click](https://github.com/LzLang/Praxisphase/blob/42be1f484280a211f4093c8382eb2e0b888e08c8/methods/main.py#L173)


#### calc_log2fc
Calculates the log2 fold change (Log2FC) between normal and tumor samples in the given DataFrame. <br />
Input: dataframe (DataFrame) - The DataFrame containing gene expression data. <br />
Input: metadata (DataFrame) - The metadata DataFrame containing sample information. <br />
Input: output (bool) - Whether to save the log2 fold change results to a CSV file (default is False). <br />
Output: A Series containing the log2 fold change (Log2FC) values <br />
Code: [Click](https://github.com/LzLang/Praxisphase/blob/42be1f484280a211f4093c8382eb2e0b888e08c8/methods/main.py#LXXX)

---

### Functions for basic analysis and visualization
#### heatmap
Creates and displays a heatmap of the given DataFrame. <br />
Input: dataframe (DataFrame) - The DataFrame to create a heatmap from. <br />
Input: output_path (str) - The path where the heatmap will be saved. <br />
Input: output_name (str) - The name of the output heatmap file (default is "heatmap"). <br />
Input: method (str) - The method to use for creating the heatmap ('seaborn' or 'plotly'). <br />
Input: top_var_genes_num (int) - The number of most variable genes to display in the heatmap (default is 250). <br/>
Input: mean_value (float) - The minimum mean value for filtering genes (default is 0).<br />
Input: show (bool) - Whether to display the heatmap (default is True).<br />
Output: None - Displays a heatmap of the DataFrame.<br />
Code: [Click](https://github.com/LzLang/Praxisphase/blob/42be1f484280a211f4093c8382eb2e0b888e08c8/methods/main.py#LXXX)


#### pca
Performs PCA on the given DataFrame and displays a scatter plot of the results. <br />
Input: dataframe (DataFrame) - The DataFrame to perform PCA on. <br />
Input: metadata (DataFrame) - The metadata DataFrame containing sample information. <br />
Input: output_path (str) - The path where the PCA plot will be saved. <br />
Input: output_name (str) - The name of the output PCA plot file (default is "pca"). <br />
Input: group (str) - The column name in metadata to use for grouping samples. <br />
Input: show (bool) - Whether to display the PCA plot (default is True). <br />
Output: None - Displays a PCA plot of the DataFrame. <br />
Code: [Click](https://github.com/LzLang/Praxisphase/blob/42be1f484280a211f4093c8382eb2e0b888e08c8/methods/main.py#LXXX)


#### distribution_plot
Creates and displays a distribution plot of the specified column in the DataFrame. <br />
Input: dataframe (DataFrame) - The DataFrame to create a distribution plot from. <br />
Input: metadata (DataFrame) - The metadata DataFrame containing sample information. <br />
Input: output_path (str) - The path where the distribution plot will be saved. <br />
Input: output_name (str) - The name of the output distribution plot file (default is "distribution_plot"). <br />
Input: show (bool) - Whether to display the distribution plot (default is True). <br />
Output: None - Displays a distribution plot of the specified column. <br />
Code: [Click](https://github.com/LzLang/Praxisphase/blob/42be1f484280a211f4093c8382eb2e0b888e08c8/methods/main.py#LXXX)


#### mean_expression
Calculates the mean expression values for each group in the specified column of the metadata. <br />
Input: dataframe (DataFrame) - The DataFrame to calculate mean expression from. <br />
Input: metadata (DataFrame) - The metadata DataFrame containing sample information. <br />
Input: output_path (str) - The path where the mean expression plot will be saved. <br />
Input: output_name (str) - The name of the output mean expression plot file (default is "mean_expression"). <br />
Input: show (bool) - Whether to display the mean expression plot (default is True). <br />
Output: None - Displays a boxen plot of mean expression values by GEO_Series. <br >/
Code: [Click](https://github.com/LzLang/Praxisphase/blob/42be1f484280a211f4093c8382eb2e0b888e08c8/methods/main.py#LXXX)

---

### Functions for machine learning
#### logistic_regression
Performs logistic regression on the given DataFrame and returns the model, training and test sets. <br />
Input: dataframe (DataFrame) - The DataFrame to perform logistic regression on. <br />
Input: metadata (DataFrame) - The metadata DataFrame containing sample information. <br />
Output: A dictionary containing the results of the logistic regression analysis. <br />
Code: [Click](https://github.com/LzLang/Praxisphase/blob/42be1f484280a211f4093c8382eb2e0b888e08c8/methods/main.py#LXXX)


#### random_forest_classification
Performs random forest classification on the given DataFrame and returns the model, training and test sets. <br />
Input: dataframe (DataFrame) - The DataFrame to perform random forest classification on. <br />
Input: metadata (DataFrame) - The metadata DataFrame containing sample information. <br />
Output: A dictionary containing the results of the random forest classification analysis. <br />
Code: [Click](https://github.com/LzLang/Praxisphase/blob/42be1f484280a211f4093c8382eb2e0b888e08c8/methods/main.py#LXXX)


#### decision_tree
NOTE: This functions was never realy used. <br />
Performs decision tree analysis on the given DataFrame and displays a plot of the results. <br />
Input: dataframe (DataFrame) - The DataFrame to perform decision tree analysis on. <br />
Input: metadata (DataFrame) - The metadata DataFrame containing sample information. <br />
Input: method (str) - The method/criterion to use for the decision tree ('gini' or 'entropy'). <br />
Output: A dictionary containing the results of the decision tree analysis. <br />
Code: [Click](https://github.com/LzLang/Praxisphase/blob/42be1f484280a211f4093c8382eb2e0b888e08c8/methods/main.py#LXXX)


#### plot_logistic_regression
Plots the confusion matrix and ROC curve for the logistic regression results. <br />
Input: logreg_result (dict) - The result dictionary from logistic regression analysis. <br />
Input: output_path (str) - The path where the logistic regression plot will be saved. <br />
Input: output_name (str) - The name of the output logistic regression plot file (default is "logistic_regression"). <br />
Input: show (bool) - Whether to display the logistic regression plot (default is False). <br />
Output: None - Displays a confusion matrix and ROC curve for the logistic regression results. <br />
Code: [Click](https://github.com/LzLang/Praxisphase/blob/42be1f484280a211f4093c8382eb2e0b888e08c8/methods/main.py#LXXX)

#### plot_decision_tree
Plots the confusion matrix and decision tree for the decision tree results. <br />
Input: tree_result (dict) - The result dictionary from decision tree analysis. <br />
Input: output_path (str - The path where the decision tree plot will be saved. <br />
Input: output_name (str) - The name of the output decision tree plot file (default is "decision_tree"). <br />
Input: show (bool) - Whether to display the decision tree plot (default is False). <br />
Output: None - Displays a confusion matrix and decision tree plot for the decision tree results. <br />
Code: [Click](https://github.com/LzLang/Praxisphase/blob/42be1f484280a211f4093c8382eb2e0b888e08c8/methods/main.py#LXXX)

#### plot_feature_importance
Calculates feature importance using the specified method and displays a plot of the results. <br />
Input: algorithm_result (dict) - The result dictionary from the algorithm analysis. <br />
Input: output_path (str) - The path where the feature importance plot will be saved. <br />
Input: output_name (str) - The name of the output feature importance plot file (default is "feature_importance"). <br />
Input: method (str) - The method to use for calculating feature importance ('logistic_regression' or 'forest_classification'). <br />
Input: show (bool) - Whether to display the feature importance plot (default is False). <br />
Output: None - Displays a plot of feature importance based on the specified method. <br />
Code: [Click](https://github.com/LzLang/Praxisphase/blob/42be1f484280a211f4093c8382eb2e0b888e08c8/methods/main.py#LXXX)


#### plot_perm_importance
Calculates permutation importance for the specified algorithm results and displays a plot of the results. <br />
Input: algorithm_result (dict) - The result dictionary from the algorithm analysis. <br />
Input: output_path (str) - The path where the permutation importance plot will be saved. <br />
Input: output_name (str) - The name of the output permutation importance plot file (default is "permutation_importance"). <br />
Input: show (bool) - Whether to display the permutation importance plot (default is False). <br />
Output: None - Displays a plot of permutation importance for the specified algorithm results. <br />
Code: [Click](https://github.com/LzLang/Praxisphase/blob/42be1f484280a211f4093c8382eb2e0b888e08c8/methods/main.py#LXXX)

---

## R
### Packages
