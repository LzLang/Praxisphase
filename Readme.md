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
  - [Functions to prepare the DataFrame](#functions_to_prepare_the_dataframe)
    - [prepare_dataframe](#prepare_dataframe)
    - [read_excel_file](#read_excel_file)
    - [merge_dfs](#merge_dfs)
    - [label_group](#label_group)
    - [calc_log2fc](#calc_log2fc)
  - [Functions for basic analysis and visualization](#functions_for_basic_analysis_and_visualization)
    - [heatmap](#heatmap)
    - [pca](#pca)
    - [distribution_plot](#distribution_plot)
    - [mean_expression](#mean_expression)
  - [Functions for machine learning methods](#functions_for_machine_learning_methods)
    - [logistic_regression](#logistic_regression)
    - [decision_tree](#decision_tree)
    - [random_forest_classification](#random_forest_classification)
  - [Functions for plotting results](#functions_for_plotting_results)
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
#### prepare

## R
### Packages
