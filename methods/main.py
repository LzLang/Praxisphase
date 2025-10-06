import os
import time
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
from functools import reduce
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, accuracy_score, classification_report, roc_curve, roc_auc_score, f1_score
from sklearn.inspection import permutation_importance
from sklearn.tree import DecisionTreeClassifier, plot_tree
from typing import Literal
from scipy.stats import zscore

# Functions to prepare the DataFrame for analysis
def prepare_dataframe(dataframe):
    """
    input: dataframe (DataFrame): The DataFrame to be transformed.
    output: expression_log (DataFrame): A DataFrame containing log2 transformed expression data.
    """
    # Set Gen Name as index
    # dataframe = dataframe.set_index("gene_id_ver")
    #print(dataframe.index.isnull().sum(), "trans_id_ver-Einträge sind NaN")
    print("Preparing DataFrame...")

    # Select columns that start with "SRR"
    srr_columns = [col for col in dataframe.columns if col.startswith("SRR")]

    dataframe[srr_columns] = dataframe[srr_columns].replace("", 0)
    dataframe[srr_columns] = dataframe[srr_columns].fillna(0)

    # Create a dictionary for aggregation
    # This will sum the values for each SRR column
    agg_dict = {
        "gene_name": lambda x: x.iloc[0],  # Keep the first value of gene_id_ver
        "trans_id_ver": lambda x: "\n".join(set(x)),   # Join unique
        "trans_name": lambda x: "\n".join(set(x))   # Join unique
    }
    agg_dict.update({col: 'sum' for col in srr_columns})

    dataframe = dataframe.groupby("gene_id_ver").agg(agg_dict).reset_index()

    #print(len(srr_columns), "SRR columns found.")
    expression_data = dataframe.loc[:, srr_columns]
    # Replace commas with dots for decimal conversion
    expression_data = expression_data.replace(",", ".", regex=True)
    # Replace empty strings with zeros
    # This will ensure that any empty cells are treated as zero, which is important for numerical operations
    expression_data = expression_data.replace("", 0)
    # Replace NaN values with zeros
    # This will ensure that any missing values are treated as zero, which is important for numerical operations
    expression_data = expression_data.fillna(0)
    # Convert the DataFrame to numeric, forcing errors to NaN
    # This will convert all values to int, and any non-numeric values will become NaN
    expression_data = expression_data.astype(float)

    #Filter out rows and columns with all zero values
    # This will remove rows and columns where all values are zero
    # This is important for visualizing the data, as we want to focus on genes with some expression
    # This will help in reducing noise in the data and focusing on genes that are expressed
    expression_filtered = expression_data[(expression_data > 0).any(axis=1)]
    expression_filtered = expression_filtered.loc[:, (expression_filtered > 0).any(axis=0)]

    # Check if the DataFrame is empty after filtering
    if expression_filtered.empty:
        print("No data available for heatmap after filtering.")
        return
    
    
    print("DataFrame prepared successfully.")

    return dataframe

def read_excel_file(directory):
    """
    input: directory (str): The path to the directory containing folders with Excel files.
    output: list: A list of DataFrames read from the Excel files.
    output: None: If an error occurs while reading the files.
    Reads Excel files from each subdirectory in the specified directory.
    """
    try:
        dfs = []  # Initialize an empty list to store DataFrames
        filename = "all.samples.tpm.xlsx"  # Name der Datei, die in jedem Ordner eingelesen werden soll
        # Go through each folder in the specified directory
        for folder in os.listdir(directory):
            folder_path = os.path.join(directory, folder)
            # Check if the path is a directory
            if os.path.isdir(folder_path):
                # Construct the full file path
                file_path = os.path.join(folder_path, filename)
                if os.path.isfile(file_path):
                    current_file = os.path.join(folder, filename)
                    print(f"Reading file: {current_file}")
                    # Use pd.read_excel to read the file and append the DataFrame to the list
                    dataframe = pd.read_excel(file_path)
                    
                    prepared_dataframe = prepare_dataframe(dataframe)
                    prepared_dataframe[f"gene_name_{folder}"] = prepared_dataframe["gene_name"]
                    prepared_dataframe[f"trans_id_ver_{folder}"] = prepared_dataframe["trans_id_ver"]
                    prepared_dataframe[f"trans_name_{folder}"] = prepared_dataframe["trans_name"]
                    prepared_dataframe = prepared_dataframe.drop(columns=["gene_name", "trans_id_ver", "trans_name"])
                    dfs.append(prepared_dataframe)
                    # print out the number of rows in the DataFrame
                    print(f"Number of rows in {current_file}: {len(prepared_dataframe)}")
                    print(f"File {current_file} read succesfully.\n~~~~")
        return dfs
    except Exception as e:
        print(f"Error reading {directory}: {e}")
        return None

def merge_dfs(dfs):
    """
    input: dfs (list): A list of DataFrames to be merged.
    output: merged_dataframes: A single DataFrame containing all merged data.
    output: None: If no DataFrames are provided or if an error occurs.
    Merges a list of DataFrames into a single DataFrame.
    """
    if not dfs:
        # If the list is empty, return None
        print("No DataFrames to merge.")
        return None
    try:
        # Merge all DataFrames in the list using reduce
        print("Merging DataFrames...")
        merged_dataframes = reduce(lambda left, right: pd.merge(left, right, on="gene_id_ver", how='outer'), dfs)

        srr_columns = [col for col in merged_dataframes.columns if col.startswith("SRR")]
        merged_dataframes[srr_columns] = merged_dataframes[srr_columns].fillna(0)
        merged_dataframes[srr_columns] = merged_dataframes[srr_columns].replace("", 0)

        
        gene_name_columns = [col for col in merged_dataframes.columns if col.startswith("gene_name")]
        gene_name_series = merged_dataframes[gene_name_columns].bfill(axis=1).iloc[:, 0]
        gene_name_series = gene_name_series.fillna(merged_dataframes["gene_id_ver"])

        trans_name_cols = [col for col in merged_dataframes.columns if col.startswith("trans_name")]
        trans_name_series = merged_dataframes[trans_name_cols].apply(
            lambda row: '\n'.join(dict.fromkeys(
                val.strip()
                for cell in row.dropna()
                for val in cell.split('\n')
                if val.strip() != ''
            )), axis=1)
        
        trans_id_cols = [col for col in merged_dataframes.columns if col.startswith("trans_id_ver")]
        trans_id_series = merged_dataframes[trans_id_cols].apply(
            lambda row: '\n'.join(dict.fromkeys(
                val.strip()
                for cell in row.dropna()
                for val in cell.split('\n')
                if val.strip() != ''
            )), axis=1)

        
        dataframe = merged_dataframes["gene_id_ver"]
        dataframe = pd.concat([dataframe, gene_name_series.rename("gene_name")], axis=1)
        dataframe = pd.concat([dataframe, trans_name_series.rename("trans_name")], axis=1)
        dataframe = pd.concat([dataframe, trans_id_series.rename("trans_id_ver")], axis=1) 
        dataframe = pd.concat([dataframe, merged_dataframes[srr_columns]], axis=1)
        print("DataFrames merged successfully.")
        return dataframe
    except Exception as e:
        # If an error occurs during merging, print the error and return None
        print(f"Error merging DataFrames: {e}")
        return None

def label_group(group):
    """
    input: group (str): The group to label.  
    output: str: The label for the group ("Normal", "Tumor", or "Other").
    Labels the group based on keywords indicating tumor or normal status.
    """
    tumor_keywords = ["tumor", "liver_cancer", "advanced_hcc"]
    normal_keywords = ["normal", "control", "non_tumor", "non-tumor", "adjacent-non-tumour"]
    if "fibrosis" in group.lower():
        print("Fibrosis detected, labeling as Tumor")

    if any(keyword == group.lower() for keyword in normal_keywords):
        return "Normal"
    elif any(keyword == group.lower() for keyword in tumor_keywords):
        return "Tumor"
    else:
        print(f"Group '{group}' does not match any known tumor or normal keywords, labeling as Other")
        return "Other"

# Functions for basic analysis and visualization
def heatmap(dataframe, output_path, output_name="heatmap",top_var_genes_num=250, expression = 0, 
            method: Literal['seaborn', 'plotly']='plotly', show: Literal[True, False]=False):
    """
    input: dataframe (DataFrame): The DataFrame to create a heatmap from.
    input: output_path (str): The path where the heatmap will be saved.
    input: output_name (str): The name of the output heatmap file (default is "heatmap").
    input: method (str): The method to use for creating the heatmap ('seaborn' or 'plotly').
    input: top_var_genes_num (int): The number of most variable genes to display in the heatmap (default is 250).
    input: mean_value (float): The minimum mean value for filtering genes (default is 0).
    input: show (bool): Whether to display the heatmap (default is True).
    output: None: Displays a heatmap of the DataFrame.
    Creates and displays a heatmap of the given DataFrame.
    """
    try:    
        # Filter out rows with a mean value less than or equal to mean_value
        # This will help in focusing on genes with significant expression levels
        expression_filtered = dataframe.where(dataframe >= expression)
        
        # Select the top x most variable genes based on variance
        # This will help in visualizing the most variable genes in the heatmap
        top_var_genes = expression_filtered.var(axis=1).sort_values(ascending=False).head(top_var_genes_num).index
        heatmap_data = expression_filtered.loc[top_var_genes]
        num_rows, num_cols = heatmap_data.shape
        print(f"Creating heatmap with {num_rows} rows and {num_cols} columns for the top {top_var_genes_num} most variable genes with an expression level >= {expression}.")

        # Creating a heatmap according to the specified method
        if method == 'plotly':
            # Use Plotly Express to create an interactive heatmap
            print("Creating heatmap using Plotly...")
        
            # Create a heatmap using Plotly Express
            fig = px.imshow(heatmap_data,
                        labels={
                            'x': 'Sample (SRR)',
                            'y': 'Genes (gene_name)',
                            'color': 'log2(TPM + 1)'},
                        color_continuous_scale='viridis',
                        title=f'Heatmap of Gene Expression Levels (Top {top_var_genes_num} Most Variable Genes with an expression level >= {expression})',
                        aspect='auto')

            # Update layout for better readability
            fig.update_layout(
                width=max(1920, num_cols * 5),          # Adjust width based on number of columns
                height=max(1080, num_rows * 10),        # Adjust height based on number of rows
                margin=dict(l=20, r=20, t=50, b=20)
            )

            # Save the heatmap without traces to an HTML file
            fig.write_html(f"{output_path}{output_name}_top{top_var_genes_num}_expression{expression}.html")

            # Update traces to adjust gaps between cells
            # This will help in better visualization of the heatmap
            fig.update_traces( 
                xgap=1,
                ygap=1,  
                selector=dict(type='heatmap')  
            )

            # Save the heatmap with traces to an HTML file
            fig.write_html(os.path.join(output_path, f"{output_name}_top{top_var_genes_num}_expression{expression}.html"))
            print("Heatmap created successfully.")
            # Show the heatmap
            if show:
                fig.show()
        elif method == 'seaborn':
            # Use seaborn to create a static heatmap
            print("Creating heatmap using seaborn...")
            # Set the figure size for better readability
            plt.figure(figsize=(32, 18))
            # Create a heatmap using seaborn
            sns.heatmap(heatmap_data,
                       cmap="viridis", 
                       annot=False, 
                       cbar_kws={"label": "log2(TPM + 1)"})#,
                      # linewidths=0.01,
                       #yticklabels=heatmap_data["gene_name"])
            
            plt.yticks(rotation=45)  # Rotate y-axis labels for better readability
            plt.title(f"Heatmap of Gene Expression Levels (Top {top_var_genes_num} Most Variable Genes with an expression level >= {expression})")
            plt.xlabel("Sample (SRR)")
            plt.ylabel("Genes (gene_name)")
            plt.tight_layout()
            # Save the heatmap to a file
            plt.savefig(os.path.join(output_path, f"{output_name}_top{top_var_genes_num}_expression{expression}.png"))
            print("Heatmap created successfully.")
            # Show the heatmap
            if show:
                plt.show()
        else:
            print(f"Unsupported method: {method}. Please use 'plotly' or 'seaborn'.")
            raise ValueError(f"Unsupported method: {method}. Please use 'plotly' or 'seaborn'.")
    except Exception as e:
        print(f"Error creating heatmap: {e}")

def pca(dataframe, metadata, output_path, output_name="pca",
        show: Literal[True, False]=False, show_geo: Literal[True, False]=False):
    """
    input: dataframe (DataFrame): The DataFrame to perform PCA on.
    input: metadata (DataFrame): The metadata DataFrame containing sample information.
    output_path (str): The path where the PCA plot will be saved.
    output_name (str): The name of the output PCA plot file (default is "pca").
    input: group (str): The column name in metadata to use for grouping samples.
    input: show (bool): Whether to display the PCA plot (default is True).
    output: None: Displays a PCA plot of the DataFrame.

    Performs PCA on the given DataFrame and displays a scatter plot of the results.
    """
    try:
        print("Performing PCA...")


        # Transpose the DataFrame to have genes as rows and samples as columns
        dataframe = dataframe.T

        # Standardize the data
        scaler = StandardScaler()
        # This will scale the data to have mean=0 and variance=1
        # This is important for PCA to ensure that all features contribute equally
        scaled_data = scaler.fit_transform(dataframe)
        # Perform PCA
        # This will reduce the dimensionality of the data to 2 components for visualization
        # PCA is a technique used to reduce the dimensionality of the data while preserving as much variance as possible
        pca = PCA(n_components=2)
        # Fit PCA on the scaled data and transform it
        # This will compute the principal components and transform the data accordingly
        pca_result = pca.fit_transform(scaled_data)
        # Create a DataFrame with PCA results
        # This will create a DataFrame with the first two principal components
        # The first principal component will be in the first column and the second in the second column
        pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
        # Add a column for sample names
        # This will help in identifying which sample corresponds to which point in the PCA plot
        pca_df["Sample"] = dataframe.index
        # Merge PCA results with metadata
        # This will help in associating PCA results with sample metadata
        # GEO_Series instead of group 
        if show_geo:
            pca_annotated = pca_df.merge(metadata[["Run", "GEO_Series"]], left_on="Sample", right_on="Run", how="left")
            plt.figure(figsize=(16, 9))
            sns.scatterplot(data=pca_annotated, x='PC1', y='PC2', hue="GEO_Series", alpha=0.7)
            plt.title("PCA of Gene Expression")
            plt.legend(title="GSE-Series")
        else:
            pca_annotated = pca_df.merge(metadata[["Run", "Tumor_Status"]], left_on="Sample", right_on="Run", how="left")
            plt.figure(figsize=(16, 9))
            sns.scatterplot(data=pca_annotated, x='PC1', y='PC2', hue="Tumor_Status", palette={"Tumor": "red", "Normal": "blue", "Other": "gray"}, alpha=0.7)
            if "Other" in metadata["Tumor_Status"].unique():
                plt.title("PCA of Gene Expression (Standardized) Normal vs Tumor vs Other")
            else:
                plt.title("PCA of Gene Expression (Standardized) Normal vs Tumor")
            plt.legend(title="Tumor Status")
        print(pca.explained_variance_ratio_)


        print("PCA performed successfully.")
        
        # Use seaborn to create a scatter plot with hue based on the specified group
        # This will color the points based on the group specified in the metadata
        plt.xlabel("PC1")
        plt.ylabel("PC2")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(output_path, f"{output_name}.png"))
        if show:
            plt.show()  
    except Exception as e:
        print(f"Error performing PCA: {e}")
        return None
    
def distribution_plot(dataframe, metadata, output_path, output_name="distribution_plot",
                      show: Literal[True, False]=False):
    """
    input: dataframe (DataFrame): The DataFrame to create a distribution plot from.
    input: metadata (DataFrame): The metadata DataFrame containing sample information.
    output_path (str): The path where the distribution plot will be saved.
    output_name (str): The name of the output distribution plot file (default is "distribution_plot").
    input: show (bool): Whether to display the distribution plot (default is True).
    output: None: Displays a distribution plot of the specified column.
    
    Creates and displays a distribution plot of the specified column in the DataFrame.
    """
    try:
        print("Performing distribution plotting...")

        expression_long = dataframe.reset_index().melt(id_vars="gene_name", 
                                                         var_name="Sample", 
                                                         value_name="Expression")
        
        expression_long = expression_long.merge(metadata[["Run", "GEO_Series"]], left_on="Sample", right_on="Run", how="left")


        plt.figure(figsize=(16, 9))
        sns.histplot(data=expression_long, x="Expression", hue="GEO_Series", kde=True, bins=50, element="step", common_norm=False)
        plt.title("Distribution of Expression Levels by 'GEO_Series'")
        plt.xlabel("Z-Score(log2(Expression + 1))")
        plt.ylabel("Density")
        plt.grid(True)
        plt.savefig(os.path.join(output_path, f"{output_name}.png"))

        print("Distribution plot created successfully.")

        if show:
            plt.show()
    except Exception as e:
        print(f"Error creating distribution plot: {e}")

def mean_expression(dataframe, metadata, output_path, output_name="mean_expression",
                      show: Literal[True, False]=False):
    """
    input: dataframe (DataFrame): The DataFrame to calculate mean expression from.
    input: metadata (DataFrame): The metadata DataFrame containing sample information.
    output_path (str): The path where the mean expression plot will be saved.
    output_name (str): The name of the output mean expression plot file (default is "mean_expression").
    input: show (bool): Whether to display the mean expression plot (default is True).
    output: None: Displays a boxen plot of mean expression values by GEO_Series.
    
    Calculates the mean expression values for each group in the specified column of the metadata.
    """
    try:
        print("Calculating mean expression...")

        # Calculate mean expression for each sample
        # This will compute the mean expression for each sample across all genes
        mean_expression = dataframe.mean(axis=0)
        mean_expression_df = mean_expression.reset_index()
        mean_expression_df.columns = ['Sample', 'Mean_Expression']


        # Merge mean expression with metadata
        mean_expression_df = mean_expression_df.merge(metadata[["Run", "GEO_Series"]], left_on="Sample", right_on="Run", how="left")

        print("Mean expression calculated successfully.")
        
        # Plotting the mean expression levels
        # This will create a boxen plot to visualize the distribution of mean expression levels for each group
        plt.figure(figsize=(16, 9))
        sns.boxenplot(data=mean_expression_df, x='Mean_Expression', y="GEO_Series", orient='h')
        sns.stripplot(data=mean_expression_df, x='Mean_Expression', y="GEO_Series", color='black', alpha=0.5, jitter=True)
        plt.title("Mean Expression Levels by GEO_Series")
        #plt.xlabel("Mean Expression (unnormalized)")
        plt.xlabel("Mean Expression (log2(Expression + 1))")
        #plt.xlabel("Mean Expression Z-Score(log2(Expression + 1))")
        plt.ylabel("GEO_Series")
        plt.tight_layout()
        plt.savefig(os.path.join(output_path, f"{output_name}.png"))

        print("Mean expression plot created successfully.")

        if show:
            plt.show()
        
    except Exception as e:
        print(f"Error calculating mean expression: {e}")
        return None

# Functions for machine learning methods
def logistic_regression(dataframe, metadata):
    """
    input: dataframe (DataFrame): The DataFrame to perform logistic regression on.
    input: metadata (DataFrame): The metadata DataFrame containing sample information.
    output: dict: A dictionary containing the results of the logistic regression analysis.

    Performs logistic regression on the given DataFrame and returns the model, training and test sets.
    """
    try:
        print("Performing Logistic Regression...")

        scaler = StandardScaler()
        X = pd.DataFrame(
            scaler.fit_transform(dataframe.T),
            index=dataframe.T.index,
            columns=dataframe.T.columns
        )
        y = metadata.set_index("Run").loc[X.index, "Tumor_Status"].map({
            "Normal": 0,
            "Tumor": 1
        })

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

        model = LogisticRegression(max_iter=10000, random_state=42)

        # fit the model with data
        model.fit(X_train, y_train)
        # predict the labels for the test set
        y_pred = model.predict(X_test)

        acc_score = accuracy_score(y_test, y_pred)
        cnf_matrix = confusion_matrix(y_test, y_pred)
        f1_per_class = f1_score(y_test, y_pred, average=None)
        f1_micro = f1_score(y_test, y_pred, average='micro')
        f1_macro = f1_score(y_test, y_pred, average='macro')
        f1_weighted = f1_score(y_test, y_pred, average='weighted')
        # F1 score
        #   F1 score is the harmonic mean of precision and recall
        #   Micro-average F1 score: Calculates metrics globally by counting the total true positives, false negatives, and false positives
        #   Macro-average F1 score: Calculates metrics for each label, and finds their unweighted mean
        #   Weighted-average F1 score: Calculates metrics for each label, and finds their average weighted by support (the number of true instances for each label)
        cft_report = classification_report(y_test, y_pred, target_names=["Normal", "Tumor"])
        # Classification Report
        # Precision: How accurate/precise is the model?
        #   When a model makes a prediction, how often is it correct
        # Recall: If there are Tumors in the test set and my logistic regression model can identifiy it X% of time
        
        # Calculate coefficients and odds ratios
        coefficients = model.coef_[0]
        # Features with higher coefficients are deemed more important/significant
        # Coefficients represents the change in the log-odds of the outcome for a one-unit change in the feature (all other features held constant)
        # Pos. coefficients indicate that as the feature increases, the odds of the outcome increase
        # Neg. coefficients indicate that as the feature increases, the odds of the outcome decrease
        odds_ratios = np.exp(coefficients)
        # Odds ratios are the exponentiation of the coefficients, which represent the multiplicative change in the odds of the outcome for a one-unit change in the feature
        # Odds ratio > 1: Feature increases the odds of the outcome
        # Odds ratio < 1: Feature decreases the odds of the outcome
        # Odds ratio = 1: Feature has no effect on the odds of the outcome

        print("Logistic Regression performed successfully.")

        print(f"Accuracy: {acc_score:.4f}")
        print("Confusion Matrix:\n", cnf_matrix)
        print(f"F1 score per class: {f1_per_class}")
        print(f"Micro-average F1 score: {f1_micro:.4f}")
        print(f"Macro-average F1 score: {f1_macro:.4f}")
        print(f"Weighted-average F1 score: {f1_weighted:.4f}")
        print("Classification Report: ", cft_report)

        return {
            'model': model,
            'X': X,
            'y': y,
            'features': X.columns,
            'X_train': X_train,
            'X_test': X_test, 
            'y_train': y_train, 
            'y_test': y_test,
            'y_pred': y_pred,
            'cnf_matrix': cnf_matrix, 
            'acc_score': acc_score,
            'f1_per_class': f1_per_class,
            'f1_micro': f1_micro,
            'f1_macro': f1_macro,
            'f1_weighted': f1_weighted, 
            'cft_report': cft_report,
            'coefficients': coefficients,
            'odds_ratios': odds_ratios
        }
        
    except Exception as e:
        print(f"Error performing Logistic Regression: {e}")
        return None

def decision_tree(dataframe, metadata, method: Literal['gini', 'entropy']='gini'):
    """
    input: dataframe (DataFrame): The DataFrame to perform decision tree analysis on.
    input: metadata (DataFrame): The metadata DataFrame containing sample information.
    input: method (str): The method/criterion to use for the decision tree ('gini' or 'entropy').
    output: dict: A dictionary containing the results of the decision tree analysis.

    Performs decision tree analysis on the given DataFrame and displays a plot of the results.
    """
 
    try:
        print(f"Preparing data for Decision Tree using {method}...")

        scaler = StandardScaler()
        X = pd.DataFrame(
            scaler.fit_transform(dataframe.T),
            index=dataframe.T.index,
            columns=dataframe.T.columns
        )
        y = metadata.set_index("Run").loc[X.index, "Tumor_Status"].map({
            "Normal": 0,
            "Tumor": 1
        })

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

        # Fit a Decision Tree Classifier with gini criterion
        model = DecisionTreeClassifier(criterion=method, random_state=42)
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        
        acc_score = accuracy_score(y_test, y_pred)
        cnf_matrix = confusion_matrix(y_test, y_pred)
        cft_report = classification_report(y_test, y_pred, target_names=["Normal", "Tumor"])

        print("Decision Tree analysis performed successfully.")

        print(f"Accuracy (gini): {acc_score:.4f}")
        print("Confusion Matrix (gini):\n", cnf_matrix)
        print("Classification Report (gini): ", cft_report)
        print(f"Accuracy (entropy): {acc_score:.4f}")
        print("Confusion Matrix (entropy):\n", cnf_matrix)
        print("Classification Report (entropy): ", cft_report)

        return {
            'model': model,
            'X': X,
            'y': y,
            'criterion': method,
            'features': X.columns,
            'X_train': X_train, 
            'X_test': X_test, 
            'y_train': y_train, 
            'y_test': y_test,
            'y_pred_gini': y_pred,
            'acc_score_gini': acc_score,
            'cnf_matrix_gini': cnf_matrix, 
            'cft_report_gini': cft_report,
        }

    except Exception as e:
        print(f"Error preparing data for Decision Tree: {e}")
        return None

def random_forest_classification(dataframe, metadata):
    """
    input: dataframe (DataFrame): The DataFrame to perform random forest classification on.
    input: metadata (DataFrame): The metadata DataFrame containing sample information.
    output: dict: A dictionary containing the results of the random forest classification analysis.

    Performs random forest classification on the given DataFrame and returns the model, training and test sets.
    """
    try:
        print("Performing Random Forest Classification...")

        scaler = StandardScaler()
        X = pd.DataFrame(
            scaler.fit_transform(dataframe.T),
            index=dataframe.T.index,
            columns=dataframe.T.columns
        )
        y = metadata.set_index("Run").loc[X.index, "Tumor_Status"].map({
            "Normal": 0,
            "Tumor": 1
        })

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

        model = RandomForestClassifier(n_estimators=100, random_state=42)
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)

        acc_score = accuracy_score(y_test, y_pred)
        cnf_matrix = confusion_matrix(y_test, y_pred)
        cft_report = classification_report(y_test, y_pred, target_names=["Normal", "Tumor"])

        print("Random Forest Classification performed successfully.")

        print(f"Accuracy: {acc_score:.4f}")
        print("Confusion Matrix:\n", cnf_matrix)
        print("Classification Report: ", cft_report)

        return {
            'model': model,
            'X': X,
            'y': y,
            'features': X.columns,
            'Feature Importance': model.feature_importances_,
            'X_train': X_train, 
            'X_test': X_test, 
            'y_train': y_train, 
            'y_test': y_test,
            'y_pred': y_pred,
            'cnf_matrix': cnf_matrix, 
            'acc_score': acc_score,
            'cft_report': cft_report
        }
    except Exception as e:
        print(f"Error performing Random Forest Classification: {e}")
        return None

# Functions for plotting results
def plot_logistic_regression(logreg_result, output_path, output_name="logistic_regression",
                        show: Literal[True, False]=False):
    """
    input: logreg_result (dict): The result dictionary from logistic regression analysis.
    input: output_path (str): The path where the logistic regression plot will be saved.
    input: output_name (str): The name of the output logistic regression plot file (default is "logistic_regression").
    input: show (bool): Whether to display the logistic regression plot (default is False).
    output: None: Displays a confusion matrix and ROC curve for the logistic regression results.

    Plots the confusion matrix and ROC curve for the logistic regression results.
    """

    try:
        print("Plotting Logistic Regression results...")
        class_names = ["Normal", "Tumor"]
        # Create the heatmap
        plt.figure(figsize=(16, 9))
        sns.heatmap(pd.DataFrame(logreg_result['cnf_matrix'], index=class_names, columns=class_names), annot=True, cmap="YlGnBu", fmt='g')
        plt.title('Confusion matrix')
        plt.ylabel('Actual label')
        plt.xlabel('Predicted label')
        plt.gca().xaxis.set_label_position('top')
        plt.gca().xaxis.tick_top()  # auch die Ticks nach oben setzen
        plt.tight_layout()
        plt.savefig(os.path.join(output_path, f"{output_name}.png"))
        print("Confusion matrix plotted successfully.")
        if show:
            plt.show()
        plt.close()

        # Plotting the ROC curve
        # ROC - Receiver Operating Characteristic
        # Curve is a plot of the true positive rate against the false positive rate.
        # It shows the tradeoff between sensitivity and specificity
        # AUC score 1 represents a perfect classifier, and 0.5 represents a worthless classifier
        print("Plotting ROC curve...")
        plt.figure(figsize=(16, 9))
        y_pred_proba = logreg_result['model'].predict_proba(logreg_result['X_test'])[::,1]
        fpr, tpr, _ = roc_curve(logreg_result['y_test'],  y_pred_proba)
        auc = roc_auc_score(logreg_result['y_test'], y_pred_proba)
        plt.plot(fpr,tpr,label="data 1, auc="+str(auc))
        plt.title('Receiver Operating Characteristic (ROC) Curve')
        plt.legend(loc=4)
        plt.savefig(os.path.join(output_path, f"{output_name}_roc_curve.png"))
        print("ROC curve plotted successfully.")
        if show:
            plt.show()
    except Exception as e:
        print(f"Error plotting Logistic Regression results: {e}")
        return None

def plot_decision_tree(tree_result, output_path, output_name="decision_tree",
                       show: Literal[True, False]=False):
    """
    input: tree_result (dict): The result dictionary from decision tree analysis.
    input: output_path (str): The path where the decision tree plot will be saved.
    input: output_name (str): The name of the output decision tree plot file (default is "decision_tree").
    input: show (bool): Whether to display the decision tree plot (default is False).
    output: None: Displays a confusion matrix and decision tree plot for the decision tree results.
    
    Plots the confusion matrix and decision tree for the decision tree results.
    """
    
    try:
        print("Plotting Decision Tree results...")
        plt.figure(figsize=(16, 9))
        plot_tree(tree_result['model'],
                  feature_names=tree_result['features'],
                  class_names=["Normal", "Tumor"],
                  filled=True, 
                  rounded=True)
        plt.title(f'Decision Tree - {tree_result["criterion"]} Criterion')
        plt.tight_layout()
        plt.savefig(os.path.join(output_path, f"{output_name}_{tree_result['criterion']}.png"))
        print("Decision Tree plotted successfully.")
        if show:
            plt.show()
    except Exception as e:
        print(f"Error plotting Decision Tree results: {e}")
        return None

def plot_feature_importance(algorithm_result, output_path, output_name="feature_importance",
                       method: Literal["logistic_regression", "random_forest_classification"]="logistic_regression",
                       show: Literal[True, False] = False):
    """
    input: algorithm_result (dict): The result dictionary from the algorithm analysis.
    output_path (str): The path where the feature importance plot will be saved.
    output_name (str): The name of the output feature importance plot file (default is "feature_importance").
    input: method (str): The method to use for calculating feature importance ('logistic_regression' or 'forest_classification').
    input: show (bool): Whether to display the feature importance plot (default is False).
    output: None: Displays a plot of feature importance based on the specified method.

    Calculates feature importance using the specified method and displays a plot of the results.
    """
    try:       
        match method:
            case "logistic_regression":
                print("Calculating feature importance using Logistic Regression...")

                # Feature importance using coefficients and odds ratios
                feature_importance_df = pd.DataFrame({
                    'Gene': algorithm_result['features'],
                    'Coefficient': algorithm_result['coefficients'],
                    'Odds Ratio': algorithm_result['odds_ratios']
                })
                
                # Sort the DataFrame by coefficients in descending order     
                feature_importance_df = feature_importance_df.sort_values(by='Coefficient', ascending=False)
                
                plt.figure(figsize=(15, 9))
                sns.scatterplot(data=feature_importance_df, x='Gene', y='Coefficient', alpha=0.7)
                plt.title("Feature Importance using Logistic Regression")
                plt.tight_layout()
                plt.savefig(os.path.join(output_path, f"{output_name}_logistic_regression.png"))
                if show: 
                    plt.show()
     
            case "random_forest_classification":
                print("Calculating feature importance using Random Forest Classification...")

                # Feature importance using feature importances from the Random Forest model
                feature_importance_df = pd.DataFrame({
                    'Gene': algorithm_result['features'],
                    'Feature Importance': algorithm_result['Feature Importance']
                })

                # Sort the DataFrame by feature_importance in descending order     
                feature_importance_df = feature_importance_df.sort_values(by='Feature Importance', ascending=False)

                plt.figure(figsize=(15, 9))
                sns.scatterplot(data=feature_importance_df, x='Gene', y='Feature Importance', alpha=0.7)
                plt.title("Feature Importance using Random Forest Classification")
                plt.tight_layout()
                plt.savefig(os.path.join(output_path, f"{output_name}_random_forest_classification.png"))
                print("Feature importance calculated successfully.")
                if show: 
                    plt.show()
            case _:
                print(f"Unsupported method: {method}. Please use 'logistic_regression', or 'random_forest_classification'.")
    except Exception as e:
        print(f"Error calculating feature importance: {e}")
        return None      

def plot_perm_importance(algorithm_result, output_path, output_name="permutation_importance",
                           show: Literal[True, False]=False):
    """
    input: algorithm_result (dict): The result dictionary from the algorithm analysis.
    output_path (str): The path where the permutation importance plot will be saved.
    output_name (str): The name of the output permutation importance plot file (default is "permutation_importance").
    input: show (bool): Whether to display the permutation importance plot (default is False).
    output: None: Displays a plot of permutation importance for the specified algorithm results.

    Calculates permutation importance for the specified algorithm results and displays a plot of the results.
    """
    
    print("Calculating permutation importance...")
    start_time = time.time()
    perm_importance = permutation_importance(algorithm_result['model'], algorithm_result['X_test'], algorithm_result['y_test'], n_repeats=30, random_state=42, n_jobs=-1)
    elapsed_time = time.time() - start_time
    perm_importance_df = pd.DataFrame({
        'Feature': algorithm_result['features'],
        'Importance Mean': perm_importance.importances_mean,
        'Importance Std': perm_importance.importances_std
    })
    print(f"Elapsed time to compute the permutation importances: {elapsed_time:.3f} seconds")
    
    perm_importance_df = perm_importance_df.sort_values(by='Importance Mean', ascending=False)
                
    plt.figure(figsize=(15, 9))
    sns.scatterplot(data=perm_importance_df, x='Gene', y='Importance Mean', alpha=0.7)
    plt.title("Permutation Importance")
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, f"{output_name}.png"))
    print("Permutation importance calculated successfully.")
    if show: 
        plt.show()

def calc_log2fc(dataframe, metadata, output=False):
    """
    input: dataframe (DataFrame): The DataFrame containing gene expression data.
    input: metadata (DataFrame): The metadata DataFrame containing sample information.
    output (bool): Whether to save the log2 fold change results to a CSV file (default is False).
    output: Series: A Series containing the log2 fold change (Log2FC) values

    Calculates the log2 fold change (Log2FC) between normal and tumor samples in the given DataFrame.
    """
    dataframe = dataframe.set_index("gene_id_ver")
    normal_runs = metadata[metadata["Tumor_Status"] == "Normal"]["Run"].values
    tumor_runs = metadata[metadata["Tumor_Status"] == "Tumor"]["Run"].values
    result_df = pd.DataFrame(index=dataframe.index, columns=["Normal", "Tumor", "Log2FC"])
    result_df["Normal"] = dataframe[normal_runs].mean(axis=1)
    result_df["Tumor"] = dataframe[tumor_runs].mean(axis=1)
    result_df["Log2FC"] = np.log2(result_df["Normal"]+1) - np.log2(result_df["Tumor"]+1)
    
    if output: result_df.to_csv("log2fc_results.csv", index=True)
    return result_df["Log2FC"]

# Main function to execute the script
def main():
    # Main function to execute the script

    # Ask if the wanna open and combine datasets or if the wanna open a premerged dataset
    merged_dataframes = None
    current_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    input_choice = '2' #input("Do you want to read and merge datasets (1) or use a premerged dataset (2)? Enter 1 or 2: ")
    if input_choice == '1':
        #directory = input("Enter the directory path: ")
        directory = os.path.join(current_path, "data")
        dfs = read_excel_file(directory)
        merged_dataframes = merge_dfs(dfs)
        if merged_dataframes is not None:
            output_file = os.path.join(directory, "merged_data.csv.gz")
            # Save the merged DataFrame to an excel file
            print(f"Saving merged DataFrame to {output_file}...")
            #merged_dataframes.to_excel(output_file, index=False)
            merged_dataframes.to_csv(output_file, index=False, compression="gzip")
            print(f"Merged DataFrame saved to {output_file}")
        else:
            print("No DataFrames to merge or an error occurred during merging.")
            return
    elif input_choice == '2':
        #dataframe_path = input("Enter the file path: ")
        dataframe_path = os.path.join(current_path,"data","merged_data.csv.gz")
        # Read the DataFrame from the specified path
        print(f"Reading DataFrame from {dataframe_path}...")
        # This assumes that the merged DataFrame has been saved as a excel file
        merged_dataframes = pd.read_csv(dataframe_path, compression="gzip")
        merged_dataframes = merged_dataframes.set_index("gene_id_ver")
        if merged_dataframes.empty:
            print("The DataFrame is empty. Please check the file path or the contents of the file.")
            return
    else:
        print("Invalid input. Please enter 1 or 2.")
        return
    
    if merged_dataframes is None or merged_dataframes.empty:
        print("No data available for further processing. Please check the input files.")
        return


    # Read the metadata from the specified path
    metadata_path = os.path.join(current_path, "data", "metadata.csv.gz")

    #metadata_path = input("Enter the metadata file path: ")
    print(f"Reading metadata from {metadata_path}...")

    # Die folgenden 3 Zeilen sind nur notwendig, wenn die metadaten nochmal neu gemerged werden müssen
    # Beachte: Weiter unten, wird die Metadatei angepasst, ebenfalls dann auskommentieren
    #Datasets = pd.read_csv(os.path.join(current_path, "Datasets.csv"))
    #metadata = pd.merge(Mega, Datasets, left_on="GEO_Series", right_on="GEO_Series", how="left")
    #metadata = metadata.rename(columns={"Run": "Run", "GEO_Series": "GEO_Series", "group": "group"})
    metadata = pd.read_csv(metadata_path, compression="gzip")

    if "GEO_Series" not in metadata.columns:
        print("Error: 'GEO_Series' column not found in metadata.")
        raise ValueError("'GEO_Series' column not found in metadata.")

    # Transform the DataFrame
    # Exclude specific GEO Series if needed
    # This will help in focusing on specific GEO Series for analysis
    excluded_series = [""] #GSE144269
    excluded_runs = metadata[metadata['GEO_Series'].isin(excluded_series)]['Run'].unique()
    print(f"Excluding runs from GEO Series: {excluded_series}")
    dataframe = merged_dataframes.drop(columns=excluded_runs, errors='ignore').reset_index()

    srr_columns = [col for col in dataframe.columns if col.startswith("SRR")]
    # Falls man sich für nur für den exludierten Run interessiert:
    #srr_columns = [col for col in excluded_runs if col.startswith("SRR")]

    # Folgende 4 bzw. 5 Zeilen erstellen die Counts Datei für R
    #dataframe["gene_id"] = dataframe["gene_id_ver"].str.split(".").str[0]
    #R_counts = dataframe[["gene_id"] + srr_columns].copy().set_index("gene_id")
    #R_counts = R_counts.apply(lambda x: np.log2(x + 1))  # Log2 transformation
    #R_counts.to_csv(os.path.join(current_path, "R-Project", "data", "merged_data_R.csv.gz"), index=True, compression="gzip")
    #breakpoint()

    dataframe_prepared = dataframe[["gene_id_ver"]+ srr_columns].copy().set_index("gene_id_ver")


    

    metadata = metadata[metadata['Run'].isin(srr_columns)].copy()
    # Anpassung der Metadaten, falls oben neu gemerged wird
    #metadata["Tumor_Status"] = metadata["group"].apply(label_group)
    #metadata = metadata[["Run", "Tumor_Status", "GEO_Series"]].copy()
    #metadata.to_csv(os.path.join(current_path, "data", "metadata.csv.gz"), index=False, compression="gzip")

    # Speichern der Metadaten für R
    #R_meta = metadata[["Run", "Tumor_Status"]].copy().set_index("Run")
    #R_meta.to_csv(os.path.join(current_path, "R-Project", "data", "merged_metadata.csv"), index=True)


    # Apply log2 transformation to the expression data
    # This transformation is commonly used in gene expression analysis to normalize the data
    dataframe_log = dataframe_prepared.apply(lambda x: np.log2(x + 1))

    # Apply z-score normalization to the log2 transformed DataFrame
    # This normalization will standardize the data to have mean=0 and variance=1
    # This is important for downstream analyses such as PCA
    dataframe_scaled = dataframe_log.apply(zscore, axis=0)
    print("DataFrame prepared successfully.")

    output_path = os.path.join(current_path, "output")
    
    # Create a heatmap from the DataFrame
    #heatmap(dataframe_log, output_path, output_name="heatmap_with_outlier", top_var_genes_num=1000, method='seaborn')
    #heatmap(dataframe_log, output_path, output_name="heatmap_without_outlier", top_var_genes_num=1000, method='seaborn')
    
    # Perform PCA on the DataFrame
    #pca(dataframe_scaled, metadata, output_path, output_name="pca_with_outlier", show_geo=True)
    #pca(dataframe_scaled, metadata, output_path, output_name="pca_without_outlier", show_geo=False)

    # Create a distribution plot from the DataFrame
    #distribution_plot(dataframe_scaled, metadata, output_path, output_name="distribution_plot_with_outlier")
    #distribution_plot(dataframe_scaled, metadata, output_path, output_name="distribution_plot_without_outlier")

    # Calculate mean expression values
    #mean_expression(dataframe_scaled, metadata, output_path, output_name="mean_expression_with_outlier")
    #mean_expression(dataframe_scaled, metadata, output_path, output_name="mean_expression_without_outlier")

    # Perform logistic regression on the DataFrame
    logreg_result = logistic_regression(dataframe_scaled, metadata)
    logreg_df = pd.DataFrame({
        'Gene': logreg_result['features'],
        'Coefficient': logreg_result['coefficients'],
        'Odds Ratio': logreg_result['odds_ratios']
    })
    logreg_df = logreg_df.sort_values(by='Coefficient', ascending=False)
    print(f"Top 10 genes with highest coefficients:\n{logreg_df.head(10)}")
    #print(f"Top 10 genes with lowest coefficients:\n{logreg_df.sort_values(by='Coefficient', ascending=True).head(10)}")
    #plot_logistic_regression(logreg_result, output_path, output_name="logreg_with_outlier")
    #plot_logistic_regression(logreg_result, output_path, output_name="logreg_without_outlier")

    # Perform decision tree analysis on the DataFrame
    # Hat in der Arbeit keiner relevanz gehabt/wurde nicht genutzt
    #decision_tree_result = decision_tree(dataframe_scaled, metadata, method='gini')
    ######
    #plot_decision_tree(decision_tree_result, output_path, output_name="decision_tree_gini_with_outlier")
    #plot_decision_tree(decision_tree_result, output_path, output_name="decision_tree_gini_without_outlier")

##############

    # Perform random forest classification on the DataFrame
    rf_result = random_forest_classification(dataframe_scaled, metadata)
    rf_df = pd.DataFrame({
        'Gene': rf_result['features'],
        'Feature Importance': rf_result['Feature Importance']
    })
    rf_df = rf_df.sort_values(by='Feature Importance', ascending=False)
    print(f"Top 10 genes with highest feature importance:\n{rf_df.head(10)}")

    # Export the results to CSV files
    log2fc = calc_log2fc(dataframe[["gene_id_ver"] + srr_columns], metadata)
    log2fc.index = log2fc.index.str.split(".").str[0]
    logreg = logreg_df.sort_values(by='Coefficient', ascending=False)['Gene'].reset_index(drop=True).str.split(".").str[0]
    randforest = rf_df.sort_values(by='Feature Importance', ascending=False)['Gene'].reset_index(drop=True).str.split(".").str[0]
    
    # ml Features für R exportieren
    #ml_features = pd.DataFrame({
    #    'logistic_regression': logreg,
    #    'logreg_log2fc': log2fc[logreg].values,
    #    'random_forest_classification': randforest,
    #    'rf_log2fc': log2fc[randforest].values
    #})
    #print(ml_features.head(10))
    #ml_features.to_csv(os.path.join(current_path, "R-Project", "data", "ml_features.csv"), index=False)

    # Perform feature importance analysis
    #plot_feature_importance(logreg_result, output_path, output_name="fi_with_outlier", method="logistic_regression")
    #plot_feature_importance(logreg_result, output_path, output_name="fi_without_outlier", method="logistic_regression")
    #plot_feature_importance(rf_result, output_path, output_name="fi_with_outlier", method="random_forest_classification")
    #plot_feature_importance(rf_result, output_path, output_name="fi_without_outlier", method="random_forest_classification")

    # Perform permutation importance analysis
    # Hat nicht funktioniert, wurde nicht genutzt
    #plot_perm_importance(logreg_result, output_path, output_name="permutation_importance_logistic_regression", show=True)



if __name__ == '__main__':
    main()


    # Sehr viel manuelles für Venn-Diagramme und Korrelation, teilweise abhängig von Dateien in R erzeugt
    """
    directory = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
    GSE = ["GSE25599", "GSE77314", "GSE82177", "GSE101432", "GSE105130", "GSE114564", "GSE144269", "GSE148355"]
    filename = "all.samples.tpm.xlsx"  # Name der Datei, die in jedem Ordner eingelesen werden soll
    Genes_datasets = pd.DataFrame()
    
    for gse in GSE:
        path = os.path.join(directory, gse, filename)
        print(f"Reading file: {path}")
        df = pd.read_excel(path)
        df["gene_id"] = df['gene_id_ver'].str.split(".").str[0]
        expr_cols = df.columns.difference(["trans_id_ver", "trans_name", "gene_id_ver", "gene_name", "gene_id"])
        df = df[["gene_id"]].copy().assign(mean_expr=df[expr_cols].mean(axis=1, skipna=True))
        df = df.groupby("gene_id", as_index=False)["mean_expr"].sum()
        df = df.set_index("gene_id")
        df = df.rename(columns={"mean_expr": gse})
        if Genes_datasets.empty:
            Genes_datasets = df
        else:
            Genes_datasets = Genes_datasets.join(df, how='outer')
        print(Genes_datasets[gse][:5])  # Print the first 5 gene IDs to verify

    print(Genes_datasets.head())
    Genes_datasets.to_csv(os.path.join(directory, "Genes_datasets.csv"), index=True)
    """


    """
    Genes = pd.read_csv("/home/laszl/Praxisphase/Phase 1/data/Genes.csv")
    Genes_datasets = pd.read_csv(os.path.join(directory, "Genes_datasets.csv"), index_col=0)
    Pathways = pd.read_csv(os.path.join(directory, "pathways.csv"))
    Quality_Paths = pd.read_csv(os.path.join(directory, "pathways_positively_correlated_to_low_quality.tsv"), sep="\t")
    Quality_Paths = Quality_Paths["pathway"].tolist()
    print(Genes_datasets.head())
    gene_mean_dict = Genes_datasets.mean(axis=1, skipna=True).to_dict()
    DGEWith = Genes["DGEWith"].tolist() 
    BatchWith = Genes["BatchWith"].tolist() 
    LogRegWith = Genes["LogRegWith"].tolist() 
    RandForestWith = Genes["RandForestWith"].tolist() 
    DGEWithout = Genes["DGEWithout"].tolist() 
    BatchWithout = Genes["BatchWithout"].tolist() 
    LogRegWithout = Genes["LogRegWithout"].tolist() 
    RandForestWithout = Genes["RandForestWithout"].tolist() 
    all_genes_with = DGEWith + BatchWith + LogRegWith + RandForestWith
    all_genes_without = DGEWithout + BatchWithout + LogRegWithout + RandForestWithout
    Genes_with = [DGEWith, BatchWith, LogRegWith, RandForestWith]
    Genes_without = [DGEWithout, BatchWithout, LogRegWithout, RandForestWithout]
    common_with = list(set.intersection(*map(set, Genes_with)))
    common_without = list(set.intersection(*map(set, Genes_without)))
    print("Common genes with outlier:", len(common_with))
    print("Common genes without outlier:", len(common_without))
    common_DGE_BatchCor_with = list(set(DGEWith).intersection(set(BatchWith)) - set(common_with))
    common_DGE_BatchCor_without = list(set(DGEWithout).intersection(set(BatchWithout)) - set(common_without))
    common_Forest_LogReg_with = list(set(RandForestWith).intersection(set(LogRegWith)) - set(common_with))
    common_Forest_LogReg_without = list(set(RandForestWithout).intersection(set(LogRegWithout)) - set(common_without))
    common_Forest_Logreg_BatchCor_with = list(set(RandForestWith).intersection(set(LogRegWith)).intersection(set(BatchWith)) - set(common_with))
    common_Forest_Logreg_BatchCor_without = list(set(RandForestWithout).intersection(set(LogRegWithout)).intersection(set(BatchWithout)) - set(common_without))
    common_DGE_Forest_Logreg_with = list(set(DGEWith).intersection(set(RandForestWith)).intersection(set(LogRegWith)) - set(common_with))
    common_DGE_Forest_Logreg_without = list(set(DGEWithout).intersection(set(RandForestWithout)).intersection(set(LogRegWithout)) - set(common_without))

    counts_with = Counter(all_genes_with)
    counts_without = Counter(all_genes_without)
    UniqueDGEWith = [x for x in DGEWith if counts_with[x] == 1]
    UniqueBatchWith = [x for x in BatchWith if counts_with[x] == 1]
    UniqueLogRegWith = [x for x in LogRegWith if counts_with[x] == 1]
    UniqueRandForestWith = [x for x in RandForestWith if counts_with[x] == 1]
    UniqueDGEWithout = [x for x in DGEWithout if counts_without[x] == 1]
    UniqueBatchWithout = [x for x in BatchWithout if counts_without[x] == 1]
    UniqueLogRegWithout = [x for x in LogRegWithout if counts_without[x] == 1]
    UniqueRandForestWithout = [x for x in RandForestWithout if counts_without[x] == 1]

    complete_overlap = list(set(common_with).intersection(set(common_without)))
    complete_difference = list(set(common_with).union(set(common_without)) - set(complete_overlap))
    all_genes_overlap = list(set(all_genes_with).intersection(set(all_genes_without)))
    all_genes_difference = list(set(all_genes_with).union(set(all_genes_without)) - set(all_genes_overlap))
    print("Complete overlap genes:", len(complete_overlap))
    print("Complete difference genes:", len(complete_difference))
    #print(len(common_with), len(common_without), len(common_DGE_BatchCor_with), len(common_DGE_BatchCor_without))
    #print(len(common_Forest_LogReg_with), len(common_Forest_LogReg_without), len(common_Forest_Logreg_BatchCor_with), len(common_Forest_Logreg_BatchCor_without))
    #print(len(common_DGE_Forest_Logreg_with), len(common_DGE_Forest_Logreg_without))

    #Batch_overlap = list(set(BatchWith).intersection(set(BatchWithout)))

    # Correlation between different GEO-Series
    #GEO = Genes_datasets.corr()
    #sns.heatmap(GEO, annot=True, cmap="coolwarm", center=0)
    #plt.title("Correlation between different GEO Series")
    #plt.tight_layout()
    #plt.show()
    # Welche Pathways sind identisch
    for path in Quality_Paths:
        if path in Pathways["complete_overlap"].values:
            print(f"Pathway {path} is in complete_overlap")
        if path in Pathways["complete_difference"].values:
            print(f"Pathway {path} is in complete_difference")
        if path in Pathways["DGE_difference"].values:
            print(f"Pathway {path} is in pathways_DGE_difference")
        if path in Pathways["BatchCor_difference"].values:
            print(f"Pathway {path} is in pathways_BatchCor_difference")
        if path in Pathways["LogReg_overlap"].values:
            print(f"Pathway {path} is in LogReg_overlap")
        if path in Pathways["LogReg_difference"].values:
            print(f"Pathway {path} is in LogReg_difference")
        if path in Pathways["RandForest_overlap"].values:
            print(f"Pathway {path} is in RandForest_overlap")
        if path in Pathways["RandForest_difference"].values:
            print(f"Pathway {path} is in RandForest_difference")

    # Correlation between different GEO-Series and variables
    DGE_overlap = list(set(UniqueDGEWith).intersection(set(UniqueDGEWithout)))
    print(len(DGE_overlap))
    DGE_difference = list(set(UniqueDGEWith).union(set(UniqueDGEWithout)) - set(set(UniqueDGEWith).intersection(set(UniqueDGEWithout))))
    Batch_overlap = list(set(UniqueBatchWith).intersection(set(UniqueBatchWithout)))
    print(len(Batch_overlap))
    Batch_difference = list(set(UniqueBatchWith).union(set(UniqueBatchWithout)) - set(set(UniqueBatchWith).intersection(set(UniqueBatchWithout))))
    LogReg_overlap = list(set(UniqueLogRegWith).intersection(set(UniqueLogRegWithout)))
    LogReg_difference = list(set(UniqueLogRegWith).union(set(UniqueLogRegWithout)) - set(set(UniqueLogRegWith).intersection(set(UniqueLogRegWithout))))
    RandForest_overlap = list(set(UniqueRandForestWith).intersection(set(UniqueRandForestWithout)))
    RandForest_difference = list(set(UniqueRandForestWith).union(set(UniqueRandForestWithout)) - set(set(UniqueRandForestWith).intersection(set(UniqueRandForestWithout))))
    #print(len(DGE_difference))
    #print(len(Batch_difference))
    #print(len(LogReg_overlap))
    #print(len(LogReg_difference))
    #print(len(RandForest_overlap))
    #print(len(RandForest_difference))
    correlation_variables_unique = {
        "DGE (A)": DGE_overlap,
        "DGE (B)": DGE_difference,
        "BatchCor (A)": Batch_overlap,
        "BatchCor (B)": Batch_difference,
        "LogReg (A)": LogReg_overlap,
        "LogReg (B)": LogReg_difference,
        "RandForest (A)": RandForest_overlap,
        "RandForest (B)": RandForest_difference,
    }
    correlation_variables_methods = {
        "DGE (A)": DGEWith,
        "DGE  (B)": DGEWithout,
        "BatchCor (A)": BatchWith,
        "BatchCor  (B)": BatchWithout,
        "LogReg (A)": LogRegWith,
        "LogReg  (B)": LogRegWithout,
        "RandForest (A)": RandForestWith,
        "RandForest  (B)": RandForestWithout,
        "Common (A∩B)": complete_overlap,
        "Uncommon (AΔB)": complete_difference
    }
    correlation_output_unique = pd.DataFrame()
    correlation_output_methods = pd.DataFrame()
    for variable_name, variable in correlation_variables_unique.items():
        for gene_id in variable:
            if gene_id not in gene_mean_dict:
                print(f"Gene ID {gene_id} not found in gene_mean_dict.")
        Gen_set = Genes_datasets.copy()
        Gen_set[variable_name] = pd.Series({gene_id: gene_mean_dict[gene_id] for gene_id in variable})
        GEO = Gen_set.corr()
        correlation_output_unique[variable_name] = GEO[variable_name].drop(variable_name)
    for variable_name, variable in correlation_variables_methods.items():
        for gene_id in variable:
            if gene_id not in gene_mean_dict:
                print(f"Gene ID {gene_id} not found in gene_mean_dict.")
        Gen_set = Genes_datasets.copy()
        Gen_set[variable_name] = pd.Series({gene_id: gene_mean_dict[gene_id] for gene_id in variable})
        GEO = Gen_set.corr()
        correlation_output_methods[variable_name] = GEO[variable_name].drop(variable_name)
    #variable = Batch_difference
    #variable_name = 'Batch_difference'
    #variable_title = 'the different Genes identified by Batch Correction with and without the outlier dataset'
    #for gene_id in variable:
    #    if gene_id not in gene_mean_dict:
    #        print(f"Gene ID {gene_id} not found in gene_mean_dict.")
    #Genes_datasets[variable_name] = pd.Series({gene_id: gene_mean_dict[gene_id] for gene_id in variable})
    #GEO = Genes_datasets.corr()
    #print(GEO[variable_name].drop(variable_name))
    #print(GEO[variable_name])
    output = pd.concat([correlation_output_unique, correlation_output_methods], axis=1)
    #sns.set(font_scale=1.4)
    fig, ax = plt.subplots(figsize=(26, 6))
    sns.heatmap(correlation_output_unique.T, annot=True, cmap="coolwarm", center=0, annot_kws={"size": 30, "weight": "bold", "color": "black"})
    plt.title("Method-Specific Feature Sets - Consistently (A) vs. Variably (B) Identified", fontsize=30, weight="bold")
    plt.xticks(rotation=0, fontsize=23, weight="bold")
    plt.yticks(fontsize=30, weight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(directory, "correlation_unique.png"))
    #plt.show()
    fig, ax = plt.subplots(figsize=(26, 7.5))
    sns.heatmap(correlation_output_methods.T, annot=True, cmap="coolwarm", center=0, annot_kws={"size": 30, "weight": "bold", "color": "black"})
    plt.title("Method-Derived Feature Sets - With (A) vs. Without (B) Outlier", fontsize=30, weight="bold")
    plt.xticks(rotation=0, fontsize=23, weight="bold")
    plt.yticks(fontsize=30, weight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(directory, "correlation_method.png"))
    #plt.show()
    """

    """
    variable_name = 'all_genes_difference'
    variable = all_genes_difference
    blob = []
    for gse in Genes_datasets:
        blob += list(Genes_datasets[gse].dropna())
    myset = set(blob)
    converter = {}
    counter = 0
    for item in myset:
        converter[item] = counter
        counter += 1
    converted_Geneset = {}
    Genes_datasets_corr = pd.DataFrame(Genes_datasets)
    Genes_datasets_corr[variable_name] = pd.Series(variable)
    for gse in Genes_datasets_corr:
        converted_Geneset[gse] = []
        for gene in Genes_datasets_corr[gse].dropna():
            converted_Geneset[gse].append(converter[gene])

    converted_Geneset = pd.DataFrame.from_dict(converted_Geneset, orient='index').transpose()
    corr = converted_Geneset.corr()
    sns.heatmap(corr, annot=True, cmap="coolwarm", center=0)
    plt.title("Correlation between different GEO Series and " + variable_name)
    plt.tight_layout()
    plt.show()
    print(corr)
    
    overlap = {}
    for gse in Genes_datasets:
        a = set(variable)
        b = set(Genes_datasets[gse].dropna())
        overlap[gse] = len(a&b)/len(a|b)
  
    sns.barplot(x=list(overlap.keys()), y=list(overlap.values()))
    plt.ylabel("Jaccard Index")
    plt.xlabel("GEO Series")
    plt.show()

    print(len(UniqueDGEWith), len(UniqueBatchWith), len(UniqueLogRegWith), len(UniqueRandForestWith))
    print(len(UniqueDGEWithout), len(UniqueBatchWithout), len(UniqueLogRegWithout), len(UniqueRandForestWithout))
    """
    