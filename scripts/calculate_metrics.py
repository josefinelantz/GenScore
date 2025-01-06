import numpy as np
import pandas as pd
from sklearn.metrics import precision_score, recall_score, f1_score

def prepare_data(filtered_df, use_col="GROUP", groups=["benign", "pathogenic"]):
    """ Retains rows that have GROUP in groups. 
    Extracts true labels for groups as binary. 
    Args:
        filtered_df (pd.DataFrame): the data to extract true_labels from. 
        use_col (str): column to use for generating true labels
        groups (list[str]): the groups to generate true labels for, where positive group is at index 1. 
    Returns:
        filtered_data (pd.DataFrame): The filtered dataframe. 
        y_scores (pd.Series): The scores to be used to classify variants, y_pred.
        y_true (pd.Series): The binary labels
    """
    # Filter the dataset for benign and pathogenic variants
    filtered_data = filtered_df[filtered_df[use_col].isin(groups)].copy()

    # Generate true labels
    true_labels = (filtered_df["GROUP"] == "pathogenic").astype(int)

    print(f"Filtered dataset has {len(filtered_data)} rows.") #Filtered dataset has 1381 rows.
    return filtered_data, true_labels

def calculate_metrics_by_threshold(scores, labels, thresholds):
    """Calculate precision, recall, and F1-score for various thresholds.
    Args:
        scores (pd.Series): Predicted scores or probabilities.
        labels (pd.Series): True binary labels (0 or 1).
        thresholds (iterable): Thresholds to evaluate.
    Returns:
        pd.DataFrame: A DataFrame containing metrics for each threshold.
    """
    
    results = []
    for threshold in thresholds:
        # Generate predictions based on the threshold
        y_pred = (scores >= threshold).astype(int)
        
        # Calculate metrics, handle undefined precision with zero_division=0
        precision = precision_score(labels, y_pred, zero_division=0)
        recall = recall_score(labels, y_pred)
        f1 = f1_score(labels, y_pred, zero_division=0)
        results.append({"threshold": threshold, "precision": precision, "recall": recall, "f1-Score": f1})
    return pd.DataFrame(results)


def calculate_optimal_threshold_metrics(df, score_column="REFINED_RANK_SCORE", optimal_threshold=20):
    """
    Calculate precision, recall, F1-score, and AUC at the optimal threshold.
    """
    precision, recall, f1, auc = calculate_performance_metrics(
        df, score_column=score_column, threshold=optimal_threshold
    )
    print(f"Optimal Threshold Metrics (Threshold = {optimal_threshold}):")
    print(f"Precision: {precision:.2f}, Recall: {recall:.2f}, F1-Score: {f1:.2f}, AUC: {auc:.2f}")
    return precision, recall, f1, auc

# 2. Analyze False Positives and Negatives
def analyze_misclassifications(df, score_column="REFINED_RANK_SCORE", optimal_threshold=20):
    """
    Identify false positives and false negatives at the optimal threshold.
    """
    df["y_pred"] = df[score_column].apply(lambda x: 1 if x >= optimal_threshold else 0)
    
    false_positives = df[(df["y_pred"] == 1) & (df["y_true"] == 0)]
    false_negatives = df[(df["y_pred"] == 0) & (df["y_true"] == 1)]

    print(f"Number of False Positives: {len(false_positives)}")
    print(f"Number of False Negatives: {len(false_negatives)}")
    
    return false_positives, false_negatives

# 3. Refine Feature Weights Based on Misclassifications
def refine_weights_based_on_misclassifications(df, false_positives, false_negatives):
    """
    Analyze the feature contributions in false positives and false negatives to guide reweighting.
    """
    print("\nFeature Contributions in False Positives:")
    print(false_positives[["AF", "PP", "CON", "VCQF", "LIN", "CLIN"]].describe())

    print("\nFeature Contributions in False Negatives:")
    print(false_negatives[["AF", "PP", "CON", "VCQF", "LIN", "CLIN"]].describe())

    # Recommendations for reweighting can be derived from these summaries

# Example Usage
def main_optimal_threshold_analysis(df):
    # Step 1: Calculate metrics at the optimal threshold
    optimal_threshold = 20
    precision, recall, f1, auc = calculate_optimal_threshold_metrics(df, optimal_threshold=optimal_threshold)

    # Step 2: Analyze false positives and false negatives
    false_positives, false_negatives = analyze_misclassifications(df, optimal_threshold=optimal_threshold)

    # Step 3: Refine weights based on misclassifications
    refine_weights_based_on_misclassifications(df, false_positives, false_negatives)

    return df, false_positives, false_negatives

# Uncomment and run when `df` and `y_true` exist
# df, false_positives, false_negatives = main_optimal_threshold_analysis(df)