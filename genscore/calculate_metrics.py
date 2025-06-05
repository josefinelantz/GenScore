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

        print(f"Optimal Threshold Metrics (Threshold = {threshold}):")
        print(f"Precision: {precision:.2f}, Recall: {recall:.2f}, F1-Score: {f1:.2f}")
    return pd.DataFrame(results)

def recalculate_rank_scores(df, weights):
    df["UPDATED_RANK_SCORE"] = (
        df["AF"] * weights["AF"] +
        df["PP"] * weights["PP"] +
        df["CON"] * weights["CON"] +
        df["VCQF"] * weights["VCQF"] +
        df["LIN"] * weights["LIN"] +
        df["CLIN"] * weights["CLIN"]
    )
    return df

