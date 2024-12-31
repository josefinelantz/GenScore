def calculate_metrics(y_true, y_pred, y_scores):
    metrics = {
        "Precision": precision_score(y_true, y_pred, average="binary"),
        "Recall": recall_score(y_true, y_pred, average="binary"),
        "F1-Score": f1_score(y_true, y_pred, average="binary"),
        "AUC": roc_auc_score(y_true, y_scores)
    }
    return metrics

def plot_roc_curve(y_true, y_scores):
    fpr, tpr, _ = roc_curve(y_true, y_scores)
    plt.plot(fpr, tpr, label="ROC Curve")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC Curve")
    plt.legend()
    plt.show()


# 1. Set the Optimal Threshold and Calculate Metrics
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