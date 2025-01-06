import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, classification_report
from scripts.process_data import classify_variants

def plot_scatter_with_controls(df, score_col="RANK_SCORE", group_col="GROUP", control_col="IS_CONTROL", output_path="scatter_by_group_with_controls.png"):
    """
    Creates a scatter plot with jitter, coloring by group, and highlighting controls.

    Args:
    - df (pd.DataFrame): The DataFrame containing the data.
    - score_col (str): Column for the score to plot on the y-axis.
    - group_col (str): Column for the groups to color points.
    - control_col (str): Column indicating control variants (boolean).
    """
    # Generate jitter for x-axis
    jitter = np.random.uniform(-1.0, 1.0, size=len(df))

    # Assign colors for groups
    group_colors = {
        "benign": "blue",
        "pathogenic": "red",
        "other": "gray",
        "uncertain": "orange",
    }

    # Create the scatter plot
    plt.figure(figsize=(10, 6))
    for group, color in group_colors.items():
        subset = df[df[group_col] == group]
        plt.scatter(
            x=jitter[subset.index],
            y=subset[score_col],
            s=60,
            label=group.capitalize(),
            color=color,
            alpha=0.6,
            edgecolor="black" if control_col in subset.columns else None,
            linewidth=0.2,
        )

    # Highlight controls
    controls = df[df[control_col]]
    plt.scatter(
        x=jitter[controls.index],
        y=controls[score_col],
        facecolors="#9AFF60",
        edgecolors="black",
        label="Control",
        s=120,
    )

    # Plot formatting
    plt.axhline(0, color="black", linewidth=0.8, linestyle="--")
    plt.xlabel("Variants (Jittered)")
    plt.ylabel(score_col)
    plt.title("Scatter Plot of Variants by Group")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)

def plot_violin_adjusted_scores_by_group(filtered_df, score_col="ADJUSTED_SCORE", group_col="GROUP", output_path="violin_adjusted_scores_by_group.png"):
    """
    Plot RankScore for benign and pathogenic variants. 
    Parameters:
        - df (pd.DataFrame): DataFrame with columns "ADJUSTED_SCORE" and "GROUP"
    """
    plt.figure(figsize=(8, 6))
    sns.violinplot(
        data=filtered_df, 
        hue=group_col,
        y=score_col
        )
    plt.title("Violin Plot of Adjusted Scores Between Benign and Pathogenic Groups")
    plt.ylabel("Adjusted Score")
    plt.grid() 
    plt.legend(loc="lower right")
    plt.savefig(output_path)

def plot_density_adjusted_scores_by_group(filtered_df, score_col="ADJUSTED_SCORE", group_col="GROUP", output_path="density_plot_adjusted_scores_by_group.png"):
    plt.figure(figsize=(10, 6))
    sns.kdeplot(
        data=filtered_df,
        x=score_col,
        hue=group_col,
        fill=True,
        common_norm=False,
        alpha=0.7
    )
    plt.title("Density Plot of Adjusted Scores - Benign and Pathogenic Groups")
    plt.xlabel("Adjusted Score")
    plt.ylabel("Density")
    plt.grid()
    plt.savefig(output_path)

def visualize_metrics(metrics_df, output_path="visualize_threshold_metrics.png"):
    """
    Visualizes precision, recall, and F1-score across thresholds.
    Args:
        metrics_df (pd.DataFrame): A DataFrame containing 'threshold', 'precision', 'recall', and 'f1_score'.
    """
    plt.figure(figsize=(10, 6))
    plt.plot(metrics_df["Threshold"], metrics_df["Precision"], label="Precision", color="blue")
    plt.plot(metrics_df["Threshold"], metrics_df["Recall"], label="Recall", color="orange")
    plt.plot(metrics_df["Threshold"], metrics_df["F1-Score"], label="F1-Score", color="green")
    plt.axvline(metrics_df.loc[metrics_df["fit_score"].idmax(), "threshold"], color="red", linestyle="--", label="Optmimal Threshold")
    plt.xlabel("Threshold")
    plt.ylabel("Metric Score")
    plt.title("Precision, Recall, and F1-Score vs Threshold")
    plt.legend()
    plt.grid(True)
    plt.savefig(output_path)

def plot_confusion_matrix(df, threshold, output_path="confusion_matrix.png"):
    # Classify variants
    df = classify_variants(df, threshold)
    
    # Extract true and predicted labels
    y_true = df["y_true"]
    y_pred = df["y_pred"]
    
    # Generate confusion matrix
    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=["Benign", "Pathogenic"])
    
    # Plot confusion matrix
    disp.plot(cmap="Blues", values_format="d")
    disp.ax_.set_title(f"Confusion Matrix (Threshold = {threshold})")
    
    # Save as PNG
    plt.title(f"Confusion Matrix (Threshold = {threshold}")
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.show()
    
    # Print classification report
    print("Classification Report:")
    print(classification_report(y_true, y_pred, target_names=["Benign", "Pathogenic"]))
    
    return cm

def plot_feature_contributions(melted_df, output_path="feature_contributions_before_reweighting.png"):
    # Plot feature contributions
    plt.figure(figsize=(12, 6))
    sns.boxplot(data=melted_df, x="Feature", y="Score", hue="GROUP", palette="Set2")
    plt.title("Feature Contributions by Group")
    plt.xlabel("Feature")
    plt.ylabel("Score")
    plt.legend(title="Group")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.show()
