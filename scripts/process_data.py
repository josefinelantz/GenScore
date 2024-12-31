import pandas as pd 

from constants import GROUPS_WITH_LABELS

def group_variants(df, from_col="CLNSIG", groups=GROUPS_WITH_LABELS, to_col="GROUP"):
    """ 
    Groups variants by 'groups' based on information in 'from_col'.
    Args:
    df (pd.DataFrame): data with variants to group
    from_col (str): the column in df to use for grouping (required in df)
    to_col (str): the name of the new column that will contain variant group labels
    Returns:
    df (pd.DataFrame): the df extended with 'to_column' values. 
    """
    def classify(clnsig):
        for group in groups.keys():
            if clnsig in groups[group]:
                return group
        return clnsig  # Return original label if no match is found

    # Apply the classification logic to the DataFrame
    df[to_col] = df[from_col].apply(classify)
    return df

def mark_controls(df, control_tsv_path):
    """
    Matches controls in VCF
    Assumes control_tsv_path has column:
    chr_pos
    Returns df with column IS_CONTROL
    """
    controls = pd.read_csv(control_tsv_path, sep="\t")
    controls["chr_pos"] = controls["Chrom"].astype(str) +"_"+ controls["Position"].astype(str)
    controls.set_index("chr_pos", inplace=True)
    positive_controls = list(controls.index)
    df["is_control"] = df["CHROM_POS"].isin(positive_controls)
    return df

def subtract_clin_score(df):
    """
    Returns df with column NO_CLIN_RANK_SCORE
    """
    df["NO_CLIN_RANK_SCORE"] = df["RANK_SCORE"] - df["CLIN"]
    return df

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