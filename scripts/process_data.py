import pandas as pd

from constants import GROUPS_WITH_LABELS, CONTROLS

def group_variants(df, from_col="CLNSIG", groups=GROUPS_WITH_LABELS, to_col="GROUP", controls=CONTROLS):
    """ 
    Groups variants by 'groups' based on information in 'from_col'.
    Adjusts RANK_SCORE by subtracting CLIN Score to mitigate bias.
    Args:
    df (pd.DataFrame): data with variants to group
    from_col (str): the column in df to use for grouping (required in df)
    to_col (str): the name of the new column that will contain variant group labels
    Returns:
    df (pd.DataFrame): the df extended with 'to_column' values and 'ADJUSTED_SCORE'. 
    """
    df[from_col] = df[from_col].str.strip()
    
    def classify(clnsig):
        for group, labels in groups.items():
            if clnsig in labels:
                return group
        print(f"Unmatched CLNSIG value: {clnsig}")  # Log unmatched values
        return "other"  # Return original label if no match is found

    # Apply the classification logic to the DataFrame
    df[to_col] = df[from_col].apply(classify)
    
    #df["ADJUSTED_SCORE"] = df["RANK_SCORE"] - df["CLIN"]
    
    return mark_controls(df, controls)

def mark_controls(df, controls):
    c = pd.read_csv(controls, sep="\t")
    control_set = set(c["chrom_pos"])
    df["IS_CONTROL"] = df["CHROM_POS"].isin(control_set)
    return df

def filter_data_and_adjust_scores(df, groups=["benign", "pathogenic"]):
    # retain rows for benign and pathogenic groups
    filtered_df = df[df["GROUP"].isin(groups)].copy()

    filtered_df["ADJUSTED_SCORE"] = filtered_df["RANK_SCORE"] - filtered_df["CLIN"]

    # Summarize the filtered dataset
    filtered_summary = {
        "Total Rows": filtered_df.shape[0],
        "Group Counts": filtered_df['GROUP'].value_counts().to_dict(),
        "Adjusted Score Range": (filtered_df['ADJUSTED_SCORE'].min(), filtered_df['ADJUSTED_SCORE'].max())
    }
    print(filtered_summary)
    return filtered_df

def classify_variants(df, threshold):
    """Generates predictions based on optimal threshold.
    """
    df = df.copy()
    df["y_pred"] = (df["ADJUSTED_SCORE"] >= threshold).astype(int)
    df["y_true"] = df["GROUP"].map({"benign": 0, "pathogenic": 1})
    return df

def melt_data(df):
    melted_data = df.melt(
        id_vars=["GROUP"], 
        value_vars=["AF", "PP", "CON", "VCQF", "LIN", "CLIN"], 
        var_name="Feature", 
        value_name="Score"
    )
    return melted_data