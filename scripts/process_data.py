import pandas as pd 

from constants import GROUPS_WITH_LABELS

def group_variants(df, from_col="CLNSIG", groups=GROUPS_WITH_LABELS, to_col="GROUP", controls="controls_match.tsv"):
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
    def mark_controls(df, controls):
        c = pd.read_csv(controls, sep="\t")
        control_set = set(c["chrom_pos"])
        df["IS_CONTROL"] = df["CHROM_POS"].isin(control_set)
        return df
    
    def classify(clnsig):
        for group, labels in groups.items():
            if clnsig in labels:
                return group
        return "other"  # Return original label if no match is found
    
    df["ADJUSTED_SCORE"] = df["RANK_SCORE"] - df["CLIN"]

    # Apply the classification logic to the DataFrame
    df[to_col] = df[from_col].apply(classify)
    
    return mark_controls(df, controls)

