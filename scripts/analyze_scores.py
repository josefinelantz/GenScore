import pandas as pd 

def analyze_scores(df):
    """
    Returns the df with column GROUP
    """
    # Group variants into benign/pathogenic
    benign_terms = {"Benign", "Likely_benign"}
    pathogenic_terms = {"Pathogenic", "Likely_pathogenic"}
    
    def classify_group(row):
        if row["CLNSIG"] in benign_terms:
            return "benign"
        elif row["CLNSIG"] in pathogenic_terms:
            return "pathogenic"
        return "control"

    df["GROUP"] = df.apply(classify_group, axis=1)
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

