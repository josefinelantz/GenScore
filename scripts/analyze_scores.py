import pandas as pd 

def analyze_scores(df):
    # Group variants into benign/pathogenic
    benign_terms = {"Benign", "Likely_benign"}
    pathogenic_terms = {"Pathogenic", "Likely_pathogenic"}
    def classify_group(row):
        if row["CLNSIG"] in benign_terms:
            return "benign"
        elif row["CLNSIG"] in pathogenic_terms:
            return "pathogenic"
        return "control"

    def assign_color(group):
        return {"benign": "blue", "pathogenic": "red", "control": "yellow"}[group]

    df["GROUP"] = df.apply(classify_group, axis=1)
    df["COLOR"] = df["GROUP"].apply(assign_color)
    return df

def mark_controls(df, control_tsv_path):
    controls = pd.read_csv(control_tsv_path, sep="\t")
    controls["chr_pos"] = controls["Chrom"].astype(str) +"_"+ controls["Position"].astype(str)
    controls.set_index("chr_pos", inplace=True)
    positive_controls = list(controls.index)
    df["is_control"] = df["CHROM_POS"].isin(positive_controls)
    #control_set = set(controls["CHROM_POS"])
    #df["IS_CONTROL"] = df["CHROM_POS"].isin(control_set)
    return df

def subtract_clin_score(df):
    df["NO_CLIN_RANK_SCORE"] = df["RANK_SCORE"] - df["CLIN"]
    return df

