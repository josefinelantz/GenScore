# scripts/analyze_scores.py

import pandas as pd


ClinVar = ["Benign", 
          "Pathogenic", 
          "Likely_benign", 
          "Likely_pathogenic", 
          "Benign/Likely_benign", 
          "Pathogenic/Likely_pathogenic", 
          "Uncertain_significance", 
          "not_reported", 
          "CONFLICTING_classifications_of_pathogenicity", 
          "drug_response"
          ]

def group_variants(df):
    benign = df[df["CLNSIG"].isin(["Benign", "Likely_benign", "Benign/Likely_benign"])]
    #benign = df[df["CLNSIG"].isin(["benign", "likely_benign"])]
    pathogenic = df[df["CLNSIG"].isin(["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"])]
    not_reported = df[df["CLNSIG"].isin(["not_reported"])]
    #unknown = df[~df["CLNSIG"].isin(["benign", "likely_benign", "pathogenic", "likely_pathogenic"])]
    

    groups = {
        "benign": benign,
        "pathogenic": pathogenic,
        "not_reported": not_reported
    }
    #group_df = pd.concat[benign, pathogenic, not_reported]
    #grouped_data = pd.concat([benign, pathogenic, not_reported]).groupby("CLNSIG")["RankScore"].agg(["mean", "std", "count"])
    grouped_data = df.groupby("CLNSIG")["RankScore"].agg(["mean", "std", "count"])

    counts = {k: v.shape[0] for k, v in groups.items()}
    return groups, counts, grouped_data

def calculate_rank_score(df):
    # Already calculated during parsing, if needed can be recalculated here
    return df