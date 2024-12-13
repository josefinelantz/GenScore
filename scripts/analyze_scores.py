# scripts/analyze_scores.py

import pandas as pd

# 'not_reported' 'Benign' 'Uncertain_significance' 'Likely_benign'
#  'Benign/Likely_benign' 'Conflicting_classifications_of_pathogenicity'
#  'Pathogenic' 'drug_response' 'Pathogenic/Likely_pathogenic'
#  'Conflicting_classifications_of_pathogenicity|drug_response|other'
#  'Benign|drug_response' 'Pathogenic|other' 'Likely_pathogenic'
#  'Likely_benign|other' 'Likely_pathogenic|association'
#  'Pathogenic/Likely_pathogenic|risk_factor'
#  'Conflicting_classifications_of_pathogenicity|association' 'not_provided'
#  'Likely_benign|drug_response|other'
BENIGN = [
    "Benign", 
    "Likely_benign", 
    "Benign/Likely_benign"
]
PATHOGENIC = [
    "Pathogenic", 
    "Likely_pathogenic", 
    "Pathogenic/Likely_pathogenic",
    "Likely_pathogenic|association"
]
UNCERTAIN = [
    "Uncertain_significance", 
    "Conflicting_classifications_of_pathogenicity",
]
def group_variants(df):
    # benign = df[df["CLNSIG"] == "Benign"]
    # pathogenic = df[df["CLNSIG"] == "Pathogenic"]
    # uncertain = df[df["CLNSIG"] == "Uncertain_significance"]
    benign = df[df["CLNSIG"].isin(BENIGN)]
    pathogenic = df[df["CLNSIG"].isin(PATHOGENIC)]
    uncertain = df[df["CLNSIG"].isin(UNCERTAIN)]

    groups = {
        "benign": benign,
        "pathogenic": pathogenic,
        "uncertain": uncertain
    }
    grouped_data_groups = pd.concat([benign, pathogenic, uncertain]).groupby("CLNSIG")["RankScore"].agg(["mean", "std", "count"])
    grouped_data_all = df.groupby("CLNSIG")["RankScore"].agg(["mean", "std", "count"])

    # Count number of variants in each group 
    counts = {k: v.shape[0] for k, v in groups.items()}
    return groups, counts, grouped_data_groups, grouped_data_all