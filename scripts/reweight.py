# Implementation to handle the low-scoring controls

# 1. Adjust Handling of `CLIN` for Controls
def adjust_clin_for_controls(df):
    """
    Neutralize negative `CLIN` contributions for controls with `Uncertain_significance`.
    """
    df.loc[(df["IS_CONTROL"]) & (df["CLNSIG"] == "Uncertain_significance"), "CLIN"] = 0
    print("Adjusted CLIN values for controls with Uncertain_significance.")
    return df

# 2. Reweight Scoring Formula
def reweight_scoring_formula(df, weights):
    """
    Reweight the scoring formula and recalculate the RANK_SCORE.
    """
    df["REWEIGHTED_RANK_SCORE"] = (
        df["AF"] * weights["AF"] +
        df["PP"] * weights["PP"] +
        df["CON"] * weights["CON"] +
        df["VCQF"] * weights["VCQF"] +
        df["LIN"] * weights["LIN"] +
        df["CLIN"] * weights["CLIN"]
    )
    print("Recalculated RANK_SCORE using new weights.")
    return df

# 3. Update Grouping for Controls
def update_control_grouping(df):
    """
    Reassign the group for controls with `Uncertain_significance` to `control_uncertain`.
    """
    df.loc[(df["IS_CONTROL"]) & (df["GROUP"] == "other"), "GROUP"] = "control_uncertain"
    print("Updated GROUP for controls with Uncertain_significance.")
    return df

# 4. Verify Changes for Specific Controls
def verify_control_updates(df, variant_id):
    """
    Print updated details for a specific control variant.
    """
    control_data = df.loc[df["VARIANT"] == variant_id, ["RANK_SCORE", "NO_CLIN_RANK_SCORE", "CLIN", "GROUP"]]
    print(f"Updated details for variant {variant_id}:")
    print(control_data)

# Example Usage
def main_handle_low_scoring_controls(df):
    # Step 1: Adjust CLIN values for controls
    df = adjust_clin_for_controls(df)

    # Step 2: Reweight scoring formula
    weights = {"AF": 1.0, "PP": 1.0, "CON": 1.5, "VCQF": 1.0, "LIN": 1.0, "CLIN": 0.5}
    df = reweight_scoring_formula(df, weights)

    # Step 3: Update control grouping
    df = update_control_grouping(df)

    # Step 4: Verify changes for the low-scoring control
    verify_control_updates(df, variant_id="21_36206711_C_T")

    return df

# Uncomment and run when `df` is defined
# df = main_handle_low_scoring_controls(df)

# Recommended adjustments based on feature overlaps 
# 1. Adjust CLIN to Avoid Negative Penalties
def neutralize_negative_clin(df):
    """
    Set negative CLIN values to 0 to avoid penalizing controls.
    """
    df["CLIN"] = df["CLIN"].apply(lambda x: 0 if x < 0 else x)
    print("Negative CLIN values have been neutralized.")
    return df

# 2. Reweight Features in the Scoring Formula
def reweight_features(df, weights):
    """
    Apply new weights to feature components and recalculate RANK_SCORE.
    """
    df["REWEIGHTED_RANK_SCORE"] = (
        df["AF"] * weights["AF"] +
        df["PP"] * weights["PP"] +
        df["CON"] * weights["CON"] +
        df["VCQF"] * weights["VCQF"] +
        df["LIN"] * weights["LIN"] +
        df["CLIN"] * weights["CLIN"]
    )
    print("Recalculated RANK_SCORE with updated weights.")
    return df

# 3. Reevaluate Overlap Region
def reevaluate_overlap_region(df, lower_threshold=5, upper_threshold=15):
    """
    Analyze and print statistics for the overlap region after reweighting.
    """
    overlap_controls = df[(df["IS_CONTROL"]) & (df["REWEIGHTED_RANK_SCORE"] > lower_threshold) & (df["REWEIGHTED_RANK_SCORE"] < upper_threshold)]
    print(f"Number of controls in overlap region ({lower_threshold} < REWEIGHTED_RANK_SCORE < {upper_threshold}): {len(overlap_controls)}")
    print("\nDescriptive statistics for feature components in the overlap region:")
    print(overlap_controls[["AF", "PP", "CON", "VCQF", "LIN", "CLIN", "REWEIGHTED_RANK_SCORE"]].describe())
    return overlap_controls

# 4. Visualize Updated Contributions
def visualize_updated_contributions(overlap_controls):
    """
    Plot boxplot of feature contributions for the overlap region after reweighting.
    """
    overlap_controls[["AF", "PP", "CON", "VCQF", "LIN", "CLIN"]].boxplot(figsize=(10, 6))
    plt.title("Feature Components in Overlap Region After Reweighting")
    plt.ylabel("Scores")
    plt.xticks(rotation=45)
    plt.savefig("updated_contributions.png")
    plt.show()

# Example Usage
def main_adjust_and_evaluate(df):
    # Step 1: Neutralize negative CLIN values
    df = neutralize_negative_clin(df)

    # Step 2: Reweight features
    weights = {"AF": 1.5, "PP": 1.2, "CON": 1.0, "VCQF": 1.0, "LIN": 1.0, "CLIN": 0.5}
    df = reweight_features(df, weights)

    # Step 3: Reevaluate the overlap region
    overlap_controls = reevaluate_overlap_region(df)

    # Step 4: Visualize updated contributions
    visualize_updated_contributions(overlap_controls)

    return df, overlap_controls

# Uncomment and run when `df` is defined
# df, updated_overlap_controls = main_adjust_and_evaluate(df)

# 1. Fine-Tune Feature Weights
def fine_tune_weights(df):
    """
    Apply fine-tuned weights to the scoring formula for further improvement.
    """
    fine_tuned_weights = {"AF": 2.0, "PP": 1.5, "CON": 0.8, "VCQF": 1.2, "LIN": 1.0, "CLIN": 0.6}
    df["FINE_TUNED_RANK_SCORE"] = (
        df["AF"] * fine_tuned_weights["AF"] +
        df["PP"] * fine_tuned_weights["PP"] +
        df["CON"] * fine_tuned_weights["CON"] +
        df["VCQF"] * fine_tuned_weights["VCQF"] +
        df["LIN"] * fine_tuned_weights["LIN"] +
        df["CLIN"] * fine_tuned_weights["CLIN"]
    )
    print("Applied fine-tuned weights and recalculated FINE_TUNED_RANK_SCORE.")
    return df

# 2. Evaluate Performance Metrics
from sklearn.metrics import precision_score, recall_score, f1_score, roc_auc_score

def calculate_performance_metrics(df, score_column="FINE_TUNED_RANK_SCORE", threshold=10):
    """
    Calculate precision, recall, F1-score, and AUC for the scoring model.
    """
    y_true = df["y_true"]
    y_pred = df[score_column].apply(lambda x: 1 if x >= threshold else 0)

    precision = precision_score(y_true, y_pred, zero_division=0)
    recall = recall_score(y_true, y_pred, zero_division=0)
    f1 = f1_score(y_true, y_pred, zero_division=0)
    auc = roc_auc_score(y_true, df[score_column])

    print(f"Performance Metrics (Threshold = {threshold}):")
    print(f"Precision: {precision:.2f}")
    print(f"Recall: {recall:.2f}")
    print(f"F1-Score: {f1:.2f}")
    print(f"AUC: {auc:.2f}")

    return precision, recall, f1, auc

# 3. Remove LIN Feature if Needed
def remove_lin_feature(df):
    """
    Remove LIN column if it consistently contributes zero to the scoring model.
    """
    if df["LIN"].sum() == 0:
        df.drop(columns=["LIN"], inplace=True)
        print("LIN feature removed from the dataset.")
    else:
        print("LIN feature retained in the dataset.")
    return df

# 4. Visualize Fine-Tuned Scores Across Entire Dataset
def visualize_fine_tuned_scores(df):
    """
    Plot the distribution of FINE_TUNED_RANK_SCORE across the entire dataset.
    """
    sns.histplot(df["FINE_TUNED_RANK_SCORE"], bins=30, kde=True, color="blue")
    plt.title("Distribution of FINE_TUNED_RANK_SCORE Across Dataset")
    plt.xlabel("FINE_TUNED_RANK_SCORE")
    plt.ylabel("Count")
    plt.savefig("reweighted.png")
    plt.show()

# 5. Apply All Steps
def main_further_refinements(df):
    # Step 1: Fine-tune weights
    df = fine_tune_weights(df)

    # Step 2: Evaluate performance metrics
    calculate_performance_metrics(df)

    # Step 3: Remove LIN feature if it contributes nothing
    df = remove_lin_feature(df)

    # Step 4: Visualize fine-tuned scores
    visualize_fine_tuned_scores(df)

    return df

# Uncomment and run when `df` and `y_true` exist
# df = main_further_refinements(df)