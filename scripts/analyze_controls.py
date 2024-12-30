# 1. Investigate low-scoring controls
def investigate_low_scoring_controls(df, score_threshold=10):
    """
    Investigate controls with low RANK_SCORE.
    """
    low_scoring_controls = df[(df["is_control"]) & (df["RANK_SCORE"] < score_threshold)]
    print(f"Number of low-scoring controls (RANK_SCORE < {score_threshold}): {len(low_scoring_controls)}")
    print("\nDescriptive statistics for low-scoring controls:")
    print(low_scoring_controls[["AF", "PP", "CON", "VCQF", "LIN", "CLIN"]].describe())
    
    # Boxplot for feature components
    low_scoring_controls[["AF", "PP", "CON", "VCQF", "LIN", "CLIN"]].boxplot(figsize=(10, 6))
    plt.title(f"Feature Components for Low-Scoring Controls (RANK_SCORE < {score_threshold})")
    plt.ylabel("Scores")
    plt.xticks(rotation=45)
    plt.show()

    return low_scoring_controls

# 2. Decompose scoring system for overlap analysis
def analyze_overlap_features(df, lower_threshold=5, upper_threshold=15):
    """
    Analyze feature components for variants in the overlap region of RANK_SCORE.
    """
    overlap_df = df[(df["RANK_SCORE"] > lower_threshold) & (df["RANK_SCORE"] < upper_threshold)]
    print(f"Number of variants in the overlap region ({lower_threshold} < RANK_SCORE < {upper_threshold}): {len(overlap_df)}")
    print("\nDescriptive statistics for feature components in the overlap region:")
    print(overlap_df[["AF", "PP", "CON", "VCQF", "LIN", "CLIN"]].describe())
    
    # Boxplot for feature components in overlap
    overlap_df[["AF", "PP", "CON", "VCQF", "LIN", "CLIN"]].boxplot(figsize=(10, 6))
    plt.title(f"Feature Components in Overlap Region ({lower_threshold} < RANK_SCORE < {upper_threshold})")
    plt.ylabel("Scores")
    plt.xticks(rotation=45)
    plt.show()

    return overlap_df

# 3. Visualize individual feature distributions for controls vs non-controls
def visualize_feature_distributions(df, feature):
    """
    Visualize distributions of a specific feature for controls and non-controls.
    """
    sns.boxplot(data=df, x="is_control", y=feature, palette="pastel")
    plt.title(f"{feature} Distribution for Controls and Non-Controls")
    plt.xlabel("Control Status")
    plt.ylabel(feature)
    plt.show()

# Example Usage
def main_investigate_controls(df):
    # Step 1: Investigate low-scoring controls
    low_scoring_controls = investigate_low_scoring_controls(df, score_threshold=10)

    # Step 2: Analyze overlap features
    overlap_df = analyze_overlap_features(df, lower_threshold=5, upper_threshold=15)

    # Step 3: Visualize individual feature distributions
    for feature in ["AF", "PP", "CON", "VCQF", "LIN", "CLIN"]:
        visualize_feature_distributions(df, feature)

    return low_scoring_controls, overlap_df

# Uncomment and run when the dataset `df` is defined
# low_scoring_controls, overlap_df = main_investigate_controls(df)