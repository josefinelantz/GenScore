# Assuming the dataset 'df' has the necessary columns: GROUP, NO_CLIN_RANK_SCORE, y_true, IS_CONTROL

# 1. Focus on CLNSIG-annotated variants (benign and pathogenic)
def filter_annotated_variants(df):
    # Retain rows where GROUP is "benign" or "pathogenic" OR where is_control is True
    filtered_df = df[(df["GROUP"].isin(["benign", "pathogenic"])) | (df["IS_CONTROL"])]
    return filtered_df

# 2. Quantify group separation
def compare_group_statistics(filtered_df):
    stats = filtered_df.groupby("GROUP")["NO_CLIN_RANK_SCORE"].describe()
    print("Statistics for (Rank Score - CLIN Score) by GROUP:\n")
    print(stats)
    print()

    # Visualize histograms
    sns.histplot(data=filtered_df, x="NO_CLIN_RANK_SCORE", hue="GROUP", kde=True, bins=30, palette="pastel")
    #plt.title("Histogram of NO_CLIN_RANK_SCORE by GROUP")
    plt.xlabel("Rank Score - CLIN Score")
    plt.ylabel("Count")
    plt.savefig("clnsig_compare_group_histogram.png")
    plt.show()

# 3. Identify and analyze overlaps
def identify_overlaps(filtered_df, threshold_low, threshold_high):
    overlap_df = filtered_df[
        (filtered_df["NO_CLIN_RANK_SCORE"] > threshold_low) &
        (filtered_df["NO_CLIN_RANK_SCORE"] < threshold_high)
    ]
    print(f"Number of overlapping variants between {threshold_low} and {threshold_high}: {len(overlap_df)}")
    return overlap_df

# 4. Fine-tune scoring thresholds
def evaluate_thresholds(filtered_df, thresholds):
    for t in thresholds:
        # Use .loc to avoid SettingWithCopyWarning
        filtered_df.loc[:, "y_pred"] = filtered_df["NO_CLIN_RANK_SCORE"].apply(lambda x: 1 if x >= t else 0)
        precision = precision_score(filtered_df["y_true"], filtered_df["y_pred"])
        recall = recall_score(filtered_df["y_true"], filtered_df["y_pred"])
        print(f"Threshold: {t}, Precision: {precision:.2f}, Recall: {recall:.2f}")

# 5. Simplify visualization with controls highlighted
def scatterplot_with_controls_clnsig(df):
    #sns.scatterplot(data=df, x="RANK_SCORE", hue="GROUP", style="is_control", palette="pastel")
    #plt.title("Scatterplot with Pathogenic Controls Highlighted")
    #plt.xlabel("Rank Score")
    #plt.ylabel("fdaf")
    #plt.legend(title="Group")
    #plt.show()
    
    plt.figure(figsize=(10, 6))
    
    benign = df[df["GROUP"] == "benign"]
    plt.scatter(
        df.index, 
        df["RANK_SCORE"], c="blue", alpha=0.7, edgecolor="black", s=30, label="Benign Variants"
    )

    pathogenic = df[df["GROUP"] == "pathogenic"]
    plt.scatter(
        df.index, 
        df["RANK_SCORE"], c="red", alpha=0.7, edgecolor="black", s=30, label="Pathogenic Variants"
    )

    # Scatter controls
    controls = df[df["IS_CONTROL"]]
    plt.scatter(
        controls.index, controls["RANK_SCORE"],
        c="gold", edgecolor="black", s=100, label="Positive Controls"
    )

    # Add threshold line (optional)
    #threshold = 10
    #plt.axhline(y=threshold, color="red", linestyle="--", label=f"Threshold = {threshold}")

    # Titles and labels
    #plt.title("Controls Relative All Variants", fontsize=14)
    #plt.xlabel("Variant Index", fontsize=12)
    plt.ylabel("Rank Score", fontsize=12)

    #plt.savefig("scatter_plot_highlighted_controls.png")
    #plt.tight_layout()
    plt.savefig("clnsig_scatter_with_controls.png")
    plt.show()