import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def violin_plot_by_group(df, column):
    sns.violinplot(data=df, hue="GROUP", y=column)
    #plt.title(f"Violin Plot by Group for {column}")
    plt.ylabel("Rank Score - Cinical Significance Score")
    plt.legend(loc="lower right", title="Groups")  # Adjust legend position
    plt.savefig("violin_plot_subtracted_clin_score.png")
    plt.show()

def add_jitter(values, jitter_amount=1.5):
    return values + np.random.uniform(-jitter_amount, jitter_amount, size=len(values))

def scatter_plot_all_by_group(df, x_col, y_col, hue="GROUP"):
    plt.figure(figsize=(12, 6))
    df["VariantIndex"] = range(len(df))
    # Adding jitter to separate overlapping points slightly
    jittered_index = add_jitter(df['VariantIndex'])
    # Create scatterplot with jittered points
    #plt.figure(figsize=(12, 6))
    plt.scatter(jittered_index, df['RANK_SCORE'], c=df['Color'], alpha=0.8, edgecolor='white', s=90)

    # Add labels and title
    plt.xlabel("Variant Index", fontsize=12)
    plt.ylabel("Rank Score", fontsize=12)
    #plt.title("Variants by Rank Score and CLNSIG Group", fontsize=14)
    plt.xticks([])
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    # Add legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label='Benign', markerfacecolor='blue', markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', label='Pathogenic', markerfacecolor='red', markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', label='Other', markerfacecolor='gray', markersize=10)
    ]
    plt.legend(handles=legend_elements, loc='lower left')

    # Show plot
    plt.tight_layout()
    plt.savefig("jitter_scatter_variants_by_rankscore_group.png")
    plt.show()

def scatter_plot_with_controls(df, x_column, y_column):
    plt.figure(figsize=(12, 6))
    #df["VariantIndex"] = range(len(df))
    
    plt.scatter(
        df.index, 
        df["RANK_SCORE"], c="gray", alpha=0.6, edgecolor="black", s=30, label="All Variants"
    )

    # Scatter controls
    controls = df[df["is_control"]]
    plt.scatter(
        controls.index, controls["RANK_SCORE"],
        c="gold", edgecolor="black", s=100, label="Positive Controls"
    )

    # Add threshold line (optional)
    #threshold = 10
    #plt.axhline(y=threshold, color="red", linestyle="--", label=f"Threshold = {threshold}")

    # Titles and labels
    #plt.title("Controls Relative All Variants", fontsize=14)
    plt.xlabel("Variant Index", fontsize=12)
    plt.ylabel("Rank Score", fontsize=12)

    plt.savefig("scatter_plot_highlighted_controls.png")
    #plt.tight_layout()
    plt.show()

    
    #sns.scatterplot(data=df, x=x_column, y=y_column, hue="GROUP", style="IS_CONTROL", palette="pastel")
   # plt.title(f"Scatter Plot of {x_column} vs {y_column}")
    #plt.show()

def stacked_barplot_for_controls(df):
    """
    Generates a stacked bar plot showing the contribution of each score category 
    to the total RANK_SCORE for the known controls.

    Parameters:
        df (pd.DataFrame): DataFrame containing the processed variant data. Must include 
                           IS_CONTROL column (True for controls) and individual score columns.

    Returns:
        None
    """
    # Filter for controls
    control_data = df[df["is_control"]]

    # Select only the score columns for the plot
    score_columns = ["AF", "PP", "CON", "VCQF", "LIN", "CLIN"]
    control_scores = control_data[score_columns]

    # Plot the stacked barplot
    control_scores.plot(
        kind="bar", 
        stacked=True, 
        figsize=(12, 8), 
        color=["blue", "green", "orange", "red", "purple", "pink"]
    )
    
    #stacked_data = control_data[["AF", "PP", "CON", "VCQF", "LIN", "CLIN"]].mean()
    #plt.title("Category Contribution for Known Controls")
    plt.savefig("stacked_barplot_controls_category_scoring.png")
    plt.xlabel("Positive Controls", fontsize=12)
    plt.ylabel("Category Score Contribution", fontsize=12)
    plt.xticks([])  # Rotate x-axis labels for clarity
    #plt.xticks(rotation=45, ha="right")  # Rotate x-axis labels for clarity
    plt.legend(title="Scoring Categories", loc="upper right", bbox_to_anchor=(0.93, 0.95))
    plt.tight_layout()  # Adjust layout for better visualization
    plt.show()