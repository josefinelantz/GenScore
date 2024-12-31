import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def violin_plot_by_group(df, column):
    sns.violinplot(data=df, hue="GROUP", y=column)
    plt.title(f"Violin Plot by Group for {column}")
    plt.ylabel("Rank Score - Cinical Significance Score")
    plt.legend(loc="lower right", title="Groups")  # Adjust legend position
    plt.savefig("violin_plot_subtracted_clin_score.png")
    plt.show()

def violin_benign_pathogenic(filtered_df):
    sns.violinplot(data=filtered_df, hue="y_true", y="NO_CLIN_RANK_SCORE", palette="pastel", legend=False)
    plt.title("Score Distributions for Benign and Pathogenic Variants")
    plt.xlabel("True Label (0=Benign, 1=Pathogenic)")
    plt.ylabel("Rank Score - CLIN Score")
    plt.savefig("violin_ben_path.png")
    plt.show()

def add_jitter(values, jitter_amount=1.5):
    return values + np.random.uniform(-jitter_amount, jitter_amount, size=len(values))

def scatter_plot_all_by_group(df, x_col, y_col, hue="GROUP"):
    plt.figure(figsize=(12, 6))
    df["VariantIndex"] = range(len(df))
    # Adding jitter to separate overlapping points slightly
    jittered_index = add_jitter(df['VariantIndex'])
    # Create scatterplot with jittered points
    plt.scatter(jittered_index, df['RANK_SCORE'], c=df['Color'], alpha=0.8, edgecolor='white', s=90)

    # Add labels and title
    plt.xlabel("Variant Index", fontsize=12)
    plt.ylabel("Rank Score", fontsize=12)
    plt.title("Variants by Rank Score and CLNSIG Group", fontsize=14)
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

def scatter_plot_with_controls(df):
    plt.figure(figsize=(12, 6))

    sns.scatterplot(data=filtered_df, x=filtered_df.index, y="NO_CLIN_RANK_SCORE", hue="GROUP")
    # Scatter controls
    controls = df[df["IS_CONTROL"]]
    plt.scatter(
        controls.index, controls["RANK_SCORE"],
        c="gold", edgecolor="black", s=50, label="Positive Controls"
    )
    plt.xticks([])
    plt.title("Scatterplot with Pathogenic Controls Highlighted")
    plt.legend(loc="lower right")
    plt.savefig("scatter_controls.png")
    plt.show()

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
    plt.title("Category Contribution for Known Controls")
    plt.savefig("stacked_barplot_controls_category_scoring.png")
    plt.xlabel("Positive Controls", fontsize=12)
    plt.ylabel("Category Score Contribution", fontsize=12)
    plt.xticks([])  # Rotate x-axis labels for clarity
    #plt.xticks(rotation=45, ha="right")  # Rotate x-axis labels for clarity
    plt.legend(title="Scoring Categories", loc="upper right", bbox_to_anchor=(0.93, 0.95))
    plt.tight_layout()  # Adjust layout for better visualization
    plt.show()