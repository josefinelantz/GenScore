import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def violin_plot_by_group_with_counts(df, column="ADJUSTED_SCORE"):
    # Group our dataset with our 'Group' variable
    grouped = df.groupby('GROUP')[column]
    medians = grouped.median().values
    nobs = df['GROUP'].value_counts().values
    nobs = [str(x) for x in nobs.tolist()]
    nobs = ["n: " + i for i in nobs]

    # Init a figure and axes
    fig, ax = plt.subplots(figsize=(10, 8))

    # Create the plot with different colors for each group using hue
    ax = sns.violinplot(x="GROUP", y=column, hue="GROUP", data=df)
    
    # Add text to the figure
    pos = range(len(nobs))
    for tick, label in zip(pos, ax.get_xticklabels()):
        ax.text(pos[tick], medians[tick] + 2, nobs[tick],
        horizontalalignment='left',
        size='medium',
        color='black')

    # Add a title and axis label
    ax.set_title('Distribution by Group for Adjusted Score')
    plt.xlabel("")
    plt.xticks([])

    # Add a legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label='Other', markerfacecolor='#5975A4', markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', label='Benign', markerfacecolor='#CC8963', markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', label='Uncertain', markerfacecolor='#609E6E', markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', label='Pathogenic', markerfacecolor='#B65D60', markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', label='Drug_Response', markerfacecolor='#857AAB', markersize=10)
    ]
    plt.legend(handles=legend_elements, loc='best', bbox_to_anchor=(0.5, 0., 0.5, 0.5))
    plt.savefig("distribution_by_group_for_adjusted_score.png")
    plt.show()

def violin_plot_by_group(df, column):
    sns.violinplot(data=df, hue="GROUP", y=column)
    plt.title(f"Violin Plot by Group for {column}")
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

    def plot_box_distribution(df, output_path="rank_score_dist_box.png"):
    """Plots the distribution of RankScores within each CLNSIG group"""
    plt.figure(figsize=(10, 6))
    sns.boxplot(x="CLNSIG", y="RankScore", data=df)
    plt.ylabel("RankScore")
    plt.title("RankScore Distribution by CLNSIG")
    plt.savefig(output_path)
    plt.close()
# def plot_box_distribution(df, output_path="rank_score_dist_box.png"):
#     """Plots the distribution of RankScores within each CLNSIG group"""
#     df.boxplot(x="CLNSIG", y="RankScore")
#     plt.ylabel("RankScore")
#     plt.title("RankScore Distribution by CLNSIG")
#     plt.savefig(output_path)
#     plt.close()

def plot_hist_distribution(df, output_path="rank_score_dist_hist.png", bins=20):
    """Plots the distribution of RankScores within each CLNSIG group"""
    df.hist("CLNSIG", "RankScore")
    plt.xlabel("RankScore")
    plt.ylabel("Frequency")
    plt.title("RankScore Distribution by CLLNSIG")
    plt.savefig(output_path)

def plot_grouped_data(grouped_data, output_path="grouped_data.png"):
    """Plot mean, std and counts for groups
    Parameters: pd.DataFrame with columns "mean", "std", "count"
    Plots a bar plot where each bar represents a CLNSIG group. 
    The height of the bar shows the mean RankScore, and the error bars indicate the standard deviation. 
    """
    grouped_data.plot(kind="bar", y="mean", yerr="std", rot=0)
    plt.xlabel("CLNSIG")
    plt.ylabel("Mean RankScore")
    plt.title("Mean RankScore By CLNSIG")
    plt.savefig(output_path)

def plot_box(groups, output_path="boxplot.png"):
    # Combine all groups into one DataFrame
    combined = pd.concat([
        group.assign(group=group_name) 
        for group_name, group in groups.items()
    ])
    
    plt.figure(figsize=(10, 6))
    sns.boxplot(x="group", y="RankScore", data=combined)
    plt.title("RankScore Distribution by CLIN_SIG Group")
    plt.xlabel("CLIN_SIG Group")
    plt.ylabel("RankScore")
    plt.savefig(output_path)
    plt.close()

def plot_stacked(groups, output_path="stackedplot.png"):
    # Calculate the sum of RankResult categories for each group
    # Assuming RankResult has 6 categories separated by "|"
    stack_data = {}
    for group_name, group in groups.items():
        rank_results = group["RankResult"].str.split("|").apply(lambda x: list(map(int, x)))
        summed = pd.DataFrame(rank_results.tolist()).sum()
        stack_data[group_name] = summed
    
    stack_df = pd.DataFrame(stack_data)
    stack_df.plot(kind="bar", stacked=True, figsize=(10, 6))
    plt.title("RankResult Category Contributions by CLIN_SIG Group")
    plt.xlabel("CLIN_SIG Group")
    plt.ylabel("Total Score per Category")
    plt.legend(title="RankResult Categories")
    plt.savefig(output_path)
    plt.close()