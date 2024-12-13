# scripts/visualize.py

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

def plot_variant_scores(df, output_path="variant_rank_scores.png"):
    """
    Plot RankScore for benign and pathogenic variants. 
    Parameters:
        - df (pd.DataFrame): DataFrame with columns "RankScore" and "Group"
    """
    plt.figure(figsize=(10, 6))
    sns.boxplot(x=df.index.values, y="RankScore", data=df)
    plt.title("Comparison of Rank Scores between Benign and Pathogenic Variants")
    plt.ylabel("RankScore")
    plt.xlabel("Group")
    plt.savefig(output_path)

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