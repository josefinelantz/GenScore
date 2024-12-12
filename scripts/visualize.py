# scripts/visualize.py

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

# def plot_scatter(df, output_path="scatter.png"):
#     """Plots distribution based on rank_score interval and clnsig group.
#     Plot rankscore on x-axis and a random y-axis. 
#     Color code the points based on their CLNSIG group. 
#     Shows the overlaps between different CLNSIG groups and how they distribute across the RankScore range"""
    
#     filtered_df = df[(df["CLNSIG"] == "Benign") | (df["CLNSIG"] == "Likely_benign") | 
#                      (df["CLNSIG"] == "Pathogenic") | (df["CLNSIG"] == "Likely_pathogenic")]
    
#     colors = {"Benign": "blue", "Likely_benign": "lightblue", "Pathogenic": "red", "Likely_pathogeic": "pink"}
#     filtered_df["color"] = filtered_df["CLNSIG"].map(colors)
#     #benign = df.loc[df["CLNSIG"].isin([["Benign", "Likely_benign", "Benign/Likely_benign"]])] 
#     #pathogenic = df.loc[df["CLNSIG"].isin(["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"])]

#     plt.scatter(filtered_df["RankScore"], np.random.rand(len(filtered_df)), c=filtered_df["color"])
#    # plt.scatter(pathogenic["RankScore"], np.random.rand(len(pathogenic)), color="red", label="Pathogenic")

#     plt.xlabel("RankScore")
#     plt.ylabel("Random Y-axis")
#     plt.title("RankScore Distribution by CLNSIG")
#     plt.legend()
#     plt.savefig(output_path)
#     plt.close()

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