# main.py
import pandas as pd
import argparse
import matplotlib.pyplot as plt
from constants import CONTROLS
from scripts.extract_data import parse_vcf
from scripts.process_data import group_variants, filter_data_and_adjust_scores, melt_data
from scripts.visualize_data import plot_scatter_with_controls, plot_violin_adjusted_scores_by_group, plot_density_adjusted_scores_by_group, visualize_metrics, plot_confusion_matrix, plot_feature_contributions
from scripts.calculate_metrics import prepare_data, calculate_metrics_by_threshold

def main(vcf_path, output_dir, controls_path=CONTROLS):
    ### PARSE ###
    print("Parsing VCF...")
    df = parse_vcf(vcf_path)

    ### GROUP VARIANTS, MARK CONTROLS ###
    print("Grouping variants...")
    df = group_variants(df.copy())

    ### SCATTER WITH CONTROLS ### 
    plot_scatter_with_controls(df, "RANK_SCORE", output_path=f"{output_dir}/scatter_plot_by_group_with_controls.png")

    ### FILTER DATA AND ADJUST SCORES ### 
    filtered_df = filter_data_and_adjust_scores(df)
   
    ### VIOLIN PLOT BY GROUP (PATHOGENIC vs. BENIGN) ###
    plot_violin_adjusted_scores_by_group(filtered_df, "ADJUSTED_SCORE", output_path=f"{output_dir}/violin_plot_by_group.png")

    ### DENSITY PLOT BY GROUP (PATHOGENIC vs. BENIGN) ###
    plot_density_adjusted_scores_by_group(filtered_df, "ADJUSTED_SCORE", output_path=f"{output_dir}/density_plot_by_group.png")
    
    ### THRESHOLD ANALYSIS
    # set thresholds and generate true_labels 
    thresholds = range(-15, 20)
    filtered_data, y_true = prepare_data(filtered_df)
    # calculate metrics for thresholds 
    threshold_metrics = calculate_metrics_by_threshold(filtered_data["ADJUSTED_SCORE"], y_true, thresholds)
    # visualize metrics 
    visualize_metrics(threshold_metrics)

    # calculate and visualize metrics for optimal threshold 
    optimal_threshold = threshold_metrics.loc[threshold_metrics["f1_score"].idmax(), "threshold"]
    y_pred = (filtered_data["ADJUSTED_SCORE"] >= optimal_threshold).astype(int)
    
    # plot the confusion matrix
    plot_confusion_matrix(y_true, y_pred, threshold=optimal_threshold)
    
    # melt the data for plotting feature contributions 
    melted_data = melt_data(filtered_data)

    # plot feature contributions
    plot_feature_contributions(melted_data)

    # reweight, recalculate, and revisualize
    weights = {"AF": 0.5, "PP": 1.0, "CON": 1.2, "VCQF": 1.2, "LIN": 1.0, "CLIN": 0.5}
    df["REWEIGHTED_RANK_SCORE"] = (
    df["AF"] * weights["AF"] +
    df["PP"] * weights["PP"] +
    df["CON"] * weights["CON"] +
    df["VCQF"] * weights["VCQF"] +
    df["LIN"] * weights["LIN"] +
    df["CLIN"] * weights["CLIN"]
)
    df[(df["IS_CONTROL"]) & (df["RANK_SCORE"] > 5) & (df["RANK_SCORE"] < 15)]

    df["CLIN"] = df["CLIN"].apply(lambda x: 0 if x < 0 else x)
    
    updated_overlap_controls = df[(df["IS_CONTROL"]) & (df["REWEIGHTED_RANK_SCORE"] > 5) & (df["REWEIGHTED_RANK_SCORE"] < 15)]
    updated_overlap_controls[["AF", "PP", "CON", "VCQF", "LIN", "CLIN"]].boxplot(figsize=(10, 6))
    plt.title("Feature Components in Overlap Region After Reweighting")
    plt.ylabel("Scores")
    plt.xticks(rotation=45)
    plt.savefig("reweighted.png")
    plt.show()
    uncertain_variants = df[df["CLNSIG"] == "Uncertain_significance"]
    uncertain_variants.head()
    low_scoring_controls[["CLIN", "RANK_SCORE"]].describe()

    metrics = calculate_metrics(y_true, y_pred, y_scores)
    print("Metrics:")
    for metric, value in metrics.items():
        print(f"{metric}: {value:.2f}")

    analyze_overlap_features(df)
    ## Investigate Feature Components for Overlapping Variants
    overlap_controls = df[(df["IS_CONTROL"]) & (df["RANK_SCORE"] > 5) & (df["RANK_SCORE"] < 15)]
    print(overlap_controls[["AF", "PP", "CON", "VCQF", "LIN", "CLIN"]].describe())
    
    df, updated_overlap_controls = main_adjust_and_evaluate(df)
    df = main_further_refinements(df)
    
    df, false_positives, false_negatives = main_optimal_threshold_analysis(df)
    #print(uncertain_variants[["VARIANT", "RANK_SCORE", "NO_CLIN_RANK_SCORE", "GROUP"]].head())
    #plot_variant_scores(df[["Group", "RankScore"]])
   # print("Grouping variants...")
   # groups, counts, grouped_data_groups, grouped_data_all = group_variants(df)
    #print("Variant counts:", counts)

    #print("Generating box plot...")
    #plot_box(groups, output_path=f"{output_dir}/boxplot.png")
    
   # print("Plotting rank_score distribution....")
    #plot_box_distribution(df, output_path=f"{output_dir}/rank_score_dist_box.png")
    # plot_hist_distribution(df, output_path=f"{output_dir}/rank_score_dist_hist.png")

    # print("Plotting grouped_data....")
    #plot_grouped_data(grouped_data_groups, output_path=f"{output_dir}/grouped_data.png")
    
    # print("Generating stacked plot...")
    # plot_stacked(groups, output_path=f"{output_dir}/stackedplot.png")
    
    # print("Plots generated successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate Genmod scores from VCF files.")
    parser.add_argument("vcf", help="Path to the VCF file.")
    parser.add_argument("--output", default="output", help="Directory to save plots.")
    parser.add_argument("controls", help="Path to the controls_match file")
    args = parser.parse_args()
    main(args.vcf, args.output, args.controls)