# main.py
import pandas as pd
import argparse
import matplotlib.pyplot as plt
from scripts.extract_data import parse_vcf
from scripts.process_data import analyze_scores, mark_controls, subtract_clin_score, calculate_metrics
from scripts.visualize_scores import violin_plot_by_group, violin_benign_pathogenic, scatter_plot_with_controls
from scripts.analyze_annotated import filter_annotated_variants, compare_group_statistics, identify_overlaps, evaluate_thresholds, main_investigate_controls, investigate_low_scoring_controls, analyze_overlap_features, main_adjust_and_evaluate, main_further_refinements
from scripts.reweight.py import main_handle_low_scoring_controls
from scripts.calculate_metrics import main_optimal_threshold_analysis

def main(vcf_path, output_dir, controls_path):
    ### PARSE ###
    print("Parsing VCF...")
    df = parse_vcf(vcf_path)

    ### GROUP VARIANTS AND MARK CONTROLS###
    print("Grouping variants...")
    df = analyze_scores(df)

    ### MARK CONTROLS ### 
    controls = pd.read_csv(controls_path, sep="\t")

    df = mark_controls(df, controls)

    ### SUBTRACT CLIN SCORE ### 
    df = subtract_clin_score(df)

    ### VIOLIN BY GROUP ###

    violin_plot_by_group(df, "NO_CLIN_RANK_SCORE", output_path=f"{output_dir}/violin_by_group.png")

    ### Convert GROUP to binary labels 
    df["y_true"] = df["GROUP"].apply(lambda x: 1 if x == "pathogenic" else 0)  # Convert GROUP to binary labels
    threshold = 10.0  # Example threshold for classification

    # Predicted labels based on the threshold
    df["y_pred"] = df["NO_CLIN_RANK_SCORE"].apply(lambda x: 1 if x >= threshold else 0)

    # Scores to use for AUC calculation
    y_true = df["y_true"].values  # Ground truth labels
    y_pred = df["y_pred"].values  # Predicted binary labels
    y_scores = df["NO_CLIN_RANK_SCORE"].values  # Continuous scores
    
    metrics = calculate_metrics(y_true, y_pred, y_scores)
    
    ### Filter annotated variants
    filtered_df = filter_annotated_variants(df).copy() # ensure independent DataFrame

    ### Violin benign and pathogenic
    violin_benign_pathogenic(filtered_df)

    ### Compare Group statistics
    compare_group_statistics(filtered_df)

    ### Identify overlaps
    threshold_low, threshold_high = 5, 15 # example thresholds
    overlap_df = identify_overlaps(filtered_df, threshold_low, threshold_high)

    ### Evaluate Thresholds
    thresholds = [2, 4, 6, 8, 10]
    evaluate_thresholds(filtered_df, threshold)

    scatter_plot_with_controls(filtered_df)

    # investigate low-scoring controls 
    low_scoring_controls, overlap_df = main_investigate_controls(df)

    low_scoring_controls = df[(df["IS_CONTROL"]) & (df["RANK_SCORE"] < 10)]
    low_scoring_controls.describe()

    # 1. Investigate low-scoring controls
    df_low_controls = investigate_low_scoring_controls(df)
    
    df.loc[(df["IS_CONTROL"]) & (df["CLNSIG"] == "Uncertain_significance"), "CLIN"] = 0
    
    df.loc[(df["IS_CONTROL"]) & (df["GROUP"] == "other"), "GROUP"] = "control_uncertain"
    print(df.loc[df["VARIANT"] == "21_36206711_C_T", ["RANK_SCORE", "NO_CLIN_RANK_SCORE", "CLIN"]])

    df = main_handle_low_scoring_controls(df)

    # reweight
    df = main_handle_low_scoring_controls(df)

    overlap_controls = df[(df["IS_CONTROL"]) & (df["RANK_SCORE"] > 5) & (df["RANK_SCORE"] < 15)]
    print(overlap_controls[["AF", "PP", "CON", "VCQF", "LIN", "CLIN"]].describe())

    updated_overlap_controls = df[(df["IS_CONTROL"]) & (df["REWEIGHTED_RANK_SCORE"] > 5) & (df["REWEIGHTED_RANK_SCORE"] < 15)]
    print(updated_overlap_controls[["AF", "PP", "CON", "VCQF", "LIN", "CLIN", "REWEIGHTED_RANK_SCORE"]].describe())
    updated_overlap_controls[["AF", "PP", "CON", "VCQF", "LIN", "CLIN"]].boxplot(figsize=(10, 6))
    plt.title("Feature Components in Overlap Region After Reweighting")
    plt.ylabel("Scores")
    plt.xticks(rotation=45)
    plt.show()

    weights = {"AF": 1.0, "PP": 1.0, "CON": 1.2, "VCQF": 1.2, "LIN": 1.0, "CLIN": 0.5}
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