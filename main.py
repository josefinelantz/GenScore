# main.py
import pandas as pd
import argparse
import matplotlib.pyplot as plt
from constants import CONTROLS
from scripts.extract_data import parse_vcf
from scripts.process_data import group_variants, filter_data_and_adjust_scores, melt_data
from scripts.visualize_data import plot_scatter_with_controls, plot_violin_adjusted_scores_by_group, plot_density_adjusted_scores_by_group, visualize_metrics, plot_confusion_matrix, plot_feature_contributions
from scripts.calculate_metrics import prepare_data, calculate_metrics_by_threshold, recalculate_rank_scores

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
    # assuming df has cols ADJUSTED_SCORE, GROUP
    filtered_data, y_true = prepare_data(filtered_df)
    # calculate metrics for thresholds 
    threshold_metrics = calculate_metrics_by_threshold(filtered_data["ADJUSTED_SCORE"], y_true, thresholds)
    # visualize metrics 
    visualize_metrics(threshold_metrics)

    # calculate and visualize metrics for optimal threshold 
    optimal_threshold = threshold_metrics.loc[threshold_metrics["f1_score"].idmax(), "threshold"]
    y_pred = (filtered_data["ADJUSTED_SCORE"] >= optimal_threshold).astype(float)
    
    ### VISUALIZE METRICS FOR OPTIMAL THRESHOLD
    # plot the confusion matrix
    plot_confusion_matrix(filtered_data, optimal_threshold)
    
    # melt the data for plotting feature contributions 
    melted_data = melt_data(filtered_data)

    # plot feature contributions
    plot_feature_contributions(melted_data, optimal_threshold)

    ### REWEIGHT, RECALCULATE AND REVISUALIZE 
    # define new weights 
    weights = {
        "AF": 0.5, 
        "PP": 1.0, 
        "CON": 1.2, 
        "VCQF": 1.2, 
        "LIN": 1.0,
        "CLIN": 0.5
    }

    # recalculate scores 
    df = recalculate_rank_scores(filtered_data.copy())

    # repeat threshold analysis
    threshold_metrics_updated = calculate_metrics_by_threshold(df["UPDATED_RANK_SCORE"], y_true, thresholds)
    
    # visualize updated metrics 
    visualize_metrics(threshold_metrics_updated)

    # plot Confusion Matrix for Updated Scores
    optimal_threshold_updated = threshold_metrics_updated.loc[threshold_metrics_updated["f1_score"].idxmax(), "threshold"]
    y_pred_updated = (df["UPDATED_RANK_SCORE"] >= optimal_threshold_updated).astype(int)
    plot_confusion_matrix(df, optimal_threshold_updated)
     # melt data
    melted_data = melt_data(df)
    plot_feature_contributions(melted_data, optimal_threshold_updated)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate Genmod scores from VCF files.")
    parser.add_argument("vcf", help="Path to the VCF file.")
    parser.add_argument("--output", default="output", help="Directory to save plots.")
    parser.add_argument("controls", help="Path to the controls_match file")
    args = parser.parse_args()
    main(args.vcf, args.output, args.controls)