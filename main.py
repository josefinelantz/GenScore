# main.py

import argparse
from scripts.parse_vcf import parse_vcf
from scripts.analyze_scores import group_variants
from scripts.visualize import plot_box, plot_stacked, plot_box_distribution, plot_hist_distribution, plot_grouped_data, plot_variant_scores

def main(vcf_path, output_dir):
    ### PARSE ###
    print("Parsing VCF...")
    df = parse_vcf(vcf_path)

    ### INSPECT DATA ###
    print(df.head())

    ### SCATTER PlOTS ### 

    ### (Variants on x-axis, RankScore on y-axis)

    ### (Variants on x-axis, AF score on y-axis)

    ### (Variants on x-axis, PP score on y-axis)

    ### Variants on x-axis, CON score on y-axis)

    ### (Variants on x-axis, VCQF score on y-axis)

    ### (Variants on x-axis, VAF score on y-axis)

    ### (Variants on x-axis, CLIN score on y-axis)


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
    
    args = parser.parse_args()
    main(args.vcf, args.output)