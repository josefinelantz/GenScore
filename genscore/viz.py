import matplotlib.pyplot as plt

def plot_pathway_impact_barplot(pathway_summary: dict, top_n: int = 15, output_file: str = None):
    """
    Plot a horizontal bar chart of pathways by total variant impact score.

    Args:
        pathway_summary (dict): Output from aggregate_by_pathway()
        top_n (int): Number of top pathways to display
        output_file (str): If given, save plot to this file (e.g., "results/pathways.png")
    """
    # Prepare data
    records = [
        {"pathway": k, "total_score": v["total_score"], "genes": ", ".join(v["genes"])}
        for k, v in pathway_summary.items()
    ]
    records = sorted(records, key=lambda x: -x["total_score"])[:top_n]

    if not records:
        print("No data to plot.")
        return

    pathways = [r["pathway"] for r in records]
    scores = [r["total_score"] for r in records]

    # Plot
    plt.figure(figsize=(10, 6))
    bars = plt.barh(pathways, scores)
    plt.xlabel("Total Variant Score")
    plt.title("Top Impacted Pathways")
    plt.gca().invert_yaxis()
    plt.tight_layout()

    # Save or show
    if output_file:
        plt.savefig(output_file)
        print(f"Saved plot to {output_file}")
    else:
        plt.show()