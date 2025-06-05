import json

def load_gene_to_pathway_map(path="data/gene_pathway_map.json"):
    with open(path, "r") as f:
        return json.load(f)

def map_genes_to_pathways(gene_scores, gene_to_pathway):
    """
    Adds pathway info to each gene in gene_scores.

    Args:
        gene_scores (dict): Output from aggregate_by_gene()
        gene_to_pathway (dict): Gene → [pathways]

    Returns:
        dict: Updated gene_scores with 'pathways' field
    """
    for gene, info in gene_scores.items():
        info["pathways"] = gene_to_pathway.get(gene, [])
    return gene_scores