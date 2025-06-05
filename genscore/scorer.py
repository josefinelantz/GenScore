from collections import defaultdict 

def aggregate_by_gene(variant_list):
    """
    Aggregate variant info per gene from parsed VCF variant list.

    Args:
        variant_list (list[dict]): Output from parse_vcf()

    Returns:
        dict[str, dict]: Aggregated gene scores and variant info
    """

    gene_data = defaultdict(lambda: {
        "total_score": 0.0, 
        "max_score": float("-inf"), 
        "n_variants": 0, 
        "impact_levels": set(), 
        "consequences": set() 
    })

    for var in variant_list:
        gene = var["gene"]
        score = var["score"]
        if score is None: 
            continue 

        try: 
            score = float(score)
        except ValueError: 
            continue 

        gene_info = gene_data[gene]
        gene_info["total_score"] += score 
        gene_info["max_score"] = max(gene_info["max_score"], score)
        gene_info["n_variants"] += 1
        gene_info["impact_levels"].add(var.get("impact", ""))
        gene_info["consequences"].add(var.get("consequence", ""))

    # Convert sets to sorted lists for JSON/export
    for gene, info in gene_data.items():
        info["impact_levels"] = sorted(info["impact_levels"])
        info["consequences"] = sorted(info["consequences"])

    return dict(gene_data)