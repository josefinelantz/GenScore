from genscore.io import parse_vcf, extract_genes_from_vcf
from genscore.scorer import aggregate_by_gene
from genscore.pathway import load_gene_to_pathway_map, map_genes_to_pathways
from genscore.viz import plot_pathway_impact_barplot
from genscore.utils import load_pathway_summary_csv_to_dict

pathway_dict = load_pathway_summary_csv_to_dict("data/Pathway_Impact_Summary.csv")

variants = parse_vcf("data/mock.vcf.gz")

genes = extract_genes_from_vcf("data/mock.vcf.gz", "data/mock_genes.json")

gene_scores = aggregate_by_gene(variants)

gene_to_pathway = load_gene_to_pathway_map("data/gene_pathway_map.json")

updated_scores = map_genes_to_pathways(gene_scores, gene_to_pathway)

#for gene in sorted(genes): 
    #print(gene)

plot_pathway_impact_barplot(pathway_dict, top_n=10, output_file="results/pathway_impact.png")