from cyvcf2 import VCF
import json
from collections import defaultdict


def extract_gene_variants(path_to_vcf, output_path):
    vcf = VCF(path_to_vcf)
    gene_variants = defaultdict(list)

    for variant in vcf:
        chrom = variant.CHROM
        pos = variant.POS
        ref = variant.REF
        alt = variant.ALT[0] if len(variant.ALT) else ""

        annotations = variant.INFO.get("CSQ").split(",")
        for ann in annotations:
            fields = ann.split("|")
            gene_name = fields[3]
            consequence = fields[1]
            impact = fields[2]


            gene_variants[gene_name].append({
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "consequence": consequence,
                "impact": impact
            })

    with open(output_path, "w") as f:
        json.dump(gene_variants, f, indent=2)
