from cyvcf2 import VCF
from collections import defaultdict
import json
import re
import gzip

def parse_vcf(file_path, max_variants=None):
    """
    Parses a VEP-annotated, Genmod-scored VCF file and returns structured variant data.

    Args:
        vcf_path (str): Path to input VCF (e.g. from Balsamic + VEP + Genmod)
        score_field (str): INFO field for pathogenicity score (e.g. RankScore)
        csq_field (str): INFO field containing VEP annotations (default: CSQ)
        max_variants (int): Optional limit for parsed variants (for testing)

    Returns:
        List[dict]: Parsed variants, each with gene, impact, and score info
    """
    # Load VCF into cyvcf2.VCF generator object
    vcf = VCF(file_path)

    # Get the headers for VEP annotations
    try: 
        desc = vcf.get_header_type("CSQ")["Description"]
        csq_format = desc.split("Format:")[-1].strip().strip('"')
        csq_fields = csq_format.split("|")
    except Exception as e: 
        raise RuntimeError(f"Failed to extract CSQ field format from VCF header: {e}")

    parsed_variants = []

    for i, record in enumerate(vcf):
        if max_variants and i >= max_variants:
            break 

        csq_entries = record.INFO.get("CSQ")
        if not csq_entries:
            continue 

        for entry in csq_entries.split(","):
            fields = entry.split("|")
            csq_dict = dict(zip(csq_fields, fields))
            gene = csq_dict.get("SYMBOL", "")
            if not gene: 
                continue 

            variant_data = {
                "chrom": record.CHROM,
                "pos": record.POS,
                "ref": record.REF,
                "alt": record.ALT[0] if record.ALT else None, 
                "gene": gene, 
                "impact": csq_dict.get("IMPACT", ""),
                "consequence": csq_dict.get("Consequence", ""), 
                "score": float(record.INFO.get("RankScore")[2:])
            }

            parsed_variants.append(variant_data) 
        
    return parsed_variants 

def extract_genes_from_vcf(path_to_vcf: str, output_path: str) -> set:
    """
    Extract unique gene symbols from the CSQ field in a VEP-annotated VCF.

    Args:
        path_to_vcf (str): Path to a VCF file (gzipped)

    Returns:
        Set of gene symbols (str)
    """
    genes = set()
    csq_fields = []
    symbol_index = None

    with gzip.open(path_to_vcf, "rt") as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                match = re.search(r'Format: (.+)">', line)
                if match:
                    csq_fields = match.group(1).split("|")
                    symbol_index = csq_fields.index("SYMBOL") if "SYMBOL" in csq_fields else None
            elif line.startswith("#"):
                continue
            else:
                if symbol_index is None:
                    continue
                parts = line.strip().split("\t")
                info = parts[7]
                match = re.search(r'CSQ=([^;]+)', info)
                if match:
                    csq_data = match.group(1).split(",")
                    for ann in csq_data:
                        fields = ann.split("|")
                        if len(fields) > symbol_index:
                            gene = fields[symbol_index]
                            if gene:
                                genes.add(gene)
    return genes