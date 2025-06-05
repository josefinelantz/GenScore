from cyvcf2 import VCF

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
