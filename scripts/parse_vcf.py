# scripts/parse_vcf.py

from cyvcf2 import VCF
import pandas as pd
import os

# All possible consequence types (SO-terms) fetched from Ensembl release 113 (November 2024).
# The consequence types are ordered by their severeness lower index more severe
# https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html

# TODO - Move this list to a separate file 

CONSEQUENCE_ORDER = [
    "transcript_ablation",
    "splice_acceptor_variant", 
    "splice_donor_variant", 
    "stop_gained",         
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "feature_elongation",
    "feature_elongation",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_donor_5th_base_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "regulatory_region_variant",
    "intergenic_variant",
    "sequence_variant",
    "not_reported"
]
BENIGN = [
    "Benign", 
    "Likely_benign", 
    "Benign/Likely_benign"
]
PATHOGENIC = [
    "Pathogenic", 
    "Likely_pathogenic", 
    "Pathogenic/Likely_pathogenic",
    "Likely_pathogenic|association"
]
UNCERTAIN = [
    "Uncertain_significance", 
    "Conflicting_classifications_of_pathogenicity",
]
def classify(clnsig):
    if clnsig.casefold() in BENIGN:
        return "BENIGN"
    elif clnsig.casefold() in PATHOGENIC:
        return "PATHOGENIC"
    elif clnsig.casefold() in UNCERTAIN:
        return "UNCERTAIN"

def parse_vcf(vcf_file):
    """ Parses vcf_file using cyvcf2 
        Parameters: (str): path to pre-processed VCF file containing somatic SNVs scored with Genmod. 
        Returns: A Pandas DataFrame with extracted data.  
    """
    # Read VCF file into a VCF-object
    vcf = VCF(vcf_file)
    
    data = []
     
    for variant in vcf:
        chrom = variant.CHROM
        pos = variant.POS
        ref = variant.REF
        alt = variant.ALT[0]
        info = dict(variant.INFO)

        variant = f"{chrom}_{pos}_{ref}_{alt}"
        rank_result = info.get("RankResult", "0|0|0|0|0|0")
        af, pp, con, vcqf, vaf, clin = [float(x) for x in rank_result.split("|")]
        rank_score = float(info["RankScore"][2:])
        clnsig = info.get("CLNSIG", "not_reported") 

        # if clnsig and rank_score:
        #     if "benign" in clnsig.lower():
        #         data["Group"].append("benign")
        #     elif "pathogenic" in clnsig.lower():
        #         data["Group"].append("pathogenic")
        #     elif "uncertain" in clnsig.lower():
        #         data["Group"].append("uncertain")
        #     else:
        #         continue 
        group = classify(clnsig)

        data.append({
            "VARIANT": variant,
            "AF": af,
            "PP": pp,
            "CON": con,
            "VCQF": vcqf,
            "VAF": vaf,
            "CLIN": clin,
            "CLNSIG": clnsig,
            "RankScore": rank_score, 
            "Group": group
        })

    df = pd.DataFrame(data)
    # benign = df[df["Group"] == "BENIGN"]
    # pathogenic = df[df["Group"] == "PATHOGENIC"]
    # uncertain = df[df["Group"] == "UNCERTAIN"]
    # grouped_df = pd.concat([benign, pathogenic, uncertain])
    return df