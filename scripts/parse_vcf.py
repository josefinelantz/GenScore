# scripts/parse_vcf.py

from cyvcf2 import VCF
import pandas as pd

def parse_vcf(file_path):
    """ 
    Reads VCF and returns df with colums:
    CHROM_POS_REF_ALT, CHROM_POS, AF, PP, CON, VCQF,
    LIN, CLIN, CLNSIG, RANK_SCORE
    """
    # load the VCF 
    vcf = VCF(file_path)
    # declare variable for extracted data
    data = []
    # iterate over variant in the VCF and extract required fields
    for variant in vcf:
        chrom = variant.CHROM
        pos = variant.POS
        ref = variant.REF
        alt = variant.ALT[0]
        info = variant.INFO
        # format strings for extracting controls 
        chrom_pos = f"{chrom}_{pos}"
        ref_alt = f"{ref}_{alt}"
        # separate category scores and convert to numerical values 
        af, pp, con, vcqf, lin, clin = [float(x) for x in info.get("RankResult", "0|0|0|0|0|0").split("|")]
        # extract rank score
        rank_score = float(info.get("RankScore")[2:])
        # extract CLNSIG, if not exist, return "no_value" to avoid missing values 
        clnsig = info.get("CLNSIG", "no_value") 
        # append the extracted data to the output list
        data.append({
            "VARIANT": f"{chrom_pos}_{ref_alt}",
            "CHROM_POS": chrom_pos,
            "AF": af,
            "PP": pp,
            "CON": con,
            "VCQF": vcqf,
            "LIN": lin,
            "CLIN": clin,
            "CLNSIG": clnsig,
            "RANK_SCORE": rank_score
        })

    return pd.DataFrame(data)
