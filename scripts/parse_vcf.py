from cyvcf2 import VCF
import pandas as pd

def parse_csq(fields, csq_str):
    csq = []
    annotation_list = csq_str.split(",")
    for annotation in annotation_list:   
        return csq.append(dict(zip(fields, annotation.split("|"))))

def parse_vcf(file_path):
    """ 
    Reads VCF, extracts and formats variant information 
    relevant to scoring.  
    Args: 
        file_path (str): path to a valid VCF. 
    
    Returns (pd.DataFrame): df with colums:
    CHROM_POS_REF_ALT, CHROM_POS, AF, PP, CON, VCQF,
    LIN, CLIN, CLNSIG, RANK_SCORE, PARSED_CSQ, AAF(aaf), FILTER, VAF(AF), COVERAGE(DP)
    """
    # Load VCF into cyvcf2.VCF generator object 
    vcf = VCF(file_path)

    # Get the headers for VEP annotations
    csq_fields = vcf.get_header_type("CSQ")["Description"][51:-1].split("|")
    
    data = []
    # Iterate the variants loaded and extract required fields
    for variant in vcf:
        chrom = variant.CHROM
        pos = variant.POS
        ref = variant.REF
        alt = variant.ALT[0]
        info = variant.INFO
        aaf = variant.aaf 
        filter_value = variant.FILTER

        # extract INFO-keys
        csq_str = info.get("CSQ", "") 
        af, pp, con, vcqf, lin, clin = [float(x) for x in info.get("RankResult", "0|0|0|0|0|0").split("|")]   
        rank_score = float(info.get("RankScore")[2:])    
        clnsig = info.get("CLNSIG", "no_value") 
        vaf = info.get("AF", "")
        coverage = info.get("DP", "")

        # Extract CSQ-keys 
        csq = parse_csq(csq_fields, csq_str)

        # format strings for extracting controls 
        chrom_pos = f"{chrom}_{pos}"
        ref_alt = f"{ref}_{alt}"
       
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
            "RANK_SCORE": rank_score,
            "PARSED_CSQ": csq,
            "AAF(aaf)": aaf,
            "FILTER": filter_value,
            "VAF(AF)": vaf,
            "COVERAGE(DP)": coverage
        })

    return pd.DataFrame(data)
