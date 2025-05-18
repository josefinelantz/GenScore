from cyvcf2 import VCF
import json
from output.constants import INFO_KEYS, CSQ_KEYS

def parse_csq(fields, csq_str):
    csq = []
    annotation_list = csq_str.split(",")
    for annotation in annotation_list:
        csq.append(dict(zip(fields, annotation.split("|"))))
        return csq

def extract_csq_keys(parsed_csq, keys=CSQ_KEYS):
    annos = []
    values = []

    for anno in parsed_csq:
        for item in anno.items():
            if item[0] in keys:
                annos.append(item[0])
                values.append(item[1])
    return dict(zip(annos, values))

def extract_info_keys(info_field, keys=INFO_KEYS):
    annos = []
    values = []

    for key in keys:
        annos.append(key)
        values.append(info_field.get(key, ""))

    return dict(zip(annos, values))

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
        filter_value = variant.FILTER

        info = variant.INFO
        #info_keys = ["CSQ", "RankResult", "RankScore", "CLNSIG", "AF", "DP", "GNOMAD_AF"]
        # extract INFO-keys
        info_keys = extract_info_keys(info, keys=["CSQ", "RankResult", "RankScore", "CLNSIG", "AF", "DP", "GNOMAD_AF"])

        #csq_str = info.get("CSQ", "")
        af, pp, con, vcqf, lin, clin = [float(x) for x in info_keys.get("RankResult", "0|0|0|0|0|0").split("|")]
        rank_score = float(info_keys.get("RankScore")[2:])
        clnsig = info_keys.get("CLNSIG", "")
        vaf = info_keys.get("AF", "")
        coverage = info_keys.get("DP", "")

        # Parse CSQ
        csq = parse_csq(csq_fields, info_keys.get("CSQ", ""))
        csq_keys = extract_csq_keys(csq, keys = ["Consequence", "SIFT", "PolyPhen", "gnomAD_AF", "COSMIC", "CLIN_SIG"])

        # Format strings for extracting controls
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
            "Consequence": csq_keys.get("Consequence", ""),
            "SIFT": csq_keys.get("SIFT", ""),
            "PolyPhen": csq_keys.get("PolyPhen", ""),
            "gnomAD_AF": csq_keys.get("gnomAD_AF", ""),
            "COSMIC": csq_keys.get("COSMIC", ""),
            "CLIN_SIG": csq_keys.get("CLIN_SIG", ""),
            "FILTER": filter_value,
            "VAF(AF)": vaf,
            "COVERAGE(DP)": coverage
        })

    return pd.DataFrame(data)
