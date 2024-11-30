# scripts/parse_vcf.py

from cyvcf2 import VCF
import pandas as pd

def parse_vcf(vcf_path):
    vcf = VCF(vcf_path)
    data = []
    
    for variant in vcf:
        info = variant.INFO
        csq = info.get('CSQ', [])
        clin_sig = info.get('CLIN_SIG', 'unknown')
        rank_result = info.get('RankResult', '0|0|0|0|0|0')
        rank_score = sum(map(int, rank_result.split('|')))
        
        # Assuming CSQ is a list with fields separated by '|'
        # Adjust the parsing based on actual CSQ format
        csq_fields = csq[0].split('|') if csq else [''] * 6
        consequence = csq_fields[0] if len(csq_fields) > 0 else ''
        cosmic = csq_fields[1] if len(csq_fields) > 1 else ''
        clin_var = clin_sig
        gnomad_af = float(csq_fields[2]) if len(csq_fields) > 2 else 0.0
        sift = csq_fields[3] if len(csq_fields) > 3 else ''
        polyphen = csq_fields[4] if len(csq_fields) > 4 else ''
        
        data.append({
            'consequence': consequence,
            'cosmic': cosmic,
            'clin_var': clin_var,
            'gnomad_af': gnomad_af,
            'sift': sift,
            'polyphen': polyphen,
            'rank_result': rank_result,
            'rank_score': rank_score
        })
    
    df = pd.DataFrame(data)
    return df
