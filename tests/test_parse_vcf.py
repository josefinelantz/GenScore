import pytest
from scripts.extract_data import parse_vcf
import pandas as pd

def test_parse_vcf():
    test_vcf = "data/SNV.somatic.sweetelf.merged.clinical.ranked.vcf.gz"
    df = parse_vcf(test_vcf)
    
    assert isinstance(df, pd.DataFrame)
    assert "VARIANT" in df.columns
    assert "CHROM_POS" in df.columns
    assert "VARIANT" in df.columns
    assert  "AF" in df.columns
    assert  "PP" in df.columns
    assert  "CON" in df.columns
    assert  "VCQF" in df.columns
    assert  "LIN" in df.columns
    assert  "CLIN" in df.columns
    assert  "CLNSIG" in df.columns
    assert  "RANK_SCORE" in df.columns
   # assert  "PARSED_CSQ" in df.columns
   # assert  "AAF(aaf)" in df.columns
   # assert  "VAF(AF)" in df.columns
   # assert  "COVERAGE(DP)" in df.columns
