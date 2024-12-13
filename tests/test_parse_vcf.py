# tests/test_parse_vcf.py

import pytest
from scripts.parse_vcf import parse_vcf
import pandas as pd

def test_parse_vcf():
    test_vcf = "data/SNV.somatic.sweetelf.merged.clinical.ranked.vcf.gz"
    df = parse_vcf(test_vcf)
    
    assert isinstance(df, pd.DataFrame)
    assert "CLNSIG" in df.columns
    assert "RankScore" in df.columns
    assert "VARIANT" in df.columns
    assert  "AF" in df.columns
    assert  "PP" in df.columns
    assert  "CON" in df.columns
    assert  "VCQF" in df.columns
    assert  "VAF" in df.columns
    assert  "CLIN" in df.columns
    assert  "CLNSIG" in df.columns
    assert  "RankScore" in df.columns
    assert  "Group" in df.columns
   