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
    # Add more assertions based on known test data