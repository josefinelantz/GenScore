import pytest
from genscore.io import parse_vcf

def test_parse_mock_vcf():
    test_vcf = "data/mock.vcf.gz"
    variants = parse_vcf(test_vcf, max_variants=10)
    
    assert isinstance(variants, list)
    assert len(variants) > 0 

    for variant in variants:
        assert "gene" in variant
        assert "score" in variant 
        assert variant["gene"] != ""