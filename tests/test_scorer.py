import pytest
from genscore.scorer import aggregate_by_gene 

def test_aggregate_by_gene():
    variant_list = [
        {"gene": "TP53", "score": 10.0, "impact": "HIGH", "consequence": "stop_gained"},
        {"gene": "TP53", "score": 7.5, "impact": "MODERATE", "consequence": "missense_variant"},
        {"gene": "BRCA1", "score": 5.0, "impact": "HIGH", "consequence": "frameshift_variant"},
        {"gene": "BRCA1", "score": None, "impact": "LOW", "consequence": "synonymous_variant"},
        {"gene": "BRCA2", "score": "invalid", "impact": "MODERATE", "consequence": "inframe_deletion"}
    ]

    result = aggregate_by_gene(variant_list)

    assert "TP53" in result
    assert result["TP53"]["total_score"] == 17.5
    assert result["TP53"]["n_variants"] == 2
    assert "HIGH" in result["TP53"]["impact_levels"]
    assert "missense_variant" in result["TP53"]["consequences"]

    assert "BRCA1" in result
    assert result["BRCA1"]["total_score"] == 5.0
    assert result["BRCA1"]["n_variants"] == 1  # The None score variant is skipped

    assert "BRCA2" not in result  # Invalid score string is skipped