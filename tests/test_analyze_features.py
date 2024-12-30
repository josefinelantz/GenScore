# tests/test_analyze_scores.py

import pytest
import pandas as pd
from scripts.analyze_scores import group_variants

def test_group_variants():
    data = {
        "CLNSIG": ["Benign", "Pathogenic", "Likely_benign", "Likely_pathogenic", "not_reported"],
        "RankScore": [1, 5, 2, 6, 3]
    }
    df = pd.DataFrame(data)
    groups, counts = group_variants(df)
    
    assert groups["benign"].shape[0] == 2
    assert groups["pathogenic"].shape[0] == 2
    assert groups["not_reported"].shape[0] == 1
    assert counts["benign"] == 2
    assert counts["pathogenic"] == 2
    assert counts["not_reported"] == 1