import pytest
import pandas as pd
from scripts.process_data import group_variants

# Mock DataFrame for testing
@pytest.fixture
def mock_data():
    return pd.DataFrame({
        "CLNSIG": [
            "Benign", 
            "Pathogenic", 
            "Uncertain_significance", 
            "not_provided", 
            "Likely_benign", 
            "random_label"
        ]
    })

# Expected results for the mock data
EXPECTED_RESULTS = [
    "benign",      # Matches "Benign" in benign group
    "pathogenic",  # Matches "Pathogenic" in pathogenic group
    "uncertain",   # Matches "Uncertain_significance" in uncertain group
    "other",       # Matches "not_provided" in other group
    "benign",      # Matches "Likely_benign" in benign group
    "random_label" # No match, returns the original label
]

# Test the group_variants function
def test_group_variants(mock_data):
    from_col = "CLNSIG"
    to_col = "GROUP"
    result = group_variants(mock_data, from_col=from_col, to_col=to_col)
    assert to_col in result.columns, "The output DataFrame should contain the new 'GROUP' column."
    assert result[to_col].tolist() == EXPECTED_RESULTS, "The grouping results are incorrect."