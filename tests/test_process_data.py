import pytest
import pandas as pd
from scripts.process_data import group_variants, mark_controls

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
        ],
        "RANK_SCORE": [10, 20, 15, 5, -9, 19],
        "CLIN": [1, -2, 1, -2, 2, 1]
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
    adjusted_score = "ADJUSTED_SCORE"
    result = group_variants(mock_data, from_col=from_col, to_col=to_col)
    assert to_col in result.columns, "The output DataFrame should contain the new 'GROUP' column."
    assert adjusted_score in result.columns, "The output DataFrame should contain the new 'ADJUSTED_SCORE' column."
    assert result[to_col].tolist() == EXPECTED_RESULTS, "The grouping results are incorrect."


# Test function
def test_mark_controls():
    # Create a sample DataFrame
    test_df = pd.DataFrame({
        "CHROM_POS": ["chr1_12345", "chr2_67890", "chr3_11111", "chr4_22222"],
        "RANK_SCORE": [10, 20, 15, 5],
    })
    
    # Create a controls DataFrame
    controls_df = pd.DataFrame({
        "chrom_pos": ["chr1_12345", "chr3_11111"]
    })
    
    # Expected output
    expected_output = pd.DataFrame({
        "CHROM_POS": ["chr1_12345", "chr2_67890", "chr3_11111", "chr4_22222"],
        "RANK_SCORE": [10, 20, 15, 5],
        "IS_CONTROL": [True, False, True, False]
    })

    # Run the function
    result = mark_controls(test_df.copy(), controls_df)

    # Assert the results
    pd.testing.assert_frame_equal(result, expected_output)