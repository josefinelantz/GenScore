import pytest
import pandas as pd
import numpy as np
from scripts.process_data import group_variants, mark_controls

# Sample GROUPS_WITH_LABELS for testing
GROUPS_WITH_LABELS = {
    "benign": ['Benign', 'Likely_benign'],
    "pathogenic": ['Pathogenic', 'Likely_pathogenic'],
    "uncertain": ['Uncertain_significance', 'Conflicting_classifications_of_pathogenicity'],
    "other": ['', 'not_provided']
}

# Sample test data
@pytest.fixture
def sample_df():
    data = {
        "VARIANT": ["1_1000_A_T", "1_2000_G_C", "1_3000_T_A", "1_4000_C_G", "1_5000_A_C"],
        "CHROM_POS": ["1_1000", "1_2000", "1_3000", "1_4000", "1_5000"],
        "RANK_SCORE": [10.0, 15.0, -5.0, 8.0, 2.0],
        "CLIN": [2.0, 5.0, 0.0, 3.0, 1.0],
        "CLNSIG": ["Benign", "Pathogenic", "Uncertain_significance", "Unknown", "not_provided"],
    }
    return pd.DataFrame(data)

# Sample controls file (mocked as a DataFrame for testing)
@pytest.fixture
def controls_df(tmp_path):
    control_data = {
        "chrom_pos": ["1_1000", "1_3000"],
    }
    file_path = tmp_path / "controls_match.tsv"
    pd.DataFrame(control_data).to_csv(file_path, sep="\t", index=False)
    return file_path


# Test for group_variants
def test_group_variants(sample_df, controls_df):
    # Act: Call group_variants
    result = group_variants(
        sample_df.copy(),
        from_col="CLNSIG",
        groups=GROUPS_WITH_LABELS,
        to_col="GROUP",
        controls=controls_df,
    )

    # Assert: Check group assignment
    assert result.loc[0, "GROUP"] == "benign"
    assert result.loc[1, "GROUP"] == "pathogenic"
    assert result.loc[2, "GROUP"] == "uncertain"
    assert result.loc[3, "GROUP"] == "other"
    assert result.loc[4, "GROUP"] == "other"

    # Assert: Check adjusted score calculation
    assert result.loc[0, "ADJUSTED_SCORE"] == 10.0 - 2.0
    assert result.loc[1, "ADJUSTED_SCORE"] == 15.0 - 5.0

    # Assert: Check control marking
    assert result.loc[0, "IS_CONTROL"] is True
    assert result.loc[1, "IS_CONTROL"] is False
    assert result.loc[2, "IS_CONTROL"] is True
    assert result.loc[4, "IS_CONTROL"] is False


# Test for unmatched CLNSIG values
def test_group_variants(sample_df, controls_df):
    # Act: Call group_variants
    result = group_variants(
        sample_df.copy(),
        from_col="CLNSIG",
        groups=GROUPS_WITH_LABELS,
        to_col="GROUP",
        controls=controls_df,
    )

    # Assert: Check group assignment
    assert result.loc[0, "GROUP"] == "benign"
    assert result.loc[1, "GROUP"] == "pathogenic"
    assert result.loc[2, "GROUP"] == "uncertain"
    assert result.loc[3, "GROUP"] == "other"
    assert result.loc[4, "GROUP"] == "other"

    # Assert: Check adjusted score calculation
    assert result.loc[0, "ADJUSTED_SCORE"] == 10.0 - 2.0
    assert result.loc[1, "ADJUSTED_SCORE"] == 15.0 - 5.0

    # Assert: Check control marking
    assert result.loc[0, "IS_CONTROL"] == True  # Use `==` instead of `is`
    assert result.loc[1, "IS_CONTROL"] == False
    assert result.loc[2, "IS_CONTROL"] == True
    assert result.loc[3, "IS_CONTROL"] == False
    assert result.loc[4, "IS_CONTROL"] == False


# Test for mark_controls
def test_mark_controls(sample_df, controls_df):
    # Act: Call mark_controls
    result = mark_controls(sample_df.copy(), controls_df)

    # Assert: Check control marking
    assert result.loc[0, "IS_CONTROL"] == True  # Matches "1_1000"
    assert result.loc[1, "IS_CONTROL"] == False
    assert result.loc[2, "IS_CONTROL"] == True  # Matches "1_3000"
    assert result.loc[3, "IS_CONTROL"] == False
    assert result.loc[4, "IS_CONTROL"] == False


# Test for handling missing or NaN values in CLNSIG
def test_group_variants_nan_values(sample_df, controls_df):
    # Introduce NaN values in CLNSIG
    sample_df.loc[0, "CLNSIG"] = np.nan

    # Act: Call group_variants
    result = group_variants(
        sample_df.copy(),
        from_col="CLNSIG",
        groups=GROUPS_WITH_LABELS,
        to_col="GROUP",
        controls=controls_df,
    )

    # Assert: Check that NaN values are classified as "other"
    assert result.loc[0, "GROUP"] == "other"