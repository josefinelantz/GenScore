# tests/test_visualize.py

import pytest
import pandas as pd
from scripts.visualize import plot_box, plot_stacked
import os

def test_plot_box(tmp_path):
    groups = {
        "benign": pd.DataFrame({"RankScore": [1, 2, 3]}),
        "pathogenic": pd.DataFrame({"RankScore": [5, 6, 7]})
    }
    output_path = os.path.join(tmp_path, "boxplot.png")
    plot_box(groups, output_path)
    assert os.path.exists(output_path)

def test_plot_stacked(tmp_path):
    groups = {
        "benign": pd.DataFrame({"RankResult": ["1|0|0|0|0|0", "2|0|0|0|0|0"]}),
        "pathogenic": pd.DataFrame({"RankResult": ["5|0|0|0|0|0", "6|0|0|0|0|0"]})
    }
    output_path = os.path.join(tmp_path, "stackedplot.png")
    plot_stacked(groups, output_path)
    assert os.path.exists(output_path)