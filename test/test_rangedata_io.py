import anndata as ad
import pandas as pd
import pyranges1 as pr
import pytest

from rangedata.rangedata import RangeData
from pathlib import Path
from pandas.testing import assert_frame_equal

base_dir = "Library/CloudStorage/OneDrive-WageningenUniversity&Research/InternshipTBB/Data"

# --- LOADING DATA ---
DATA_FILE = (Path.home() / base_dir / "gast_H3K4me3_2000_bins_ctype.h5ad")

EXPECTED_N_OBS = 18_279
EXPECTED_N_VARS = 1_364_141

@pytest.fixture(scope="module")
def range_data() -> RangeData:
    if not DATA_FILE.is_file():
        pytest.skip(f"Test data file not found: {DATA_FILE}")

    return RangeData.read_h5ad(DATA_FILE)

def test_loads_as_rangedata_and_anndata(range_data: RangeData) -> None:
    assert isinstance(range_data, RangeData)
    assert isinstance(range_data, ad.AnnData)

def test_dimensions(range_data: RangeData) -> None:
    assert range_data.n_obs == EXPECTED_N_OBS
    assert range_data.n_vars == EXPECTED_N_VARS
    assert range_data.X.shape == (EXPECTED_N_OBS, EXPECTED_N_VARS)

def test_var_remains_a_dataframe(range_data: RangeData) -> None:
    assert isinstance(range_data.var, pd.DataFrame)
    assert {"Chromosome", "Start", "End"}.issubset(range_data.var.columns)

def test_var_ranges_is_stored_and_has_expected_columns(
    range_data: RangeData,
) -> None:
    ranges = range_data.var_ranges

    assert isinstance(ranges, pr.PyRanges)
    assert {"Chromosome", "Start", "End"}.issubset(ranges.columns)

    # Confirms this is the stored object, not a newly constructed one.
    assert ranges is range_data.var_ranges

def test_var_ranges_has_one_entry_per_variable(
    range_data: RangeData,
) -> None:
    ranges = range_data.var_ranges

    assert len(ranges) == range_data.n_vars
    assert "_range_position" in ranges.columns
    assert ranges["_range_position"].is_unique

def test_var_ranges_matches_var_coordinates(
    range_data: RangeData,
) -> None:
    expected = (
        range_data.var
        .rename(columns={
            "chrom": "Chromosome",
            "start": "Start",
            "end": "End",
        })
        [["Chromosome", "Start", "End"]]
        .reset_index(drop=True)
    )
    expected.insert(0, "_range_position", range(len(expected)))

    actual = (
        range_data.var_ranges.df
        .sort_values("_range_position")
        .reset_index(drop=True)
    )
    actual = actual[expected.columns]

    assert_frame_equal(actual, expected)