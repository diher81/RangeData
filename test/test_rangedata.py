import anndata as ad
import numpy as np
import pandas as pd
import pyranges1 as pr
import pytest

from rangedata import RangeData


def make_adata(var: pd.DataFrame) -> ad.AnnData:
    n_obs = 3
    n_vars = len(var)

    obs = pd.DataFrame(
        {"cell_type": ["A", "B", "A"]},
        index=["cell_1", "cell_2", "cell_3"],
    )
    X = np.arange(n_obs * n_vars).reshape(n_obs, n_vars)

    result = ad.AnnData(X=X, obs=obs, var=var)
    result.layers["counts"] = X.copy()
    return result


@pytest.fixture
def source_adata() -> ad.AnnData:
    var = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2", "chr2"],
            "start": [0, 10, 20, 30],
            "end": [5, 15, 25, 35],
            "name": ["gene_a", "gene_b", "gene_c", "gene_d"],
        },
        index=["id_a", "id_b", "id_c", "id_d"],
    )
    return make_adata(var)


@pytest.fixture
def range_data(source_adata: ad.AnnData) -> RangeData:
    return RangeData(source_adata)


def test_initializes_from_anndata(
    source_adata: ad.AnnData,
    range_data: RangeData,
) -> None:
    assert isinstance(range_data, RangeData)
    assert isinstance(range_data, ad.AnnData)

    assert range_data.shape == source_adata.shape
    np.testing.assert_array_equal(range_data.X, source_adata.X)
    assert range_data.obs.equals(source_adata.obs)
    assert range_data.layers["counts"].shape == source_adata.layers["counts"].shape


def test_renames_anndata_range_columns(range_data: RangeData) -> None:
    assert {"Chromosome", "Start", "End"}.issubset(range_data.var.columns)
    assert "chrom" not in range_data.var.columns
    assert "start" not in range_data.var.columns
    assert "end" not in range_data.var.columns

    assert range_data.var.index.tolist() == ["id_a", "id_b", "id_c", "id_d"]


def test_accepts_var_with_pyranges_column_names() -> None:
    var = pd.DataFrame(
        {
            "Chromosome": ["chr1", "chr2"],
            "Start": [0, 10],
            "End": [5, 15],
        },
        index=["id_a", "id_b"],
    )
    source = make_adata(var)

    result = RangeData(source)

    assert {"Chromosome", "Start", "End"}.issubset(result.var.columns)
    assert result.var.index.tolist() == ["id_a", "id_b"]


def test_rejects_var_without_required_range_columns() -> None:
    var = pd.DataFrame(
        {
            "chromo": ["chr1", "chr1", "chr2", "chr2"],
            "start": [0, 10, 20, 30],
            "end": [5, 15, 25, 35],
            "gene_name": ["gene_a", "gene_b", "gene_c", "gene_d"]},
        index=["id_a", "id_b", "id_c", "id_d"],
    )
    source = make_adata(var)

    with pytest.raises(ValueError, match="Input must contain either"):
        RangeData(source)


def test_var_ranges_is_stored_and_aligned(range_data: RangeData) -> None:
    ranges = range_data.var_ranges

    assert isinstance(ranges, pr.PyRanges)
    assert {"Chromosome", "Start", "End"}.issubset(ranges.columns)
    assert len(ranges) == range_data.n_vars

    assert "_range_position" in ranges.columns
    assert ranges["_range_position"].is_unique
    assert sorted(ranges["_range_position"].tolist()) == list(
        range(range_data.n_vars)
    )


def test_sync_var_from_var_ranges_updates_var_and_slices_x(
    range_data: RangeData,
) -> None:
    transformed = range_data.var_ranges
    mask = transformed["_range_position"] == 1
    transformed.loc[mask, "Start"] = 11

    updated = range_data.sync_var_from_var_ranges(transformed)

    assert updated.var.loc["id_b", "Start"] == 11
    assert updated.n_vars == range_data.n_vars
    assert updated.X.shape == (range_data.n_obs, range_data.n_vars)

    # The original object remains unchanged.
    assert range_data.var.loc["id_b", "Start"] == 10


def test_var_ranges_returns_a_defensive_copy(
        range_data: RangeData,
) -> None:
    original_starts = range_data.var_ranges["Start"].tolist()

    exposed_ranges = range_data.var_ranges
    exposed_ranges["Start"] = [-1] * len(exposed_ranges)

    # Editing the returned copy must not change the stored ranges.
    assert range_data.var_ranges["Start"].tolist() == original_starts


def test_string_and_repr(range_data: RangeData) -> None:
    expected = (
        "rangedata object with n_obs × n_vars = 3 × 4\n"
        "    obs: ['cell_type']\n"
        "    var: ['Chromosome', 'Start', 'End', 'name']\n"
    )

    assert str(range_data) == expected
    assert repr(range_data) == expected