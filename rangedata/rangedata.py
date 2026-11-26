import pyranges1 as pr
import numpy as np

from anndata import AnnData

from .constants import ANNDATA_COLUMNS, PYRANGES_COLUMNS, COLUMN_MAPPING

class RangeData(AnnData):

    @property
    def var_ranges(self) -> pr.PyRanges:
        return self._var_ranges.copy()

    def __init__(
            self,
            adata: AnnData
    ):
        if ANNDATA_COLUMNS.issubset(adata.var.columns):
            var = adata.var.rename(columns=COLUMN_MAPPING)
        elif PYRANGES_COLUMNS.issubset(adata.var.columns):
            var = adata.var.copy()
        else:
            raise ValueError(
                f"Input must contain either {ANNDATA_COLUMNS} "
                f"or {PYRANGES_COLUMNS}."
            )

        super().__init__(
            X=adata.X,
            obs=adata.obs,
            var=var,
            uns=adata.uns,
            obsm=adata.obsm,
            varm=adata.varm,
            layers=adata.layers,
            obsp=adata.obsp,
            varp=adata.varp,
            raw=adata.raw,
        )

        self._var_ranges = self._build_var_ranges(self.var)

    def _build_var_ranges(self, var) -> pr.PyRanges:
        ranges_df = var.copy()

        # This ID maps each range back to its row in the current .var.
        position_column = "_range_position"
        while position_column in ranges_df.columns:
            position_column = f"_{position_column}"

        self._range_position_column = position_column
        ranges_df.insert(0, position_column, range(len(ranges_df)))

        return pr.PyRanges(ranges_df.reset_index(drop=True))

    @property
    def var_ranges(self) -> pr.PyRanges:
        # Return a copy, so callers cannot mutate the stored object directly.
        return self._var_ranges.copy()


    # TODO this method is AI generated and still needs to be checked
    def sync_var_from_var_ranges(
            self,
            ranges: pr.PyRanges,
    ) -> "RangeData":
        """Return a new RangeData aligned to the supplied range results."""
        ranges_df = ranges.copy().reset_index(drop=True)
        position_column = self._range_position_column

        if position_column not in ranges_df.columns:
            raise ValueError(
                f"Range results must retain the {position_column!r} column."
            )

        raw_positions = ranges_df[position_column].to_numpy()

        try:
            positions = raw_positions.astype(np.intp)
        except (TypeError, ValueError) as exc:
            raise ValueError("Range positions must be integers.") from exc

        if not np.array_equal(raw_positions, positions):
            raise ValueError("Range positions must be integers.")

        if len(np.unique(positions)) != len(positions):
            raise ValueError(
                "Each range must map to a unique variable position."
            )

        if np.any(positions < 0) or np.any(positions >= self.n_vars):
            raise ValueError("A range position is outside the current .var.")

        # Slice AnnData in the same order as the range results.
        subset = self[:, positions].copy()
        result = type(self)(subset)

        # Carry range-column edits into the corresponding .var rows.
        updated_var = result.var.copy()
        for column in ranges_df.columns:
            if column != position_column:
                updated_var[column] = ranges_df[column].to_numpy()

        result.var = updated_var

        # Build ranges for the new object and assign fresh position IDs.
        result._var_ranges = result._build_var_ranges(result.var)
        return result


    # TODO
    # create x.original
    

    # TODO range methods
    # Structure
    # def some_range_operation(self, ...):
    #     transformed_ranges = self._var_ranges.some_pyranges_operation(...)
    #     return self.sync_var_from_var_ranges(transformed_ranges)

    def __str__(self) -> str:
        return (
            f"rangedata object with n_obs × n_vars = "
            f"{self.n_obs} × {self.n_vars}\n"
            f"    obs: {list(self.obs.columns)}\n"
            f"    var: {list(self.var.columns)}\n"
        )

    def __repr__(self) -> str:
        return self.__str__()