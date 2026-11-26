

import pyranges1 as pr
from anndata import AnnData

class RangeData(AnnData):

    def __init__(
        self,
        adata: AnnData,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

        # todo extract these C's outside the class
        required_columns1 = set("chrom","start", "end")
        required_columns2 = set("Chromosome","Start", "End")

        col1_ok = required_columns1.issubset(self.var.columns)
        col2_ok = required_columns2.issubset(self.var.columns)

        if not col1_ok and not col2_ok:
            raise ValueError(f"Missing required columns: input must contain {required_columns1} or {required_columns2}")

# Figure this out: if col2 ok, then you don't want to do this:
        var = adata.var.rename(
            columns={
                "chrom": "Chromosome",
                "start": "Start",
                "end": "End",
            }
        )
        self.var = var

        self._var_ranges = self._build_var_ranges(self.var)


    @staticmethod
    def _build_var_ranges(var) -> pr.PyRanges:
        ranges_df = var.copy()

        # Keep each range linked to its row in AnnData.var.
        position_column = "_range_position"
        while position_column in ranges_df.columns:
            position_column = f"_{position_column}"

        ranges_df.insert(0, position_column, range(len(ranges_df)))
        return pr.PyRanges(ranges_df.reset_index(drop=True))




    @property
    def var_ranges(self) -> pr.PyRanges:
        return self._var_ranges

    def refresh_var_ranges(self) -> None:
        """Rebuild var_ranges after changing the range columns in .var."""
        self._var_ranges = self._build_var_ranges(self.var)

    # TODO
    # create x.bins

    def __str__(self) -> str:
        return (
            f"rangedata object with n_obs × n_vars = "
            f"{self.n_obs} × {self.n_vars}\n"
            f"    obs: {list(self.obs.columns)}\n"
            f"    var: {list(self.var.columns)}\n"
        )

    def __repr__(self) -> str:
        return self.__str__()