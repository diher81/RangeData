import pyranges1 as pr
from anndata import AnnData

from .constants import ANNDATA_COLUMNS, PYRANGES_COLUMNS, COLUMN_MAPPING

class RangeData(AnnData):

    def __init__(
            self,
            adata: AnnData,
            *args,
            **kwargs,
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