from anndata import AnnData
from .constants import REQUIRED_COLUMNS_1, REQUIRED_COLUMNS_2

def __init__(
    self,
    adata: AnnData,
    *args,
    **kwargs,
):
    super().__init__(*args, **kwargs)

    col1_ok = REQUIRED_COLUMNS_1.issubset(self.var.columns)
    col2_ok = REQUIRED_COLUMNS_2.issubset(self.var.columns)

    if not col1_ok and not col2_ok:
        raise ValueError(f"Missing required columns: input must contain {REQUIRED_COLUMNS_1} or {REQUIRED_COLUMNS_2}")

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