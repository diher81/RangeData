import anndata as ad
from pathlib import Path
from rangedata import RangeData


def read_h5ad(cls, filename: str | Path) -> RangeData:
    adata = ad.read_h5ad(filename)

    rdata = RangeData(adata)

    return rdata

