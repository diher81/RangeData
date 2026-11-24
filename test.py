import scanpy as sc

from main import RangeData
from pathlib import Path

base_dir = Path.home() / "Library/CloudStorage/OneDrive-WageningenUniversity&Research/InternshipTBB/Data"
output_dir = base_dir / "output"

# --- LOADING DATA ---
adata_H3K4me3 = sc.read_h5ad(
    f"{base_dir}/gast_H3K4me3_2000_bins_ctype.h5ad"
)

# Create a RangedData object
rangeData = RangeData(adata_H3K4me3)

# Test the dimensions
assert rangeData.n_obs == 18279
assert rangeData.n_vars == 1364141

# Test your custom string representation
print(rangeData)

# Test repr
print(repr(rangeData))

# Test inheritance
assert isinstance(rangeData, RangeData)
assert isinstance(rangeData, __import__("anndata").AnnData)

print("All tests passed!")