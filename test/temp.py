from rangedata.rangedata import RangeData
from pathlib import Path

base_dir = "Library/CloudStorage/OneDrive-WageningenUniversity&Research/InternshipTBB/Data"

# --- LOADING DATA ---
DATA_FILE = (
    Path.home()
    / base_dir
    / "gast_H3K4me3_2000_bins_ctype.h5ad"
)

range_data = RangeData.from_h5ad(DATA_FILE)

# Test custom string representation
print(range_data)

# Test repr
print(repr(range_data))