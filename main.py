from anndata import AnnData

class RangeData(AnnData):

    # This method defines the behavior of `print(RangedData)`
    def __str__(self) -> str:
        return (
            f"RangeData object with n_obs × n_vars = "
            f"{self.n_obs} × {self.n_vars}\n"
            f"    obs: {list(self.obs.columns)}\n"
            f"    var: {list(self.var.columns)}\n"
        )

    def __repr__(self) -> str:
        return self.__str__()
