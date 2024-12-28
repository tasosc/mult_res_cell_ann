from pydantic import BaseModel

from model.enums import ReportingOptions, SvdSolverOptions


class Settings(BaseModel):
    cluster_resolution: float = 0.4
    svd_solver: SvdSolverOptions = SvdSolverOptions.arpack
    leiden_key: str = "leiden"
    n_genes_min: int = 1000
    n_genes_max: int = 10000
    min_genes: int = 100
    min_cells: int = 3
    n_counts_max: int = 30000
    pc_mito: int = 20
    pc_rib: int = 25
    n_neigh: int = 10
    n_pcs: int = 40
    csv_delimiter: str = ","
    normalize_total_counts: bool = False
    only_highly_significant_genes: bool = False
    verbosity: int = 1
    output: ReportingOptions = ReportingOptions.as_progress | ReportingOptions.pdf
# TODO init from Configuration