#
# This file is part of the mult_res_cell_ann distribution (https://github.com/tasosc/mult_res_cell_ann).
# Copyright (c) 2024 Anastasios Chronis.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
from io import StringIO
import logging
from os import PathLike
from typing import Iterator
import scanpy as sc
import anndata
import warnings
import gzip
from pathlib import Path

warnings.filterwarnings("ignore")
from utilities import Config, Render

class PreProcessing:
    """
    Contains the methods for pre-processing

    The workflow is based on  https://github.com/kostaslazaros/cell_annotation_web_app/blob/main/decoupler_cell_annotation.ipynb

    How to construct an instance of this class
    -----------

    1. If you have already `anndata` then simply use the constructor of this class: `pp = PreProcessing(anndata, config)`
    2. Else if you have a csv (gzip'ed or not) or a h5ad file, use the class method `pp = PreProcessing.build_from(uploaded_file, config)`

    How to use
    ----------
    Ideally all commands must be run in following order.
    ```python
    pp.qc()
    pp.normalization()
    pp.feature_selection()
    pp.dimensionality_reduction()
    pp.visualization()
    ```

    Note the `pp.adata` is modified through out the steps and will not be the same object

    """
    logger = logging.getLogger("PreProcessing")
    def __init__(self, adata: anndata.AnnData, config: Config) -> None:
        """
        Create a new instance of PreProcessing class

        Parameters
        ----------
        adata : anndata.AnnData
            The single cell RNA sequencing dataset
        confg: Config
            Configuration options for various preprocessing commands
        """
        self.config = config
        self.adata = adata
        self.render = Render()

    def qc(self) -> None:
        """
        Run the Quality Control step
        """
        # Based on https://github.com/kostaslazaros/cell_annotation_web_app/blob/main/adata_preprocessor.py#L10
        self.render.render_text(self.adata,3)
        fadata = self.adata
        n_genes_min=self.config.defaults["n_genes_min"]
        n_genes_max=self.config.defaults["n_genes_max"]
        min_genes=self.config.defaults["min_genes"]
        min_cells=self.config.defaults["min_cells"]
        n_counts_max=self.config.defaults["n_counts_max"]
        pc_mito = self.config.defaults["pc_mito"]
        pc_rib = self.config.defaults["pc_rib"]
        # Pre-filtering
        sc.pp.filter_cells(fadata, min_genes=min_genes)  # Equivalent to min.features in Seurat.
        self.render.render_text(f"Filtering cells with number of genes < {min_genes}: {fadata.shape}",2)

        sc.pp.filter_genes(fadata, min_cells=min_cells)  # Equivalent to min.cells in Seurat.
        self.render.render_text(f"Filtering genes expressed in < {min_cells} cells: {fadata.shape}",2)

        # Calculate the percentage of mitochondrial genes.
        mito_genes = fadata.var_names.str.startswith(tuple(['MT-', 'mt-', 'MT.', "mt."]))
        fadata.obs['prc_mt'] = (fadata[:, mito_genes].X.sum(axis=1) / fadata.X.sum(axis=1)) * 100
        self.render.render_text("Mitochondrial gene percentage calculated and annotated in the prc_mt observation",2)

        # Calculate the percentage of ribosomal genes.
        ribo_genes = fadata.var_names.str.startswith('RPS')
        fadata.obs['prc_rb'] = (fadata[:, ribo_genes].X.sum(axis=1) / fadata.X.sum(axis=1)) * 100
        self.render.render_text("Ribosomal gene percentage calculated and annotated in the prc_rb observation",2)

        # Calculate number of genes and counts for each cell.
        fadata.obs['n_genes'] = (fadata.X > 0).sum(axis=1)
        self.render.render_text("Calculate number of genes with non-zero counts",2)

        fadata.obs['n_counts'] = fadata.X.sum(axis=1)
        self.render.render_text("Calculate total number of counts for each cell",2)

        # Subsetting the data based on the calculated values.
        fadata = fadata[fadata.obs['n_genes'] > n_genes_min, :]
        self.render.render_text(f"Filter cells with too few genes detected: {fadata.shape}",2)

        fadata = fadata[fadata.obs['n_genes'] < n_genes_max, :]
        self.render.render_text(f"Filter cells with too many genes detected: {fadata.shape}",2)

        fadata = fadata[fadata.obs['n_counts'] < n_counts_max, :]
        self.render.render_text(f"Filter cells with too many counts detected: {fadata.shape}",2)

        fadata = fadata[fadata.obs['prc_mt'] < pc_mito, :]
        self.render.render_text(f"Filter cells with too many mitochondrial genes expressed: {fadata.shape}",2)

        fadata = fadata[fadata.obs['prc_rb'] < pc_rib, :]
        self.render.render_text(f"Filter cells with too many ribosomal genes expressed: {fadata.shape}",2)
        self.adata = fadata

    def feature_selection(self) -> None:
        """
        Run the feature selection step, only if `config.get_bool("only_highly_significant_genes")` returns true
        """
        if self.config.get_bool("only_highly_significant_genes"):
            self.render.render_text("Identify the most highly variable genes")
            # Identify the most highly variable genes
            sc.pp.highly_variable_genes(self.adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
            # Filter higly variable genes
            self.render.render_text("Filter high variable genes")
            self.adata.raw = self.adata
            self.adata = self.adata[:, self.adata.var.highly_variable]
            self.render.render_text(self.adata,3)
        else:
            self.render.render_text("Skipped")

    def normalization(self) -> None:
        """
        Run the normalization step
        """
        # normalization
        if self.config.get_bool("normalize_total_counts"):
            sc.pp.normalize_total(self.adata, target_sum=1e4)
        sc.pp.log1p(self.adata)

    def dimensionality_reduction(self) -> None:
        """
        Run the dimension reduction step
        """
        sc.tl.pca(self.adata, svd_solver=self.config.defaults["svd_solver"])

    def visualization(self):
        """
        Run the visualization step
        """
        # Compute distances in the PCA space, and find cell neighbors
        sc.pp.neighbors(
            self.adata, n_neighbors=self.config.defaults["n_neigh"], n_pcs=self.config.defaults["n_pcs"]
        )

        # Perform leiden clustering
        sc.tl.leiden(
            self.adata,
            resolution=self.config.defaults["cluster_resolution"],
            key_added=self.config.defaults["leiden_key"],
        )

        # visualization
        # Calculate UMAP embeddings
        sc.tl.umap(self.adata)
        self.render.render_fig(sc.pl.umap(
            self.adata,
            color="leiden",
            title=f'Leiden clustering (Resolution: {self.config.defaults["cluster_resolution"]})',
            frameon=True,
            legend_fontweight="normal",
            legend_fontsize=10,
            return_fig=True,
        ), 1)

    @classmethod
    def build_from_csv(cls, path: PathLike|Iterator[str], config: Config):
        """
        Build a class instance from a csv file.

        Parameters
        ----------
        path : PathLike|Iterator[str]
            The single cell RNA sequencing dataset as a file or something read-able
        confg: Config
            Configuration options for various preprocessing commands
        """
        delimiter = '\t' if config.defaults["csv_delimiter"] == 'T' else config.defaults["csv_delimiter"]
        adata = sc.read_csv(path, delimiter=delimiter, first_column_names=True).T
        return PreProcessing(adata, config)
    @classmethod
    def build_from_hdf5(cls, path: PathLike|Iterator[str], config: Config):
        """
        Build a class instance from a HDF5 file.

        Parameters
        ----------
        path : PathLike|Iterator[str]
            The single cell RNA sequencing dataset as a file or something read-able
        confg: Config
            Configuration options for various preprocessing commands
        """
        adata = sc.read_h5ad(path)
        return PreProcessing(adata=adata, config=config)
    @classmethod
    def build_from(cls, uploaded_file, config :Config):
        """
         Build a class instance from a HDF5 or CSV (gzip'ed or not) file.

        Parameters
        ----------
        uploaded_file : PathLike|Iterator[str] with .name and .type fields
            The single cell RNA sequencing dataset as a file or something read-able
        confg: Config
            Configuration options for various preprocessing commands
        """
        if config is None:
            raise ValueError("Internal error, configuration not found")
        if uploaded_file is None:
            config.defaults["data"]=Path(config.defaults["h5ad_path"]).stem
            return cls.build_from_hdf5(path=config.defaults["h5ad_path"], config=config)
        cls.logger.info(f"filename:{uploaded_file.name} and type {uploaded_file.type}")
        config.defaults["data"]=Path(uploaded_file.name).stem
        if uploaded_file.type != 'application/gzip':
            if uploaded_file.name.endswith(".csv"):
                return cls.build_from_csv(path=StringIO(uploaded_file.getvalue().decode("utf-8")), config=config)
            return cls.build_from_hdf5(path=uploaded_file, config=config)

        with gzip.open(uploaded_file, mode='rt') as gcsv:
            return cls.build_from_csv(path=gcsv, config=config)
