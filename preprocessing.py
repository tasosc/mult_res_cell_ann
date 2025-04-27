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
import matplotlib

from render_to_pdf import RenderToPdf
matplotlib.use('Agg')
from io import StringIO
import logging
from os import PathLike
from typing import Iterator
import warnings
import gzip
from pathlib import Path
import scanpy as sc
import anndata
from model.settings import Settings
warnings.filterwarnings("ignore")

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

    Note the `pp.adata` is modified throughout the steps and will not be the same object

    """
    logger = logging.getLogger("PreProcessing")

    def __init__(self, adata: anndata.AnnData, config: Settings, render: RenderToPdf = None) -> None:
        """
        Create a new instance of PreProcessing class

        Parameters
        ----------
        adata : anndata.AnnData
            The single cell RNA sequencing dataset
        config: Settings
            Configuration options for various preprocessing commands
        """
        self.config = config
        self.adata = adata
        self.render = render

    def set_render(self, render: RenderToPdf):
        self.render = render
    
    def qc(self) -> None:
        """
        Run the Quality Control step
        """
        # Based on https://github.com/kostaslazaros/cell_annotation_web_app/blob/main/adata_preprocessor.py#L10
        self.render.write_text(str(self.adata), 3)
        fadata = self.adata
        n_genes_min = self.config.n_genes_min
        n_genes_max = self.config.n_genes_max
        min_genes = self.config.min_genes
        min_cells = self.config.min_cells
        n_counts_max = self.config.n_counts_max
        pc_mito = self.config.pc_mito
        pc_rib = self.config.pc_rib
        # Pre-filtering
        sc.pp.filter_cells(fadata, min_genes=min_genes)  # Equivalent to min.features in Seurat.
        self.render.write_text(f"Filtering cells with number of genes < {min_genes}: {fadata.shape}", 1)

        sc.pp.filter_genes(fadata, min_cells=min_cells)  # Equivalent to min.cells in Seurat.
        self.render.write_text(f"Filtering genes expressed in < {min_cells} cells: {fadata.shape}", 1)

        # Calculate the percentage of mitochondrial genes.
        mito_genes = fadata.var_names.str.startswith(tuple(['MT-', 'mt-', 'MT.', "mt."]))
        fadata.obs['prc_mt'] = (fadata[:, mito_genes].X.sum(axis=1) / fadata.X.sum(axis=1)) * 100
        self.render.write_text("Mitochondrial gene percentage calculated and annotated in the prc_mt observation", 1)

        # Calculate the percentage of ribosomal genes.
        ribo_genes = fadata.var_names.str.startswith('RPS')
        fadata.obs['prc_rb'] = (fadata[:, ribo_genes].X.sum(axis=1) / fadata.X.sum(axis=1)) * 100
        self.render.write_text("Ribosomal gene percentage calculated and annotated in the prc_rb observation", 1)

        # Calculate number of genes and counts for each cell.
        fadata.obs['n_genes'] = (fadata.X > 0).sum(axis=1)
        self.render.write_text("Calculate number of genes with non-zero counts", 1)

        fadata.obs['n_counts'] = fadata.X.sum(axis=1)
        self.render.write_text("Calculate total number of counts for each cell", 1)

        # Subsetting the data based on the calculated values.
        fadata = fadata[fadata.obs['n_genes'] > n_genes_min, :]
        self.render.write_text(f"Filter cells with too few genes detected: {fadata.shape}", 1)

        fadata = fadata[fadata.obs['n_genes'] < n_genes_max, :]
        self.render.write_text(f"Filter cells with too many genes detected: {fadata.shape}", 1)

        fadata = fadata[fadata.obs['n_counts'] < n_counts_max, :]
        self.render.write_text(f"Filter cells with too many counts detected: {fadata.shape}", 1)

        fadata = fadata[fadata.obs['prc_mt'] < pc_mito, :]
        self.render.write_text(f"Filter cells with too many mitochondrial genes expressed: {fadata.shape}", 1)

        fadata = fadata[fadata.obs['prc_rb'] < pc_rib, :]
        self.render.write_text(f"Filter cells with too many ribosomal genes expressed: {fadata.shape}", 1)
        self.adata = fadata

    def feature_selection(self) -> None:
        """
        Run the feature selection step, only if `config.only_highly_significant_genes` returns true
        """
        if self.config.only_highly_significant_genes:
            self.render.write_text("Identify the most highly variable genes")
            # Identify the most highly variable genes
            sc.pp.highly_variable_genes(self.adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
            # Filter higly variable genes
            self.render.write_text("Filter high variable genes")
            self.adata.raw = self.adata
            self.adata = self.adata[:, self.adata.var.highly_variable]
            self.render.write_text(str(self.adata), 3)
        else:
            self.render.write_text("Skipped")

    def normalization(self) -> None:
        """
        Run the normalization step
        """
        # normalization
        if self.config.normalize_total_counts:
            self.render.write_text("Normalizing counts per cell with target sum 1e4")
            sc.pp.normalize_total(self.adata, target_sum=1e4, inplace=True)
        
        self.render.write_text("Logarithmizing the data matrix.")
        sc.pp.log1p(self.adata)

    def dimensionality_reduction(self) -> None:
        """
        Run the dimension reduction step
        """
        self.render.write_text("Running principal component analysis")
        sc.tl.pca(self.adata, svd_solver=self.config.svd_solver.name)

    def visualization(self):
        """
        Run the visualization step
        """
        # Compute distances in the PCA space, and find cell neighbors
        self.render.write_text("Computing distances in the PCA space, and finding cell neighbors")
        sc.pp.neighbors(
            self.adata, n_neighbors=self.config.n_neigh, n_pcs=self.config.n_pcs
        )

        # Perform leiden clustering
        self.render.write_text("Performing leiden clustering")
        sc.tl.leiden(
            self.adata,
            resolution=self.config.cluster_resolution,
            key_added=self.config.leiden_key,
        )

        # visualization
        # Calculate UMAP embeddings
        self.render.write_text("Calculating leiden clustering")
        sc.tl.umap(self.adata)
        return sc.pl.umap(
            self.adata,
            color="leiden",
            title=f'Leiden clustering (Resolution: {self.config.cluster_resolution})',
            frameon=True,
            legend_fontweight="normal",
            legend_fontsize=10,
            return_fig=True,
        )

    @classmethod
    def build_from_txt(cls, path: PathLike | Iterator[str], config: Settings):
        """
        Build a class instance from a txt  file.

        e.g. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4773521

        Parameters
        ----------
        path : PathLike|Iterator[str]
            The single cell RNA sequencing dataset as a file or something read-able
        config: Settings
            Configuration options for various preprocessing commands
        """
        adata = sc.read_text(path, first_column_names=True).T
        return PreProcessing(adata, config)

    @classmethod
    def build_from_csv(cls, path: PathLike | Iterator[str], config: Settings):
        """
        Build a class instance from a csv file.

        e.g. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5226584

        Parameters
        ----------
        path : PathLike|Iterator[str]
            The single cell RNA sequencing dataset as a file or something read-able
        config: Settings
            Configuration options for various preprocessing commands
        """
        delimiter = '\t' if config.csv_delimiter == 'T' else config.csv_delimiter
        adata = sc.read_csv(path, delimiter=delimiter, first_column_names=True).T
        return PreProcessing(adata, config)

    @classmethod
    def build_from_hdf5(cls, path: str | Path, config: Settings):
        """
        Build a class instance from a HDF5 file.

        Parameters
        ----------
        path : PathLike|Iterator[str]
            The single cell RNA sequencing dataset as a file or something read-able
        config: Settings
            Configuration options for various preprocessing commands
        """
        adata = sc.read_h5ad(path)
        return PreProcessing(adata=adata, config=config)

    @classmethod
    def build_from(cls, uploaded_file: Path, config: Settings):
        """
         Build a class instance from a HDF5 or CSV (gzip'ed or not) file.

        Parameters
        ----------
        uploaded_file : PathLike|Iterator[str] with .name and .type fields
            The single cell RNA sequencing dataset as a file or something read-able
        config: Settings
            Configuration options for various preprocessing commands
        """
        if config is None:
            raise ValueError("Internal error, configuration not found")
        if uploaded_file is None:
            return cls.build_from_hdf5(path="./data/GSM4089151_P1.h5ad", config=config)
        extension= "".join(uploaded_file.suffixes) if uploaded_file.suffix == ".gz" else uploaded_file.suffix
        cls.logger.info("filename: %s and type %s", uploaded_file.name, extension)
        if extension == '.csv':
            with StringIO(uploaded_file.getvalue().decode("utf-8")) as csv:
                return cls.build_from_csv(path=csv, config=config)
        if extension == '.txt':
            with StringIO(uploaded_file.getvalue().decode("utf-8")) as txt:
                return cls.build_from_txt(path=txt, config=config)
        if extension == '.csv.gz':
            with gzip.open(uploaded_file, mode='rt') as gcsv:
                return cls.build_from_csv(path=gcsv, config=config)
        if extension == '.txt.gz':
            with gzip.open(uploaded_file, mode='rt') as gcsv:
                return cls.build_from_txt(path=gcsv, config=config)
        return cls.build_from_hdf5(path=uploaded_file, config=config)
