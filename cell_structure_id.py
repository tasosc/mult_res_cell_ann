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
import logging
from os import PathLike
import scanpy as sc
import hdf5plugin
import numpy as np
import pandas as pd
import decoupler as dc
import anndata
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

from utilities import CellType, Config, Render

warnings.filterwarnings("ignore")
logger = logging.getLogger("Store")


class StructureIdentification:
    """
    Contains the methods for Structure Identification ( clustering, cell annotation)

    The workflow is based on  https://github.com/kostaslazaros/cell_annotation_web_app/blob/main/decoupler_cell_annotation.ipynb

    How to construct an instance of this class
    -----------

    `si = StructureIdentification(adata, config)`


    How to use
    ----------
    Ideally all commands must be run in following order.
    ```python
    si.clustering(cell_types)
    si.annotation()
    si.write_ann_ds(output)
    ```

    Note the `si.adata` is modified throughout the steps and may not be the same object

    """

    def __init__(self, adata: anndata.AnnData, config: Config) -> None:
        """
        Create a new instance of StructureIdentification class

        Parameters
        ----------
        adata : anndata.AnnData
            The single cell RNA sequencing dataset, it expects to use the :attr:`~preprocessing.PreProcessing.adata` after all Pre Processing steps have run
        config: Config
            Configuration options for various preprocessing commands
        """
        self.config = config
        self.adata = adata
        self.render = Render()
        self.acts = None

    def clustering(self, cell_types: list[CellType]):
        """
        Perform clustering

        Parameters
        ----------
        cell_types: list[CellType]
               The list of cell types to use an input in ORA after they are filtered. Only CellType with selected genes, new genes or all selected genes will be used.
        """
        # Filter cell marker dataframe to obtain markers related to give cell types
        self.render.render_text(
            "Filter cell marker dataframe to obtain markers related to give cell types",
            2,
        )
        filtered_marker_df = self.filtered_marker(cell_types)
        filtered_marker_df = filtered_marker_df[
            ~filtered_marker_df.duplicated(["cell_name", "Symbol"])
        ]
        filtered_marker_df = filtered_marker_df.dropna()
        self.render.render_text(filtered_marker_df, 3)
        # Enrichment with Over Representation Analysis (ORA)
        self.render.render_text("Enrichment with Over Representation Analysis (ORA)", 2)
        dc.run_ora(
            mat=self.adata,
            net=filtered_marker_df,
            source="cell_name",
            target="Symbol",
            min_n=3,
            verbose=True,
            use_raw=self.config.get_bool("only_highly_significant_genes"),
        )
        # The obtained scores (-log10(p-value))(ora_estimate) and p-values (ora_pvals) are stored in the .obsm key
        self.render.render_text(
            "The obtained scores (-log10(p-value))(ora_estimate) and p-values (ora_pvals) are stored in the .obsm key",
            2,
        )
        acts: anndata.AnnData = dc.get_acts(self.adata, obsm_key="ora_estimate")
        # We need to remove inf and set them to the maximum value observed for pvals=0
        acts_v = acts.X.ravel()
        max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
        acts.X[~np.isfinite(acts.X)] = max_e
        self.render.render_text(acts, 3)
        self.acts = acts
        # Create cell-type list and ORA-score dataframe
        self.render.render_text("Create cell-type list and ORA-score dataframe", 2)
        score_df = self.adata.obsm["ora_estimate"]
        ctype_lst = list(score_df.columns)
        score_df["cluster"] = self.adata.obs["leiden"]
        self.render.render_text(ctype_lst, 3)
        melted_df = self.__create_melted_df(score_df, ctype_lst)
        self.render.render_text(melted_df, 3)
        self.render.render_text(
            "Create ORA-score violin plots for all leiden clusters and cell-types", 2
        )
        self.render.render_fig(self.cluster_vln_plot(melted_df))

    @staticmethod
    def __create_melted_df(score_df, ctype_lst):
        """
        Create the melted dataframe, unpivoting (columns to values) the dataframe but leaving some IDs intact
        Copied as it is from https://github.com/kostaslazaros/cell_annotation_web_app/blob/main/adata_preprocessor.py#L96

        """
        # Step 1: Creating separate DataFrames for each cluster
        cluster_dfs = [
            score_df[score_df["cluster"] == str(i)]
            for i in range(len(score_df["cluster"].value_counts()))
        ]

        # Step 2: Adding a 'dataset' column
        dfs = [df.assign(dataset=f"dataset_{i}") for i, df in enumerate(cluster_dfs)]

        # Checking the assignment of the 'dataset' column

        # Step 3: Concatenating DataFrames into one
        final_df = pd.concat(dfs, ignore_index=True)

        # Confirming that the 'dataset' column is present in the concatenated DataFrame

        # Step 4: Melting the DataFrame for seaborn plotting
        final_df_melted = final_df.melt(
            id_vars=["dataset", "cluster"],
            value_vars=ctype_lst,
            var_name="cell_type",
            value_name="score",
        )

        # Ensuring the melted DataFrame looks as expected
        return final_df_melted

    def annotation(self):
        """
        Annotate the :attr:`~adata` and visualize the result
        """
        # Perform statistical test to annotate cell clusters automatically
        self.render.render_text(
            "Perform statistical test to annotate cell clusters automatically", 2
        )
        df = dc.rank_sources_groups(
            self.acts, groupby="leiden", reference="rest", method="t-test_overestim_var"
        )
        self.render.render_text(df, 3)
        n_ctypes = 3
        ctypes_dict = (
            df.groupby("group")
            .head(n_ctypes)
            .groupby("group")["names"]
            .apply(lambda x: list(x))
            .to_dict()
        )
        self.render.render_text(ctypes_dict, 3)
        annotation_dict = (
            df.groupby("group").head(1).set_index("group")["names"].to_dict()
        )
        self.render.render_text(annotation_dict, 3)
        self.adata.obs["cell_type"] = [
            annotation_dict[clust] for clust in self.adata.obs["leiden"]
        ]
        # Visualize final cell-type annotation result
        self.render.render_text("Visualize final cell-type annotation result", 2)
        umap_plot = sc.pl.umap(
            self.adata,
            color="cell_type",
            title="decoupleR cell annotation",
            frameon=True,
            legend_fontweight="normal",
            legend_fontsize=10,
            return_fig=True,
        )
        self.render.render_fig(umap_plot, 1)

    def write_ann_ds(self, output: PathLike):
        """
        Write the annotated :attr:`~adata` to an HDF 5 (h5ad) file.
        """
        self.adata.write_h5ad(output, compression=hdf5plugin.FILTERS["zstd"])

    def cluster_vln_plot(self, melted_df):
        """
        Plot violin plot for clusters in provided `melted_df`

        Copied from https://github.com/kostaslazaros/cell_annotation_web_app/blob/main/decoupler_cell_annotation.ipynb
        Parameters
        ----------
        melted_df: pd.DataFrame
           The melted dataframe
        """
        # Plotting
        plt.figure(figsize=(16, 6))
        ax = sns.violinplot(
            x="cluster",
            y="score",
            hue="cell_type",
            data=melted_df,
            inner="quartile",
            palette="muted",
            scale="count",
        )

        # Ensuring gridlines are below plot elements
        ax.set_axisbelow(True)

        # Adding gridlines
        ax.yaxis.grid(
            True, color="gainsboro", linestyle="dashed"
        )  # Horizontal gridlines

        # Compute midpoints between x-tick labels for grid placement
        xticks = ax.get_xticks()
        apoints = xticks[:-1] + xticks[1:]
        midpoints = [i / 2 for i in apoints]

        # Add vertical gridlines at midpoints
        ax.vlines(
            midpoints,
            ymin=ax.get_ylim()[0],
            ymax=ax.get_ylim()[1],
            colors="lightgray",
            linestyles="solid",
            linewidth=0.5,
        )

        # Customizing other plot elements
        # ax.set_facecolor("lightgray")
        plt.gca().set_facecolor("whitesmoke")

        plt.title(
            f'Violin plots of cell-type ORA scores across leiden clusters (Data: {self.config.defaults["data"]}, Leiden_res = {self.config.defaults["cluster_resolution"]})'
        )
        plt.ylabel("ORA score")
        plt.xlabel("Leiden cluster")

        plt.legend(title="Cell Type", bbox_to_anchor=(1.01, 1), loc="upper left")

        plt.tight_layout()
        return plt

    @staticmethod
    def filtered_marker(cell_types: list[CellType]) -> pd.DataFrame:
        return pd.concat([c.get_selected() for c in cell_types], ignore_index=True)
