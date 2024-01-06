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
 
import streamlit as st


def on_config_change(option: str):
    key = f"config_{option}"
    if key in st.session_state and "config" in st.session_state:
        st.session_state.config.defaults[option] = st.session_state[key]


def show_pp_settings(config):
    if not config:
        st.error("Cannot view settings as configuration was not loaded")
    if st.button("Reset defaults"):
        config.reload_from_disk()
    st.markdown("### General\n")
    st.slider(
        "Verbosity of analysis.",
        help="0 - minimum output, 1 - more text and figures, 2- more text,  3- include also tables (can consume a lot of resources on your PC)",
        min_value=0,
        max_value=3,
        on_change=on_config_change,
        key="config_verbosity",
        value=config.defaults["verbosity"],
        args=["verbosity"],
    )
    st.markdown("### Read dataset\n")
    st.text_input(
        "Reading CSV, field delimiter",
        key="config_csv_delimiter",
        max_chars=1,
        on_change=on_config_change,
        args=["csv_delimiter"],
        value=config.defaults["csv_delimiter"],
    )
    st.markdown("### Quality Control\n")
    st.number_input(
        "Lower cut-off of gene number for filtering cells",
        key="config_n_genes_min",
        min_value=1,
        on_change=on_config_change,
        args=["n_genes_min"],
        value=config.defaults["n_genes_min"],
    )
    st.number_input(
        "Upper cut-off of gene number for filtering cells",
        key="config_n_genes_max",
        min_value=1,
        on_change=on_config_change,
        args=["n_genes_max"],
        value=config.defaults["n_genes_max"],
    )
    st.number_input(
        "Minimum genes that should be expressed in each cell",
        key="config_min_genes",
        min_value=1,
        on_change=on_config_change,
        args=["min_genes"],
        value=config.defaults["min_genes"],
    )
    st.number_input(
        "Minimum number of cells expected to express each gene",
        key="config_min_cells",
        min_value=1,
        on_change=on_config_change,
        args=["min_cells"],
        value=config.defaults["min_cells"],
    )
    st.number_input(
        "Upper cut-off of counts for filtering cells",
        key="config_n_counts_max",
        min_value=1,
        on_change=on_config_change,
        args=["n_counts_max"],
        value=config.defaults["n_counts_max"],
    )
    st.number_input(
        "Percentage of mitochondrial gene expression (cut-off)",
        key="config_pc_mito",
        min_value=1,
        on_change=on_config_change,
        args=["pc_mito"],
        value=config.defaults["pc_mito"],
    )
    st.number_input(
        "Percentage of ribosomal gene expression (cut-off)",
        key="config_pc_rib",
        min_value=1,
        on_change=on_config_change,
        args=["pc_rib"],
        value=config.defaults["pc_rib"],
    )
    st.markdown("### Normalization\n")
    st.toggle(
        "Normalize total counts",
        key="config_normalize_total_counts",
        on_change=on_config_change,
        args=["normalize_total_counts"],
        value=config.defaults["normalize_total_counts"],
    )
    st.markdown("### Feature Selection\n")
    st.toggle(
        "Only highy variable genes",
        key="config_only_highly_significant_genes",
        on_change=on_config_change,
        args=["only_highly_significant_genes"],
        value=config.defaults["only_highly_significant_genes"],
    )
    st.markdown("### Dimensionality reduction\n")
    svd_solver_options: list[str] = config.defaults["svd_solver_options"]
    st.selectbox(
        label="SVD Solver",
        options=svd_solver_options,
        key="config_svd_solver",
        on_change=on_config_change,
        args=["svd_solver"],
        index=svd_solver_options.index(config.defaults["svd_solver"]),
    )
    st.markdown("### Visualization\n")
    st.number_input(
        "Number of Neighbours",
        key="config_n_neigh",
        help="See n_neighbors in https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.neighbors.html#n_neighbors",
        min_value=2,
        max_value=100,
        on_change=on_config_change,
        args=["n_neigh"],
        value=config.defaults["n_neigh"],
    )
    st.number_input(
        "Number of principal components (PCs)",
        help="See n_pcs in https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.neighbors.html#n_pcs",
        key="config_n_pcs",
        min_value=0,
        on_change=on_config_change,
        args=["n_pcs"],
        value=config.defaults["n_pcs"],
    )

    st.number_input(
        "Cluster resolution. Higher values lead to more clusters",
        key="config_cluster_resolution",
        min_value=0.1,
        max_value=1.0,
        step=0.1,
        on_change=on_config_change,
        args=["cluster_resolution"],
        value=config.defaults["cluster_resolution"],
    )
