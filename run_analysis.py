#
# This file is part of the mult_res_cell_ann distribution
# (https://github.com/tasosc/mult_res_cell_ann).
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
from pathlib import Path
from tempfile import NamedTemporaryFile

import streamlit as st

from cell_structure_id import StructureIdentification
from preprocessing import PreProcessing
from utilities import Config

logger = logging.getLogger("run_analysis")


def run_analysis(config: Config, uploaded_file, tissue_type: str, status):
    st.write("Scroll down to view more results")
    verbosity = config.defaults["verbosity"]

    def render_fig(fig, expected_verbosity=1):
        if verbosity < expected_verbosity:
            return
        st.pyplot(fig=fig)

    def render_text(something, expected_verbosity=1):
        if verbosity < expected_verbosity:
            return
        st.write(something)

    st.subheader("Preprocessing")
    render_text("Read & preprocess dataset (cell qc & filtering)", 2)
    render_text("Reading dataset", 2)
    status.update(label="Reading dataset", state="running", expanded=True)
    try:
        pp = PreProcessing.build_from(uploaded_file, config)
        pp.render.set_render_fig_lambda(render_fig)
        pp.render.set_render_text_lambda(render_text)
    except Exception as e:
        logger.info("Could not parse file. Reason %s", repr(e), exc_info=True, stack_info=True)
        status.update(label="Could not parse file", state="error", expanded=True)
        st.warning(
                f"Could not parse file. Please check if the file is not corrupted and supported. Reason: {str(e)}",
                icon="⚠️",
                )
        st.stop()

    if verbosity > 1:
        st.write("First names for cells and genes from the given dataset")
        col1, col2 = st.columns(2)
        with col1:
            st.write("Cells")
            st.write(pp.adata.obs_names[0])
        with col2:
            st.write("Genes")
            st.write(pp.adata.var_names[0])
    status.update(label="Running quality control", state="running", expanded=True)
    st.subheader("Quality Control")
    pp.qc()
    render_text(pp.adata, 3)
    status.update(label="Running normalization", state="running", expanded=True)
    st.subheader("Normalization")
    pp.normalization()

    status.update(label="Running feature selection", state="running", expanded=True)
    st.subheader("Feature selection")
    pp.feature_selection()
    st.subheader("Dimensionality reduction")
    status.update(
            label="Running Dimensionality reduction", state="running", expanded=True
            )
    pp.dimensionality_reduction()
    st.subheader("Visualize UMAP colored by leiden clusters")
    status.update(label="Running Visualization", state="running", expanded=True)
    pp.visualization()

    st.subheader("Structure identification")
    si = StructureIdentification(pp.adata, config)

    si.render.set_render_fig_lambda(render_fig)
    si.render.set_render_text_lambda(render_text)
    st.subheader("Clustering")
    render_text(f"Filter cell marker dataframe to obtain markers related to {tissue_type}", 2)
    try:
        si.clustering(st.session_state.cell_types)
    except ValueError as e:
        logger.info("Could not perform clustering Reason %s", repr(e))
        status.update(label="Could not perform clustering ", state="error", expanded=True)
        st.warning(
                f"Could not perform clustering. Please check if the tissue selected matches the dataset or try with different cells. Reason: {str(e)}",
                icon="⚠️",
                )
        st.stop()

    st.subheader("Annotation")
    si.annotation()

    st.subheader("Preparing annotated dataset")
    with NamedTemporaryFile() as tmp:
        # Write the annotated dataset to download_file
        tmp_path = Path(tmp.name)
        si.write_ann_ds(tmp_path)
        with open(tmp_path, "rb") as download_file:
            st.sidebar.download_button(
                    label="8️⃣Download annotated dataset and clear analysis",
                    data=download_file,
                    help="Download the filtered enriched with cell annotation scRNA-seq dataset in HDF5 format, ready to be used in downstream analysis",
                    file_name=f'{config.defaults["data"]}-annotated.h5ad',
                    mime="application/x-hdf5"
                    )
    status.update(label="Analysis completed.", state="complete", expanded=True)
    st.success("Automatic analysis successfully completed. You can optionally ⋮ -> Print and then Download the results.", icon="✅")
