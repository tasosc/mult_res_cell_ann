#!/usr/bin/env python
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
import sys
import streamlit as st
import logging
from cell_selection import cell_selection
from common_handlers import Buttons, clear_clicked, was_clicked
from run_analysis import run_analysis
from selection_manager import restore_selection, save_current_selection, select_tissue_repo
from store import Store
from utilities import CellType, Config
from preprocessing_settings import show_pp_settings

st.set_page_config(
    page_title="Multi Resource Cell Annotation (with decoupleR)",
    page_icon="🦠",
    layout="wide",
    initial_sidebar_state="auto",
    menu_items={
        "Get Help": "https://github.com/tasosc/mult_res_cell_ann/blob/main/README.md",
        "Report a bug": "https://github.com/tasosc/mult_res_cell_ann/issues",
        "About": """# This is a web app that can be used to automatically annotate single-cell RNA-seq dataset 
        by cell-type using prior knowledge databases and the decoupleR's python implementation. 
        It was based on https://github.com/kostaslazaros/cell_annotation_web_app""",
    },
)

logging.basicConfig(stream=sys.stdout, level=logging.INFO)

logger = logging.getLogger("main_web_app")
logger.info("Starting")
st.header("Multi Resource Cell Annotation (with decoupleR)")
# initialize


if "clicked" not in st.session_state:
    st.session_state.clicked = {Buttons.Ensemble.value: False, Buttons.Add.value: False}

if "repos" not in st.session_state:
    st.session_state.repos = set()


def clear(key: str):
    if key in st.session_state:
        st.session_state[key] = ""


if was_clicked(Buttons.Add):
    clear("new_cell_type")
    clear("new_cell_genes")
    clear_clicked(Buttons.Add)


# TODO this is left from old code
def on_upload():
    # user has not set yet the configuration settings we cannot build pre-processing object
    if "preprocessing" in st.session_state:
        del st.session_state["preprocessing"]
    if "structure_id" in st.session_state:
        del st.session_state["structure_id"]


if "config" not in st.session_state:
    st.session_state.config = Config()

store: Store = Store()
config: Config = st.session_state.config

with st.sidebar:
    H5AD_UPLOAD = st.file_uploader(
        "1️⃣Select the file containing the single cell RNA sequencing dataset.",
        help="The application supports either .h5ad or .csv files. The csv files can be compressed with gzip.By default it will use predefined Prostate test data",
        key="data_upload",
        on_change=on_upload,
        type=["h5ad", "csv", "csv.gz", "txt", "txt.gz"],
    )
    tissue_options = ["Database", "Saved selection"]
    tissue_option = st.radio(
        label="Select the source of tissues, cells and genes",
        options=tissue_options,
        captions=[
            "If not sure select this. You will need to choose the tissue/organ type, the cell types and their genes",
            "Select this if you wish to restore a previous selection of cell types and genes. You need to provide a .json file that was created by this application",
        ],
        index=0,
    )
    if tissue_option == tissue_options[0]:
        select_tissue_repo(store)

    elif tissue_option == tissue_options[1]:
        restore_selection()

    tissue_type = st.session_state.get("tissue_type", None)

    # Advanced settings
    if st.toggle("5️⃣Parameter Adjustment and advanced settings",
                 help="Parameter Adjustment: Customize general and pre-processing parameters for single-cell RNA sequencing analysis. These controls enable fine-tuning of various settings such as filtering criteria, normalization methods, dimensionality reduction techniques, and other essential parameters. Adjust these settings to tailor data processing according to specific analysis requirements."):
        show_pp_settings(config)

if was_clicked(Buttons.Ensemble):
    if "cell_types" not in st.session_state:
        st.session_state.cell_types = list(
            CellType.parse_json(
                store.cell_type_of(
                    tissue_type=tissue_type, repos=st.session_state.repos
                )
            )
        )

    if st.sidebar.button("7️⃣Run Analysis",
                         help="Run Automated Cell Type Annotation: Activate automated analysis to enrich the uploaded single-cell RNA sequencing dataset with cell type annotations. Prior selection of specific cell types and their expressing genes is required. This analysis augments the dataset with refined cell type information, leveraging chosen markers to classify and annotate cells, enhancing the understanding of cellular heterogeneity and functional characterization"):
        can_run = any(map(lambda x: x.has_selected(), st.session_state.cell_types))
        if not can_run:
            st.warning("Select some cells and genes first", icon="🚨")
        else:
            with st.status("Running analysis..", expanded=True) as status:
                run_analysis(config, H5AD_UPLOAD, tissue_type, status=status)

    cell_selection()

    save_current_selection()
