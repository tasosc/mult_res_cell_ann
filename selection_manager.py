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
import streamlit as st
from cell_selection import set_is_selected_value_of
from common_handlers import Buttons, button_clicked, clear_clicked

from datetime import datetime
from utilities import ImportExport


logger = logging.getLogger("CellType")

def ensemble_clicked():
    clear_cell_types()
    button_clicked(Buttons.Ensemble)

def clear_cell_types():
    if "cell_types" in st.session_state:
        del st.session_state["cell_types"]



def tissue_changed():
    clear_cell_types()
    clear_clicked(Buttons.Ensemble)
    st.session_state.repos = set()
    tissue_that_was_selected = st.session_state.get("tissue_selection", None)
    logger.info("Tissue selection to %s", tissue_that_was_selected)
    st.session_state.tissue_type = tissue_that_was_selected

def on_repo_checked(repo):
    logger.info("repo %s, st.session_state %s", repo, st.session_state.get(repo))
    if st.session_state.get(repo):
        st.session_state.repos.add(repo)
    else:
        st.session_state.repos.discard(repo)

def import_settings(key: str):
    if key not in st.session_state:
        logger.warning(f"Import settings clicked but {key} not in session ")
        return
    try:
        config = st.session_state.config
        nconfig, ntissue, ncells = ImportExport.import_selections(st.session_state[key])
        config.update(nconfig)
        if ntissue:
            tissue_changed()
            st.session_state.tissue_type = ntissue
            button_clicked(Buttons.Ensemble)
            st.session_state.cell_types = ncells
            for c in ncells:
                set_is_selected_value_of(c, c.is_selected)

    except Exception as e:
        logger.error(
            "While parsing the selection settins, reason: %s",
            e,
            exc_info=1,
            stack_info=True,
            stacklevel=2,
        )
        st.warning(f"Could parse the settings file. Reason {e}", icon="🤖")


def restore_selection():
    st.file_uploader(
            "2️⃣ - 4️⃣Restore a previously selected tissue, preprocessing settings, cell types and genes",
            key="upload_settings",
            help="Restore Previous Selection: Click to revert to the previously saved tissue, cell type, and gene selections. This functionality enables you to quickly retrieve and apply previously stored configurations, facilitating seamless navigation back to earlier settings without manual reselection",
            on_change=import_settings,
            args=["upload_settings"],
            type=["json"],
        )

def select_tissue_repo(store):
    st.selectbox(
            "2️⃣Select tissue type related to the uploaded file:",
            store.get_tissue_types(),
            key="tissue_selection",
            help="Select the tissue or organ type associated with the uploaded dataset. Choose from a range of options representing different anatomical structures or body systems. This selection helps categorize and filter cell types based on specific tissue or organ characteristics, aiding in targeted analysis.",
            on_change=tissue_changed,
            index=None,
            placeholder="Select tissue type",
        )
    tissue_type = st.session_state.get("tissue_type", None)
    if tissue_type:
        st.write("3️⃣Select one or more resources")
        help_resource={
                "PangalaoDB":"""Oscar Franzén, Li-Ming Gan, Johan L M Björkegren, PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data, Database, Volume 2019, 2019, baz046, doi:10.1093/database/baz046""",
                "CellMarker 2.0":"""Congxue Hu, Tengyue Li, Yingqi Xu, Xinxin Zhang, Feng Li, Jing Bai, Jing Chen, Wenqi Jiang, Kaiyue Yang, Qi Ou, Xia Li, Peng Wang, Yunpeng Zhang, CellMarker 2.0: an updated database of manually curated cell markers in human/mouse and web tools based on scRNA-seq data, Nucleic Acids Research, Volume 51, Issue D1, 6 January 2023, Pages D870–D876, https://doi.org/10.1093/nar/gkac947""",
                "7k":"""Quan, F., Liang, X., Cheng, M. et al. Annotation of cell types (ACT): a convenient web server for cell type annotation. Genome Med 15, 91 (2023). https://doi.org/10.1186/s13073-023-01249-5"""
        }
        repos=store.get_repos_for_tissue(tissue_type)
        for r in repos:
            st.checkbox(r, key=r, on_change=on_repo_checked, help=help_resource[r], args=[r])
        st.button(
                "4️⃣Query selected resources",
                help="Query Cell Types: View cell types derived from chosen prior knowledge resources corresponding to the selected tissue or organ. Explore each cell type along with its associated set of genes. Click to access a comprehensive overview of identified cell types and their associated genetic signatures, aiding in comparative analysis and understanding of cellular compositions.",
                disabled=not (st.session_state.repos),
                on_click=ensemble_clicked,
            )
    else:
        st.write("3️⃣Select a tissue to be able to select resources")

def save_current_selection():
    tissue_type = st.session_state.get("tissue_type", None)
    config = st.session_state.get("config", None)
    cells = st.session_state.get("cell_types", None)
    if not tissue_type or not config:
        return
    st.sidebar.markdown(ImportExport.export_link(config=config, tissue=tissue_type, cells=cells), unsafe_allow_html=True, help="Save Current Selection: Click to preserve the current choices of tissue, cell type, and gene selections. This feature allows you to store and recall specific configurations for future reference or analysis, enabling quick access to preferred settings without the need for reselection")

