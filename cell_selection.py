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
from common_handlers import Buttons, button_clicked 
from utilities import CellType



import streamlit as st

logger = logging.getLogger("run_analysis")


def add_clicked(new_cell_type: str):
    logger.info("Adding %s", new_cell_type)
    c = CellType(cell_name=new_cell_type)
    st.session_state.cell_types.append(c)
    button_clicked(Buttons.Add)


def genes_input_key(cell: CellType):
    return f"text_input_{cell}"

def add_gene(cell: CellType):
    key = genes_input_key(cell)
    new_genes = st.session_state[key] if key in st.session_state else None
    if new_genes and len(new_genes) > 0:
        as_list = new_genes.split(" ")
        logger.info(f"Adding genes: {as_list}")
        cell.set_new_genes(as_list)
    else:
        cell.set_new_genes(None)

def set_default_all_gene(cells: list[CellType], value: bool):
    for cell in cells:
        set_is_selected_value_of(cell, value)


def get_is_selected_value_of(cell: CellType) -> bool:
    key = cell_all_gene_toggle_key(cell=cell)
    return key in st.session_state and st.session_state[key]


def set_is_selected_value_of(cell: CellType, value: bool):
    st.session_state[cell_all_gene_toggle_key(cell=cell)] = value


def on_cell_select_all(cells: list[CellType], value: bool):
    set_default_all_gene(cells, value)


def cell_all_gene_toggle_key(cell: CellType):
    return f"all_gene_toggle_{cell.cell_type}"

def cell_selection():
    # Check if we have selected a tissue as we need it
    tissue_type = st.session_state.get("tissue_type", None)
    if not tissue_type:
        return
    st.info(f"Select cell types and genes for tissue type {tissue_type}", icon="6️⃣")
    
    # Create the search box
    search = st.text_input(
        "Search cells or add new",
        placeholder="type a cell name, for searching it can be partial, and press Enter",
        autocomplete="on",
    )

    # Get the cell types optionally filtered by the text in the searchbox
    cell_types: list[CellType] = [
        c for c in st.session_state.cell_types if not (search) or search in c.cell_type
    ]

    # Organize buttons under search box in three columns
    col1, col2, col3 = st.columns(3)
    
    # Select all cell types and all their genes button
    col1.button("Select all", help="Select all visible cell types and their genes" , on_click=on_cell_select_all, args=[cell_types, True])
    # Unselect all cell types and all their genes button
    col2.button("Deselect all", help="Unselect all visible cell types and their genes", on_click=on_cell_select_all, args=[cell_types, False])
    # Add what it was written in the search box as a new cell type
    col3.button("Add as new cell type", help="Add the typed text in the search box as new cell type", disabled=not search, on_click=add_clicked, args=[search])

    # Render the cell type and their genes in a box for each cell type
    # The columns here is to try to put two cell type boxes in one line
    columns = st.columns(2)
    for index, cell_type in enumerate(cell_types):
        with columns[index % 2]:
            logger.debug(cell_type.genes)
            # for each cell type create a box
            with st.container(border=True):
                st.write(f"{cell_type.cell_type}")
                cell_type.is_selected = st.toggle(
                    f"Select all genes for {cell_type.cell_type}",
                    value=get_is_selected_value_of(cell_type),
                )
                set_is_selected_value_of(cell_type, cell_type.is_selected)
                cell_type.gene_selection = st.multiselect(
                    "Select specific genes",
                    options=cell_type.genes,
                    placeholder=f"select genes for [{cell_type.cell_type}]",
                    disabled=get_is_selected_value_of(cell_type),
                )
                # text box to add new gene symbols. They are separated by space. 
                st.text_input(
                    label="Add new genes",
                    placeholder="Type new gene symbols separated by space and press Enter",
                    key=genes_input_key(cell_type),
                    on_change=add_gene,
                    args=[cell_type],
                )