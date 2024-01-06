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
from enum import Enum
import streamlit as st

class Buttons(Enum):
    Ensemble = 1
    Add = 2
    Run = 4
    

def button_clicked(button_name: Buttons):
    st.session_state.clicked[button_name.value] = True

def was_clicked(button_name: Buttons) -> bool:
    return (
        st.session_state.clicked
        and button_name.value in st.session_state.clicked
        and st.session_state.clicked[button_name.value]
    )

def clear_clicked(button_name: Buttons):
    st.session_state.clicked[button_name.value] = False

