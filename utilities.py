#!/usr/bin/env -- python3
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
import re
import json
import logging

class CellType:
    logger = logging.getLogger("CellType")
    cell_combi = re.compile(r"cells?$", flags=re.I)

    def __init__(self, cell_name: str) -> None:
        self.cell_type = cell_name
        self.is_selected = False
        self.genes = set()
        self.gene_selection = None
        self.new_genes = None

    def add_genes(self, genes: list[str]):
        self.genes.update(genes)
        self.logger.debug("Added genes %s to %s", genes, self.genes)

    @classmethod
    def parse_json(cls, tissue_cells: list[str]):
        cells_dict = dict()
        for json_cells in tissue_cells:
            j = json.loads(json_cells)
            for c in j:
                key = c["cell_type"]
                key = cls.cell_combi.sub("cell", key)
                if key not in cells_dict:
                    cells_dict[key] = CellType(cell_name=key)
                if "genes" in c:
                    cells_dict[key].add_genes(c["genes"])
        return sorted(cells_dict.values(), key=lambda k: k.cell_type)

