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
import base64
import io
import pandas as pd
import json
import logging
from datetime import datetime
class CellType:
    logger = logging.getLogger("CellType")

    def __init__(self, json_cell=None, cell_name: str = None) -> None:
        j = json_cell
        self.cell_type = j["cell_type"] if j and "cell_type" in j else cell_name
        self.is_selected = False
        self.genes = set(j["genes"] if j and "genes" in j else [])
        self.gene_selection = None
        self.new_genes = None

    def add_genes(self, genes: list[str]):
        self.genes.update(genes)
        self.logger.info(f"Added genes {genes} to {self.genes}")

    def set_new_genes(self, genes: list[str]):
        if genes:
            self.new_genes = set(genes)
        else:
            self.new_genes = None

    def get_selected(self) -> pd.DataFrame:
        selection = set()
        if self.is_selected:
            selection.update(self.genes)
        elif self.gene_selection:
            selection.update(self.gene_selection)
        if self.new_genes:
            selection.update(self.new_genes)
        if len(selection) == 0:
            return pd.DataFrame({"cell_name": [], "Symbol": []})
        selection.discard(None)
        df = pd.DataFrame(data=selection, columns=["Symbol"])
        df.insert(0, "cell_name", self.cell_type)
        return df
    def has_selected(self) -> bool:
        selection = set()
        if self.is_selected:
            selection.update(self.genes)
        elif self.gene_selection:
            selection.update(self.gene_selection)
        if self.new_genes:
            selection.update(self.new_genes)
        selection.discard(None)
        return len(selection) > 0 

    @classmethod
    def parse_json(cls, tissue_cells: list[str]):
        cells_dict = dict()
        for json_cells in tissue_cells:
            j = json.loads(json_cells)
            for c in j:
                key = c["cell_type"]
                if key not in cells_dict:
                    cells_dict[key] = CellType(json_cell=c)
                elif "genes" in c:
                    cells_dict[key].add_genes(c["genes"])
        return cells_dict.values()

class Config:
    def __init__(self) -> None:
        self.defaults = self.__config_load()
    
    def update(self, new_defaults: dict|None=None) -> None:
        if new_defaults and len(new_defaults) > 0:
            # first sanitze input to avoid including config values not
            d = {key: value for key, value in new_defaults.items() if key in self.defaults.keys()}
            if "h5ad_path" in d:
                del d["h5ad_path"]
            self.defaults.update(d)

    def reload_from_disk(self):
        self.defaults = self.__config_load()

    def __config_load(self) -> dict:
        with open("config.json") as f:
            return json.load(f)
    
    def get_bool(self, key):
        return key in self.defaults and self.defaults[key]
    
class ImportExport:
    logger = logging.getLogger("ImportExport")
    @classmethod
    def export(cls, config : Config, tissue: str, cells : list[CellType]|None) -> str:
        cells = cells if cells else []
        cls.logger.info("Cells given with selected value %s", any(map(lambda x: x.has_selected(), cells)))
        cells = [c.__dict__ for c in cells]
        # sanitize config.defaults
        c=dict(config.defaults)
        if "h5ad_path" in c:
            del c["h5ad_path"]
        obj={"config":c, "tissue": tissue, "cells": cells}
        return json.dumps(obj=obj, default=list, indent=4)
#        with io.BytesIO() as buf:
#            json.dump(obj=obj,fp=buf)
#            b64=base64.b64encode(buf)
#            return f'<a href="data:application/json;base64,{b64.decode()}" download={file}>Export selections</a>'
    @classmethod
    def export_link(cls, config : Config, tissue: str, cells : list[CellType]|None) -> str:
        cells = cells if cells else []
        cls.logger.info("Cells given with selected value %s", any(map(lambda x: x.has_selected(), cells)))
        cells = [c.__dict__ for c in cells]
        # sanitize config.defaults
        c=dict(config.defaults)
        if "h5ad_path" in c:
            del c["h5ad_path"]
        obj={"config":c, "tissue": tissue, "cells": cells}
        file=f'last-run-{datetime.now().isoformat().replace(":","_")}.json'
        js = json.dumps(obj=obj, default=list, indent=4)
        b64=base64.b64encode(js.encode())
        return f'<a href="data:application/json;base64,{b64.decode()}" download="{file}">Save current selection</a>'

    @classmethod
    def import_selections(cls ,uploaded_file)-> tuple[dict, str, list[CellType]]:
        obj=json.load(uploaded_file)
        cls.logger.info(obj)
        config = obj.get('config', None)
        
        tissue = obj.get('tissue', None )
        cells_json = obj.get('cells', {})
        
        cells=list()
        for cell_dict in cells_json:
            if 'cell_type' in cell_dict:
                cell = CellType(cell_name=cell_dict['cell_type'])
                cells.append(cell)
                cell.is_selected = cell_dict.get('is_selected', False)
                cell.add_genes(cell_dict.get('genes', []))
                cell.gene_selection = cell_dict.get('gene_selection', [])
                cell.new_genes = cell_dict.get('new_genes', None)
        return (config, tissue, cells)

                
            

        

        
        




