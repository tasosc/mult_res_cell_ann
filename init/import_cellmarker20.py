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
import os
import http.client
import shutil
from contextlib import closing
import logging
import sys
import urllib.parse
import json
import pandas as pd

logger = logging.getLogger(__file__)

"""
Download from http://117.50.127.228/CellMarker/CellMarker_download_files/file/Cell_marker_Human.xlsx
"""


def build_filename(cache_dir):
    return os.path.join(cache_dir, 'Cell_marker_Human.xlsx')


def build_url():
    return urllib.parse.quote('/CellMarker/CellMarker_download_files/file/Cell_marker_Human.xlsx')


def build_import_object(file_path):
    repo = list()
    logger.info("Reading excel file %s", file_path)
    df = pd.read_excel(file_path, dtype=str, usecols=[1, 6, 9])
    logger.info("Done reading excel file %s", file_path)
    df = df.dropna()
    #    df = df[['tissue_type', 'cell_name', 'Symbol']]
    tissue_dict = dict()
    logger.info("Start parsing excel in memory")
    for row in df.itertuples():
        tissue = row[1]
        cell_type = row[2]
        gene = row[3]
        # include entries only when there is a gene
        if pd.isna(gene):
            raise ValueError("Gene is NA")
        if tissue not in tissue_dict:
            tissue_dict[tissue] = dict()
        cell_dict = tissue_dict[tissue]
        if cell_type not in cell_dict:
            cell_dict[cell_type] = list()
        cell_dict[cell_type].append(gene)
    logger.info("Second pass")
    for tissue, cell_dict in tissue_dict.items():
        cells = list()
        for cell_type, genes in cell_dict.items():
            cell = {"cell_type": cell_type, "genes": genes}
            cells.append(cell)
        repo.append({"tissue": tissue, "repo": "CellMarker 2.0", "cells": cells})
    return repo


def download(http_conn: http.client.HTTPConnection, cache_file):
    url = build_url()
    logger.info("Downloading %s", url)
    http_conn.request("GET", url)
    with http_conn.getresponse() as response:
        if response.status == http.HTTPStatus.OK:
            with open(cache_file, mode='wb') as file:
                shutil.copyfileobj(response, file)
                file.flush()
                logger.info("Downloaded %s to %s", url, cache_file)
        else:
            logger.error("Cannot GET %s , server responded with %d reason %s", url, response.status, response.reason)
            raise ImportError("couldn't download")


def download_if_not_cached(cache_dir):
    with closing(http.client.HTTPConnection("117.50.127.228")) as http_conn:
        file_path = build_filename(cache_dir)
        if not os.path.exists(file_path):
            download(http_conn, file_path)
        if os.path.exists(file_path) and os.stat(file_path).st_size > 0:
            return file_path
    raise IOError("Cannot download data")


def init_cellmarker20(base_cache_dir):
    cache_dir = os.path.join(base_cache_dir, 'cellmarker')
    os.makedirs(cache_dir, exist_ok=True)
    cache_file = download_if_not_cached(cache_dir)
    return build_import_object(cache_file)


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stdout, level=logging.INFO)
    print(json.dumps(init_cellmarker20('cache'), indent=2))
