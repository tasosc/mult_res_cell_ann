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
import csv
import os
import http.client
import shutil
from contextlib import closing
import pandas as pd
import decoupler as dc
import logging
import urllib.parse
import json
import gzip

logger = logging.getLogger(__file__)

"""
Download from https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz
"""


def build_filename(cache_dir):
    return os.path.join(cache_dir, 'PanglaoDB_markers_27_Mar_2020.tsv.gz')


def build_url():
    return urllib.parse.quote('/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz')


def build_import_object(file_path):
    repo = list()
    with gzip.open(file_path, mode="rt", newline='') as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader)
        tissue_dict = dict()
        for row in reader:
            # species	official gene symbol	cell type	nicknames	ubiquitousness index	product description	gene type	canonical marker	germ layer	organ	sensitivity_human	sensitivity_mouse	specificity_human	specificity_mouse
            if 'Hs' in row[0]:  # Only Human
                cell_type = row[2]
                gene = row[1]
                tissue = row[9]
                if tissue not in tissue_dict:
                    tissue_dict[tissue] = dict()
                cell_dict = tissue_dict[tissue]
                if cell_type not in cell_dict:
                    cell_dict[cell_type] = list()
                cell_dict[cell_type].append(gene)
    return build_common_model(repo, tissue_dict)


def build_import_object_omni():
    repo = list()
    df: pd.DataFrame = dc.get_resource("PanglaoDB")
    df = df[df["human"] == 'True']
    df = df[["genesymbol", "cell_type", "organ"]]
    df = df.dropna()
    tissue_dict = dict()
    for row in df.to_numpy():
        cell_type = row[1]
        gene = row[0]
        tissue = row[2]
        if tissue not in tissue_dict:
            tissue_dict[tissue] = dict()
        cell_dict = tissue_dict[tissue]
        if cell_type not in cell_dict:
            cell_dict[cell_type] = list()
        cell_dict[cell_type].append(gene)
    return build_common_model(repo, tissue_dict)


def build_common_model(repo, tissue_dict):
    for tissue, cell_dict in tissue_dict.items():
        cells = list()
        for cell_type, genes in cell_dict.items():
            cell = {"cell_type": cell_type, "genes": genes}
            cells.append(cell)
        repo.append({"tissue": tissue, "repo": "PanglaoDB", "cells": cells})
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
    with closing(http.client.HTTPSConnection("panglaodb.se")) as http_conn:
        file_path = build_filename(cache_dir)
        if not os.path.exists(file_path):
            download(http_conn, file_path)
        if os.path.exists(file_path) and os.stat(file_path).st_size > 0:
            return file_path
    raise IOError("Cannot download data")


def init_pangalaodb(base_cache_dir):
    try:
        return build_import_object_omni()
    except Exception as e:
        logger.error("While trying to get data from omnipath %s", e, exc_info=1, stack_info=True)
        cache_dir = os.path.join(base_cache_dir, 'pangalaodb')
        os.makedirs(cache_dir, exist_ok=True)
        cache_file = download_if_not_cached(cache_dir)
        return build_import_object(cache_file)


if __name__ == '__main__':
    print(json.dumps(init_pangalaodb('cache')))
