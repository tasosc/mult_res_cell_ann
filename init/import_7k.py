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
import logging
import urllib.parse
import json

logger=logging.getLogger(__file__)

def build_filename(cache_dir, tissue_name: str):
    return os.path.join(cache_dir, f'Human_{tissue_name}')

def build_url(tissue_name: str):
    return urllib.parse.quote(f'/ACT/download/Human/Human_{tissue_name}.txt')

def extract_tissue_name(row: list[str])-> str:
    return row[1]

def download(http_conn: http.client.HTTPConnection, tissue_name: str, cache_file):
    url = build_url(tissue_name)
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

def init_7k(base_cache_dir):
    cache_dir=os.path.join(base_cache_dir, 'act')
    os.makedirs(cache_dir, exist_ok=True)
    tissue_cache=download_if_not_cached(cache_dir)
    return build_import_object(tissue_cache)

def build_import_object(tissue_cache):
    repo=list()
    for name, file_path in tissue_cache:
        with open(file_path, newline='') as file:
            reader=csv.reader(file, delimiter='\t')
            next(reader)
            cells=list()
            for row in reader:
                # Species	Tissue	CellType	Marker	Resource
                cell={"cell_type": row[2], "genes":[g.strip() for g in row[3].split(',')]}
                cells.append(cell)
            repo.append({"tissue": name, "repo": "7k", "cells":cells})
    return repo

def download_if_not_cached(cache_dir):
    tissue_cache=list()
    with closing(http.client.HTTPConnection("xteam.xbio.top")) as http_conn:
        # 7k.txt is a file containing the organs 
        with open(os.path.join('init','7k.txt'), mode='r', newline='') as csvfile:
            reader=csv.reader(csvfile, delimiter='\t')
            next(reader)
            for row in reader:
                tissue_name=extract_tissue_name(row)
                file_path=build_filename(cache_dir, tissue_name)
                if not os.path.exists(file_path):
                    download(http_conn, tissue_name, file_path)
                if os.path.exists(file_path) and os.stat(file_path).st_size > 0:
                    tissue_cache.append((tissue_name, file_path))
    return tissue_cache

if __name__ == '__main__':
    print(json.dumps(init_7k('cache')))
