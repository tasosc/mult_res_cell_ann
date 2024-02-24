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
import sqlite3
from contextlib import closing
import logging
import json
import os
import sys

from import_7k import init_7k
from import_pangalaodb import init_pangalaodb
from import_cellmarker20 import init_cellmarker20

base_folder = 'init'
base_cache_dir = os.path.join('init', 'cache')
logger = logging.getLogger("init_db")


def init_db():
    logger.info("Starting initialization")
    conn = sqlite3.connect("data/cells.db")
    logger.info("connected")
    with open(os.path.join(base_folder, 'init_sqlite.sql')) as sql_file:
        sql_script = sql_file.read()
    conn.executescript(sql_script)
    logger.info("executed script")

    with open(os.path.join(base_folder, 'repos.json')) as repos_file:
        repos = json.load(repos_file)

    logger.info("read repos from json %s", repos)
    with conn:
        conn.executemany("insert into repos (repo_id, repo_name) values (:id, :name)", repos)

    tissue_ids = dict()
    import_data = list()
    logger.info("parsing PanglaoDB")
    import_data.extend(init_pangalaodb(base_cache_dir))
    logger.info("parsing Cell Marker 2.0")
    import_data.extend(init_cellmarker20(base_cache_dir))
    logger.info("parsing 7k")
    import_data.extend(init_7k(base_cache_dir))
    logger.info("Importing data")
    import_to_db(conn, import_data, repos, tissue_ids)


def import_to_db(conn: sqlite3.Connection, import_data: list, repos, tissue_ids: dict):
    repo_ids = dict([(r['name'], r['id']) for r in repos])
    with conn:  # start transaction
        for row in import_data:
            tissue, repo, cells = row["tissue"], row["repo"], row["cells"]
            logger.info("Importing %s %s", repo, tissue)
            if tissue not in tissue_ids:
                cur = conn.execute("insert into tissues (tissue_name) values (?)", (tissue,))
                tissue_ids[tissue] = cur.lastrowid
            tissue_id = tissue_ids[tissue]
            repo_id = repo_ids[repo]
            conn.execute("insert into cells (tissue_id, repo_id, cell_data) values (?,?,?)",
                         (tissue_id, repo_id, json.dumps(cells)))


if __name__ == '__main__':
    logging.basicConfig(stream=sys.stdout, level=logging.INFO, datefmt="%Y-%m-%d'T'%H:%M:%S")
    init_db()
