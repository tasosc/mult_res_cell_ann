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


class Store:
    logger = logging.getLogger("Store")

    def __init__(self) -> None:
        self.con = sqlite3.connect("file:data/cells.db?mode=ro", uri=True)

    def __query(self, sqlquery: str, parameters=()):
        self.logger.debug("Running %s", sqlquery)
        cur: sqlite3.Cursor
        with closing(self.con.cursor()) as cur:
            res = cur.execute(sqlquery, parameters)
            self.logger.debug("done running %s", sqlquery)
            return res.fetchall()

    def __query_one_column(self, sqlquery: str, parameters=()):
        return [res[0] for res in self.__query(sqlquery, parameters)]

    def get_tissue_types(self):
        return self.__query_one_column("select tissue_name from tissues order by tissue_name")

    def cell_type_of(self, tissue_type: str, repos: set[str]):
        self.logger.debug(repos)
        parameters = ", ".join("?" for _ in repos)
        return self.__query_one_column(
            f'select cell_data from cells c where exists (select 1 from tissues t where tissue_name=? and c.tissue_id = t.rowid) and exists (select 1 from repos r where repo_name in ({parameters}) and c.repo_id = r.repo_id)',
            (tissue_type,) + tuple(repos))

    def get_repos_for_tissue(self, tissue_type: str):
        return self.__query_one_column('''
                    select r.repo_name 
                    from repos r 
                    where exists 
                    (select 1 
                        from cells c
                        inner join tissues t on c.tissue_id = t.rowid 
                        where t.tissue_name = ? and r.repo_id  = c.repo_id)
            ''', (tissue_type,))
