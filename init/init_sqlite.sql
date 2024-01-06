begin;
-- clean up db
drop table if exists cells;
drop table if exists repos;
DROP table if exists tissues;

-- Table to store the tissues
create table tissues ( 
    tissue_id INTEGER PRIMARY KEY,
    tissue_name TEXT NOT NULL,
    UNIQUE (tissue_name)
);
create table repos ( 
    repo_id INTEGER PRIMARY KEY,
    repo_name TEXT NOT NULL,
    UNIQUE (repo_id)
);
create table cells ( 
    tissue_id INTEGER NOT NULL,
    repo_id INTEGER NOT NULL, 
    cell_data TEXT NOT NULL, 
    PRIMARY KEY (tissue_id, repo_id), 
    FOREIGN KEY(tissue_id) REFERENCES tissues(tissue_id),
    FOREIGN KEY(repo_id) REFERENCES repos(repo_id)
);
commit;