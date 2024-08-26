from typing import List
from fastapi import FastAPI, UploadFile, Query

from model.enums import SvdSolverOptions
from model.settings import Settings
from store import Store
from utilities import CellType

app = FastAPI()
default_settings: Settings = Settings()


@app.get("/metadata/settings/defaults")
def read_default_settings():
    return default_settings


@app.get("/metadata/svd_solver_options")
def read_svd_solver_options():
    return SvdSolverOptions.values()


@app.get("/tissues/{tissue}/sources")
def read_sources(tissue: str):
    store: Store = Store()
    return store.get_repos_for_tissue(tissue)


@app.get("/tissues")
def read_tissues():
    store: Store = Store()
    return store.get_tissue_types()


@app.get("/tissues/{tissue}/cells")
def read_cells(tissue: str, sources: List[str] = Query([])):
    store: Store = Store()
    return list(
    CellType.parse_json(store.cell_type_of(tissue, set(sources))))


    # TODO upload scRNAseq https://fastapi.tiangolo.com/reference/uploadfile/#fastapi.UploadFile
@app.post("/session/")
def analyze_file(file: UploadFile):
    return {"filename": file.filename}

# TODO https://fastapi.tiangolo.com/tutorial/static-files/
