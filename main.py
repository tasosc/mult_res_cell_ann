from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Annotated, List
from fastapi import FastAPI, UploadFile, Query
from fastapi.responses import FileResponse
from fastapi.middleware.cors import CORSMiddleware

from model.enums import SvdSolverOptions
from model.session import Cell, SessionData, SessionManager
from model.settings import Settings
from preprocessing import PreProcessing
from store import Store
from utilities import CellType
from cell_structure_id import StructureIdentification

app = FastAPI()
default_settings: Settings = Settings()

origins = [
    "http://localhost",
    "http://localhost:8080",
    "http://localhost:5173"
]
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

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

@app.post("/session")
def create_session(settings: Settings, cells: list[Cell]):
    """
    Create a session and return the session id
    """
    return {"session": SessionManager.create_session(cells=cells, settings=settings)}

@app.post("/analyze/{session}")
async def analyze_file(session: str, file: UploadFile):
    """
    analyse uploaded file use
    """

    # upload scRNAseq https://fastapi.tiangolo.com/reference/uploadfile/#fastapi.UploadFile
    session_data : SessionData = SessionManager.get_session(session)
    session_data.file = file
    def render_fig(fig, expected_verbosity=1):
        if session_data.settings.verbosity < expected_verbosity:
            return

    def render_text(something, expected_verbosity=1):
        if session_data.settings.verbosity < expected_verbosity:
            return
    pp = PreProcessing.build_from(file, session_data.settings)
    pp.render.set_render_fig_lambda(render_fig)
    pp.render.set_render_text_lambda(render_text)
    pp.qc()
    pp.normalization()
    pp.feature_selection()
    pp.dimensionality_reduction()
    pp.visualization()
    si = StructureIdentification(pp.adata, session_data.settings)
    si.render.set_render_fig_lambda(render_fig)
    si.render.set_render_text_lambda(render_text)
    si.clustering(session_data.cells)
    si.annotation()
    with NamedTemporaryFile(delete=False) as tmp:
        # Write the annotated dataset to download_file
        tmp_path = Path(tmp.name)
        si.write_ann_ds(tmp_path)
        return FileResponse(tmp_path)
    # Save to PDF pages (report0
    # https://stackoverflow.com/questions/11328958/save-multiple-plots-in-a-single-pdf-file
    # and/or using websockets ? return json with fig & text
    # https://stackoverflow.com/questions/71936110/correct-way-of-connecting-websocket-events-to-update-my-react-component

    # https://stackoverflow.com/questions/73550398/how-to-download-a-large-file-using-fastapi

# TODO https://fastapi.tiangolo.com/tutorial/static-files/
