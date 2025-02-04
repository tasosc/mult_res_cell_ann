"""
REST API Controller for multiple resource cell annotation
"""
import datetime
import io
from pathlib import Path
import shutil
from tempfile import NamedTemporaryFile
from typing import List
from fastapi import FastAPI, UploadFile, Query, WebSocket, HTTPException
from fastapi.responses import FileResponse
from fastapi.middleware.cors import CORSMiddleware
from matplotlib.figure import Figure

from model.enums import Activity, SvdSolverOptions
from model.session import Cell, SessionData, SessionManager
from model.settings import Settings
from preprocessing import PreProcessing
from store import Store
from utilities import CellType
from cell_structure_id import StructureIdentification
from utils.feedback_socket import FeedbackSocket

app = FastAPI()
default_settings: Settings = Settings()
BUF_SIZE=128*1024
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

@app.post("/session", status_code=201)
def create_session(settings: Settings, cells: list[Cell]):
    """
    Create a session and return the session id
    """
    return {"session": SessionManager.create_session(cells=cells, settings=settings)}

@app.post("/dataset/{session}", status_code=201)
async def analyze_file(session: str, file: UploadFile):
    """
    analyse uploaded file use
    """
    if not session:
        raise HTTPException(status_code=422, detail="Session not specified")

    # upload scRNAseq https://fastapi.tiangolo.com/reference/uploadfile/#fastapi.UploadFile
    session_data : SessionData = SessionManager.get_session(session)

    if not session_data:
        raise HTTPException(status_code=404, detail="Session not found")
    filename = file.filename
    with NamedTemporaryFile(delete=False) as tmp:
        tmp_path = Path(tmp.name)
        shutil.copyfileobj(file.file, tmp)
        session_data.file = tmp_path

    file.close()
    # TODO background task for deleting session in X time

    return {
        "activity": Activity.UPLOAD_DATASET,
        "finished": datetime.datetime.now().isoformat(sep="T"),
        "duration": -1,
        "message": f"Uploaded scRNA-seq dataset {filename}",
        "link": None,
    }

@app.get("/annotated/{session}/{filename}")
async def get_annotated_dataset(session: str, filename: str):
    session_data : SessionData = SessionManager.get_session(session)
    if not session_data or not session_data.annotated or session_data.download_filename != filename:
        raise HTTPException(status_code=404, detail="Session or file not found")
    return FileResponse(path = session_data.annotated, media_type="application/octet-stream", filename=session_data.download_filename, content_disposition_type="attachment")

@app.websocket("/ws/{session}")
async def get_feedback(socket: WebSocket, session : str):
    session_data : SessionData = SessionManager.get_session(session)
    if not session_data:
        raise HTTPException(status_code=404, detail="Session not found")
    await socket.accept()
    await run_analysis(session_data=session_data, socket=socket)
    await socket.close()


async def run_analysis(session_data: SessionData, socket: WebSocket):
    feedback = FeedbackSocket(socket)
    async def render_fig(fig: Figure, expected_verbosity=1):
        if session_data.settings.verbosity < expected_verbosity:
            return
        with io.BytesIO() as buf:
            fig.savefig(buf, format='png')
            buf.seek(0)
            await socket.send_bytes(buf)

    async def render_text(something: str, current_activity: Activity = Activity.NONE, expected_verbosity=1):
        if session_data.settings.verbosity < expected_verbosity:
            return
        await feedback.completed(current_activity=current_activity, message=something)

    pp = PreProcessing.build_from(session_data.file, session_data.settings)
    pp.render.set_render_fig_lambda(render_fig)
    pp.render.set_render_text_lambda(render_text)
    feedback.start()
    pp.qc()
    await feedback.completed(Activity.PP_QC)
    pp.normalization()
    await feedback.completed(Activity.PP_NORM)
    pp.feature_selection()
    await feedback.completed(Activity.PP_FEATURE)
    pp.dimensionality_reduction()
    await feedback.completed(Activity.PP_REDUCTION)
    pp.visualization()
    await feedback.completed(Activity.PP_VISUALIAZTION)
    si = StructureIdentification(pp.adata, session_data.settings)
    si.render.set_render_fig_lambda(render_fig)
    si.render.set_render_text_lambda(render_text)
    si.clustering(session_data.cells)
    await feedback.completed(Activity.SI_CLUSTERING)
    si.annotation()
    await feedback.completed(Activity.SI_ANNOTATION)
    with NamedTemporaryFile(delete=False) as tmp:
        # Write the annotated dataset to download_file
        tmp_path = Path(tmp.name)
        si.write_ann_ds(tmp_path)
        session_data.annotated = tmp_path
        filename = session_data.file.filename if session_data.file is not None and session_data.file.filename is not None else "annotated"
        filename = filename + ".h5ad"
        session_data.download_filename = filename
        await feedback.completed(current_activity= Activity.SI_FILE, link=f"/{SessionManager.get_session_id(session_data)}/{filename}")
        # Save to PDF pages (report0
    # https://stackoverflow.com/questions/11328958/save-multiple-plots-in-a-single-pdf-file
    # and/or using websockets ? return json with fig & text
    # https://stackoverflow.com/questions/71936110/correct-way-of-connecting-websocket-events-to-update-my-react-component

    # https://stackoverflow.com/questions/73550398/how-to-download-a-large-file-using-fastapi

# TODO https://fastapi.tiangolo.com/tutorial/static-files/
