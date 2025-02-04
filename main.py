"""
REST API Controller for multiple resource cell annotation
"""
import datetime
from pathlib import Path
import shutil
from tempfile import NamedTemporaryFile
from typing import List
from fastapi import FastAPI, UploadFile, Query, WebSocket, HTTPException
from fastapi.responses import FileResponse
from fastapi.middleware.cors import CORSMiddleware

from model.enums import Activity, SvdSolverOptions
from model.session import Cell, SessionData, SessionManager
from model.settings import Settings
from worfklow import run
from store import Store
from utilities import CellType
from utils.feedback_socket import FeedbackModel, FeedbackSocket

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
    await monitor_analysis(session_data=session_data, socket=socket)
    await socket.close()


async def monitor_analysis(session_data: SessionData, socket: WebSocket):
    feedback = FeedbackSocket(socket)

    async def render_text(something: str, expected_verbosity=1):
        if session_data.settings.verbosity < expected_verbosity:
            return
        await feedback.send(FeedbackModel(activity=Activity.NONE, message=something))

    for model in run(session_data=session_data, render_text=render_text):
        await feedback.send(model)

        # Save to PDF pages (report0
    # https://stackoverflow.com/questions/11328958/save-multiple-plots-in-a-single-pdf-file
    # and/or using websockets ? return json with fig & text
    # https://stackoverflow.com/questions/71936110/correct-way-of-connecting-websocket-events-to-update-my-react-component

    # https://stackoverflow.com/questions/73550398/how-to-download-a-large-file-using-fastapi

# TODO https://fastapi.tiangolo.com/tutorial/static-files/
