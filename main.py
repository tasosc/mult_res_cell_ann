"""
REST API Controller for multiple resource cell annotation
"""
import asyncio
import datetime
import logging
from pathlib import Path
from queue import Empty
import shutil
from tempfile import NamedTemporaryFile
from typing import List
from fastapi import BackgroundTasks, FastAPI, UploadFile, Query, WebSocket, HTTPException
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
logger = logging.getLogger("uvicorn.error")
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
async def analyze_file(session: str, file: UploadFile,  background_tasks: BackgroundTasks):
    """
    analyse uploaded file use
    """
    if not session:
        raise HTTPException(status_code=422, detail="Session not specified")

    # upload scRNAseq https://fastapi.tiangolo.com/reference/uploadfile/#fastapi.UploadFile
    session_data : SessionData = SessionManager.get_session(session)

    if not session_data:
        raise HTTPException(status_code=404, detail="Session not found")
    if session_data.file:
        raise HTTPException(status_code=409, detail="file already uploaded")
    filename = file.filename
    with NamedTemporaryFile(delete=False) as tmp:
        tmp_path = Path(tmp.name)
        session_data.file = tmp_path
        session_data.download_filename = filename
        shutil.copyfileobj(file.file, tmp)

    file.close()
    background_tasks.add_task(run, session_data)
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
        await socket.close(code=4004, reason="Session not found")
        return

    if session_data.has_finished:
        await socket.close(code=4022, reason="analysis has finished")
        return
    if session_data.has_started:
        await socket.close(code=4009, reason="monitoring is active")
        return

    session_data.has_started = True
    await socket.accept()
    logger.info("socket accepted")
    try:
        if not session_data.file:
            logger.info("file not uploaded")
            await asyncio.sleep(10)
        if session_data.file:
            logger.info("file uploaded")
            await monitor_analysis(session_data=session_data, socket=socket)
    except Empty:
        pass
    await socket.close()
    session_data.has_started = False


async def monitor_analysis(session_data: SessionData, socket: WebSocket):
    feedback = FeedbackSocket(socket)
    queue = session_data.message_queue

    if queue.empty():
        logger.info("queue empty")
        await asyncio.sleep(10)
    while not queue.empty():
        logger.info("queue not empty")
        current : FeedbackModel = queue.get(block=False)
        await feedback.send(current)
        queue.task_done()

        if current.activity == Activity.END:
            logger.info("end reached")
            session_data.has_finished = True
            break

        if queue.empty():
            await asyncio.sleep(10)
        # Save to PDF pages (report0
    # https://stackoverflow.com/questions/11328958/save-multiple-plots-in-a-single-pdf-file
    # and/or using websockets ? return json with fig & text
    # https://stackoverflow.com/questions/71936110/correct-way-of-connecting-websocket-events-to-update-my-react-component

    # https://stackoverflow.com/questions/73550398/how-to-download-a-large-file-using-fastapi

# TODO https://fastapi.tiangolo.com/tutorial/static-files/
