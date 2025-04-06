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
from typing import Annotated, List
from fastapi import BackgroundTasks, Depends, FastAPI, UploadFile, Query, WebSocket, HTTPException, WebSocketDisconnect, WebSocketException, status
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

async def get_session_data(session: str):
    """
    Get session_data if it exists; else throw an 422 error 
    """
    session_data : SessionData = SessionManager.get_session(session)
    if not session_data:
        raise HTTPException(status_code=422, detail="Session or file not found")
    return session_data


@app.post("/session", status_code=201)
def create_session(settings: Settings, cells: list[Cell]):
    """
    Create a session and return the session id
    """
    return {"session": SessionManager.create_session(cells=cells, settings=settings)}

@app.post("/dataset/{session}", status_code=201)
async def analyze_file(session_data : Annotated[SessionData, Depends(get_session_data)],
                       file: UploadFile,
                       background_tasks: BackgroundTasks):
    """
    uploaded scRNA-seq dataset and create background task that starts the analysis workflow
    """
    if session_data.file:
        raise HTTPException(status_code=409, detail="file already uploaded")
    filename = file.filename
    with NamedTemporaryFile(delete=False) as tmp:
        tmp_path = Path(tmp.name)
        session_data.file = tmp_path
        session_data.download_filename = filename
        shutil.copyfileobj(file.file, tmp)

    await file.close()
    logger.info("Starting background task")
    background_tasks.add_task(run, session_data)
    logger.info("Started background  task")
    # TODO background task for deleting session in X time

    return {
        "activity": Activity.UPLOAD_DATASET,
        "finished": datetime.datetime.now().isoformat(sep="T"),
        "duration": -1,
        "message": f"Uploaded scRNA-seq dataset {filename}",
        "link": None,
    }

@app.get("/annotated/{session}/{filename}")
async def get_annotated_dataset(session_data : Annotated[SessionData, Depends(get_session_data)], filename: str):
    if not session_data.annotated or session_data.download_filename != filename:
        raise HTTPException(status_code=404, detail="file not found")
    return FileResponse(path = session_data.annotated, media_type="application/octet-stream", filename=session_data.download_filename, content_disposition_type="attachment")

@app.get("/report/{session}/{filename}")
async def get_report(session_data : Annotated[SessionData, Depends(get_session_data)], filename: str):
    if not session_data.annotated or session_data.download_report != filename:
        raise HTTPException(status_code=404, detail="file not found")
    return FileResponse(path = session_data.annotated.with_suffix(".pdf"), media_type="application/octet-stream", filename=session_data.download_report, content_disposition_type="attachment")



async def get_session_data_ws(session: str):
    session_data : SessionData = SessionManager.get_session(session)
    if not session_data:
        raise WebSocketException(code=status.WS_1008_POLICY_VIOLATION, reason="cannot find session id")

    if session_data.has_finished:
        raise WebSocketException(code=4022, reason="analysis for this has finised")

    return session_data


@app.websocket("/ws/{session}")
async def get_feedback(socket: WebSocket, session_data : Annotated[SessionData, Depends(get_session_data_ws)]):
    session_id = session_data.uuid.hex
    await socket.accept()

    logger.info("socket accepted for %s", session_id)
    try:
        if not session_data.file:
            logger.info("file not uploaded")
           # await asyncio.sleep(10)
        if session_data.file:
            logger.info("file uploaded")
            await monitor_analysis(session_data=session_data, socket=socket)
    except WebSocketDisconnect:
        logger.warning("Client disconnected, %s", session_id)
        return
    await socket.close(code=1000, reason="End of line")

async def monitor_analysis(session_data: SessionData, socket: WebSocket):
    feedback = FeedbackSocket(socket)
    queue = session_data.message_queue

    logger.info("Wait for message to start")
    is_ready = await socket.receive_text()
    logger.info("Got '%s' to start", is_ready)

    while not session_data.has_finished:
        logger.info("Getting next item in the queue")
        current : FeedbackModel = await queue.get()
        logger.info("Got... %s. Trying to send", current)
        try:
            await feedback.send(current)
        except TypeError as e:
            logger.error(e)
            logger.error("Feedback message type: %s", type(current.message))
            logger.error("Feedback message : %s", current.message)
        queue.task_done()

        if current.activity == Activity.END:
            logger.info("end reached")
            session_data.has_finished = True
            queue.shutdown()

        # Save to PDF pages (report0
    # https://stackoverflow.com/questions/11328958/save-multiple-plots-in-a-single-pdf-file
    # and/or using websockets ? return json with fig & text
    # https://stackoverflow.com/questions/71936110/correct-way-of-connecting-websocket-events-to-update-my-react-component

    # https://stackoverflow.com/questions/73550398/how-to-download-a-large-file-using-fastapi

# TODO https://fastapi.tiangolo.com/tutorial/static-files/
