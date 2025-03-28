""" Session model """
import dataclasses
from pathlib import Path
from queue import Queue
from typing import Optional
import uuid
import pandas as pd
from pydantic import BaseModel

from model.settings import Settings


class Cell(BaseModel):
    """model for a selected cell with the selected genes"""
    cell_type: str
    genes: list[str]
    def get_selected(self) -> pd.DataFrame:
        selection = set(self.genes)
        if len(selection) == 0:
            return pd.DataFrame({"cell_name": [], "Symbol": []})
        selection.discard(None)
        df = pd.DataFrame(data=selection, columns=["Symbol"])
        df.insert(0, "cell_name", self.cell_type)
        return df

@dataclasses.dataclass
class SessionData:
    """ Model for session """
    uuid: uuid.UUID
    cells: list[Cell]
    settings: Settings
    message_queue: Queue
    file: Optional[Path] = None
    annotated: Optional[Path] = None
    download_filename: Optional[str] = None
    download_report: Optional[str] = None
    has_finished: bool = False
    has_started: bool = False


class SessionManager:
    """
    Manages the lifecycle of a session
    """
    sessions: dict = {}
    
    @staticmethod
    def create_session(cells: list[Cell], settings: Settings) -> str:
        """
        Create a new session
        """
        session_id= uuid.uuid4()
        data : SessionData = SessionData(uuid=session_id, cells=cells, settings=settings, message_queue=Queue())
        SessionManager.sessions[session_id.hex]=data
        return session_id.hex
    
    @staticmethod
    def get_session(sessiod_id: str) -> SessionData:
        """
        Get session withh the specified uuid
        """
        return SessionManager.sessions[sessiod_id] if sessiod_id in SessionManager.sessions else None
    
    @staticmethod
    def delete_session(session_id: str) -> None:
        """
        Delete session
        """
        del SessionManager.sessions[session_id]
    @staticmethod
    def get_session_id(session: SessionData) -> Optional[str]:
        if session is None:
            return None
        return session.uuid.hex
