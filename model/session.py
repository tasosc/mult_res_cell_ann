""" Session model """
import uuid
import pandas as pd
from pydantic import BaseModel
from fastapi import  UploadFile

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

class SessionData(BaseModel):
    """ Model for session """
    file: UploadFile = None
    uuid: uuid.UUID
    cells: list[Cell]
    settings: Settings

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
        data : SessionData = SessionData(uuid=session_id, cells=cells, settings=settings)
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
