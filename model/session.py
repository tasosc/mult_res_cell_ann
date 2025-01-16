""" Session model """
import uuid
from pydantic import BaseModel
from fastapi import  UploadFile

from model.settings import Settings


class Cell(BaseModel):
    """model for a selected cell with the selected genes"""
    cell_type: str
    genes: list[str]

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
    
    @classmethod
    def get_session(cls, sessiod_id: str) -> SessionData:
        """
        Get session withh the specified uuid
        """
        return cls.sessions[sessiod_id] if sessiod_id in cls.sessions else None
    
    @classmethod
    def delete_session(cls, session_id: str) -> None:
        """
        Delete session
        """
        del cls.sessions[session_id]
