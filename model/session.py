""" Session model """
from uuid import UUID
from pydantic import BaseModel
from fastapi import  UploadFile

from model.settings import Settings


class Cell(BaseModel):
    """model for a selected cell with the selected genes"""
    name: str
    genes: list[str]

class SessionData(BaseModel):
    """ Model for session """
    file: UploadFile = None
    uuid: UUID
    cells: list[Cell]
    settings: Settings

class SessionManager:
    """
    Manages the lifecycle of a session
    """
    sessions: dict
    
    @classmethod
    def create_session(cls, cells: list[Cell], settings: Settings) -> str:
        """
        Create a new session
        """
        uuid : UUID = UUID()
        data : SessionData = SessionData(cells=cells, settings=settings)
        cls.sessions[uuid.hex]=data
        return uuid.hex
    
    @classmethod
    def get_session(cls, uuid: str) -> SessionData:
        """
        Get session withh the specified uuid
        """
        return cls.sessions[uuid] if uuid in cls.sessions else None
    
    @classmethod
    def delete_session(cls, uuid: str) -> None:
        """
        Delete session
        """
        del cls.sessions[uuid]
