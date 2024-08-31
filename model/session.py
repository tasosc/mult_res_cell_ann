from uuid import UUID
from pydantic import BaseModel

from utilities import CellType
from model.settings import Settings
from fastapi import  UploadFile
class SessionData(BaseModel):
        file : UploadFile
        uuid : UUID
        cells : list[CellType]= list()
        settings : Settings = Settings()
        
    
