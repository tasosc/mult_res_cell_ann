from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional
from fastapi import WebSocket
from matplotlib.figure import Figure
from pydantic import BaseModel

from model.enums import Activity

class FeedbackModel(BaseModel):
    activity: Activity
    message: Optional[str] = None
    link: Optional[str] = None
    image: Optional[Figure] = None
    outputPath : Optional[Path] = None

class FeedbackSocket:
    """
    A stopwatch that sends a message to a socket at the end
    """
    def __init__(self, socket: WebSocket):
        self.socket = socket
        self.start_time = datetime.now()

    def start(self):
        """
        Start timer
        """
        self.start_time = datetime.now()
    async def completed(self, current_activity: Activity, message: str|None = None, link: str|None = None):
        """
        Specify that the current activity has completed
        """
        end_time = datetime.now()
        duration : timedelta = end_time - self.start_time
        await self.socket.send_json({'activity': str(current_activity),
                                    'finished': end_time.isoformat(sep='T'),
                                    'duration': duration.total_seconds(), 
                                    'message': message,
                                    'link': link})
        self.start()
    
    async def send(self, feedback: FeedbackModel):
        """
        Specify that the current activity has completed
        """
        end_time = datetime.now()
        duration : timedelta = end_time - self.start_time
        await self.socket.send_json({'activity': feedback.activity,
                                    'finished': end_time.isoformat(sep='T'),
                                    'duration': duration.total_seconds(), 
                                    'message': feedback.message,
                                    'link': feedback.link})
        self.start()


