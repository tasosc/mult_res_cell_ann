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

"""
Send feedback via websockets
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import dataclasses
from datetime import datetime, timedelta
import io
import logging
from pathlib import Path
from queue import Queue
from typing import Optional
from fastapi import WebSocket
from fastapi.websockets import WebSocketState
from matplotlib.figure import Figure

from model.enums import Activity

logger = logging.getLogger("workflow")

@dataclasses.dataclass
class FeedbackModel:
    """Model for sending feedback from workflow to the FeedbackSocket """ 
    activity: Activity
    message: Optional[str] = None
    link: Optional[str] = None
    image: Optional[Figure] = None
    output_path: Optional[Path] = None
    end : datetime = dataclasses.field(default_factory=datetime.now)
    duration: float = -1

    def set_duration(self, start: datetime):
        if (start > self.end):
            logger.error("Start '%s' is after end '%s'", start, self.end)
        duration : timedelta = self.end - start
        self.duration = duration.total_seconds()


class FeedbackSocket:
    """
    A stopwatch that sends a message to a socket at the end
    """
    def __init__(self, socket: WebSocket):
        self.socket = socket

    async def send(self, feedback: FeedbackModel):
        """
        Specify that the current activity has completed
        """
        if self.socket.client_state != WebSocketState.CONNECTED:
            logger.error("Socket state is not open %s", self.socket.state)
            return

        if (feedback.image):
            with io.BytesIO() as buf:
                fig = feedback.image
                fig.savefig(buf, format='png')
                buf.seek(0)
                await self.socket.send_bytes(buf)
                return

        await self.socket.send_json({'activity': feedback.activity,
                                    'finished': feedback.end.isoformat(sep='T'),
                                    'duration': feedback.duration, 
                                    'message': feedback.message,
                                    'link': feedback.link})
class FeedbackQueue:
    """
    A stopwatch that sends a message to a socket at the end
    """
    def __init__(self, queue: Queue):
        self.queue = queue
        self.start_time = datetime.now()

    def put(self, feedback: FeedbackModel):
        """
        Specify that the current activity has completed
        """
        # Simple message no duration
        if (feedback.activity == Activity.NONE):
            self.queue.put(feedback)
            return

        feedback.set_duration(self.start_time)
        self.queue.put(feedback)
        self.start_time = datetime.now()
