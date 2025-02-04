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
from datetime import datetime, timedelta
import io
import logging
from pathlib import Path
from typing import Callable, Optional
from fastapi import WebSocket
from fastapi.websockets import WebSocketState
from matplotlib.figure import Figure

from model.enums import Activity

logger = logging.getLogger("workflow")

# @dataclasses.dataclass
class FeedbackModel:
    """Model for sending feedback from workflow to the FeedbackSocket """ 
    def __init__(
        self,
        activity: Activity,
        message: Optional[str] = None,
        link: Optional[str] = None,
        image: Optional[Figure] = None,
        output_path: Optional[Path] = None,
    ):
        self.activity: Activity = activity
        self.message: Optional[str] = message
        self.link: Optional[str] = link
        self.image: Optional[Figure] = image
        self.output_path: Optional[Path] = output_path
        self.end = datetime.now()
        self.duration = timedelta(seconds=0)

    def stop_watch(self, command: Callable[[], None]):
        start = datetime.now()
        command()
        self.end = datetime.now()
        self.duration = self.end - start

class FeedbackSocket:
    """
    A stopwatch that sends a message to a socket at the end
    """
    def __init__(self, socket: WebSocket):
        self.socket = socket
        self.start_time = datetime.now()

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

        end_time = datetime.now()
        duration : timedelta = end_time - self.start_time
        await self.socket.send_json({'activity': feedback.activity,
                                    'finished': end_time.isoformat(sep='T'),
                                    'duration': duration.total_seconds(), 
                                    'message': feedback.message,
                                    'link': feedback.link})
        self.start_time = datetime.now()
