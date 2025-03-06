#
# This file is part of the mult_res_cell_ann distribution
# (https://github.com/tasosc/mult_res_cell_ann).
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
import json
import logging
from pathlib import Path
from tempfile import NamedTemporaryFile

import pandas as pd

from cell_structure_id import StructureIdentification
from model.enums import Activity
from model.session import SessionData, SessionManager
from preprocessing import PreProcessing
from utils.feedback_socket import FeedbackModel, FeedbackQueue

logger = logging.getLogger("uvicorn.error")


def run(session_data : SessionData):
    logger.info("Started analysis for session %s", session_data.uuid.hex)
    feedback_queue = FeedbackQueue(session_data.message_queue)
    def render_text(something, expected_verbosity=1):
        if session_data.settings.verbosity < expected_verbosity:
            return
        if isinstance(something, pd.DataFrame):
            something = something.to_json()
            feedback_queue.put(FeedbackModel(activity=Activity.NONE, message=json.loads(something)))
            return
        elif not isinstance(something, str):
            logger.error("Message is not str %s", type(something))
            logger.error("as string %s", something)
        feedback_queue.put(FeedbackModel(activity=Activity.NONE, message=something))

    verbosity = session_data.settings.verbosity
    pp = PreProcessing.build_from(session_data.file, session_data.settings)
    pp.render.set_render_text_lambda(render_text)
    feedback_queue.put(FeedbackModel(activity=Activity.PARSE_DATASET, message="Parse dataset"))
 
    pp.qc()
    feedback_queue.put(FeedbackModel(activity=Activity.PP_QC, message="Pre-processing Quality Control"))
    pp.normalization()
    feedback_queue.put(FeedbackModel(activity=Activity.PP_NORM, message="Pre-Processing, Normalization"))
    pp.feature_selection()
    feedback_queue.put(FeedbackModel(activity=Activity.PP_FEATURE, message="Pre-Processing, Feature selection"))
    pp.dimensionality_reduction()
    feedback_queue.put(FeedbackModel(activity=Activity.PP_REDUCTION, message="Pre-Processing, Dimensionality Reduction"))
    fig = pp.visualization()
    if (verbosity >= 1):
        feedback_queue.put(FeedbackModel(activity=Activity.NONE, image=fig))
    feedback_queue.put(FeedbackModel(activity=Activity.PP_VISUALIAZTION))
    si = StructureIdentification(pp.adata, session_data.settings)
    si.render.set_render_text_lambda(render_text)
    fig2=si.clustering(session_data.cells)
    if (verbosity >= 1):
        feedback_queue.put(FeedbackModel(activity=Activity.NONE, image=fig2))
    feedback_queue.put(FeedbackModel(activity=Activity.SI_CLUSTERING))
    fig3=si.annotation()
    if (verbosity >= 1):
        feedback_queue.put(FeedbackModel(activity=Activity.NONE, image=fig3))
    feedback_queue.put(FeedbackModel(activity=Activity.SI_ANNOTATION))
    with NamedTemporaryFile(delete=False) as tmp:
        # Write the annotated dataset to download_file
        tmp_path = Path(tmp.name)
        si.write_ann_ds(tmp_path)
        session_data.annotated = tmp_path
        filename = session_data.file.name if session_data.file is not None and session_data.file.name is not None else "annotated"
        filename = filename + ".h5ad"
        session_data.download_filename = filename
        feedback_queue.put(FeedbackModel(activity=Activity.SI_FILE, link=f"/annotated/{SessionManager.get_session_id(session_data)}/{filename}", message="Download annotated dataset"))
    feedback_queue.put(FeedbackModel(activity=Activity.END))
    
