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
from render_to_pdf import RenderToPdf
from utils.feedback_socket import FeedbackModel, FeedbackQueue

logger = logging.getLogger("uvicorn.error")


def run(session_data : SessionData):
    logger.info("Started analysis for session %s", session_data.uuid.hex)
    feedback_queue = FeedbackQueue(session_data.message_queue)
    pdf = RenderToPdf(session_data.settings.verbosity)
    pdf.set_frontpage(session_data.file.name if session_data.file else "annotated")

    pdf.start_title("Parse dataset")
    pdf.write_text(session_data.file)
    pdf.write_text(str(session_data.settings))
    pp = PreProcessing.build_from(session_data.file, session_data.settings)
    pp.set_render(pdf)
    feedback_queue.put(FeedbackModel(activity=Activity.PARSE_DATASET, message="Parse dataset"))
 
    pdf.start_title("Pre-processing Quality Control")
    pp.qc()

    feedback_queue.put(FeedbackModel(activity=Activity.PP_QC, message="Pre-processing Quality Control"))
    pdf.start_title("Pre-processing Normalization")
    pp.normalization()
    feedback_queue.put(FeedbackModel(activity=Activity.PP_NORM, message="Pre-Processing, Normalization"))
    pdf.start_title("Pre-processing Feature Selection")
    pp.feature_selection()
    feedback_queue.put(FeedbackModel(activity=Activity.PP_FEATURE, message="Pre-Processing, Feature selection"))
    pdf.start_title("Pre-processing Dimensionlity Reduction")
    pp.dimensionality_reduction()
    feedback_queue.put(FeedbackModel(activity=Activity.PP_REDUCTION, message="Pre-Processing, Dimensionality Reduction"))
    pdf.start_title("Pre-processing Visualization")
    fig = pp.visualization()
    pdf.write_fig(fig, "Pre-processing visualization")
    feedback_queue.put(FeedbackModel(activity=Activity.PP_VISUALIAZTION, message="Pre-processing, Visualization"))
    pdf.start_title("Structure Identification")
    si = StructureIdentification(pp.adata, session_data.settings, pdf)
    pdf.start_title("Clustering")
    fig2=si.clustering(session_data.cells)
    pdf.write_fig(fig2)
    feedback_queue.put(FeedbackModel(activity=Activity.SI_CLUSTERING,message="Clustering"))
    pdf.start_title("Annotation")
    fig3=si.annotation()
    pdf.add_page()
    pdf.write_fig(fig3)
    feedback_queue.put(FeedbackModel(activity=Activity.SI_ANNOTATION, message="Annotation"))
    with NamedTemporaryFile(delete=False) as tmp:
        # Write the annotated dataset to download_file
        tmp_path = Path(tmp.name)
        si.write_ann_ds(tmp_path)
        session_data.annotated = tmp_path
        filename=session_data.download_filename if session_data.download_filename else "annotated"
        logger.info("Output filename %s", filename)
        session_data.download_report = filename + ".pdf"
        pdf.set_filename(session_data.annotated.with_suffix(".pdf"))
        pdf.close()
        feedback_queue.put(
            FeedbackModel(activity=Activity.SI_FILE, 
                          link=f"/annotated/{SessionManager.get_session_id(session_data)}/{session_data.download_filename}",
                          report_link=f"/report/{SessionManager.get_session_id(session_data)}/{session_data.download_report}",
                          message="Download annotated dataset"))
    feedback_queue.put(FeedbackModel(activity=Activity.END))
    
