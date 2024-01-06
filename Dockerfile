# 
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

FROM python:3.11-slim AS build
WORKDIR /usr/src/app
COPY requirements.init_db.txt ./requirements.txt
RUN --mount=type=cache,target=/root/.cache/pip pip install -r requirements.txt
COPY init/*.py init/*.sql init/*.json init/7k.txt ./init/
RUN --mount=type=cache,target=/usr/src/app/init/cache \
	mkdir data \
	&& python init/init_db.py

FROM python:3.11-slim
# see https://docs.streamlit.io/knowledge-base/tutorials/deploy/docker
RUN apt-get update \
	&& apt-get install -y curl \
	&& rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/* \
	&& adduser webapp
WORKDIR /home/webapp
USER webapp
ENV PATH="/home/webapp/.local/bin:${PATH}" 
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY *.py config.json ./
COPY data/*P1.h5ad ./data/
COPY --from=build /usr/src/app/data/cells.db ./data/
EXPOSE 8501

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

ENTRYPOINT ["streamlit", "run", "main_web_app.py", "--server.port=8501", "--server.address=0.0.0.0"]
