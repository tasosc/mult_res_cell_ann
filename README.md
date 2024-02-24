# Multi Resource Cell Annotation Web application

Multi Resource Cell Annotation Web application is a web application that 
uses multiple resources to help users perform automatic cell annotation.

It is based on [Kostas Lazaros decoupleR example](https://github.com/kostaslazaros/cell_annotation_web_app)

## Requirements 

- Python 3.11
or
- Docker or Podman

## Running with Docker

```sh
docker run --rm --name multi_res_cell_all -p 8501:8501  tasosc/mult_res_cell_ann:1.0
```

## Running with Podman

```sh
podman run --rm --name multi_res_cell_all -p 8501:8501  tasosc/mult_res_cell_ann:1.0
```

## Running directly with Python (Linux)

First you need to set up the virtual environment, [more information](https://docs.python.org/3/library/venv.html)
Note in some Linux distributions instead of `python` you need to call `python3.11` or `python3`

```sh
python -m venv venv
. venv/bin/activate
pip install -r requirements.init_id.txt
```
The activate the environment. 

On Linux:

```sh
. venv/bin/activate
```

On Windows with powershell:

```powershell
venv\Scripts\Activate.ps1
```

On Windows with command line (cmd.exe):

```cmd
venv\Scripts\Activate.bat
```

In the next steps we need to initialize the database, which consists of the following commands.

```sh
# Install dependencies
pip install -r requirements.init_id.txt
# Init database ( On Windows change the / to \ )
python init/init_db.py
```

Now we need to install the dependencies for web application:

```sh
pip install -r requirements.txt
```

Finally, run the application:

```sh
streamlit run main_web_app.py
```



