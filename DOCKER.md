# Docker

## Web application

### Build the image

To build the image with docker run the following command:

```sh
docker build . -t tasosc/mult_res_cell_ann:latest
```

Or to build with podman run the following

```sh
podman build . -t tasosc/mult_res_cell_ann:latest --format docker
```

### Running the container

To run the container run the following command with docker:

```sh
docker run --rm --name mult_res_cell_ann -p 8501:8501  tasosc/mult_res_cell_ann:latest
```

Or to run it with podman run:

```sh
podman run --rm --name mult_res_cell_ann -p 8501:8501  tasosc/mult_res_cell_ann:latest
```

### Accessing the application

You can then access the web application by visiting the URL

[http://localhost:8501](http://localhost:8501)