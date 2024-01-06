# Docker

## Web application

### Build the image

To build the image with docker run the following command:

```sh
docker build . -t eap/cellann_webapp
```

Or to build with podman run the following

```sh
podman build . -t eap/cellann_webapp --format docker
```

### Running the container

To run the container run the following command with docker:

```sh
docker run --rm --name cellann_webapp -p 8501:8501  eap/cellann_webapp:latest
```

Or to run it with podman run:

```sh
podman run --rm --name cellann_webapp -p 8501:8501  eap/cellann_webapp:latest
```

### Accessing the application

You can then access the web application by visiting the URL

[http://localhost:8501](http://localhost:8501)