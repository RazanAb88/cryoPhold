FROM python:3.10-slim

RUN apt-get update && apt-get install -y git build-essential wget
WORKDIR /app
COPY . /app

RUN pip install --upgrade pip
RUN pip install --extra-index-url https://storage.googleapis.com/jax-releases/jax_releases.html \
    "jax==0.4.13" \
    "jaxlib==0.4.13" \
    "numpy==1.26.4" \
    jupyterlab \
    jupyter


EXPOSE 8888
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root"]
