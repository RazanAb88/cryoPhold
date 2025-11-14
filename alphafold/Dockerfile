FROM python:3.10-slim

# Install system dependencies
RUN apt-get update && apt-get install -y git build-essential wget

# Set working directory
WORKDIR /app
COPY . /app

# Install JAX 0.4.34 and matching jaxlib
RUN pip install --upgrade pip
RUN pip install --extra-index-url https://storage.googleapis.com/jax-releases/jax_releases.html \
    "jax==0.4.34" \
    "jaxlib==0.4.34" \
    "numpy==1.26.4"
# Install JupyterLab and CLI
RUN pip install jupyterlab jupyter

# Expose Jupyter port
EXPOSE 8888

# Launch JupyterLab
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root"]