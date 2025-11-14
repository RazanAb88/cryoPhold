# Base image pinned to Python 3.10
FROM python:3.10-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    build-essential \
    wget \
    bzip2 \
    cmake \
    libomp-dev \
 && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app
COPY . /app

# Upgrade pip tooling
RUN pip install --upgrade pip setuptools wheel

# Install JAX and NumPy
RUN pip install --extra-index-url https://storage.googleapis.com/jax-releases/jax_releases.html \
    "jax==0.4.34" \
    "jaxlib==0.4.34" \
    "numpy==1.26.4"

# Install JupyterLab
RUN pip install jupyterlab jupyter

# Expose Jupyter port
EXPOSE 8888

# Install Miniconda pinned to Python 3.10
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py310_24.9.2-0-Linux-x86_64.sh && \
    bash Miniconda3-py310_24.9.2-0-Linux-x86_64.sh -b -p /opt/miniconda && \
    /opt/miniconda/bin/conda install -n base -y python=3.10

# Configure Conda channels and install HH-suite
RUN /opt/miniconda/bin/conda config --add channels defaults && \
    /opt/miniconda/bin/conda config --add channels conda-forge && \
    /opt/miniconda/bin/conda config --add channels bioconda && \
    /opt/miniconda/bin/conda config --set channel_priority strict && \
    /opt/miniconda/bin/conda install -n base hhsuite=3.3.0 -y

# Install matplotlib and pandas with versions compatible with Python 3.10
RUN pip install "matplotlib>=3.5,<4.0" "pandas==1.5.3" --prefer-binary

# Install ColabFold with AlphaFold support
RUN pip install "colabfold[alphafold]" --prefer-binary

# Launch JupyterLab
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root"]
