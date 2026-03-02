# Dockerfile for Binder (SageMath)
# Reference: https://mybinder.readthedocs.io/en/latest/tutorials/dockerfile.html

FROM ghcr.io/sagemath/sage-binder-env:10.7

USER root

# Create user with uid 1000 (Binder default)
ARG NB_USER=user
ARG NB_UID=1000
ENV NB_USER=${NB_USER}
ENV NB_UID=${NB_UID}
ENV HOME=/home/${NB_USER}
RUN adduser --disabled-password --gecos "Default user" --uid ${NB_UID} ${NB_USER}

# Copy only the module and demo notebook into the user's home
RUN mkdir -p ${HOME}/repo/notebooks
COPY parabolic_arrangements.sage ${HOME}/repo/
COPY notebooks/demo.ipynb ${HOME}/repo/notebooks/
RUN chown -R ${NB_USER}:${NB_USER} ${HOME}/repo

# Switch to the user
USER ${NB_USER}

# Install Sage kernel to Jupyter
RUN mkdir -p $(jupyter --data-dir)/kernels
RUN ln -s /sage/venv/share/jupyter/kernels/sagemath $(jupyter --data-dir)/kernels

# Make Sage accessible from anywhere
ENV PATH="/sage:$PATH"

# Start in the repository root
WORKDIR ${HOME}/repo

# Suppress the perpetual nodejs warning
RUN mkdir -p ${HOME}/.jupyter
RUN echo "\
import logging\n\
\n\
class NoNodeJSWarningFilter(logging.Filter):\n\
    def filter(self, record):\n\
        return 'Could not determine jupyterlab build status without nodejs' not in record.getMessage()\n\
\n\
logging.getLogger('LabApp').addFilter(NoNodeJSWarningFilter())\n\
" > ${HOME}/.jupyter/jupyter_lab_config.py

# Open the demo notebook by default
RUN echo \"c.ServerApp.default_url = '/lab/tree/repo/notebooks/demo.ipynb'\" \
  > ${HOME}/.jupyter/jupyter_server_config.py
