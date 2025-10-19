FROM ghcr.io/sagemath/sage-binder-env:10.7

# Crea usuario 'user' (como en el Dockerfile oficial del proyecto)
USER root
ARG NB_USER=user
ARG NB_UID=1000
ENV NB_USER=${NB_USER}
ENV NB_UID=${NB_UID}
ENV HOME=/home/${NB_USER}
RUN adduser --disabled-password --gecos "Default user" --uid ${NB_UID} ${NB_USER}

# Copia notebooks (si tienes) y TU módulo .sage a una carpeta del home
COPY notebooks/* ${HOME}/
COPY parabolic_dga.sage ${HOME}/lib/parabolic_dga.sage

# *** Solo cambia dueño de lo que acabas de copiar ***
RUN chown -R ${NB_USER}:${NB_USER} ${HOME}

# Kernel de Sage para Jupyter
USER ${NB_USER}
RUN mkdir -p $(jupyter --data-dir)/kernels && \
    ln -s /sage/venv/share/jupyter/kernels/sagemath $(jupyter --data-dir)/kernels

# Haz que Sage encuentre tu módulo sin rutas absolutas
ENV SAGE_PATH=${HOME}/lib:${SAGE_PATH}
ENV PATH="/sage:${PATH}"

WORKDIR ${HOME}

# Silenciar el warning de nodejs de JupyterLab (opcional)
RUN mkdir -p ${HOME}/.jupyter && \
    printf "import logging\n\nclass NoNodeJSWarningFilter(logging.Filter):\n    def filter(self, record):\n        return 'Could not determine jupyterlab build status without nodejs' not in record.getMessage()\n\nlogging.getLogger('LabApp').addFilter(NoNodeJSWarningFilter())\n" > ${HOME}/.jupyter/jupyter_lab_config.py
