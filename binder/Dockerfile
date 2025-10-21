FROM ghcr.io/sagemath/sage-binder-env:10.7

USER root
ARG NB_USER=user
ARG NB_UID=1000
ENV NB_USER=${NB_USER}
ENV NB_UID=${NB_UID}
ENV HOME=/home/${NB_USER}
RUN adduser --disabled-password --gecos "Default user" --uid ${NB_UID} ${NB_USER}

COPY notebooks/* ${HOME}/
COPY parabolic_dga.sage ${HOME}/parabolic_dga.sage

RUN chown -R ${NB_USER}:${NB_USER} ${HOME}

USER ${NB_USER}
RUN mkdir -p $(jupyter --data-dir)/kernels && \
    ln -s /sage/venv/share/jupyter/kernels/sagemath $(jupyter --data-dir)/kernels

ENV SAGE_PATH=${HOME}/lib:${SAGE_PATH}
ENV PATH="/sage:${PATH}"

WORKDIR ${HOME}

RUN mkdir -p ${HOME}/.jupyter && \
    printf "import logging\n\nclass NoNodeJSWarningFilter(logging.Filter):\n    def filter(self, record):\n        return 'Could not determine jupyterlab build status without nodejs' not in record.getMessage()\n\nlogging.getLogger('LabApp').addFilter(NoNodeJSWarningFilter())\n" > ${HOME}/.jupyter/jupyter_lab_config.py
