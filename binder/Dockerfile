# Dockerfile para Binder con SageMath + tu módulo parabolic_dga.sage

FROM ghcr.io/sagemath/sage-binder-env:10.7

USER root

# Crear usuario (igual que en el Dockerfile oficial)
ARG NB_USER=user
ARG NB_UID=1000
ENV NB_USER user
ENV NB_UID 1000
ENV HOME /home/${NB_USER}
RUN adduser --disabled-password --gecos "Default user" --uid ${NB_UID} ${NB_USER}

# Copiar los notebooks al home
COPY notebooks/* ${HOME}/

# 🔹 Copiar también tu módulo .sage a la instalación del sistema Sage
#   para que se pueda usar con `load("parabolic_dga.sage")` desde cualquier lugar
COPY parabolic_dga.sage /sage/

# 🔹 Opcional: si tienes más código fuente, puedes copiar todo un directorio:
# COPY src/ /sage/src/

# Asegurar permisos correctos
RUN chown -R ${NB_USER}:${NB_USER} ${HOME} /sage/

# Instalar kernel de Sage para Jupyter
RUN mkdir -p $(jupyter --data-dir)/kernels && \
    ln -s /sage/venv/share/jupyter/kernels/sagemath $(jupyter --data-dir)/kernels

# Añadir Sage al PATH
ENV PATH="/sage:$PATH"

# Cambiar a usuario normal
USER ${NB_USER}
WORKDIR /home/${NB_USER}

# Configurar Jupyter para silenciar warnings de nodejs
RUN mkdir -p /home/${NB_USER}/.jupyter && \
    echo "\
import logging\n\
\n\
class NoNodeJSWarningFilter(logging.Filter):\n\
    def filter(self, record):\n\
        return 'Could not determine jupyterlab build status without nodejs' not in record.getMessage()\n\
\n\
logging.getLogger('LabApp').addFilter(NoNodeJSWarningFilter())\n\
" > /home/${NB_USER}/.jupyter/jupyter_lab_config.py
