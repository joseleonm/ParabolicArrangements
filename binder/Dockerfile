# Usa la imagen oficial de SageMath para Binder
# Basada en SageMath 10.7 (alineado con el requisito del proyecto)
FROM ghcr.io/sagemath/sage-binder-env:10.7

# Define HOME explícitamente para evitar variables no definidas en build
ENV HOME=/home/sage

# Copia el contenido del repositorio al home del usuario sage
COPY --chown=sage:sage . ${HOME}

# Establece el directorio de trabajo
WORKDIR ${HOME}

# Abre el demo directamente al iniciar JupyterLab
RUN mkdir -p ${HOME}/.jupyter && \
    echo "c.ServerApp.default_url = '/lab/tree/notebooks/demo.ipynb'" \
    > ${HOME}/.jupyter/jupyter_server_config.py
