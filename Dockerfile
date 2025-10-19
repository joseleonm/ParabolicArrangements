FROM ghcr.io/sagemath/sage-binder-env:10.7

# Copia repo y fija permisos
COPY --chown=jovyan:jovyan . /home/jovyan/
WORKDIR /home/jovyan/

# (NUEVO) Si existe binder/postBuild, dale permiso y ejecútalo en build
# para “precocinar” Sage y acelerar el arranque.
RUN if [ -f postBuild ]; then chmod +x postBuild && /bin/bash postBuild; fi

# Asegura que el script de arranque sea ejecutable
RUN chmod +x /home/jovyan/start

RUN echo "Loaded ParabolicArrangements with SageMath Binder base (10.7)"
