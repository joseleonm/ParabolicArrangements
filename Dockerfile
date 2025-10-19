FROM ghcr.io/sagemath/sage-binder-env:10.7

# Copia el repo y fija permisos para el usuario jovyan
COPY --chown=jovyan:jovyan . /home/jovyan/
WORKDIR /home/jovyan/

# Asegura que el script de arranque (start) sea ejecutable si existe
RUN chmod +x /home/jovyan/start || true
RUN chmod +x binder/postBuild

# Mensaje de build
RUN echo "Loaded ParabolicArrangements with SageMath Binder base (10.7)"
