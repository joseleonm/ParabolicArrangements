FROM ghcr.io/sagemath/sage-binder-env:10.7

# Copia el repo (todo) al home del usuario jovyan con permisos correctos
COPY --chown=jovyan:jovyan . /home/jovyan/
WORKDIR /home/jovyan/

# Si tienes un script start, déjalo ejecutable (si no existe, no falla)
RUN chmod +x /home/jovyan/start || true

# Log útil en build
RUN echo "Loaded ParabolicArrangements with SageMath Binder base (10.7)"
