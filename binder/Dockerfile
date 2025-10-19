FROM ghcr.io/sagemath/sage-binder-env:10.7

# Copia el repo al HOME del usuario estándar de Jupyter (jovyan)
COPY --chown=jovyan:jovyan . /home/jovyan/
WORKDIR /home/jovyan/

# Warm-up opcional: verifica que Sage está accesible durante build
# (NO afecta al arranque de Binder)
RUN sage -c 'print("Sage OK:", version())'
