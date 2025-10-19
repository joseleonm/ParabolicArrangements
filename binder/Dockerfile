FROM ghcr.io/sagemath/sage-binder-env:10.7

# Copia el repo al HOME del usuario de Binder (jovyan) con permisos correctos
COPY --chown=jovyan:jovyan . /home/jovyan/
WORKDIR /home/jovyan/

# Asegura que 'start' sea ejecutable si existe
RUN chmod +x /home/jovyan/start || true

# NO pongas ENTRYPOINT; la base ya inicia bien sin él
# Lanza el script de arranque con bash -lc (hereda entorno de JupyterHub)
CMD ["bash", "-lc", "/home/jovyan/start"]
