FROM ghcr.io/sagemath/sage-binder-env:10.7

# Copia tu repo al HOME del usuario de Binder (jovyan) con permisos correctos
COPY --chown=jovyan:jovyan . /home/jovyan/
WORKDIR /home/jovyan/

# Asegura que 'start' sea ejecutable si existe
RUN chmod +x /home/jovyan/start || true

# Usa tini como entrypoint (estándar de Jupyter)
ENTRYPOINT ["tini", "-g", "--"]

# Lanza tu script de arranque (ver siguiente paso)
CMD ["/home/jovyan/start"]
