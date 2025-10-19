FROM ghcr.io/sagemath/sage-binder-env:10.7
COPY --chown=jovyan:jovyan . /home/jovyan/
WORKDIR /home/jovyan/
RUN chmod +x /home/jovyan/start
RUN echo "Loaded ParabolicArrangements with SageMath Binder base (10.7)"

