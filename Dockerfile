# Sagemath version for binder
FROM ghcr.io/sagemath/sage-binder-env:10.7
COPY --chown=jovyan:jovyan . /home/jovyan/
WORKDIR /home/jovyan/
RUN echo "Loaded ParabolicArrangements with SageMath Binder base (10.7)"
