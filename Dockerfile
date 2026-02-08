# Use SageMath Binder environment
FROM ghcr.io/sagemath/sage-binder-env:10.4

# Copy repository contents
COPY --chown=sage:sage . /home/sage
