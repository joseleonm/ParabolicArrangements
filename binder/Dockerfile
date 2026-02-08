# Usa la imagen oficial de SageMath 10.7+ para Binder
FROM ghcr.io/sagemath/sage-binder-env:10.7

# Copia el contenido del repositorio al home del usuario sage
COPY --chown=sage:sage . ${HOME}

# Establece el directorio de trabajo
WORKDIR ${HOME}