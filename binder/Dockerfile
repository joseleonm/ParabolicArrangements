# Usa la imagen oficial de SageMath para Binder
# Basada en SageMath 10.4 (imagen estable más reciente disponible en ghcr.io)
FROM ghcr.io/sagemath/sage-binder-env:10.4

# Copia el contenido del repositorio al home del usuario sage
COPY --chown=sage:sage . ${HOME}

# Establece el directorio de trabajo
WORKDIR ${HOME}