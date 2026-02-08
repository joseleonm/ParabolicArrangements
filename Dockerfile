# Use SageMath Binder environment
FROM ghcr.io/sagemath/sage-binder-env:10.4

# Copy only the necessary files to the working directory
COPY --chown=sage:sage parabolic_arrangement_cohomology.sage /home/sage/
COPY --chown=sage:sage notebooks/ /home/sage/notebooks/
