# syntax=docker/dockerfile:1

# https://hub.docker.com/r/condaforge/miniforge3
FROM condaforge/miniforge3

SHELL [ "/bin/bash", "-c" ]
ARG DEBIAN_FRONTEND=noninteractive

ARG STH_SCH_AMIS_DIR=/ntdmc/sth-sch-amis-integration
ARG STH_SCH_MODEL_DIR=${STH_SCH_AMIS_DIR}/model/ntd-model-sch
ARG FITTING_PREP_DIR=${STH_SCH_AMIS_DIR}/fitting-prep
ARG FITTING_DIR=${STH_SCH_AMIS_DIR}/fitting
ARG PROJECTIONS_PREP_DIR=${STH_SCH_AMIS_DIR}/projections-prep
ARG PROJECTIONS_DIR=${STH_SCH_AMIS_DIR}/projections

RUN apt update && apt install -y \
    build-essential \
    cmake \
    vim \
    curl \
    git \
    unzip \
    openssh-client \
    libssl-dev \
    libgdal-dev \
    libudunits2-dev

RUN conda install --override-channels -c conda-forge -c r --yes --name base \
    python=3.10 \
    pandas \
    joblib \
    r-base \
    r-reticulate \
    r-truncnorm \
    r-dplyr \
    r-magrittr \
    r-invgamma \
    r-tidyr \
    r-devtools \
    r-openxlsx \
    r-hmisc \
    r-mclust \
    r-mnormt \
    r-optparse \
    r-sf \
    r-renv \
    r-readxl \
    r-writexl \
    r-pracma \
    r-mvtnorm

# Cannot activate the conda environment easily
# So instead adjust shell to run everything inside Conda
# from here on out
# https://pythonspeed.com/articles/activate-conda-dockerfile/
SHELL ["conda", "run", "--no-capture-output", "/bin/bash", "-c"]

RUN Rscript -e "install.packages('weights', repos='https://cran.r-project.org/', lib='/opt/conda/lib/R/library')"
RUN Rscript -e "install.packages('AMISforInfectiousDiseases', repos='https://cran.r-project.org/', lib='/opt/conda/lib/R/library')"
RUN conda clean -a -y

# Verify installations
RUN R --version && python --version

# https://medium.com/datamindedbe/how-to-access-private-data-and-git-repositories-with-ssh-while-building-a-docker-image-a-ea283c0b4272
RUN mkdir -p -m 0600 ~/.ssh && \
    ssh-keyscan github.com >> ~/.ssh/known_hosts

# Copy the scripts into the container
# Order of the scripts follows the stages of the pipeline
ADD fitting-prep ${FITTING_PREP_DIR}
ADD fitting ${FITTING_DIR}
ADD projections-prep ${PROJECTIONS_PREP_DIR}
ADD projections ${PROJECTIONS_DIR}
ADD post_AMIS_analysis ${STH_SCH_AMIS_DIR}/post_AMIS_analysis

# Copy the sth_amis package and pyproject.toml
ADD sth_amis ${STH_SCH_AMIS_DIR}/sth_amis
ADD pyproject.toml ${STH_SCH_AMIS_DIR}

# Copy input data directories
ADD Maps-STH ${FITTING_PREP_DIR}/inputs/Maps-STH
ADD Maps-SCH ${FITTING_PREP_DIR}/inputs/Maps-SCH
ADD ESPEN_IU_2021 ${FITTING_PREP_DIR}/inputs/ESPEN_IU_2021

# Get STH/SCH model
ADD --keep-git-dir git@github.com:NTD-Modelling-Consortium/ntd-model-sch.git#updateImportation ${STH_SCH_MODEL_DIR}
RUN cd ${STH_SCH_MODEL_DIR}

WORKDIR ${STH_SCH_AMIS_DIR}

# Install the STH/SCH model
RUN --mount=type=cache,target=/root/.cache/pip cd ${STH_SCH_MODEL_DIR} && pip install .

# Install the sth_amis package
RUN --mount=type=cache,target=/root/.cache/pip cd ${STH_SCH_AMIS_DIR} && pip install -e .

# Set environment variables for all paths
ENV STH_SCH_AMIS_DIR=${STH_SCH_AMIS_DIR}
ENV STH_SCH_MODEL_DIR=${STH_SCH_MODEL_DIR}
ENV PATH_TO_FITTING_PREP_ARTEFACTS="$STH_SCH_AMIS_DIR/fitting-prep/artefacts"
ENV PATH_TO_FITTING_PREP_INPUTS="$STH_SCH_AMIS_DIR/fitting-prep/inputs"
ENV PATH_TO_FITTING_PREP_SCRIPTS="$STH_SCH_AMIS_DIR/fitting-prep/scripts"

ENV PATH_TO_FITTING_ARTEFACTS="$STH_SCH_AMIS_DIR/fitting/artefacts"
ENV PATH_TO_FITTING_INPUTS="$STH_SCH_AMIS_DIR/fitting/inputs"
ENV PATH_TO_FITTING_SCRIPTS="$STH_SCH_AMIS_DIR/fitting/scripts"

ENV PATH_TO_PROJECTIONS_PREP_ARTEFACTS="$STH_SCH_AMIS_DIR/projections-prep/artefacts"
ENV PATH_TO_PROJECTIONS_PREP_INPUTS="$STH_SCH_AMIS_DIR/projections-prep/inputs"
ENV PATH_TO_PROJECTIONS_PREP_SCRIPTS="$STH_SCH_AMIS_DIR/projections-prep/scripts"

ENV PATH_TO_PROJECTIONS_ARTEFACTS="$STH_SCH_AMIS_DIR/projections/artefacts"
ENV PATH_TO_PROJECTIONS_INPUTS="$STH_SCH_AMIS_DIR/projections/inputs"
ENV PATH_TO_PROJECTIONS_SCRIPTS="$STH_SCH_AMIS_DIR/projections/scripts"

ENV RETICULATE_PYTHON=/opt/conda/bin/python
ENV RETICULATE_PYTHON_FALLBACK=FALSE

# Add pipeline runner and test script
ADD run_pipeline.py ${STH_SCH_AMIS_DIR}

# Set entrypoint to run python directly (conda environment is already activated via SHELL)
ENTRYPOINT ["python", "run_pipeline.py"]

# Default command shows help
CMD ["--help"]
