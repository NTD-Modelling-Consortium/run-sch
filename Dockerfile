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
    python=3.10.4 \
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

# Verify Python version and force reinstall if needed
RUN python --version && python -c "import sys; assert sys.version_info[:2] == (3, 10) and sys.version_info[2] == 4, f'Wrong Python version: {sys.version}'"

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
# STH AMIS prior input (also kept at repo root for convenience)
ADD RawDataForPrior.csv ${FITTING_DIR}/inputs/RawDataForPrior.csv
ADD projections-prep ${PROJECTIONS_PREP_DIR}
ADD projections ${PROJECTIONS_DIR}
ADD post_AMIS_analysis ${STH_SCH_AMIS_DIR}/post_AMIS_analysis

# Copy input data directories from Google Cloud Storage
ADD https://storage.googleapis.com/ntd-data-storage/pipeline/sth/Maps-STH.tar.gz ${FITTING_PREP_DIR}/inputs/Maps-STH.tar.gz
ADD https://storage.googleapis.com/ntd-data-storage/pipeline/sch/ESPEN_IU_2021.tar.gz ${FITTING_PREP_DIR}/inputs/ESPEN_IU_2021.tar.gz
ADD https://storage.googleapis.com/ntd-data-storage/pipeline/sch/Maps-SCH.tar.gz ${FITTING_PREP_DIR}/inputs/Maps-SCH.tar.gz

# Extract input data archives
RUN mkdir -p ${FITTING_PREP_DIR}/inputs ${FITTING_PREP_DIR}/artefacts && \
    cd ${FITTING_PREP_DIR}/inputs && \
    tar -xzf Maps-STH.tar.gz && \
    tar -xzf ESPEN_IU_2021.tar.gz && \
    tar -xzf Maps-SCH.tar.gz && \
    rm Maps-STH.tar.gz ESPEN_IU_2021.tar.gz Maps-SCH.tar.gz

# Add fitting-prep artifacts from Google Cloud Storage
ADD https://storage.googleapis.com/ntd-data-storage/pipeline/sth/fitting-prep-artefacts-sth.tar.gz ${FITTING_PREP_DIR}/fitting-prep-artefacts-sth.tar.gz
ADD https://storage.googleapis.com/ntd-data-storage/pipeline/sch/fitting-prep-artefacts-sch.tar.gz ${FITTING_PREP_DIR}/fitting-prep-artefacts-sch.tar.gz

# Extract fitting-prep artifacts (strip top-level directory)
RUN mkdir -p ${FITTING_PREP_DIR}/artefacts && \
    cd ${FITTING_PREP_DIR} && \
    tar -xzf fitting-prep-artefacts-sth.tar.gz --strip-components=1 -C artefacts && \
    tar -xzf fitting-prep-artefacts-sch.tar.gz --strip-components=1 -C artefacts && \
    rm fitting-prep-artefacts-sth.tar.gz fitting-prep-artefacts-sch.tar.gz

# Get STH/SCH model (master includes *_params_projections.txt for SCH/STH)
ADD --keep-git-dir git@github.com:NTD-Modelling-Consortium/ntd-model-sch.git#master ${STH_SCH_MODEL_DIR}
RUN cd ${STH_SCH_MODEL_DIR}

WORKDIR ${STH_SCH_AMIS_DIR}

# Install the STH/SCH model
RUN --mount=type=cache,target=/root/.cache/pip cd ${STH_SCH_MODEL_DIR} && pip install .

# Create symlink to use source model data files instead of installed package files
# This ensures the model uses the data files from the source repository
RUN rm -rf /opt/conda/lib/python3.10/site-packages/sch_simulation/data && \
    ln -s ${STH_SCH_MODEL_DIR}/sch_simulation/data /opt/conda/lib/python3.10/site-packages/sch_simulation/data

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
