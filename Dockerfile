FROM ubuntu:24.04

SHELL ["/bin/bash", "-c"]

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN apt-get install -y \
    build-essential \
    cmake \
    vim curl git python3 \
    r-base-dev \
    r-cran-devtools

# Install R packages using Ubuntu package manager to ensure we get a valid version
RUN apt-get install -y r-cran-matrix r-cran-testthat

RUN Rscript -e 'devtools::install_github("evandrokonzen/AMISforInfectiousDiseases-dev")'

RUN mkdir -p /sch-model/
WORKDIR /sch-model2
RUN git clone https://github.com/NTD-Modelling-Consortium/ntd-model-sch.git
WORKDIR /sch-model2/ntd-model-sch

RUN apt-get install -y pip python3-venv
# Pip install fails if you don't set up a virtual enviroment
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Temp branch for combining the fitting team work with Thibaults work
RUN git switch cloud-fitting

RUN pip install .

ENV RETICULATE_PYTHON_ENV="/opt/venv/"
WORKDIR /sch-model2/ntd-model-sch/sch_simulation/amis_integration
# RUN source ../../.venv/bin/activate && Rscript tests/testthat.R


