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
RUN apt-get install -y r-cran-renv

RUN apt-get install -y pip python3-venv
# Pip install fails if you don't set up a virtual enviroment
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

ADD . /run-sch
WORKDIR /run-sch
RUN pip install .[dev]

ENV RETICULATE_PYTHON_ENV="/opt/venv/"

RUN Rscript setup_r_env.r
# RUN source ../../.venv/bin/activate && Rscript tests/testthat.R
