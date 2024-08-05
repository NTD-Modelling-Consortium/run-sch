# Run STH & SCH

A project for running the STH model through AMIS in a Docker container. 

## Pre-requistites

Either install [Docker](https://www.docker.com/), or refer to the Dockerfile to see
what needs to be installed and setup on your machine. 

See [Renv lock files](#Renv-lock-files) for details on getting the right set of R dependencies. 

## Running the AMIS algorithm using Docker

1. Build the Dockerfile:

```bash
$ docker build -f Dockerfile . -t sth-amis-fitting-environment
```

2. Run the appropriate fitting script in the Docker image:

```bash
$ docker run sth-amis-fitting-environment Rscript <NAME OF FITTING SCRIPT>
```

## Running the AMIS algorithm without Docker

Having setup your local environment in an equivalent way to the Docker file, 
you should be able to run `Rscript <NAME OF FITTING SCRIPT>`. 

## Renv lock files

Due to different versions of R requiring different packages, the lock file depends on the version of R.

The script [setup_r_env.R](setup_r_env.R) will automatically select the correct one if it exists and
then installs the correct set of dependencies into the R virtual environment. 

If there is no lock file for the version of R you are using, you can probably create it by doing the following:

```R
renv::init()
renv::install()
renv::snapshot()
```

The snapshot will create a renv.lock file locally with the right versions - please check this in
with the name renv_4_x_y.lock (where x and y are the minor and point version of your version of R). 
