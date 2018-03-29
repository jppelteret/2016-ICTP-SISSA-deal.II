FROM dealii/dealii:v8.5.0-gcc-mpi-fulldepscandi-debugrelease

MAINTAINER luca.heltai@gmail.com

USER root

RUN apt-get update && apt-get -yq install \
    vim \
    emacs24-nox \
    gdb

USER dealii
WORKDIR /home/dealii
