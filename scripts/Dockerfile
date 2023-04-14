# syntax=docker/dockerfile:1

FROM ubuntu:latest
RUN apt-get update
RUN apt-get install python3-pip -y
RUN python3 -m pip install astropy
RUN python3 -m pip install matplotlib
RUN python3 -m pip install astroquery
RUN python3 -m pip install scipy
RUN python3 -m pip install gaiadr3-zeropoint

ADD julia-1.6.7-linux-x86_64.tar  /soft
COPY ./run /run
ENV GAIA_ROOT /
WORKDIR /data

RUN /soft/julia-1.6.7/bin/julia /run/src/install.jl