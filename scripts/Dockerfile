# syntax=docker/dockerfile:1
# sudo chmod 666 /var/run/docker.sock

FROM ubuntu:latest
RUN apt-get update
RUN apt-get install python3-pip -y
RUN python3 -m pip install astropy
RUN python3 -m pip install matplotlib
RUN python3 -m pip install astroquery
RUN python3 -m pip install scipy
RUN python3 -m pip install gaiadr3-zeropoint

RUN apt-get install wget -y
WORKDIR /soft
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.7-linux-x86_64.tar.gz
RUN tar -zxvf julia-1.6.7-linux-x86_64.tar.gz

WORKDIR /
RUN apt-get install git -y
RUN rm -Rf run
RUN git clone https://github.com/bosscha/gaia-shock run
RUN /soft/julia-1.6.7/bin/julia /run/src/install.jl

ENV GAIA_ROOT /
WORKDIR /data