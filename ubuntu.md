---
title: Ubuntu
layout: default
---

# Installation Instructions for Ubuntu


The main software is labelled as MCMC. The other two are necessary
only for model handling and creating the hdf5 data files. You can
create models and data files via other methods.

Use a command similar to this:

~~~ bash
sudo apt install libmxml-dev libginac-dev libatlas-base-dev libgsl-dev libhdf5-dev openmpi-bin openmpi-common libopenmpi-dev libglib2.0-0 
~~~

This should work for all [ubuntu derivatives](http://www.ubuntu.com/about/about-ubuntu/derivatives)
and ubuntu based distributions like [Linux Mint](http://linuxmint.com/ and [Elementary OS](http://elementaryos.org/).

|Component|Package|URL|
|--------:|------:|:--|
|vfgen|libmxml-dev|https://packages.ubuntu.com/disco/libmxml-dev|
||libginac-dev|https://packages.ubuntu.com/disco/ginac-tools|
|MCMC|libgsl-dev|https://packages.ubuntu.com/disco/libgsl-dev|
||libatlas-base-dev|https://packages.ubuntu.com/disco/libatlas-base-dev|
||libhdf5-dev|https://packages.ubuntu.com/disco/libhdf5-dev|
||libsundials-dev (2.7.0+dfsg)|https://packages.ubuntu.com/bionic/libsundials-dev|
|OpenMPI|libopenmpi-dev|https://packages.ubuntu.com/disco/libopenmpi3|
|SBtab import|libglib2.0-0|https://packages.ubuntu.com/disco/libglib2.0-0|



Package names and versions can be searched at
[Ubuntu packages](https://packages.ubuntu.com/)

On the ubuntu versions `>bionic` the
[sundials](http://computation.llnl.gov/casc/sundials/download/download.html)
solvers are packaged in a newer version than 2.7.0. Unfortunately, we
don't yet use the newer API, and you should install the older `2.7.0`
release. If necessary, pull the [sundials git
repository](https://github.com/LLNL/sundials/tree/v2.7.0) from GitHub
and install manually as described in the [installation
instructions](https://github.com/LLNL/sundials/blob/v2.7.0/INSTALL_GUIDE.pdf).


