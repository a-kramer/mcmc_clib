---
title: windows
layout: default
---

# Installation Instructions for Microsoft Windows systems:

It is perhaps best to run a virtual machine and follow the GNU/Linux
instructions. But, since the software is supposed to run on several
nodes of a computing cluster, MS Windows is perhaps a bad choice in
the first place.

It will almost certainly be better to use the approach of
`Matlab+PESTO+AMICI` as mentioned earlier.

## Cygwin

It may be possible to use [Cygwin](https://www.cygwin.com/) to run the
software, which largely depends on the availability of the libraries
we use. Cygwin is a collection of tools, some of which are usually
preinstalled on a GNU/Linux system. Most importantly a C compiler like
[gcc](http://gcc.gnu.org/) and [GNU
make](http://www.gnu.org/software/make/).

During installation you will have a choice of which components of
cygwin you want to install, pay special attention to the categories
`Math,` `Science,` and `libraries`. They contain libraries
such as the [GNU scientific library](http://www.gnu.org/software/gsl/)
and the basic linear algebra subroutines:
[BLAS](https://www.gnu.org/software/gsl/manual/html_node/BLAS-Support.html). If
in doubt install everything in those categories. Make sure that `gnu
make` is selected for installation as well.

If you want to use [vfgen](http://www.warrenweckesser.net/vfgen/) to
create the C model files, you will need to compile it first, this is
the first thing the `make` command will do. To successfully compile
[vfgen](https://github.com/WarrenWeckesser/vfgen) you will need
[ginac](http://www.ginac.de/) and its dependeny
[CLN](http://www.ginac.de/CLN/), and also
[minixml](http://sourceforge.net/projects/minixml/). For
parallelization, we also need any OpenMPI packages and for data
storage: hdf5.

1. [CLN](http://www.ginac.de/CLN/)
2. [minixml](http://sourceforge.net/projects/minixml/)
3. [ginac](http://www.ginac.de/)
4. [hdf5](https://support.hdfgroup.org/HDF5/)
5. [OpenMPI](https://www.open-mpi.org/)


