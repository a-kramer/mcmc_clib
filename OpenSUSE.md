---
title: openSUSE
layout: default
---

# Installation Instructions for OpenSUSE

Use either zypper or YaST to install the following packages, if you
don't have any equivalents of those already:

~~~ bash
# vfgen
zypper install mxml ginac

# MCMC (if necessary)
zypper install openblas_serial hdf5 

# SBtab to hdf5
zypper install glib2
~~~

The MCMC software also requires the
[sundials](https://github.com/LLNL/sundials) solvers, currently
version 2.7.0. You can checkout that version in the repository and
compile it using the installation guide included in the repository.

There is also the [openSUSE Build
Service](https://build.opensuse.org/package/show/science/sundials)
package of the sundials suite, in the [science
project](https://build.opensuse.org/project/show/science). However,
due to API changes that are not yet reflected on our end, we cannot
use this version.
