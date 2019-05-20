---
title: Arch
layout: default
---

# Installation Instructions for Arch Linux

The main software is labelled as MCMC. The other two are necessary
only for model handling and creating the hdf5 data files. You can
create models and data files via other methods.

|Component|Package|URL|
|--------:|------:|:--|
|vfgen|mxml|https://www.archlinux.org/packages/community/x86_64/mxml/|
||ginac|https://www.archlinux.org/packages/community/x86_64/ginac/|
|MCMC|gsl|https://www.archlinux.org/packages/extra/x86_64/gsl/|
||lapack|https://www.archlinux.org/packages/extra/x86_64/lapack/|
||hdf5|https://www.archlinux.org/packages/community/x86_64/hdf5/|
|SBtab import|glib2|https://www.archlinux.org/packages/core/x86_64/glib2/|
|hdf5|https://www.archlinux.org/packages/community/x86_64/hdf5/|


Package names and versions can be searched at the 
[arch linux package database](https://www.archlinux.org/packages/)

You should also install the [sundials](http://computation.llnl.gov/casc/sundials/download/download.html)
as well.

solvers. However, the code of this software package is not yet
compatible with the current version of the [sundials solvers on
arch](https://aur.archlinux.org/packages/sundials/). So, as a
workaround until this is changed you should checkout the 2.7.0 version
from this repository from [github](https://github.com/LLNL/sundials).


