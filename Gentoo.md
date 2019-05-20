---
title: gentoo
layout: default
---


Gentoo Installation Instructions
================================

Most of the dependencies are available in Portage. You might want to
install better performing linear algebra packages from the [science
overlay](https://wiki.gentoo.org/wiki/Project:Science/Overlay;
e.g. atlas. The [sundials](https://github.com/LLNL/sundials) suite has
changed its API somewhat. Until our mcmc code is changed to adapt to
these changes, you should install version `2.7.0`.

~~~ bash
  # vfgen
  emerge dev-libs/mini-xml \
         sci-mathematics/ginac

  # MCMC
  USE="hl" emerge sci-libs/gsl \
         sci-libs/cblas-reference \
	 =sci-libs/sundials-2.7.0
	 
  # SBtab to hdf5
  emerge dev-libs/glib
~~~

The high level API of hdf5 is used by the sampler. So, it's probably best
to add it to the hdf5 [USE flags](https://wiki.gentoo.org/wiki/USE_flag).

Some of these packages might not have been declared stable. In this
case you will have to add them to the
[package.accept_keywords](https://www.gentoo.org/doc/en/handbook/handbook-x86.xml?part=3&chap=1)
file.

The MCMC solver uses the
[CVODES](http://computation.llnl.gov/casc/sundials/download/download.html).
solver for ODE model integration with forward sensitivity analysis.

As mentioned before, you might want to install
[atlas](http://math-atlas.sourceforge.net/) from the science
overlay. If you do, you have to migrate to the science overlay eselect
modules. Another alternative is `sci-libs/gotoblas2` or
`sci-libs/openblas/`, also from the [science
overlay](https://github.com/gentoo/sci/tree/master/sci-libs/openblas).