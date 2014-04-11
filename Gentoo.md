Gentoo Installation Instructions
================================

some of the dependencies are available in Portage (though you might
want to install packages from overlays) 
``` emerge dev-libs/mini-xml \
sci-libs/gsl \ sci-libs/cblas-reference \ sci-libs/cln \
sci-mathematics/ginac ``` 

Some of these packages might not have been declared stable in the
required version yet. In this case you will have to add them to the
[package.accept_keywords](https://www.gentoo.org/doc/en/handbook/handbook-x86.xml?part=3&chap=1) file.

Unfortunately the sundials ode solvers are not available through
portage. You can manually install them after downloading
[CVODES](http://computation.llnl.gov/casc/sundials/download/download.html). Gain
```root``` priviliges, extract and follow the instructions. It should
be very easy to install cvodes.

You might want to install [ATLAS](http://math-atlas.sourceforge.net/)
from the science overlay. If you do, you have to migrate to the
science overlay eselect modules.