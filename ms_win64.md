---
title: windows
layout: default
---

# Installation Instructions for Microsoft Windows systems:

If you already know how to obtain the required libraries and know how
to compile a C project, don't read further. 

Until we can provide a binary release, you will have to compile the
sources and dependencies yourself. Keep in mind that this is more
complicated than on a Linux system, because there is no package
manager to do it for you.

Disclaimer: We do not officially recommend or support a windows
installation. So far, we have no windows developer and all information
below is tentative.

To emulate the steps of an installation on a Linux system, we
recommend the use of [Cygwin](https://www.cygwin.com/). Cygwin is a
collection of handy tools which are usually preinstalled on a Linux
system. Most importantly a c compiler [gcc](http://gcc.gnu.org/) and
[GNU make](http://www.gnu.org/software/make/).

During installation you will have a choice of which components of
cygwin you want to install, pay special attention to the categories
```Math,``` ```Science,``` and ```libraries```. They contain libraries
such as the [GNU scientific library](http://www.gnu.org/software/gsl/)
and the basic linear algebra subroutines:
[BLAS](https://www.gnu.org/software/gsl/manual/html_node/BLAS-Support.html). If
in doubt install everything in those categories. Make sure that gnu
make is selected for installation as well.

If you want to use [vfgen](http://www.warrenweckesser.net/vfgen/) to
create the C model files, you will need to compile it first, this is
the first thing the ```make``` command will do. To successfully
compile [vfgen](http://www.warrenweckesser.net/vfgen/) you will need
[ginac](http://www.ginac.de/) and its dependencies:
[CLN](http://www.ginac.de/CLN/) and
[minixml](http://sourceforge.net/projects/minixml/). Since you don't
have a package manager, you will have to do it manually. To do that,
follow the instructions on the respective webpages:

1. [CLN](http://www.ginac.de/CLN/)
2. [minixml](http://sourceforge.net/projects/minixml/)
3. [ginac](http://www.ginac.de/)

in that order. download the appropriate source archives and extract
them in your home folder inside the cygwin installation. Then start
the cygwin terminal and read the instructions in each of the INSTALL
files (one for each of the three above). To change folders, use the
```cd «folder»``` command. To show the contents of a folder use ```ls
-l```. To read the contents of a text file: ```less «file»``` (quit
with ```q```). If you need more commands, it is highly recommended
that you read the manual of the shell you are using (most likely
[bash](http://www.gnu.org/software/bash/)).

When you have compiled those libraries following the instructions in
the INSTALL files, change into the man folder of our software and run
the ```make -B``` command to recompile everything.

If the compilation process returns no errors you can now run the test
shell script: ```$ ./run_test.sh```

