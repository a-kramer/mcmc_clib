mcmc_clib
=========

A program written in C for simplified manifold Metropolis adjusted
Lengevin algorithm (SMMALA) sampling of the parameters of ordinary
differential euqations. The main purpose is to fit models to data in
systems biology.

This project is supported by the [Human Brain
Project](https://www.humanbrainproject.eu/en/) (EU) and
[EBRAINS](https://ebrains.eu/) (EU), details in
[ACKNOWLEDGMENTS.md](ACKNOWLEDGMENTS.md).

The main program relies on the [GNU Scientific
Library](http://www.gnu.org/software/gsl/doc/html/index.html), and the
[SUNDIALS solver](https://computation.llnl.gov/projects/sundials)
cvodes for initial value problem integration with (forward) sensitiviy
analysis.

## Target Users

Researchers in the natural sciences studying ordinary differential
equation models. This software can be used to obtain a parameter
sample for the model.

We have been specifically motivated by examples in systems
biology. But, the sampler is not specific to biological models.

## Models

Systems biology models are often quite large, so handling
them should be automatic as much as possible. If a user already has a
good way to do that and can create `cvodes` compatible C files and
headers, then the tools mentioned in this section are not needed.

Internally, we have chosen an [SBtab](https://sbtab.net) to vfgen workflow for
autmoation of source creation. We rely on the
[vfgen](https://github.com/WarrenWeckesser/vfgen) tool to create C
code for the model, including a Jacobian and sensitivity
equations. VFGEN depends [CLN](https://www.ginac.de/CLN/),
[Ginac](https://www.ginac.de/), and
[mini-xml](https://github.com/michaelrsweet/mxml)
(https://www.msweet.org/mxml/).

[SBtab](https://sbtab.net) is a format in systems biology and is
suited for writing/drafting biological models. It is an external
project and not covered by documentation here and is not mandatory.
However, our [SBtabVFGEN](https://github.com/a-kramer/SBtabVFGEN)
repository has some documentation on how to write SBtab files.

Functions in the [SBtabVFGEN](https://github.com/a-kramer/SBtabVFGEN)
package print vfgen files from source files written using the SBtab
format.

We don't recommend that you start by writing the model in some
other environment and then export it. SBtab's role here is only to
write the model from scratch. It's a spreadsheet format and allows
fairly easy model editing. 

If you already have a model, e.g. as a
[COPASI](https://github.com/copasi/COPASI) project, then it may be
more convenient to create `.vf` files or C code directly from COPASI's
export functions. 

The workflow could be:


```
  User written:                           generated          MCMC Sampling
  +-----------+      +---(use)---+      +------------+      +==============+
  | SBtab (M) |      |           |      |  (CVODES)  |      I              I
  |           +--+-->+   VFGEN   +----->+  ODE code  +--+-->+ ./ode_smmala I
  |*.{tsv,ods}|  |   |           |      |  model.vf  |  |   I              I
  +----+------+  |   +-----------+      +------------+  |   +==========+===+
       |         |                                      |              ^
       |    +----+-----(use)----+                 +-----+-(create)-+   |
       |    | sbtab_to_vfgen()  |                 |[gcc -shared]   |   |
       |    +-------------------+                 |                |   |
       |                                          |    model.so    |   |
       |                                          +----------------+   |
       |        +-----(use)------+         +-(created)-+               |
       |        |                |         |           |               |
       +------->+ ./sbtab_import +-------->+  data.h5  +---------------+
                |                |         |           |
                +----------------+         +-----------+
```


But, this is also possible:

```
  User written:          generated          MCMC Sampling
  +-----------+       +-------------+      +==============+
  |           |       |[gcc -shared]|      I              I
  | model.vf  +---+-->+  ODE code   +----->+ ./ode_smmala I
  |           |   |   |  model.so   |      I              I
  +-----------+   |   +-------------+      +========+=====+
                  |                                 ^
            +-----+-(create)------------------+     |
            | vfgen cvodes:func=yes model.vf  |     | 
            +---------------------------------+     |
                                                    |
                                                    |
                +----(use any)---+         +--------+--+
                | R (hdf5r),     |         |           |
                | MATLAB,        +-------->+  data.h5  |
                | h5import       |         |           |
                +----------------+         +-----------+

``` 

where
[h5import](https://support.hdfgroup.org/HDF5/doc/RM/Tools.html#Tools-Import) is one of the many 
[command line tools](https://support.hdfgroup.org/products/hdf5_tools/#cmd) from the HDF5 project.

## Slightly more Detailed Installation Guide

This project has a gh page:
[http://a-kramer.github.io/mcmc_clib/](http://a-kramer.github.io/mcmc_clib/),
with more detailed instructions.

Quick Installation Instructions
=========================

Dependencies (libraries)

|Component|library|version|
|--------:|:-----:|:------|
|data storage| hdf5| |
|Parallelization|OpenMPI| |
| |OpenMP| |
|ODEs|	Sundials CVODES |=2.7.0|
|MCMC|	GSL |>=2.3|
| |CBLAS| |
|VFGEN|	CLN |>=1.3.2|
| |GiNaC |>=1.6.2|
||mini XML| >=2.6|

for example, on ubuntu, you can install the following packages, or similar:

    libmxml-dev 
    libatlas-base-dev 
    libginac-dev 
    libcln-dev
    libsundials-dev
    libgsl-dev
    libhdf5-dev
    libopenmpi-dev
    libgomp1

to satisfy these dependencies. 

There is no build system yet. To compile the code you will need to edit the Makefile.
Makefile.ubuntu is agood place to start. Once you are satisfied:
```
make
```

Usage
=====

```
mpirun -N 12 ./ode_smmala -l ./model.so -d ./data.h5 -s ${sample_size} 

-a $ACCEPTANCE_RATE
			Target acceptance value (all markov chains will be 
			tuned for this acceptance).

-d, --hdf5 ./data.h5
			data.h5 is a file that contains the data points and 
			the conditions of measurement in hdf5 format. A suitable 
			h5 file is produced by the hdf5_import program 
			bundled with ode_smmala.

-g $G
			This will define how the inverse MCMC temperatures β are chosen: 
			β = (1-rank/R)^G, 
			where R is MPI_Comm_Size.

-i $STEP_SIZE
			The initial step size of each markov chain, this will 
			usually be tuned later to get the desired acceptance rate, 
			see $A (-a $A).

-l ./ode_model.so
			ode_model.so is a shared library containing the 
			CVODE functions of the model.

-m $M
			If this number is larger than 1.0, each MPI rank will get 
			a different initial step size s: step_size(rank)=STEP_SIZE*M^(rank).

-o ./output_file.h5
			Filename for hdf5 output. This file will contain the 
			log-parameter sample and log-posterior values. The samples 
			will have attributes that reflect the markov chain setup.

-p, --prior-start
			Start the markov chain at the center of the prior. Otherwise 
			it will be started from the DefaultParameters in the vfgen file.
-r, --resume
			Resume from last sampled MCMC point. Only the last MCMC position 
			is read from the resume file. Everything else about the problem 
			can be changed.
-s $N
			$N sample size. default N=10.

-t,--init-at-t $T_INITIAL
			Specifies the initial time «t0» of the model integration 


--seed $SEED
		Set the gsl pseudo random number generator seed to $SEED. 
		(e.g. --seed $RANDOM)
```

The allowed contents of the data file are described in [documentation.pdf](./doc/documentation.pdf)

The software is parallelized, see [parallelization.md](./parallelization.md).

Result
======

If the run is successful, the result is an hdf5 file (see `-o` option) containing the
MCMC sample in log-space (natural logarithm), log-posterior
probability values with some annotation saved as hdf5 attributes.

In [GNU Octave](https://www.gnu.org/software/octave/) the file can be
loaded as is via `load sample.h5`. The `load()` function reads
everything and disregards all attributes.

In [R](https://www.r-project.org), the file can be loaded using the
[hdf5r package](https://www.bioconductor.org/packages/release/bioc/html/rhdf5.html),
which on [ubuntu](https://ubuntu.com/) systems is also available in
the package manager as `r-bioc-rhdf5` (e.g. `r-bioc-rhdf5/focal,focal
2.30.1+dfsg-1build1`). The hdf5r package allows fine grained control
and can import attributes.

MATLAB has import
[functions](https://www.mathworks.com/help/matlab/ref/h5read.html) for
hdf5 natively, e.g.: `data = h5read(filename,ds)`.
