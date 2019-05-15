mcmc_clib
=========

A program written in C for simplified manifold Metropolis adjusted
Lengevin algorithm (SMMALA) sampling of the parameters of ordinary
differential euqations. The main purpose is to fit models to data in
systems biology.

The main program relies on the [GNU Scientific
Library](http://www.gnu.org/software/gsl/doc/html/index.html), and the
[SUNDIALS solver](https://computation.llnl.gov/projects/sundials)
cvodes for initial value problem integration with (forward) sensitiviy
analysis.

Systens biology models are quite large, so handling them should be
automatic. For this purpose we heavily rely on the
[vfgen](https://github.com/WarrenWeckesser/vfgen) tool, which creates
C code for the model. VFGEN depends [CLN](https://www.ginac.de/CLN/),
[Ginac](https://www.ginac.de/), and
[mini-xml](https://github.com/michaelrsweet/mxml)
(https://www.msweet.org/mxml/).

The code is intended to run on a computing cluster for big problems
(~100 parameters) but can also run on a workstation for problens with
~25 parameters. The appropriate number of nodes is ~16 for bigger
problems. We use a parallel tempering scheme, so each node processes
the problem using a different thermodynamic temperature beta, with

    gamma = 2 # by default
    beta=(1-MPI_rank/MPI_Comm_size)^gamma
    Posterior(theta|data,beta) = Likelihood(D|theta)^beta * prior(theta)


Experimental data is compared

This project has a gh page:
[http://a-kramer.github.io/mcmc_clib/](http://a-kramer.github.io/mcmc_clib/),
with more detailed instructions.

Installation Instructions
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
|VFGEN|	CLN |>=1.3.2 or later|
| |GiNaC |>=1.6.2 or later|
||mini XML| >=2.6 or later|

for example, on ubuntu, you can install the following packages, or similar:

    libmxml-dev 
    libmxml1 
    libatlas-base-dev 
    libatlas3-base
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
	Examples:
	./bin/ode_smmala --help
			 to get a list of command line options
			 
	./ode_smmala -l ./model.so -d ./data.h5 \
                     -s ${sample_size} > mcmc_run_1.log 2> mcmc_run_1.err

	important options
	-d data.h5
	      HDF5 file with: data, reference data, inputs, and
              prior

	-l model.so
	      shared library file

	-s $N
              sample size


The contents of the data file are described in doc/documentation.pdf

Result
======

If the run is successful, the result is an hdf5 file containing the
MCMC sample in log-space (natural logarithm), log-posterior
probability values with som eannotation saved as hdf5 attributes and
some description of the contents.

In [GNU Octave](https://www.gnu.org/software/octave/) the file can be
loaded as is via `load sample.h5`.
