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

## Models

Systems biology models are often quite large, so handling
them should be automatic as much as possible. If a user already has a
good way to do that and can create `cvodes` compatible C files and
headers, then the tools mentioned in this section are not needed.

We have chosen an Sbtab to vfgen workflow for autmoation of source
creation. We rely on the
[vfgen](https://github.com/WarrenWeckesser/vfgen) tool to create C
code for the model, including a Jacobian and sensitivity
equations. VFGEN depends [CLN](https://www.ginac.de/CLN/),
[Ginac](https://www.ginac.de/), and
[mini-xml](https://github.com/michaelrsweet/mxml)
(https://www.msweet.org/mxml/).

SBtab is used to write the model down. It's a spreadsheet format and
allows fairly easy model editing.

## Parallelization


The code is intended to run on a computing cluster for big problems
(~100 parameters) but can also run on a workstation for problems with
~25 parameters `theta`. The appropriate number of nodes is ~16 for bigger
problems.

### MPI

We use a parallel tempering scheme, so each node processes the problem
using a different (thermodynamic) inverse temperature `beta`, with

    gamma = 2 # by default
    beta=(1-MPI_rank/MPI_Comm_size)^gamma
    posterior(theta|data,beta) = likelihood(data|theta)^beta * prior(theta)
    # or equivalently
    LogPosterior(theta|data,beta) = beta*LogLikelihood(data|theta) + LogPrior(theta)

The different MPI workers sometimes exchange positions, if it is a
benefitial switch. So, lower temperature chains can be seens as
exploring the space for chains with higher inverse temperature. This
reduces the risk of chains being stuck in locally isolated modes.


### OpenMP

Experimental data is compared to simulations (one per experiment). The
simulation are done in parallel via OpenMP (parallel for loop), so one MPI instance
per Node is ok if the number of experiments is higher than the number
of cores/node.

For optimal load balancing these two properties have to map to one another nicely:

1. Problem size
   - number of sensible temperature regime
   - number of experiments to process
2. Number of workers
   - MPI nodes
   - number of cores per MPI worker

### Sensible HPC splits

These are guidelines: each parameterization `theta` of a given model
requires a different number of solver steps, so the workers do not return
simultaneously.

Some Examples:

1. 4 different temperatures, say `[1.0 0.5625 0.25 0.0625]`,
   are enough to allow swaps between the regimes and 6 experiments
   have to be simulated, then on a workstation with 8 cores
   this would make sense
   - 4 MPI workers, with 2 OpenMP cores each (with 3 simulations each).
2. 16 temperatures needed (many parameters make swaps more
   difficult), 24 experiments are available, and the machine
   has 32 cores per node.
   - 16 MPI instances on 2 nodes (8 per node),
     but now there are 4 cores per worker on each node, who
     will work on the 24 experiments (6 each).
   - 16 MPI instances on 4 nodes (4 per node), with 8 cores per MPI worker.
     Each OpenMP thread works on 3 experiments.

If the numbers don't split up well, it's probably still ok as the
simulation times don't align well in any case.

It's probably a lot better to have more experiments per node than
choosing a bad number of temperatures.

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
