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
  +-----------+      +---(use)---+      +------------+      +--------------+
  |           |      |           |      |  (CVODES)  |      |              |
  | SBtab (M) +--+-->+   VFGEN   +----->+  ODE code  +--+-->+ ./ode_smmala |
  |*.{tsv,ods}|  |   |           |      |  model.vf  |  |   |              |
  +----+------+  |   +-----------+      +------------+  |   +--------------+
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
  +-----------+       +-------------+      +--------------+
  |           |       |[gcc -shared]|      |              |
  | model.vf  +---+-->+  ODE code   +----->+ ./ode_smmala |
  |           |   |   |  model.so   |      |              |
  +----+------+   |   +-------------+      +--------+-----+
       |          |                                 ^
       |    +-----+-(create)------------------+     |
       |    | vfgen cvodes:func=yes model.vf  |     | 
       |    +---------------------------------+     |
       |                                            |
       |                                            |
       |        +-----(use)------+         +--------+--+
	   |        |                |         |           |
       +------->+ ./sbtab_import +-------->+  data.h5  |
	            |                |         |           |
				+----------------+         +-----------+

```


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
	Examples:
```bash
mpirun -N 12 ./ode_smmala -l ./model.so -d ./data.h5 -s ${sample_size} 

-a $ACCEPTANCE_RATE
			Target acceptance value (all markov chains will be tuned for this acceptance).

-d, --hdf5 ./data.h5
			data.h5 is a file that contains the data points and the conditions of measurement in hdf5 format. A suitable h5 file is produced by the hdf5_import program bundled with ode_smmala.

-g $G
			This will define how the inverse MCMC temperatures β are chosen: β = (1-rank/R)^G, where R is MPI_Comm_Size.

-i $STEP_SIZE
			The initial step size of each markov chain, this will usually be tuned later to get the desired acceptance rate $A (-a $A).

-l ./ode_model.so
			ode_model.so is a shared library containing the CVODE functions of the model.

-m $M
			If this number is larger than 1.0, each MPI rank will get a different initial step size s: step_size(rank)=STEP_SIZE*M^(rank).

-o ./output_file.h5
			Filename for hdf5 output. This file will contain the log-parameter sample and log-posterior values. The samples will have attributes that reflect the markov chain setup.

-p, --prior-start
			Start the markov chain at the center of the prior. Otherwise it will be started from the DefaultParameters in the vfgen file.
-r, --resume
			Resume from last sampled MCMC point. Only the last MCMC position is read from the resume file. Everything else about the problem can be changed.
-s $N
			$N sample size. default N=10.

-t,--init-at-t $T_INITIAL
			Specifies the initial time «t0» of the model integration [initial value problem for the ordinary differential equation in x; x(t0)=x0]

--seed $SEED
		Set the gsl pseudo random number generator seed to $SEED. (perhaps --seed $RANDOM)
```

The contents of the data file are described in [documentation.pdf](./doc/documentation.pdf)

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
