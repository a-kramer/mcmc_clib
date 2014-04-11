mcmc_clib
=========

is a C program for simplified manifold MALA sampling of ODE model parameters.


Installation Instructions
=========================

install the following dependencies (libraries)
```
External library dependencies
	For ODE: 
		Sundials 2.4.0 or later
	For MCMC:
		GSL	1.15 or later
		CBLAS

For VFGEN:
	CLN 1.3.2 or later
	GiNaC 1.6.2 or later
	mini XML 2.6 or later
```
for example, on ubuntu, you can install the following packages:

    libmxml-dev 
    libmxml1 
    libatlas-base-dev 
    libatlas3gf-base
    libginac-dev 
    libcln-dev 
    libsundials-serial-dev 
    libsundials-cvode1
    libsundials-cvodes2
    libatlas-dev
    libgsl0-dev

to satisfy these dependencies. 

inspect the Makefile.

type
```
make -B
```
note that the option "-B" will compile everything, even if the
binaries appear up to date. It is not necessary once you have compiled
the binaries once for your architecture.

Usage
=====

	./ode_rmhmc_binary -l ./model_esens.so -c ./data.cfg -b \
                           -o sampleND.double -s ${sample_size} > ttr.out

	-b 
	      switches output mode to binary, otherwise text mode (printf)
	      Burn-In Sample will always be printed to stdout

	-c data.cfg
	      configuration (data, reference data, inputs, output function,
              prior, etc.)

	-l ode_model.so
	      shared library file

	-s N
              sample size, Burn-In will be $((sample_size/10))

	-h 
              prints help


The configuration file includes the measurement time specifications,
the data, the standard deviation of observations and the
hyperparameters \mu and \Sigma^{-1} of the Gaussian prior.  Optionally
some of the sampling parameters can be set there for convenience.

The following Parameters can be set in the cfg file:

Property     |  Setting
-----------: | :------------
sample size  |  ```sample_size=[integer]```
step size    |  ```step_size=[double]```
target acceptance | ```acceptance=[0, 1]```
sample file name  | ```output=[string]```
initial condition time | ```t0=[double]```

These definitions should not have spaces before the '=' sign since
everything before = is checked for matches with option names (i.e.  
```« t0 »=0.0;``` is not the same as ```«t0»=0.0;``` ).


Hints
=====

it might be helpful to specify the sample size like this in bash: 
```
   -s $((2**14))
   -s $((10**6))
```
because ```1e6``` or ```1E6``` will not work. (the sample size is an integer).


====================================================================================

Output of parameters sample {p} will have the following structure (in
binary and text mode; though, in binary mode there are no newlines):

sampling is done in logarithmic space, so use exp(p[i]) if you want to simulate trajectories

columns: parameters and log-posterior
   rows: Markov Chain members (i.e. the sample members)
```
   Line 
   L1   p[0] p[1] p[2] ... p[n-1] log-Posterior {\n}
   L2   p[0] p[1] p[2] ... p[n-1] log-Posterior {\n}
   L3   p[0] p[1] p[2] ... p[n-1] log-Posterior {\n}
   L4   p[0] p[1] p[2] ... p[n-1] log-Posterior {\n}
   ...
   L{sample_size}      ...
```
you can load the binary sample using an fread type function.

e.g. in octave: ```[VAL, COUNT] = fread (FID, SIZE, PRECISION,SKIP, ARCH)```

so:
```
     FID=fopen("MySample.double","r");
     sample=fread(FID,[n+1,sample_size],"double");
     fclose(FID);
```
should work fine. In fact, the sources include octave scripts that read and process the sample. 

