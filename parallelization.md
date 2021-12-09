# Parallelization


The code is intended to run on a computing cluster for big problems
(~100 parameters) but can also run on a workstation for problems with
~25 parameters `theta`. The appropriate number of nodes is ~16 for bigger
problems.

## MPI

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


## OpenMP

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

## Sensible HPC splits

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

