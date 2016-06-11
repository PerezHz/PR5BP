# PR5BP: A Taylor integrator for a minimalistic model of the epsilon ring of Uranus

The model includes Uranus (with an oblate-panet potential), Cordelia and Ophelia (the shepherd
moons of the epsilon ring), Ariel, and an ensemble of non-interacting particles.

## Usage @ Miztli

First, load the `mpi` and `blitz` modules

`module load mpi`
`module load blitz<TAB>`

To compile at  Miztli, use

`mpicxx -I $TMPU/<blitz-dir> -O2 -o <executable-file> PRNplus1BP_blitz.cpp`

To submit a job at Miztli, use

`bsub -oo <output-file> -eo <error-file> -q <queue-name> -n <number-of-processes> mpirun <executable-file>`