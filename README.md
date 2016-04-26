# PR5BP: A Taylor integrator for a minimalistic model of the epsilon ring of Uranus

The model includes Uranus (with an oblate-panet potential), Cordelia and Ophelia (the shepherd
moons of the epsilon ring), Ariel, and an ensemble of non-interacting particles.

## Usage @ Miztli

To compile at  Miztli, use

`mpicxx -I $TMPU/<my-working-dir>/PR5BP/blitz-0.10 -O2 -o SF15T.o PRNplus1BP_blitz.cpp`

To submit a job at Miztli, use

`bsub -oo salida -eo error -q q_8p_1h -n 8 mpirun ./SF15T.o`