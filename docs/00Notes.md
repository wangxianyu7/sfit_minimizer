#sfit_minimizer

A new minimization routine based on the algorithm A. Gould built into sfit.

## Goals

- It should be called like and return output following the conventions of scipy.optimize.minimize()
- Demonstrations:
    - Works like sfit
    - Works better than "Newton-CG" or "BFGS" (or whatever that was called).

## Unit Tests

- Unit tests based on:
    - simple gradient surface
    - microlening event (parallax and PSPL?)

# To Do List

- Architecture:
    1. write high-level use case
        - Review scipy.optiimize.minimize() architecture
    2. Setup high-level sfit_minimizer architecture 
    3. draft internal architecture
- Tests:    
    - Write simple gradient unit tests
    - Find a test PSPL event that doesn't minimize well with scipy.
        - Create demonstrations comparing:
            - original sfit
            - Newton-CG
            - BFGS
        - Extract unit tests based on sfit
        - Create demonstration for python sfit