#sfit_minimizer

A new minimization routine based on the algorithm A. Gould built into sfit.

## Goals

- It should be called like and return output following the conventions of 
scipy.optimize.minimize()
- Demonstrations:
    - Works like sfit
    - Works better than "Newton-CG" or "BFGS"

## Unit Tests

- Unit tests based on:
    - simple gradient surface
    - microlening event (parallax and PSPL?)

# To Do List

- Architecture:
    1. [DONE] write high-level use case
    2. [Working] Review scipy.optiimize.minimize() architecture
    3. [Started] Setup high-level sfit_minimizer architecture
    4. Think about algorithm components and how it interacts with jac, etc.
    (necessary to determine inputs) 
    5. draft internal architecture
- Tests:    
    - Write simple gradient unit tests
    - Find a test PSPL event that doesn't minimize well with scipy.
        - Create demonstrations comparing:
            - original sfit
            - Newton-CG
            - BFGS
        - Extract unit tests based on sfit
        - Create demonstration for python sfit
        
# Notes

From the scipy.optimize.minimize documentation:

**Custom minimizers**
    It may be useful to pass a custom minimization method, for example
    when using a frontend to this method such as `scipy.optimize.basinhopping`
    or a different library.  You can simply pass a callable as the ``method``
    parameter.
    The callable is called as ``method(fun, x0, args, **kwargs, **options)``
    where ``kwargs`` corresponds to any other parameters passed to `minimize`
    (such as `callback`, `hess`, etc.), except the `options` dict, which has
    its contents also passed as `method` parameters pair by pair.  Also, if
    `jac` has been passed as a bool type, `jac` and `fun` are mangled so that
    `fun` returns just the function values and `jac` is converted to a function
    returning the Jacobian.  The method shall return an `OptimizeResult`
    object.
    The provided `method` callable must be able to accept (and possibly ignore)
    arbitrary parameters; the set of parameters accepted by `minimize` may
    expand in future versions and then these parameters will be passed to
    the method.  You can find an example in the scipy.optimize tutorial.
    
So this is what I'm working for (=method, a callable):

    if meth == '_custom':
        # custom method called before bounds and constraints are 'standardised'
        # custom method should be able to accept whatever bounds/constraints
        # are provided to it.
        return method(fun, x0, args=args, jac=jac, hess=hess, hessp=hessp,
                      bounds=bounds, constraints=constraints,
                      callback=callback, **options)
    
    
Newton-CG source code:

https://github.com/scipy/scipy/blob/9b0c66f564d85d0c3c6eea2a2704df1e3b09b48f/scipy/optimize/optimize.py#L1523