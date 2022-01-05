"""
Tests for the MulensModel functions = direct comparisons to Andy's fortran sfit.
"""

"""
# Thoughts on test data:

Should start by generating two perfect test datasets with only a few (but different numbers of) points each and with 
different values of fs, fb. Then, run through sfit.f to retrieve bmat, dvec, cmat, at several successive steps and then
final values and sigmas.

# Tests:

1. single dataset
2. both datasets
3. fixing fb = 0 for one datasets
4. fixing fs = X for one dataset
5. Different step sizes: fixed, adaptive

"""