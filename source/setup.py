from setuptools import setup, find_packages

setup(
    name='sfit_minimizer',
    version='1.0.2',
    packages=find_packages(where='source'),
    package_dir={'': 'source'},
)
