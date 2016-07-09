#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

# Dynamically figure out the version
version = __import__('ribolands').__version__

setup(
    name='ribolands',
    description='energy landscapes and folding kinetics of nucleic acids',
    long_description=readme,
    version=version,
    license=license,
    author='Stefan Badelt',
    author_email='stef@tbi.univie.ac.at',
    # TODO: url='http://rna.tbi.univie.ac.at/ribolands/', 
    # TODO: #install_requires=['matplotlib','numpy','pandas'],
    test_suite='tests',
    packages=['ribolands'],
    scripts=['scripts/pylands_spatch.py',
      'scripts/cofold_kinetics.py', 
      'scripts/BarMap.py',
      'scripts/DrTransformer.py']
)

