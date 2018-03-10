#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Python 3 compatibility
from __future__ import absolute_import

from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

version = 0.6

setup(
    name='ribolands',
    description='energy landscapes and folding kinetics of nucleic acids',
    long_description=readme,
    version=version,
    license=license,
    author='Stefan Badelt',
    author_email='stef@tbi.univie.ac.at',
    url='https://github.com/bad-ants-fleet/ribolands',
    install_requires=['matplotlib', 'networkx', 'pandas', 'crnsimulator>=0.1'],
    dependency_links=[
        'https://github.com/bad-ants-fleet/crnsimulator/tarball/master#egg=crnsimulator-0.1'],
    test_suite='tests',
    packages=['ribolands'],
    scripts=['scripts/BarMap.py',
             'scripts/DrTransformer.py']
)
