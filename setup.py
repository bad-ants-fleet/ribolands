#!/usr/bin/env python

# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

setup(
    name='ribolands',
    description='nucleic acid energy landscapes and kinetic folding',
    long_description=LONG_DESCRIPTION,
    version='0.8.2',
    license='MIT',
    author='Stefan Badelt',
    author_email='stef@tbi.univie.ac.at',
    url='https://github.com/bad-ants-fleet/ribolands',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3',
        ],
    install_requires=[
        'networkx>=2.4', 
        'matplotlib', 
        'seaborn', 
        'crnsimulator>=0.8'],
    test_suite='tests',
    packages=['ribolands'],
    scripts=['scripts/BarMap.py',
             'scripts/DrTransformer.py']
)
