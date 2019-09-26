#!/usr/bin/env python

# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

LONG_DESCRIPTION="""
A package to analyze energy landscapes and folding kinetics of nucleic acids.
It provides wrapper functions for the programs `RNAsubopt`, `barriers` and
`treekin`, which have to installed separately. 

Two scripts for cotranscriptional folding are part of `ribolands`: 

  * **DrTransformer**: Short for "DNA-to-RNA Transformer", the program
    computes cotranscriptional folding of larger RNAs by generating a
    heuristic energy landscape at every transcription step.  `DrTransformer`
    uses the `ViennaRNA package` to calculate transition rates and `treekin` to
    simulate folding kinetics.
    ```sh
    echo "CUGCGGCUUUGGCUCUAGCC" | DrTransformer.py --visualize pdf
    ```

  * **BarMap**: folding kinetics on dynamic energy landscapes. For each
    sequence length, the coarse-grained barriers landscape is computed. During
    kinetic simulations, a mapping between subsequent landscapes is used to
    transfer occupancy from one landscape to the next. This is mostly a
    reimplementation of `BarMap` by [Hofacker et al. (2010)], but it makes use
    of more recent functionality of `barriers` and `treekin`.
    ```sh
    echo "CUGCGGCUUUGGCUCUAGCC" | BarMap.py --pyplot
    ```
"""

setup(
    name='ribolands',
    description='nucleic acid energy landscapes and kinetic folding',
    long_description=LONG_DESCRIPTION,
    version='0.8',
    license='MIT',
    author='Stefan Badelt',
    author_email='stef@tbi.univie.ac.at',
    url='https://github.com/bad-ants-fleet/ribolands',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        ],
    install_requires=[
        'networkx', 
        'matplotlib', 
        'future', 
        'seaborn', 'pandas', 
        'crnsimulator>=0.4'],
    test_suite='tests',
    packages=['ribolands'],
    scripts=['scripts/BarMap.py',
             'scripts/DrTransformer.py']
)
