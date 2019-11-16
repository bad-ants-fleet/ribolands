#!/usr/bin/env python3

from setuptools import setup, find_packages

LONG_DESCRIPTION="""
A package to analyze energy landscapes and folding kinetics of nucleic acids.
It provides wrapper functions for the programs `RNAsubopt`, `Kinfold`,
`barriers` and `treekin`, which have to installed separately. 

Two scripts for cotranscriptional folding are part of `ribolands`: 

  * **BarMap**: folding kinetics on dynamic energy landscapes. First, the
    thermodynamic energy landscape is coarse-grained using barriers for each
    sequence length. During kinetic simulations, occupancy is transferred from
    one landscape to the next. This script is mostly a reimplementation of
    `BarMap` by [Hofacker et al. (2010)], but it makes use of more recent
    functionality of `barriers` and `treekin`.
    ```sh
    echo "CUGCGGCUUUGGCUCUAGCC" | BarMap.py --pyplot
    ```

  * **DrTransformer**: Short for "DNA-to-RNA Transformer", the program computes
    cotranscriptional folding of RNA molecules by generating a heuristic energy
    landscape at every transcription step.  `DrTransformer` uses the `ViennaRNA
    package` to calculate free energies and transition rates between structures
    and `treekin` to simulate folding kinetics.
    ```sh
    echo "CUGCGGCUUUGGCUCUAGCC" | DrTransformer.py --visualize pdf
    ```
"""

setup(
    name = 'ribolands',
    description = 'Nucleic acid energy landscapes and folding kinetics.',
    long_description = LONG_DESCRIPTION,
    version = '0.9',
    license = 'MIT',
    author = 'Stefan Badelt',
    author_email = 'stef@tbi.univie.ac.at',
    url = 'https://github.com/bad-ants-fleet/ribolands',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3',
        ],
    install_requires = [
        'networkx>=2.4', 
        'matplotlib', 
        'seaborn', 
        'pandas', 
        'crnsimulator>=0.7'],
    test_suite = 'tests',
    packages = ['ribolands', 'ribolands.parser'],
    scripts = [
        'scripts/BarMap.py',
        'scripts/DrTransformer.py']
    )

