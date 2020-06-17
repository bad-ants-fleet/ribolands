#!/usr/bin/env python3

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

setup(
    name = 'ribolands',
    description = 'Nucleic acid energy landscapes and folding kinetics.',
    long_description = LONG_DESCRIPTION,
    version = '0.9',
    license = 'MIT',
    author = 'Stefan Badelt',
    author_email = 'bad-ants-fleet@posteo.eu',
    url = 'https://github.com/bad-ants-fleet/ribolands',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        ],
    install_requires = [
        'networkx>=2.4', 
        'matplotlib', 
        'seaborn', 
        'pandas', 
        'crnsimulator>=0.7.1'],
    test_suite = 'tests',
    packages = ['ribolands', 'ribolands.parser'],
    entry_points = {
        'console_scripts': [
            'findpath=ribolands.pathfinder:main',
            'DrTransformer=ribolands.trafo:main'
            ],
        }
)

