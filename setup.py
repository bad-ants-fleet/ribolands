#!/usr/bin/env python

from setuptools import setup

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

setup(
    name = 'ribolands',
    version = '0.9',
    description = 'Nucleic acid energy landscapes and folding kinetics.',
    long_description = LONG_DESCRIPTION,
    long_description_content_type = 'text/markdown',
    author = 'Stefan Badelt',
    author_email = 'bad-ants-fleet@posteo.eu',
    maintainer = 'Stefan Badelt',
    maintainer_email = 'bad-ants-fleet@posteo.eu',
    url = 'https://github.com/bad-ants-fleet/ribolands',
    license = 'MIT',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Intended Audience :: Science/Research',
        ],
    python_requires = '>=3.8',
    install_requires = [
        'networkx>=2.4', 
        'matplotlib', 
        'seaborn', 
        'natsort',
        'packaging',
        'pandas', 
        'crnsimulator>=0.7.1'],
    packages = ['ribolands', 'ribolands.parser'],
    test_suite = 'tests'
)

