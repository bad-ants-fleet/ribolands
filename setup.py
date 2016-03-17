# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='rnaworld',
    version='0.0.1',
    description='RNA folding kinetics using ViennaRNA programs',
    long_description=readme,
    author='Stefan Badelt',
    author_email='stef@tbi.univie.ac.at',
    #url='http://www.tbi.univie.ac.at/~stef/rnaworld/',
    license=license,
    #packages=find_packages(exclude=('tests', 'docs'))
    packages=['rnaworld']

    #scripts=['scripts/spatch.py'] #, 'scripts/barmap.py', 'scripts/drtransformer.py']
)

#   install_requires=['matplotlib','numpy','pandas'],
#   classifiers=['Development Status :: 4 - Beta',\
#                    'Programming Language :: Python :: 2.7',\
#                    'License :: OSI Approved :: MIT License',\
#                    'Operating System :: Linux',\
#                    'Intended Audience :: Science/Research',\
#                    'Topic :: Scientific/Engineering :: Visualization']

