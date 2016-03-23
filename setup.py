# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

# Dynamically figure out the version
version = __import__('rnaworld').__version__

setup(
    name='rnaworld',
    version=version,
    description='RNA folding kinetics using ViennaRNA programs',
    long_description=readme,
    author='Stefan Badelt',
    author_email='stef@tbi.univie.ac.at',
    #url='http://www.tbi.univie.ac.at/~stef/rnaworld/',
    license=license,
    #packages=find_packages(exclude=('tests', 'docs'))
    packages=['rnaworld'],
    scripts=['scripts/pylands_spatch.py',
      'scripts/cofold_kinetics.py', 
      'scripts/BarMap.py',
      'scripts/DrTransformer.py']
)

#   install_requires=['matplotlib','numpy','pandas'],
#   classifiers=['Development Status :: 4 - Beta',\
#                    'Programming Language :: Python :: 2.7',\
#                    'License :: OSI Approved :: MIT License',\
#                    'Operating System :: Linux',\
#                    'Intended Audience :: Science/Research',\
#                    'Topic :: Scientific/Engineering :: Visualization']

