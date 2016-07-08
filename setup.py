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
    version=version,
    description='RNA folding kinetics using ViennaRNA',
    long_description=readme,
    author='Stefan Badelt',
    author_email='stef@tbi.univie.ac.at',
    #url='http://www.tbi.univie.ac.at/~stef/ribolands/',
    license=license,
    #packages=find_packages(exclude=('tests', 'docs'))
    packages=['ribolands'],
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

