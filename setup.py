#! /usr/bin/env python

from setuptools import setup

setup(name='sbmlodex',
      version='1.0.2',
      description='SBML ODE extractor',
      author='J Kyle Medley',
      url='https://github.com/0u812/sbmlodex',
      packages=['sbmlodex'],
      install_requires=['tesbml>=5.15.0'])
