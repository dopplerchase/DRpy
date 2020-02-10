#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup

long_description = """A Python code for loading GPM-DPR files into xarray datasets
"""

setup(name='DRpy',
      description='Dual-frequency precipitation Radar PYthon package',
      author='Randy J. Chase',
      author_email='randy.chase12@gmail.com',
      url='https://github.com/dopplerchase/DRpy',
      packages=['drpy','drpy.core','drpy.graph','drpy.util'],
      long_description = long_description,
      license = 'MIT',
      requires = ['h5py', 'xarray']
     )
