
.. -*- mode: rst -*-
Dual-frequency precipitation Radar PYthon package (DRpy)
=====================================
|Tweet|

.. |Tweet| image:: https://img.shields.io/twitter/url/http/shields.io.svg?style=social
    :target: https://twitter.com/dopplerchase


(pronounced derpy)

This package is designed for reading hdf5 files from the NASA's Global Preicpitation Measurement mission Dual-Frequency Precipitation Radar (GPM-DPR) into xarray datasets. 

The reason for creating this package is not to reinvent the wheel (i.e. h5py works just fine), but to allow users to access benefical functions in xarray. 

To get specific, the datafiles currently supported are the level 2 DPR files (2A.DPR*). You can download them from here for free once you have an account: ftp://arthurhou.pps.eosdis.nasa.gov__ 

__ ftp://arthurhou.pps.eosdis.nasa.gov 
