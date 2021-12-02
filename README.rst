
.. -*- mode: rst -*-
Dual-frequency precipitation Radar PYthon package (DRpy)
=====================================
|Tweet|

.. |Tweet| image:: https://img.shields.io/twitter/url/http/shields.io.svg?style=social
    :target: https://twitter.com/dopplerchase


(pronounced derpy)

**NOTE: Version 7 is coming!

Hello all, there is a new version of GPM-DPR data in the pipeline see the below messages:


===


Dear Science User,

As of November 30, 2021, PPS has ceased production of  GPM V06 DPR and Combined products. Starting on December 6, 2021 we will start GPM V07 processing beginning with Data date Dec 1, 2021 for GPM DPR ONLY. On December 10, 2021 we will start GPM (reprocessing) RP for V07 back to  the beginning of the mission for GPM DPR ONLY. At approximately 10:30 UTC 5 December 2021 until about 00:30 UTC 6 December. The GPM near realtime will be transitioning its software to product V07 radar products.  The V07 radar products have a different format from V06. 

===

With this new version there is some major changes to the file structure. But the codebase here currently does support auto determining what version you have. That being said expect bugs and issues.... 

############################################

drpy.core.GPMDPR():

This class is designed for reading hdf5 files from the NASA's Global Precipitation Measurement mission Dual-Frequency Precipitation Radar (GPM-DPR) into xarray datasets. 

drpy.core.APR():

This class is designed for reading hdf5 files from the NASA's Airborne Precipitation Radar (mainly developed in OLYMPEX, so if you are using a different campaign beware, some dataset names might be different). 

############################################

The reason for creating this package is not to reinvent the wheel (i.e. h5py works just fine), but to allow users to access beneficial functions in xarray. 

To get specific, the datafiles currently supported are the level 2 DPR files (2A.DPR*). You can download them from here for free once you have an account: ftp://arthurhou.pps.eosdis.nasa.gov__ 

__ ftp://arthurhou.pps.eosdis.nasa.gov 
