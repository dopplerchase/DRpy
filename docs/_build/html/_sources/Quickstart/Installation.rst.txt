============
Installation 
============

The easiest way to install DRpy is to install it within a new anaconda environment. If you don't already have conda, please
see go install it from here: `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_

Now that you have conda, either download the zipped file from github, or use the following (if you have git)

.. code-block:: bash

    $ git clone https://github.com/dopplerchase/DRpy.git

cd into the DRpy director 

.. code-block:: bash

    $ cd ./DRpy

use conda to make the new enviroment and install the dependances 

.. code-block:: bash

    $ conda env create -f environment.yml

finally, install the package 

.. code-block:: bash

    $ pip install . 