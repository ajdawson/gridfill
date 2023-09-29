gridfill
========

Fill missing values in a grid by solving Poisson's equation via an iterative
relaxation scheme.


Installation
------------

The easiest way to install is with conda_::

    conda install -c conda-forge gridfill

You can always build from source yourself once you have the necessary
dependencies installed (see the requirements section below). Checkout the
gridfill repository from Github_ and::

    cd gridfill/
    python setup.py install


Requirements
------------

For installation gridfill requires:

* numpy
* Cython
* setuptools

gridfill can also operate on `iris` cubes, which requires the iris_
package to be installed.


Developers
----------

For development the pytest package is required in order to run tests.


.. _conda: https://docs.conda.io/en/latest/
.. _Github: https://github.com/ajdawson/gridfill
.. _iris: https://scitools-iris.readthedocs.io/en/stable/
