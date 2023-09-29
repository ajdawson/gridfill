.. default-role:: py:obj

.. toctree::
   :maxdepth: 2
   :hidden:

   api


Gridfill
========

Gridfill is a Python package for filling in missing values in gridded data.
Missing values are filled by solving Poisson's equation using an iterative relaxation scheme at the missing values, resulting in a smooth filling.

The core package runs on Python 3, on Linux and OSX (maybe on Windows too), and requires numpy_.
Gridfill has an optional interface for working with iris_ cubes.


How to Install
--------------

The easiest way to install is using the conda_ package manager::

    conda install -c conda-forge gridfill

If you want to install from source you will need to have setuptools_, cython_, and numpy_ installed before you start.
Once you have the dependencies installation should be as simple as entering the gridfill source directory and doing::

    python setup.py install


Getting Started
---------------

Gridfill has one core function `gridfill.fill`.
This function accepts a masked array array as input, as well as some parameters that control the numerical scheme, and returns a 2-tuple consisting of the filled in array and some convergence information:

.. code-block:: python

   import gridfill

   # Fill an array where the y-dimension is first and the x-dimension is
   # second, with the x-dimension being cyclic, using a tolerance of 0.01
   # (tolerance depends on your data values, choose appropriately):
   filled_array, converged = gridfill.fill(masked_array, 1, 0, 0.01, cyclic=True)


If you are working with `iris.cube.Cube` objects then the convenient `gridfill.fill_cube` function is probably the tool for the job:

.. code-block:: python

   import gridfill

   # With a Cube as input we don't need to give dimension numbers:
   filled_cube = gridfill.fill_cube(cube_with_missing, 0.01, cyclic=True)


Development and Contributing
----------------------------

Development is done through Github_.
Bug reports and feature requests should be made using the repository issues_ page.
Pull requests for new features and bug fixes are welcome!


.. _numpy: https://numpy.org
.. _iris: https://scitools-iris.readthedocs.io/en/stable/
.. _conda: https://docs.conda.io/en/latest/
.. _setuptools: https://setuptools.pypa.io/en/latest/
.. _cython: https://cython.org
.. _Github: https://github.com/ajdawson/gridfill
.. _issues: https://github.com/ajdawson/gridfill/issues
