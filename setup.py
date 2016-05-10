"""Build and install gridfill."""
# Copyright (c) 2012-2016 Andrew Dawson
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import setuptools
from setuptools import setup, Extension


try:
    from Cython.Distutils import build_ext
except ImportError:
    raise ImportError('Cython 0.15.1+ is required to install gridfill')
try:
    import numpy as np
except ImportError:
    raise ImportError('NumPy 1.6+ is required to install gridfill')

# Define the required dependencies:
install_requires = [
    'numpy>=1.6',
    'Cython>=0.15.1',
    'setuptools>=0.7.2',
]

# Get the library version:
for line in open('gridfill/__init__.py').readlines():
    if line.startswith('__version__'):
        exec(line.strip())

# Define packages and package data:
packages = [
    'gridfill',
    'gridfill.tests',
]
package_data = {
    'gridfill.tests': ['data/*.npy'],
}

# Define extension modules:
ext_modules = [
    #The core implemented as a Cython extension:
    Extension(
        'gridfill._gridfill',
        ['gridfill/_gridfill.pyx'],
        include_dirs=[np.get_include()],
    ),
]

classifiers = [
    'License :: OSI Approved :: MIT License',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: POSIX :: Linux',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Mathematics',
    'Topic :: Scientific/Engineering :: Physics',
    'Topic :: Scientific/Engineering :: Atmospheric Science',
]

setup(
    name='gridfill',
    version=__version__,
    author="Andrew Dawson",
    url='https://github.com/ajdawson/gridfill',
    description='Fill missing values in grids using iterative relaxation',
    license='MIT',
    install_requires=install_requires,
    packages=packages,
    package_data=package_data,
    ext_modules=ext_modules,
    cmdclass={'build_ext': build_ext},
    classifiers=classifiers,
)
