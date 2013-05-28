"""Build and install script."""
import os
import sys

from numpy.distutils.core import setup, Extension


if __name__ == '__main__':
    sources = 'src/gridfill.f90'.split()
    fext = Extension(name='gridfill_f',
            sources=['src/gridfill_f.pyf']+sources)
    setup(name='gridfill',
          version='0.1',
          description='Fill missing values.',
          author='Andrew Dawson',
          author_email='dawson _at_ atm.ox.ac.uk',
          ext_modules=[fext],
          packages=['gridfill'],
          package_dir={'': 'lib'},)

