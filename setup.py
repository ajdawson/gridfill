"""Build and install script."""
# Copyright (c) 2012-2013 Andrew Dawson
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
from __future__ import print_function
import sys
from numpy.distutils.core import setup, Extension, Command


for line in open('lib/gridfill/__init__.py').readlines():
    if line.startswith('__version__'):
        exec(line.strip())


class TestRunner(Command):

    description = "run the test suite"
    user_options = []
    default_test_modules = ['gridfill.tests.test_fill']

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            import nose
            nose.main(defaultTest=self.default_test_modules)
        except ImportError:
            print('the test suite requires nose... skipping', file=sys.stderr)


packages = ['gridfill']
package_data = {}
extension = Extension(name='gridfill_f',
                      sources=['src/gridfill_f.pyf', 'src/gridfill.f90'],)

setup(name='gridfill',
      version=__version__,
      description='fill missing values in grids',
      author='Andrew Dawson',
      author_email='dawson@atm.ox.ac.uk',
      url='https://github.com/ajdawson/gridfill',
      ext_modules=[extension],
      packages=packages,
      package_dir={'': 'lib'},
      package_data=package_data,
      cmdclass={'test': TestRunner},)
