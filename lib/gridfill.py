"""Define a procedure for filling missing values."""
# Copyright (c) 2012 Andrew Dawson
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
from gridfill_f import poisson_fill as poisson_fill_f


def poisson_fill(grid, eps, relax=.6, itermax=100, initzonal=False,
        cyclic=False, verbose=False):
    """
    Fill missing values in a grid with values derived by solving
    Poisson's equation using a relaxation scheme.

    **Arguments:**

    *grid*
        An array with missing values to fill. This should be a
        masked array.

    *eps*
        Tolerance for ending the relaxation scheme before the maximum
        number of iterations.

    **Optional arguments:**

    *relax*
        Relaxation constant. Usually 0.45 <= *relc* <= 0.6. Defaults to
        0.6.

    *itermax*
        Maximum number of iterations of the relaxation scheme. Defaults
        to 100 iterations.

    *initzonal*
        If *True* missing values will be initialized with the zonal mean
        prior to the relaxation scheme. If *False* missing values will
        be initialized to zero. Defaults to *False*.

    *cyclic*
        Set to *True* if the rightmost dimension of *grid* is cyclic.
        Set to *False otherwise. Defaults to *False*.

    *verbose*
        If *True* information about convergence will be printed. If
        *False* nothing is printed. Defaults to *False*.

    """
    # Fill the missing values in the grid with 1.e20.
    missing = 1.e20
    grid = grid.filled(fill_value=missing)
    # Call the Fortran subroutine to fill the missing values.
    fgrid, resmax, niter = poisson_fill_f(grid, missing, itermax, eps, relax,
            initzonal, cyclic)
    converged = not (niter == itermax and resmax > eps)
    # Optionally print convergence information.
    if verbose:
        if converged:
            converged_str = 'converged'
        else:
            converged_str = 'did not converge'
        message = 'relaxation {} ({} iterations with maximum residual {})'
        print message.format(converged_str, niter, resmax)
    return fgrid, converged

