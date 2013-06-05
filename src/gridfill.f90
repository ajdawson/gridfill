!
! Contents
! ========
!
! Provided is a subroutine for filling missing values in a grid using an
! iterative relaxation scheme, and supporting subroutines.
!
! License
! =======
!
! Copyright (c) 2012-2013 Andrew Dawson
! 
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.


subroutine initialize_missing(nlat, nlon, grid, missing, initzonal, mask)
    implicit none
!
! Calling arguments:
    integer,                             intent(in)    :: nlat, nlon
    real(kind=8), dimension(nlat, nlon), intent(inout) :: grid
    real(kind=8),                        intent(in)    :: missing
    logical,                             intent(in)    :: initzonal
    integer,      dimension(nlat, nlon), intent(out)   :: mask
!
! Purpose
! =======
!
! Locate missing values in a grid and replace them with an initial
! value.
!
! Author
! ======
!
! Andrew Dawson <dawson _at_ atm.ox.ac.uk>
!

    ! Local variable declarations.
    integer :: i, j, n
    real(kind=8) :: zonal_mean, init

    ! Create a mask to indicate missing values and compute zonal sums.
    ! mask=1 where a point is missing and mask=0 where a point is not
    ! missing.
    misval_lat: do i = 1, nlat

        n = 0
        zonal_mean = 0.0_8
        misval_lon1: do j = 1, nlon
            if (grid(i, j) .eq. missing) then
                mask(i, j) = 1
            else
                n = n + 1
                zonal_mean = zonal_mean + grid(i, j)
                mask(i, j) = 0
            end if
        end do misval_lon1

        ! Compute the zonal mean by dividing the zonal sum by the number of
        ! non-missing points at this latitude.
        if (n .gt. 0) then
            zonal_mean = zonal_mean / real(n)
        end if

        ! Set the initial value of missing points depending on the user
        ! specified guess type. If guess is 0 use 0 as the initial value, if
        ! guess is 1 use the zonal mean.
        if (initzonal) then
            init = zonal_mean
        else
            init = 0.0_8
        end if
        misval_lon2: do j = 1, nlon
            if (grid(i, j) .eq. missing) then
                grid(i, j) = init
            end if
        end do misval_lon2

    end do misval_lat

    return
end subroutine initialize_missing


subroutine latitude_indices(i, nlat, im1, ip1)
    implicit none
!
! Calling arguments:
    integer, intent(in) :: i, nlat
    integer, intent(out) :: im1, ip1
!
! Purpose
! =======
!
! Determine indices of the latitude points adjacent to a given point
! taking into account boundaries.
!
!
! Author
! ======
!
! Andrew Dawson <dawson _at_ atm.ox.ac.uk>
!
    ! Set initial values, valid for all interior points.
    im1 = i - 1
    ip1 = i + 1

    ! Bottom boundary points.
    if (i .eq. 1) then
        im1 = 2
    end if

    ! Top boundary points.
    if (i .eq. nlat) then
        ip1 = nlat - 1
    end if

    return
end subroutine latitude_indices


subroutine longitude_indices(j, nlon, cyclic, jm1, jp1)
    implicit none
!
! Calling arguments:
    integer, intent(in)  :: j, nlon
    logical, intent(in)  :: cyclic
    integer, intent(out) :: jm1, jp1
!
! Purpose
! =======
!
! Determine indices of the longitude points adjacent to a given point
! taking into account boundaries and cyclicity.
!
!
! Author
! ======
!
! Andrew Dawson <dawson _at_ atm.ox.ac.uk>
!
    ! Set initial values, valid for all interior points.
    jm1 = j - 1
    jp1 = j + 1

    ! Left boundary points.
    if (j .eq. 1) then
        if (.not. cyclic) then
            jm1 = 2
        else
            jm1 = nlon
        end if
    end if

    ! Right boundary points.
    if (j .eq. nlon) then
        if (.not. cyclic) then
            jp1 = nlon - 1
        else
            jp1 = 1
        end if
    end if

    return
end subroutine longitude_indices


subroutine poisson_fill (nlat, nlon, grid, missing, itermax, tolerance, &
        relaxc, initzonal, cyclic, resmax, numiter)
    implicit none
!
! Calling arguments:
    integer,                             intent(in)    :: nlat
    integer,                             intent(in)    :: nlon  
    real(kind=8), dimension(nlat, nlon), intent(inout) :: grid            
    real(kind=8),                        intent(in)    :: missing
    integer,                             intent(in)    :: itermax
    real(kind=8),                        intent(in)    :: tolerance
    real(kind=8),                        intent(in)    :: relaxc
    logical,                             intent(in)    :: initzonal
    logical,                             intent(in)    :: cyclic
    real(kind=8),                        intent(out)   :: resmax
    integer,                             intent(out)   :: numiter
!
! f2py directives
    !f2py intent(in,out) :: grid
    !f2py intent(hide)   :: nlat, nlon
!
! Purpose
! =======
!
! Fill missing values in a grid by iteratively solving Poisson's
! equation.
!
!
! Arguments
! =========
!
! nlat         (input) integer
!              Number of latitudes in the input grid.
!
! nlon         (input) integer
!              Number of longitudes in the input grid.
!
! grid         (input,output) real(kind=8) dimension(nlat, nlon)
!              A grid with missing values.
!
! missing      (input) real(kind=8)
!              The missing value in the input grid.
!
! itermax      (input) integer
!              The maximum number of iterations allowed for the
!              relaxation scheme.
!
! tolerance    (input) real(kind=8)
!              Tolerance for convergence of the relaxation scheme.
!
! relaxc       (input) real(kind=8)
!              Relaxation constant, typically 0.45 <= relaxc <= 0.6.
!
! initzonal    (input) logical
!              Specifies the initial guess used for the missing values:
!                  = .TRUE. : use the zonal mean as the initial guess
!                  = .FALSE.: use zero as the initial guess
!
! cyclic       (input) logical
!              Specifies whether or not the input grid is cyclic in its
!              second dimension:
!                  = .FALSE.: grid is not cyclic
!                  = .TRUE. : grid is cyclic
!
! resmax       (output) real(kind=8)
!              Maximum residual at the end of the iterative scheme.
!
! numiter      (output) integer
!              The number of iterations used to achieve the solution.
!
!
! Author
! ======
!
! Andrew Dawson <dawson _at_ atm.ox.ac.uk>
!
    ! Local variable declarations
    integer                             :: i, j, im1, ip1, jm1, jp1
    real(kind=8)                        :: dp25, residual
    integer,      dimension(nlat, nlon) :: mask
    
    ! Define double precision quarter for finite differencing.
    dp25 = 0.25_8

    ! Create a mask to indicate missing values and compute zonal sums.
    ! mask=1 where a point is missing and mask=0 where a point is not
    ! missing.
    call initialize_missing(nlat, nlon, grid, missing, initzonal, mask)
    if (sum(mask) .eq. 0) then
        ! No missing elements so we can return now.
        resmax = 0.0
        numiter = 0
        return
    endif
    
    ! Initialize the iteration counter to zero.
    numiter = 0
    
    ! Iterate until the iteration limit is met.
    iterloop: do while (numiter .lt. itermax)
    
        ! Initialize the current maximum residual to 0 and increment the
        ! iteration counter.
        resmax = 0.0_8
        numiter  = numiter + 1
    
        ! Loop over latitude.
        iterloop_lat: do i = 1, nlat

            ! Create indices for the latitudes north and south of the current
            ! latitude.
            call latitude_indices(i, nlat, im1, ip1)
            
            ! Loop over longitude.
            iterloop_lon: do j = 1, nlon

                ! Check that this point is a missing value in the original
                ! input. If it is not missing then nothing needs to be done.
                if (mask(i, j) .eq. 1) then

                    ! Create indices for the longitudes east and west of the
                    ! current longitude.
                    call longitude_indices(j, nlon, cyclic, jm1, jp1)

                    ! Compute the Laplacian residual at the current point.
                    residual = dp25 * (grid(im1, j) + grid(ip1, j) + &
&                           grid(i, jm1) + grid(i, jp1)) - grid(i, j)

                    ! Multiply the residual by the relaxation constant.
                    residual = residual * relaxc

                    ! Add the residual at the missing point to the
                    ! input/output array.
                    grid(i, j) = grid(i, j) + residual

                    ! Set this value of the residual to be the maximum if it
                    ! is bigger than all the other computed so far.
                    resmax = max(abs(residual), resmax)
                end if
            end do iterloop_lon
        end do iterloop_lat

        ! If the maximum Laplacian value at this iteration is less than the
        ! error tolerance then stop the iteration process.
        if (resmax .le. tolerance) then
            exit
        end if

    end do iterloop

    return
end subroutine poisson_fill


subroutine poisson_fill_grids(nlat, nlon, ng, grid, missing, itermax, &
                              tolerance, relaxc, initzonal, cyclic, resmax, &
                              numiter)
    implicit none
!
! Calling arguments:
    integer,                                 intent(in)    :: nlat
    integer,                                 intent(in)    :: nlon
    integer,                                 intent(in)    :: ng
    real(kind=8), dimension(nlat, nlon, ng), intent(inout) :: grid
    real(kind=8),                            intent(in)    :: missing
    integer,                                 intent(in)    :: itermax
    real(kind=8),                            intent(in)    :: tolerance
    real(kind=8),                            intent(in)    :: relaxc
    logical,                                 intent(in)    :: initzonal
    logical,                                 intent(in)    :: cyclic
    real(kind=8), dimension(ng),             intent(out)   :: resmax
    integer, dimension(ng),                  intent(out)   :: numiter
!
! f2py directives
    !f2py intent(in,out) :: grid
    !f2py intent(hide)   :: nlat, nlon, ng
!
! Purpose
! =======
!
! Fill missing values in multiple grids by iteratively solving Poisson's
! equation.
!
!
! Arguments
! =========
!
! nlat         (input) integer
!              Number of latitudes in the input grid.
!
! nlon         (input) integer
!              Number of longitudes in the input grid.
!
! ng           (input) integer
!              Number of latitude-longitude grids in the input.
!
! grid         (input,output) real(kind=8) dimension(nlat, nlon, ng)
!              A grid with missing values.
!
! missing      (input) real(kind=8)
!              The missing value in the input grid.
!
! itermax      (input) integer
!              The maximum number of iterations allowed for the
!              relaxation scheme.
!
! tolerance    (input) real(kind=8)
!              Tolerance for convergence of the relaxation scheme.
!
! relaxc       (input) real(kind=8)
!              Relaxation constant, typically 0.45 <= relaxc <= 0.6.
!
! initzonal    (input) logical
!              Specifies the initial guess used for the missing values:
!                  = .TRUE. : use the zonal mean as the initial guess
!                  = .FALSE.: use zero as the initial guess
!
! cyclic       (input) logical
!              Specifies whether or not the input grid is cyclic in its
!              second dimension:
!                  = .FALSE.: grid is not cyclic
!                  = .TRUE. : grid is cyclic
!
! resmax       (output) real(kind=8) dimension(ng)
!              Maximum residual at the end of the iterative scheme for
!              each grid.
!
! numiter      (output) integer dimension(ng)
!              The number of iterations used to achieve the solution for
!              each grid.
!
!
! Author
! ======
!
! Andrew Dawson <dawson _at_ atm.ox.ac.uk>
!
    ! Local variable declarations
    integer :: i

    ! Fill each grid
    grid_loop: do i = 1, ng
        call poisson_fill (nlat, nlon, grid(:, :, i), missing, itermax, &
                           tolerance, relaxc, initzonal, cyclic, resmax(i), &
                           numiter(i))
    end do grid_loop

    return
end subroutine poisson_fill_grids
