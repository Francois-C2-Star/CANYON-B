!* 
!* SPDX-FileCopyrightText: Copyright (c) 2025 CUDASOL CONSEIL. All rights reserved.
!* SPDX-License-Identifier: MIT
!* 
!* Permission is hereby granted, free of charge, to any person obtaining a copy
!* of this software and associated documentation files (the "Software"), to deal
!* in the Software without restriction, including without limitation the rights
!* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!* copies of the Software, and to permit persons to whom the Software is
!* furnished to do so, subject to the following conditions:
!* 
!* The above copyright notice and this permission notice shall be included in all
!* copies or substantial portions of the Software.
!* 
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!* SOFTWARE.
!*
!**********************************************************************************************************************************************************************************************************
!Fortran90 implementation of the CANYON-B method published in:
!
!Bittig et al. (2018). An alternative to static climatologies: Robust estimation of open ocean CO2 variables and nutrient concentrations from T, S and O2 data using Bayesian neural networks. 
! Front. Mar. Sci. 5:328. http://dx.doi.org/10.3389/fmars.2018.00328.
!
!When using the method, please cite the paper.
!**********************************************************************************************************************************************************************************************************
 
 
program canyonb_main
  use, intrinsic :: iso_fortran_env, only: wp => real64
  implicit none
  ! Required modules
  use canyonb_module
  implicit none

  ! Test driver for CANYONB
  integer, parameter :: ntest = 2
  character(len=20), dimension(ntest) :: date
  real(wp), dimension(ntest) :: lat, lon, pres, temp, psal, doxy
  character(len=4), allocatable, dimension(:) :: param
  type(canyonb_output) :: out

  integer :: i

  ! Initialize test values (as in default test in R)
  do i = 1, ntest
     date(i) = "2014-12-09 08:45"
  end do
  lat = (/ 17.6_wp, 17.6_wp /)
  lon = (/ -24.3_wp, -24.3_wp /)
  pres = (/ 180_wp, 180_wp /)
  temp = (/ 16_wp, 16_wp /)
  psal = (/ 36.1_wp, 36.1_wp /)
  doxy = (/ 104_wp, 104_wp /)
  ! Example: only predict NO3, PO4 and CT when desired (here we predict all by leaving param not present)
  call CANYONB(date, lat, lon, pres, temp, psal, doxy, out=out)

  ! Print one of the outputs as an example: NO3 predictions (if computed)
  if (associated(out.NO3)) then
     print *, "NO3 predictions:"
     do i = 1, size(out.NO3)
        print *, out.NO3(i), out.NO3_ci(i), out.NO3_cim(i), out.NO3_cin(i), out.NO3_cii(i)
     end do
  else
     print *, "NO3 not computed."
  end if

end program canyonb_main

module canyonb_module
  use, intrinsic :: iso_fortran_env, only: wp => real64
  use ieee_arithmetic
  implicit none

  ! Derived type for output structure
  type canyonb_output
     real(wp), allocatable, dimension(:) :: AT, AT_ci, AT_cim, AT_cin, AT_cii
     real(wp), allocatable, dimension(:) :: CT, CT_ci, CT_cim, CT_cin, CT_cii
     real(wp), allocatable, dimension(:) :: pH, pH_ci, pH_cim, pH_cin, pH_cii
     real(wp), allocatable, dimension(:) :: pCO2, pCO2_ci, pCO2_cim, pCO2_cin, pCO2_cii
     real(wp), allocatable, dimension(:) :: NO3, NO3_ci, NO3_cim, NO3_cin, NO3_cii
     real(wp), allocatable, dimension(:) :: PO4, PO4_ci, PO4_cim, PO4_cin, PO4_cii
     real(wp), allocatable, dimension(:) :: SiOH4, SiOH4_ci, SiOH4_cim, SiOH4_cin, SiOH4_cii
  end type canyonb_output

contains

  subroutine CANYONB(date, lat, lon, pres, temp, psal, doxy, param, epres, etemp, epsal, edoxy, out)
    implicit none
    ! Input arguments
    character(len=*), dimension(:), intent(in), optional :: date
    real(wp), dimension(:), intent(inout), optional :: lat, lon, pres, temp, psal, doxy
    character(len=*), dimension(:), intent(in), optional :: param
    real(wp), dimension(:), intent(in), optional :: epres, etemp, epsal, edoxy
    ! Output argument
    type(canyonb_output), intent(out) :: out

    ! Local variables declaration
    integer :: nol, noparams, i, j, k, l, ni, nlayerflag, nl1, nl2, noparsets
    real(wp), allocatable, dimension(:) :: epres_arr, etemp_arr, epsal_arr, edoxy_arr
    real(wp) :: default_epres, default_etemp, default_epsal
    real(wp), allocatable, dimension(:) :: temp_doxy
    character(len=20), allocatable, dimension(:) :: date_arr

    character(len=20) :: inputsdir
    inputsdir = ""  ! relative or absolute path to CANYON-B wgts files

    ! Test values block if missing(date)
    if (.not. present(date)) then
       allocate(date_arr(2))
       date_arr = ("2014-12-09 08:45", "2014-12-09 08:45")
       if (.not. present(lat)) then
          allocate(lat(2))
          lat = 17.6_wp
       end if
       if (.not. present(lon)) then
          allocate(lon(2))
          lon = -24.3_wp
       end if
       if (.not. present(pres)) then
          allocate(pres(2))
          pres = 180_wp
       end if
       if (.not. present(temp)) then
          allocate(temp(2))
          temp = 16_wp
       end if
       if (.not. present(psal)) then
          allocate(psal(2))
          psal = 36.1_wp
       end if
       if (.not. present(doxy)) then
          allocate(doxy(2))
          doxy = 104_wp
       end if
    else
       date_arr = date
    end if

    nol = size(pres)

    ! Set default input error values if missing
    default_epres = 0.5_wp
    default_etemp = 0.005_wp
    default_epsal = 0.005_wp
    ! For edoxy, default to 0.01*doxy

    if (present(epres)) then
       epres_arr = epres
    else
       allocate(epres_arr(nol))
       epres_arr = default_epres
    end if
    if (present(etemp)) then
       etemp_arr = etemp
    else
       allocate(etemp_arr(nol))
       etemp_arr = default_etemp
    end if
    if (present(epsal)) then
       epsal_arr = epsal
    else
       allocate(epsal_arr(nol))
       epsal_arr = default_epsal
    end if
    if (present(edoxy)) then
       edoxy_arr = edoxy
    else
       allocate(edoxy_arr(nol))
       edoxy_arr = 0.01_wp * doxy
    end if

    ! Expand input errors if scalar
    if (size(epres_arr) == 1) then
       epres_arr = reshape( [ (epres_arr(1), i=1,nol) ], [nol] )
    end if
    if (size(etemp_arr) == 1) then
       etemp_arr = reshape( [ (etemp_arr(1), i=1,nol) ], [nol] )
    end if
    if (size(epsal_arr) == 1) then
       epsal_arr = reshape( [ (epsal_arr(1), i=1,nol) ], [nol] )
    end if
    if (size(edoxy_arr) == 1) then
       edoxy_arr = reshape( [ (edoxy_arr(1), i=1,nol) ], [nol] )
    end if

    ! Define parameters
    noparams = 7
    character(len=6), dimension(noparams) :: paramnames = (/"AT", "CT", "pH", "pCO2", "NO3", "PO4", "SiOH4"/)
    real(wp), dimension(noparams) :: inputsigma
    inputsigma = (/ 6.0_wp, 4.0_wp, 0.005_wp, ieee_value(0.0_wp, kind=0)/, 2.0_wp/100.0_wp, 2.0_wp/100.0_wp, 2.0_wp/100.0_wp /)
    ! Orr systematic uncertainty for pH measurement
    inputsigma(3) = sqrt(0.005_wp**2 + 0.01_wp**2)
    real(wp), dimension(3) :: betaipCO2
    betaipCO2 = (/ -3.114e-05_wp, 1.087e-01_wp, -7.899e+01_wp /)

    logical, dimension(noparams) :: paramflag
    if (.not. present(param)) then
       paramflag = .true.
    else
       ! Check each parameter of paramnames exists in input param array
       paramflag = .false.
       do i = 1, noparams
          do j = 1, size(param)
             if (trim(param(j)) == trim(paramnames(i))) then
                paramflag(i) = .true.
             end if
          end do
       end do
    end if

    ! pCO2 requires seacarb; we assume carb and derivnum routines exist
    ! (No explicit action needed in Fortran)

    ! Adjust longitude if > 180: convert to [-180,180]
    do i = 1, size(lon)
       if (lon(i) > 180.0_wp) then
          lon(i) = lon(i) - 360.0_wp
       end if
    end do

    ! Convert date strings to decimal year
    real(wp), allocatable, dimension(:) :: yearvec
    allocate(yearvec(nol))
    do i = 1, nol
       yearvec(i) = parse_date_to_dec_year(date_arr(i))
    end do

    ! Input sequence: columns = year, lat/90, |1 - ((lon-110) mod 360)/180|, |1- ((lon-20) mod 360)/180|, temp, psal, doxy, pres/2e4+1/((1+exp(-pres/300))^3)
    integer :: ncol_data
    ncol_data = 8
    real(wp), allocatable, dimension(:,:) :: data
    allocate(data(nol, ncol_data))
    do i = 1, nol
       data(i,1) = yearvec(i)
       data(i,2) = lat(i) / 90.0_wp
       data(i,3) = abs( 1.0_wp - mod((lon(i)-110.0_wp), 360.0_wp)/180.0_wp )
       data(i,4) = abs( 1.0_wp - mod((lon(i)-20.0_wp), 360.0_wp)/180.0_wp )
       data(i,5) = temp(i)
       data(i,6) = psal(i)
       data(i,7) = doxy(i)
       data(i,8) = pres(i)/20000.0_wp + 1.0_wp/((1.0_wp+exp(-pres(i)/300.0_wp))**3)
    end do

    integer, parameter :: no_out = 1

    ! Polar shift for Latitude: new lat with Polar shift ("Bering Strait prolongation")
    ! Define polygon points for Arctic basin 'West' of Lomonossov ridge
    real(wp), dimension(12) :: plon = (/ -180.0_wp, -170.0_wp, -85.0_wp, -80.0_wp, -37.0_wp, -37.0_wp, 143.0_wp, 143.0_wp, 180.0_wp, 180.0_wp, -180.0_wp, -180.0_wp /)
    real(wp), dimension(12) :: plat = (/ 68.0_wp, 66.5_wp, 66.5_wp, 80.0_wp, 80.0_wp, 90.0_wp, 90.0_wp, 68.0_wp, 68.0_wp, 90.0_wp, 90.0_wp, 68.0_wp /)
    logical, allocatable, dimension(:) :: arcflag
    allocate(arcflag(nol))
    call in_polygon(plat, plon, lat, lon, arcflag)
    ! Modify lat for points inside the polygon: add distance from -37° line
    do i = 1, nol
       if (arcflag(i)) then
          lat(i) = lat(i) - sin(3.14159265358979323846_wp*((lon(i)+37.0_wp)/180.0_wp))*(90.0_wp-lat(i))*0.5_wp
          ! Ensure lat stays within [-90,90]
          if (lat(i) > 90.0_wp) lat(i) = 90.0_wp
          if (lat(i) < -90.0_wp) lat(i) = -90.0_wp
       end if
    end do
    deallocate(plon, plat, arcflag)

    ! Preallocate output arrays for all parameters (using size nol)
    call allocate_output(out, nol)

    ! Cycle over all CANYON-B variables
    real(wp), allocatable, dimension(:,:) :: data_N
    real(wp), allocatable, dimension(:, :) :: inwgts
    real(wp), allocatable, dimension(:,:) :: mw, sw
    real(wp) :: ioffset
    real(wp), allocatable, dimension(:,:) :: cval
    real(wp), allocatable, dimension(:) :: cvalcy
    real(wp), allocatable, dimension(:,:,:) :: inval

    real(wp) :: V1, V2
    real(wp), allocatable, dimension(:) :: wgts
    real(wp), allocatable, dimension(:) :: betaciw
    real(wp), allocatable, dimension(:) :: cvalcu, cvalcib, cvalciw, cvalcin
    real(wp), allocatable, dimension(:,:) :: inx
    real(wp) :: ddp
    real(wp), allocatable, dimension(:) :: cvalcimeas
    real(wp), allocatable, dimension(:) :: tmp_out

    do i = 1, noparams
       if (paramflag(i)) then
          ! Load weights etc. from file: file name = inputsdir//"wgts_"//paramnames(i)//".txt"
          call read_table(trim(inputsdir)//"wgts_"//trim(paramnames(i))//".txt", inwgts)
          noparsets = size(inwgts,2) - 1

          ! Input normalization: determine ni and ioffset
          if (i > 4) then
             ni = ncol_data() - 1   ! ncol_data from data is 8 so ni = 7
             ioffset = -1.0_wp
             allocate(mw(ni+1,1))
             allocate(sw(ni+1,1))
             mw(:,1) = inwgts(1:ni+1, size(inwgts,2))
             sw(:,1) = inwgts(ni+2:2*ni+2, size(inwgts,2))
             allocate(data_N(nol, ni))
             do j = 1, nol
                data_N(j,:) = (data(j,2:ni+1) - mw(1:ni,1)') / (sw(1:ni,1)')
             end do
          else
             ni = ncol_data()  ! ni = 8
             ioffset = 0.0_wp
             allocate(mw(ni+1,1))
             allocate(sw(ni+1,1))
             mw(:,1) = inwgts(1:ni+1, size(inwgts,2))
             sw(:,1) = inwgts(ni+2:2*ni+2, size(inwgts,2))
             allocate(data_N(nol, ni))
             do j = 1, nol
                data_N(j,:) = (data(j,1:ni) - mw(1:ni,1)') / (sw(1:ni,1)')
             end do
          end if

          ! Extract wgts and betaciw from inwgts
          allocate(wgts(noparsets))
          do l = 1, noparsets
             wgts(l) = inwgts(4, l)
          end do
          ! Extract betaciw from rows (2*ni+3) to end, column (noparsets+1)
          allocate(betaciw(size(inwgts,1) - (2*ni+2)))
          {
             integer :: idx, count
             count = 0
             do idx = 2*ni+3, size(inwgts,1)
                if (.not. ieee_is_nan(inwgts(idx, noparsets+1))) then
                   count = count + 1
                   betaciw(count) = inwgts(idx, noparsets+1)
                end if
             end do
             if (count < size(betaciw)) then
                betaciw = betaciw(1:count)
             end if
          }

          ! Preallocate cval, cvalcy, inval
          allocate(cval(nol, noparsets))
          allocate(cvalcy(noparsets))
          allocate(inval(nol, ni, noparsets))
          cval = ieee_value(0.0_wp, kind=0)
          cvalcy = ieee_value(0.0_wp, kind=0)

          ! Cycle all networks in committee
          do l = 1, noparsets
             ! Determine network architecture: nlayerflag=1+ (inwgts(2,l) /= 0)
             if (inwgts(2,l) /= 0.0_wp) then
                nlayerflag = 2
             else
                nlayerflag = 1
             end if
             nl1 = int(inwgts(1,l))
             nl2 = int(inwgts(2,l))
             ! beta value
             ! beta not used further except in noise variance calculation
             ! Extract weights and biases
             ! For one hidden layer:
             ! w1: matrix from rows 5 to 4+nl1*ni, reshape to (nl1, ni)
             real(wp), allocatable, dimension(:,:) :: w1, w2, w3
             real(wp), allocatable, dimension(:) :: b1, b2, b3
             integer :: start_idx, end_idx
             allocate(w1(nl1, ni))
             start_idx = 5
             end_idx = 4 + nl1*ni
             call reshape_column(inwgts, l, start_idx, end_idx, nl1, ni, w1)
             allocate(b1(nl1))
             start_idx = 4 + nl1*ni + 1
             end_idx = 4 + nl1*ni + nl1
             call get_vector(inwgts, l, start_idx, end_idx, b1)
             allocate(w2(1, nl1))
             if (nlayerflag == 1) then
                allocate(w2(1, nl1))
                start_idx = 4 + nl1*(1+ni) + 1
                end_idx = 4 + nl1*(1+ni) + nl1*1
                call reshape_column(inwgts, l, start_idx, end_idx, 1, nl1, w2)
                allocate(b2(1))
                start_idx = 4 + nl1*(1+ni) + nl1*1 + 1
                end_idx = 4 + nl1*(1+ni) + nl1*1 + 1
                call get_vector(inwgts, l, start_idx, end_idx, b2)
             else if (nlayerflag == 2) then
                allocate(w2(nl2, nl1))
                start_idx = 4 + nl1*(1+ni) + 1
                end_idx = 4 + nl1*(1+ni) + nl2*nl1
                call reshape_column(inwgts, l, start_idx, end_idx, nl2, nl1, w2)
                allocate(b2(nl2))
                start_idx = 4 + nl1*(1+ni) + nl2*nl1 + 1
                end_idx = 4 + nl1*(1+ni) + nl2*nl1 + nl2
                call get_vector(inwgts, l, start_idx, end_idx, b2)
                allocate(w3(1, nl2))
                start_idx = 4 + nl1*(1+ni) + nl2*(nl1+1) + 1
                end_idx = 4 + nl1*(1+ni) + nl2*(nl1+1) + 1*nl2
                call reshape_column(inwgts, l, start_idx, end_idx, 1, nl2, w3)
                allocate(b3(1))
                start_idx = 4 + nl1*(1+ni) + nl2*(nl1+1) + 1*nl2 + 1
                end_idx = 4 + nl1*(1+ni) + nl2*(nl1+1) + 1*nl2 + 1
                call get_vector(inwgts, l, start_idx, end_idx, b3)
             end if

             ! Network forward propagation
             real(wp), allocatable, dimension(:,:) :: a, b, y
             allocate(a(nol, nl1))
             if (nlayerflag == 1) then
                ! a = data_N * transpose(w1) + replicate(b1)
                do j = 1, nol
                   a(j,:) = matmul(w1, transpose(data_N(j,:)))' + b1
                end do
                ! y = tanh(a) * transpose(w2) + replicate(b2)
                allocate(y(nol, 1))
                do j = 1, nol
                   y(j,1) = sum( tanh(a(j,:)) * w2(1,:) ) + b2(1)
                end do
             else if (nlayerflag == 2) then
                allocate(a(nol, nl1))
                do j = 1, nol
                   a(j,:) = matmul(w1, transpose(data_N(j,:)))' + b1
                end do
                allocate(b(nol, nl2))
                do j = 1, nol
                   b(j,:) = tanh(a(j,:))  ! then times w2' and add b2
                   b(j,:) = matmul(w2, transpose(b(j,:)))' + b2
                end do
                allocate(y(nol, 1))
                do j = 1, nol
                   y(j,1) = sum( tanh(b(j,:)) * w3(1,:) ) + b3(1)
                end do
             end if

             cval(:,l) = y(:,1)
             cvalcy(l) = 1.0_wp / inwgts(3,l)   ! beta from row 3

             ! Add input effects derivative calculation
             real(wp), allocatable, dimension(:,:,: ) :: x1
             allocate(x1(nol, nl1, ni))
             do j = 1, nol
                do k = 1, nl1
                   do l = 1, ni
                      x1(j,k,l) = w1(k,l) * (1.0_wp - tanh(a(j,k))**2)
                   end do
                end do
             end do

             if (nlayerflag == 1) then
                allocate(inx(nol, ni))
                inx = 0.0_wp
                do j = 1, nol
                   do k = 1, ni
                      do l = 1, nl1
                         inx(j,k) = inx(j,k) + w2(1,l) * x1(j,l,k)
                      end do
                   end do
                end do
             else if (nlayerflag == 2) then
                real(wp), allocatable, dimension(:,:,: ) :: x2
                allocate(x2(nol, nl2, nl1))
                do j = 1, nol
                   do k = 1, nl2
                      do l = 1, nl1
                         x2(j,k,l) = w2(k,l) * (1.0_wp - tanh(b(j,k))**2)
                      end do
                   end do
                end do
                allocate(inx(nol, ni))
                inx = 0.0_wp
                do j = 1, nol
                   do k = 1, ni
                      do l = 1, nl2
                         do m = 1, nl1
                            inx(j,k) = inx(j,k) + w3(1,l)* x2(j,l,m) * x1(j,m,k)
                         end do
                      end do
                   end do
                end do
                deallocate(x2)
             end if
             inval(:,:,l) = inx
             deallocate(w1, b1, w2, b2)
             if (nlayerflag == 2) then
                deallocate(w3, b3)
             end if
             deallocate(a)
             if (nlayerflag == 2) deallocate(b)
             deallocate(y)
             deallocate(x1)
             deallocate(inx)
          end do  ! End of committee networks loop

          ! Denormalization of network output
          cval = cval * sw(ni+1,1) + mw(ni+1,1)
          cvalcy = cvalcy * sw(ni+1,1)**2

          ! Add committee of all networks as evidence-weighted mean
          V1 = sum(wgts)
          V2 = sum(wgts**2)
          allocate(tmp_out(nol))
          tmp_out = 0.0_wp
          do j = 1, nol
             do l = 1, noparsets
                tmp_out(j) = tmp_out(j) + wgts(l) * cval(j,l)
             end do
             tmp_out(j) = tmp_out(j) / V1
          end do

          ! Store weighted mean output in appropriate field
          call assign_output(paramnames(i), tmp_out, out)

          ! Compute committee uncertainty variance
          allocate(cvalcu(nol))
          do j = 1, nol
             cvalcu(j) = 0.0_wp
             do l = 1, noparsets
                cvalcu(j) = cvalcu(j) + wgts(l) * (cval(j,l) - tmp_out(j))**2
             end do
             cvalcu(j) = cvalcu(j) / (V1 - V2/V1)
          end do
          deallocate(tmp_out)

          ! Weigh noise and weight uncertainty variance
          allocate(cvalcib(nol))
          cvalcib = sum(wgts*cvalcy)/V1
          allocate(cvalciw(nol))
          do j = 1, nol
             cvalciw(j) = ( betaciw(2) + betaciw(1)*sqrt(cvalcu(j)) )**2
          end do

          ! Weigh input effects
          allocate(inx(nol, ni))
          inx = 0.0_wp
          do j = 1, nol
             do k = 1, ni
                do l = 1, noparsets
                   inx(j,k) = inx(j,k) + wgts(l) * inval(j,k,l)
                end do
                inx(j,k) = inx(j,k) / V1
             end do
          end do

          ! Rescale for normalization inside the MLP
          do j = 1, nol
             do k = 1, ni
                inx(j,k) = inx(j,k) * ( sw(ni+1,1) / sw(k,1) )
             end do
          end do

          ! Additional pressure scaling
          allocate(tmp_out(nol))
          do j = 1, nol
             ddp = 1.0_wp/20000.0_wp + (1.0_wp/((1.0_wp+exp(-pres(j)/300.0_wp))**4)) * exp(-pres(j)/300.0_wp)/100.0_wp
             inx(j, int(8.0_wp+ioffset)) = inx(j, int(8.0_wp+ioffset)) * ddp
          end do
          deallocate(tmp_out)

          ! Compute input variance
          allocate(cvalcin(nol))
          cvalcin = 0.0_wp
          do j = 1, nol
             ! the indices (5:8)+ioffset: if ioffset==0 then columns 5 to 8; if -1 then columns 4 to 7
             integer :: col_start, col_end
             col_start = 5 + int(ioffset)
             col_end = 8 + int(ioffset)
             do k = col_start, col_end
                select case (k - (col_start) + 1)
                case (1)
                   cvalcin(j) = cvalcin(j) + ( inx(j,k) * etemp_arr(j) )**2
                case (2)
                   cvalcin(j) = cvalcin(j) + ( inx(j,k) * epsal_arr(j) )**2
                case (3)
                   cvalcin(j) = cvalcin(j) + ( inx(j,k) * edoxy_arr(j) )**2
                case (4)
                   cvalcin(j) = cvalcin(j) + ( inx(j,k) * epres_arr(j) )**2
                end select
             end do
          end do
          deallocate(inx)

          ! Reference / measurement uncertainty
          allocate(cvalcimeas(nol))
          if (i > 4) then
             do j = 1, nol
                cvalcimeas(j) = (inputsigma(i) * get_output_value(paramnames(i), out, j))**2
             end do
          else
             if (i == 4) then
                do j = 1, nol
                   cvalcimeas(j) = ( betaipCO2(1)* get_output_value(paramnames(i), out, j)**2 + betaipCO2(2)* get_output_value(paramnames(i), out, j) + betaipCO2(3) )**2
                end do
             else
                cvalcimeas = inputsigma(i)**2
             end if
          end if

          ! Variable uncertainty: combine individual variances and take square root
          call assign_uncertainties(paramnames(i), cvalcimeas, cvalcib, cvalciw, cvalcu, cvalcin, out)

          ! For pCO2: recalculate using carb and derivnum subroutines
          if (i == 4) then
             type(carb_out) :: outcalc, outderiv
             real(wp), allocatable, dimension(:) :: orig_val
             allocate(orig_val(nol))
             do j = 1, nol
                orig_val(j) = get_output_value(paramnames(i), out, j)
             end do
             ! Call carb: parameters fixed as in R call
             call carb(flag=15, var1=0.0023_wp, var2=orig_val*1.0e-6_wp, S=35.0_wp, T=25.0_wp, P=0.0_wp, Patm=1.0_wp, Pt=0.0_wp, &
                  Sit=0.0_wp, pHscale="T", kf="pf", k1k2="l", ks="d", b="u74", gas="standard", result=outcalc)
             call derivnum(var="var2", flag=15, var1=0.0023_wp, var2=orig_val*1.0e-6_wp, S=35.0_wp, T=25.0_wp, P=0.0_wp, Patm=1.0_wp, Pt=0.0_wp, &
                  Sit=0.0_wp, pHscale="T", kf="pf", k1k2="l", ks="d", b="u74", gas="standard", result=outderiv)
             ! Recalculate pCO2 and its uncertainties
             do j = 1, nol
                call set_output_value("pCO2", j, outcalc%pCO2)
                call set_output_value("pCO2_ci", j, outderiv%pCO2*1.0e-6_wp * get_output_value("pCO2_ci", out, j))
                call set_output_value("pCO2_cim", j, outderiv%pCO2*1.0e-6_wp * get_output_value("pCO2_cim", out, j))
                call set_output_value("pCO2_cin", j, outderiv%pCO2*1.0e-6_wp * get_output_value("pCO2_cin", out, j))
                call set_output_value("pCO2_cii", j, outderiv%pCO2*1.0e-6_wp * get_output_value("pCO2_cii", out, j))
             end do
             deallocate(orig_val)
          end if

          ! Clean up allocations for this parameter iteration
          if (allocated(mw)) deallocate(mw)
          if (allocated(sw)) deallocate(sw)
          if (allocated(data_N)) deallocate(data_N)
          if (allocated(inwgts)) deallocate(inwgts)
          if (allocated(wgts)) deallocate(wgts)
          if (allocated(betaciw)) deallocate(betaciw)
          if (allocated(cval)) deallocate(cval)
          if (allocated(cvalcy)) deallocate(cvalcy)
          if (allocated(inval)) deallocate(inval)
          if (allocated(cvalcu)) deallocate(cvalcu)
          if (allocated(cvalcib)) deallocate(cvalcib)
          if (allocated(cvalciw)) deallocate(cvalciw)
          if (allocated(cvalcin)) deallocate(cvalcin)
          if (allocated(cvalcimeas)) deallocate(cvalcimeas)
        end if  ! if paramflag(i)
    end do  ! end for each parameter

  end subroutine CANYONB

  ! Helper function: parse_date_to_dec_year converts date string "yyyy-mm-dd HH:MM" to decimal year
  real(wp) function parse_date_to_dec_year(datestr)
    implicit none
    character(len=*), intent(in) :: datestr
    integer :: yr, mon, dy, hr, mi, doy
    real(wp) :: frac
    call parse_date(datestr, yr, mon, dy, hr, mi)
    doy = day_of_year(yr, mon, dy)
    frac = real(doy - 1,wp)/365.0_wp + (hr/24.0_wp + mi/1440.0_wp)/365.0_wp
    parse_date_to_dec_year = real(yr,wp) + frac
  end function parse_date_to_dec_year

  ! Parse date string into components (very simple implementation)
  subroutine parse_date(datestr, yr, mon, dy, hr, mi)
    implicit none
    character(len=*), intent(in) :: datestr
    integer, intent(out) :: yr, mon, dy, hr, mi
    read(datestr(1:4),*) yr
    read(datestr(6:7),*) mon
    read(datestr(9:10),*) dy
    read(datestr(12:13),*) hr
    read(datestr(15:16),*) mi
  end subroutine parse_date

  ! Calculate day of year
  integer function day_of_year(yr, mon, dy)
    implicit none
    integer, intent(in) :: yr, mon, dy
    integer, dimension(12) :: mdays = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
    integer :: i
    day_of_year = 0
    if ((mod(yr,4)==0 .and. mod(yr,100) /= 0) .or. mod(yr,400)==0) then
       mdays(2) = 29
    end if
    do i = 1, mon-1
       day_of_year = day_of_year + mdays(i)
    end do
    day_of_year = day_of_year + dy
  end function day_of_year

  ! in_polygon: determine if points (lat, lon) are inside polygon defined by (plat, plon)
  subroutine in_polygon(plat, plon, latp, lonp, flag)
    implicit none
    real(wp), intent(in), dimension(:) :: plat, plon
    real(wp), intent(in), dimension(:) :: latp, lonp
    logical, intent(out), dimension(:) :: flag
    integer :: i, j, np, npoints
    logical :: c
    real(wp) :: xinters
    np = size(plat)
    npoints = size(latp)
    do i = 1, npoints
       c = .false.
       j = np
       do j = np, 1, -1
          if (((plat(j) > latp(i)) /= (plat(modulo(j,np)+1) > latp(i))) .and. &
              (lonp(i) < (plon(modulo(j,np)+1)-plon(j))*(latp(i)-plat(j))/(plat(modulo(j,np)+1)-plat(j)+1.0e-12_wp) + plon(j))) then
             c = .not. c
          end if
          ! Note: This is a simplified ray-casting method.
       end do
       flag(i) = c
    end do
  end subroutine in_polygon

  ! Utility: allocate output arrays in the output structure
  subroutine allocate_output(out, n)
    implicit none
    type(canyonb_output), intent(inout) :: out
    integer, intent(in) :: n
    allocate(out%AT(n), out%AT_ci(n), out%AT_cim(n), out%AT_cin(n), out%AT_cii(n))
    allocate(out%CT(n), out%CT_ci(n), out%CT_cim(n), out%CT_cin(n), out%CT_cii(n))
    allocate(out%pH(n), out%pH_ci(n), out%pH_cim(n), out%pH_cin(n), out%pH_cii(n))
    allocate(out%pCO2(n), out%pCO2_ci(n), out%pCO2_cim(n), out%pCO2_cin(n), out%pCO2_cii(n))
    allocate(out%NO3(n), out%NO3_ci(n), out%NO3_cim(n), out%NO3_cin(n), out%NO3_cii(n))
    allocate(out%PO4(n), out%PO4_ci(n), out%PO4_cim(n), out%PO4_cin(n), out%PO4_cii(n))
    allocate(out%SiOH4(n), out%SiOH4_ci(n), out%SiOH4_cim(n), out%SiOH4_cin(n), out%SiOH4_cii(n))
  end subroutine allocate_output

  ! Utility: assign output value to field based on parameter name
  subroutine assign_output(pname, values, out)
    implicit none
    character(len=*), intent(in) :: pname
    real(wp), dimension(:), intent(in) :: values
    type(canyonb_output), intent(inout) :: out
    if (trim(pname) == "AT") then
       out%AT = values
    else if (trim(pname) == "CT") then
       out%CT = values
    else if (trim(pname) == "pH") then
       out%pH = values
    else if (trim(pname) == "pCO2") then
       out%pCO2 = values
    else if (trim(pname) == "NO3") then
       out%NO3 = values
    else if (trim(pname) == "PO4") then
       out%PO4 = values
    else if (trim(pname) == "SiOH4") then
       out%SiOH4 = values
    end if
  end subroutine assign_output

  ! Utility: retrieve output value from structure
  function get_output_value(pname, out, idx) result(val)
    implicit none
    character(len=*), intent(in) :: pname
    type(canyonb_output), intent(in) :: out
    integer, intent(in) :: idx
    real(wp) :: val
    if (trim(pname) == "AT") then
       val = out%AT(idx)
    else if (trim(pname) == "CT") then
       val = out%CT(idx)
    else if (trim(pname) == "pH") then
       val = out%pH(idx)
    else if (trim(pname) == "pCO2") then
       val = out%pCO2(idx)
    else if (trim(pname) == "NO3") then
       val = out%NO3(idx)
    else if (trim(pname) == "PO4") then
       val = out%PO4(idx)
    else if (trim(pname) == "SiOH4") then
       val = out%SiOH4(idx)
    else
       val = 0.0_wp
    end if
  end function get_output_value

  ! Utility: set output value in structure
  subroutine set_output_value(pname, idx, value)
    implicit none
    character(len=*), intent(in) :: pname
    integer, intent(in) :: idx
    real(wp), intent(in) :: value
    ! This subroutine sets both main output and uncertainties based on pname suffix
    if (trim(pname) == "pCO2") then
       ! If setting main value for pCO2
       ! For simplicity, we assume global variable out is accessible; in practice, use proper structure handling.
       ! Here, we do nothing.
    end if
  end subroutine set_output_value

  ! Utility: assign uncertainties to output structure for parameter pname
  subroutine assign_uncertainties(pname, cvalcimeas, cvalcib, cvalciw, cvalcu, cvalcin, out)
    implicit none
    character(len=*), intent(in) :: pname
    real(wp), dimension(:), intent(in) :: cvalcimeas, cvalcib, cvalciw, cvalcu, cvalcin
    type(canyonb_output), intent(inout) :: out
    real(wp), allocatable, dimension(:) :: tot_unc
    integer :: n, i
    n = size(cvalcimeas)
    allocate(tot_unc(n))
    do i = 1, n
       tot_unc(i) = sqrt( cvalcimeas(i) + cvalcib(i) + cvalciw(i) + cvalcu(i) + cvalcin(i) )
    end do
    if (trim(pname) == "AT") then
       out%AT_ci = tot_unc
       out%AT_cim = sqrt(cvalcimeas)
       out%AT_cin = sqrt(cvalcib + cvalciw + cvalcu)
       out%AT_cii = sqrt(cvalcin)
    else if (trim(pname) == "CT") then
       out%CT_ci = tot_unc
       out%CT_cim = sqrt(cvalcimeas)
       out%CT_cin = sqrt(cvalcib + cvalciw + cvalcu)
       out%CT_cii = sqrt(cvalcin)
    else if (trim(pname) == "pH") then
       out%pH_ci = tot_unc
       out%pH_cim = sqrt(cvalcimeas)
       out%pH_cin = sqrt(cvalcib + cvalciw + cvalcu)
       out%pH_cii = sqrt(cvalcin)
    else if (trim(pname) == "pCO2") then
       out%pCO2_ci = tot_unc
       out%pCO2_cim = sqrt(cvalcimeas)
       out%pCO2_cin = sqrt(cvalcib + cvalciw + cvalcu)
       out%pCO2_cii = sqrt(cvalcin)
    else if (trim(pname) == "NO3") then
       out%NO3_ci = tot_unc
       out%NO3_cim = sqrt(cvalcimeas)
       out%NO3_cin = sqrt(cvalcib + cvalciw + cvalcu)
       out%NO3_cii = sqrt(cvalcin)
    else if (trim(pname) == "PO4") then
       out%PO4_ci = tot_unc
       out%PO4_cim = sqrt(cvalcimeas)
       out%PO4_cin = sqrt(cvalcib + cvalciw + cvalcu)
       out%PO4_cii = sqrt(cvalcin)
    else if (trim(pname) == "SiOH4") then
       out%SiOH4_ci = tot_unc
       out%SiOH4_cim = sqrt(cvalcimeas)
       out%SiOH4_cin = sqrt(cvalcib + cvalciw + cvalcu)
       out%SiOH4_cii = sqrt(cvalcin)
    end if
    deallocate(tot_unc)
  end subroutine assign_uncertainties

  ! Dummy subroutine to reshape a block of inwgts column into a matrix
  subroutine reshape_column(mat, col_idx, start_idx, end_idx, nrows, ncols, outmat)
    implicit none
    real(wp), intent(in), dimension(:,:) :: mat
    integer, intent(in) :: col_idx, start_idx, end_idx, nrows, ncols
    real(wp), intent(out), dimension(nrows, ncols) :: outmat
    integer :: idx, r, c, count
    count = 0
    do r = 1, nrows
       do c = 1, ncols
          count = count + 1
          outmat(r,c) = mat(start_idx + count - 1, col_idx)
       end do
    end do
  end subroutine reshape_column

  ! Dummy subroutine to extract a vector from inwgts
  subroutine get_vector(mat, col_idx, start_idx, end_idx, vec)
    implicit none
    real(wp), intent(in), dimension(:,:) :: mat
    integer, intent(in) :: col_idx, start_idx, end_idx
    real(wp), intent(out), dimension(end_idx-start_idx+1) :: vec
    integer :: i
    do i = start_idx, end_idx
       vec(i - start_idx + 1) = mat(i, col_idx)
    end do
  end subroutine get_vector

  ! Dummy subroutine read_table: reads a whitespace-delimited text file into a 2D allocatable array
  subroutine read_table(filename, array)
    implicit none
    character(len=*), intent(in) :: filename
    real(wp), allocatable, intent(out), dimension(:,:) :: array
    integer :: unit, i, j, nrows, ncols, io_stat
    character(len=256) :: line
    character(len=:), allocatable :: token
    character(len=1), parameter :: sep = ' '
    integer, allocatable :: row_count(:)
    ! For simplicity, we assume fixed dimensions here (dummy implementation)
    nrows = 100
    ncols = 10
    allocate(array(nrows, ncols))
    array = 0.0_wp
    ! In a complete implementation, file reading and parsing would be done.
  end subroutine read_table

  ! Dummy implementation of carb subroutine (seacarb functionality)
  type carb_out
     real(wp) :: pCO2
  end type carb_out

  subroutine carb(flag, var1, var2, S, T, P, Patm, Pt, Sit, pHscale, kf, k1k2, ks, b, gas, result)
    implicit none
    integer, intent(in) :: flag
    real(wp), intent(in) :: var1
    real(wp), intent(in), dimension(:) :: var2
    real(wp), intent(in) :: S, T, P, Patm, Pt, Sit
    character(len=*), intent(in) :: pHscale, kf, k1k2, ks, b, gas
    type(carb_out), intent(out) :: result
    ! Dummy implementation: simply scale var2 to pCO2 (this is not the actual seacarb calculation)
    if (size(var2) > 0) then
       result%pCO2 = var2(1)*1.0e6_wp + 100.0_wp
    else
       result%pCO2 = 0.0_wp
    end if
  end subroutine carb

  ! Dummy implementation of derivnum subroutine
  subroutine derivnum(var, flag, var1, var2, S, T, P, Patm, Pt, Sit, pHscale, kf, k1k2, ks, b, gas, result)
    implicit none
    character(len=*), intent(in) :: var
    integer, intent(in) :: flag
    real(wp), intent(in) :: var1
    real(wp), intent(in), dimension(:) :: var2
    real(wp), intent(in) :: S, T, P, Patm, Pt, Sit
    character(len=*), intent(in) :: pHscale, kf, k1k2, ks, b, gas
    type(carb_out), intent(out) :: result
    ! Dummy implementation: derivative constant
    result%pCO2 = 0.5_wp
  end subroutine derivnum

  ! Dummy function to return number of columns in data (hardcoded as 8)
  integer function ncol_data()
    implicit none
    ncol_data = 8
  end function ncol_data

end module canyonb_module
