!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP CHyM.
!
!    ICTP CHyM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP CHyM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP CHyM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
module mod_runparams

      real , parameter :: chym_version = 5.0
      character(len=50) , parameter :: &
        CHYMNAME = 'CHyM - CETEMPS Hydrological Model'


      integer :: nlon
      integer :: nlat

      integer :: chym_sdate
      integer :: chym_edate
      integer :: chym_steps

      integer :: chym_tempfl = 1
      integer :: chym_savet = 0
      integer :: chym_savep = 0
      integer :: chym_restart = 0
      integer :: chym_river
      integer :: chym_grads = 0
      integer :: chym_netcdf = 1
      integer :: chym_modis = 0
      integer :: chym_tplot = 1
      integer :: chym_iplot = 0
      integer :: chym_zoom = 0

      logical :: chym_verbose = .false.

      real :: slon
      real :: slat
      real :: dij
      real :: chym_radius = 10.5
      real :: chym_regcm_rad = -1.0

      integer :: nsli = 10
      integer :: nsave = 1
      integer :: rsave = 12
      integer :: angiocycle = 1
      integer :: demf = 20

      character(len=256) :: chym_mfile = ' '
      character(len=256) :: chym_ofile = 'tmp/test1.chym'
      character(len=256) :: chym_sfile = 'tmp/test1.chym'
      character(len=256) :: chym_rfile = 'tmp/test1.chym'
      character(len=256) :: chym_pfile = 'tmp/rainfall.chym'
      character(len=256) :: chym_tfile = 'tmp/temperature.chym'
      character(len=256) :: chym_manning = 'doc/manning.coeff'

      character(len=256) :: chym_dsource = 'era'
      character(len=256) :: chym_ifile1 = ' '
      character(len=256) :: chym_ifile2 = ' '
      character(len=256) :: chem_symtype = 'Exercise with CHyM'
      character(len=256) :: chym_savefld = 'por,rai,evp,rsr,tem'

      integer :: dsfreq = 1

      integer :: ilon1 = 1
      integer :: ilon2 = 150
      integer :: ilat1 = 1
      integer :: ilat2 = 150

      integer :: ncyc1 = 100
      integer :: ncyc2 = 1000
      integer :: ncyc3 = 21
      integer :: angionp = 40

      integer :: integrflag = 0

      integer :: deflate_level = 0

      integer :: num_chunky = 6
      integer :: num_chunkx = 6

      integer :: packnetcdf = 0

      real ::  cpar1=4e-07! Return flow factor (4.8e-07)
      real ::  cpar2=0.0015 ! Alpha coefficients for hydraulic radius (0.0015)
      real ::  cpar3=0.050  ! Beta coefficients for hydraulic radius (0.050)
      real ::  cpar4=0.050  ! Melting temperature factor (0.050)
      real ::  cpar5=0.0094 ! Melting shortwave rad. factor (0.0094)
      real ::  cpar6=500.0  ! River/land threshold (Km2) (500.0)
      real ::  cpar7=90.0   ! Number of days to consider for return flow (90)
      real ::  cpar8=4.5    ! Reduction of land/channel manning coefficient
      real ::  cpar9=200.0  ! River/land threshold (Km2) for returnflow
      real ::  cpar10=0.0    ! Not yet used
      real ::  infiltr = 40.0
      real ::  infi_lago=0.0
      real ::  infi_fiume=0.0
      real ::  infi_ice=0.0

end module mod_runparams
