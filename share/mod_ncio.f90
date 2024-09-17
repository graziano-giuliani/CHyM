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

module mod_ncio

  use netcdf
  use mod_internal
  use mod_statparams
  use mod_runparams
  use mod_libmv
  use mod_time
  use mod_mssg
  use mod_museo
  use mod_varandtypes
  use mod_mpimess

  private

  public :: createfile , closefile , closerstfile
  public :: createfile_rst
  public :: write_statvar , add_timestep , write_dynvar
  public :: writerst_dynvar , addrst_timestep
  public :: read_restart_dyn_NC, getasterdata
  public :: read_restart_stat_NC, gethydrodata
  public :: gethydrodata_real,gethydrodata_int                        !interface between gethydrodata gethydrodata_real ?
  public :: ein75readnc , persiannreadncll
  public :: trimmreadncll, erareadncll, erareadncll2
  public :: erareadncll3
  public :: writerst_dynvar_real
  public :: hourlyrainnc

  interface chymreadrec
     module procedure chymreadirec
     module procedure chymreadrrec
  end interface

  interface write_statvar
    module procedure write_statvar_real
    module procedure write_statvar_integer
  end interface write_statvar

  interface write_dynvar
    module procedure write_dynvar_real
    module procedure write_dynvar_integer2
    module procedure write_dynvar_integer4
  end interface write_dynvar

  character(len=64) , parameter :: institution = 'ICTP'

  type chym_variable
    character(len=8) :: vname
    character(len=8) :: vunit
    character(len=32) :: lname
    character(len=32) :: sname
    integer :: itype
    integer :: varid
  end type chym_variable

  type global_domain
    integer :: global_ni
    integer :: global_nj
    integer :: ntiles
    integer , dimension(2) :: ni
    integer , dimension(2) :: igstart
    integer , dimension(2) :: igstop
    integer :: nj
    integer :: jgstart
    integer :: jgstop
  end type global_domain

  integer , parameter :: npar = 50
  integer , parameter :: slen = 50
  integer , parameter :: maxf = 33
  integer , parameter :: start_dyn = 13
  integer , parameter :: start_rst = 29
  integer , parameter :: maxdyn = 16
  integer , parameter :: maxrst = 5

  character(len=slen) , dimension(maxdyn+maxrst) :: schymes
!  character(len=slen) , dimension(maxrst) :: rchymes

  type(chym_variable) , dimension(maxf) :: var =  (/ &
   chym_variable('dem','m','DEM Elevation','elevation',1,-1), &
   chym_variable('lat','degrees_north','Grid Latitude','grid_latitude',1,-1), &
   chym_variable('lon','degrees_east','Grid Longitude','grid_longitude',1,-1), &
   chym_variable('fdm','1','Flow direction matrix','flow_direction',0,-1), &
   chym_variable('acc','1','Tangent angle for flow direction','slope',1,-1), &
   chym_variable('lus','1','Land Use Category','soil_category',0,-1), &
   chym_variable('aer','km^2','Cell area','cell_area',1,-1), &
   chym_variable('dra','km^2','Drainage area','drainage_area',1,-1), &
   chym_variable('run','hours','Runoff time','runoff_time',1,-1), &
   chym_variable('alf','km/h','Flow velocity','alfa',1,-1), &
   chym_variable('ctr','m','Flow control','flow control',1,-1), &
   chym_variable('bas','1','Basins','basins',1,-1), &
   chym_variable('mchym','1','mchym array','mchym',0,-1), &
   chym_variable('rchym','1','rchym array','rchym',1,-1), &
   chym_variable('schym','1','schym array','schym',1,-1), &
   chym_variable('rai','mm/h','Input rain','precipitation_flux',1,-1), &
   chym_variable('rsr','1','Rain Source','rain_source_class',1,-1), &
   chym_variable('ara','mm','Water available for runoff','pot_runoff',1,-1), &
   chym_variable('por','m^3 s-1','Flow discharge','flow_discharge',1,-1), &
   chym_variable('wet','m^2','Wetted area x cell-channell lenght','wet',1,-1), &
   chym_variable('gwt','mm','Ground and Vegetation water content','gwt',1,-1), &
   chym_variable('h2o','mm','Total water content','h2o',1,-1), &
   chym_variable('evp','mm','Evaporation','evaporation',1,-1), &
   chym_variable('sno','mm','Snow cover','snow_amount',1,-1), &
   chym_variable('tem','C','Temperature','temperature',1,-1), &
   chym_variable('dgw','mm','Deep soil Water content','dswc',1,-1), &
   chym_variable('ddw','mm','Acc Matr Deep soil Water content','dswd',1,-1), &
   chym_variable('unk','unknown','unknown','unknown',1,-1), &
   chym_variable('rgwt','mm','Ground and Vegetation water content','gwt',1,-1), &
   chym_variable('rh2o','mm','Total water content','h2o',1,-1), &
   chym_variable('rsno','mm','Snow cover','snow_amount',1,-1), &
   chym_variable('rpor','m^3 s-1','Flow discharge','flow_discharge',1,-1), &
   chym_variable('rdgw','mm','Deep soil Water content','dswc',1,-1) /)
  character(len=3) :: fld
  character(len=2) :: ftype
  integer :: isize , ildate

  logical , parameter :: static = .false.
  logical , parameter :: dynamic = .true.
  character(len=32) :: isodate0 , isodate1 , source
  character(len=64) :: timeunit
  integer :: nstep
  real :: lat0 , lon0 , dlat , dlon
  integer :: nstatic , ndynamic , nrestart
  integer , dimension(4) :: id_dims
  integer , dimension(4) :: rd_dims
  integer :: lon_varid , lat_varid , time_varid , timerst_varid
  integer :: mchym_varid , rchym_varid, schym_varid
  real , dimension(:) , allocatable :: lats , lons
  real , dimension(:,:) , allocatable :: vals
  integer , dimension(:,:) , allocatable :: ivals
  integer :: year0, month0, day0, hour0
  integer :: year1, month1, day1, hour1
  integer :: ncid
  integer :: ncidrst
  integer , dimension(3) :: istart , icount
  integer , dimension(:) , allocatable :: map_ids , maprst_ids
  integer :: istep , rstep , ivar , rvar
  logical :: bisest
  character(len=6) , dimension(12) :: xmese
  character(len=10) :: amese
  character(len=4) :: xyear
  character(len=2) :: xmonth


  contains

!  subroutine createfile(outname,mchym,rchym,schym,ncdata)
  subroutine createfile(outname,ncdata)
    implicit none
    include '../doc/settings.inc'
!    integer , dimension(npar) , intent(in) :: mchym
!    real , dimension(npar) , intent(in) :: rchym
!    character(len=slen) , dimension(npar) , intent(in) :: schym
    character(len=*) , intent(in) :: outname
    integer, intent(in) :: ncdata
    logical :: fcentury
    integer :: i,j

    print*,''
    fcentury = .true.
    schymes(1:maxdyn) = schym(39:54)

    nlon = mchym(2)
    nlat = mchym(3)

    lon0 = rchym(1)
    lat0 = rchym(2)
    dlon = rchym(8)
    dlat = rchym(9)

    nstatic  = mchym(7)+3
    ndynamic = mchym(8)
    allocate(map_ids(ndynamic))

    nstep = mchym(6)

    write(source,'(a,a)') 'CHYM model version ', VERSION
    !iretval = nf90_create(outname,nf90_clobber,ncid)
    if (deflate_level>0) then
      iretval = nf90_create(outname,nf90_netcdf4,ncid)
    else if (deflate_level==0) then
      iretval = nf90_create(outname,nf90_clobber,ncid)
    else
      print*,"deflate_level not valid EXIT"
    end if

    call checknetcdfresult('create output file')

    iretval = nf90_def_dim(ncid,'lon',nlon,id_dims(1))
    call checknetcdfresult('Add lon dimension')
    iretval = nf90_def_dim(ncid,'lat',nlat,id_dims(2))
    call checknetcdfresult('Add lat dimension')
    iretval = nf90_def_dim(ncid,'nchymp',nchymp,id_dims(4))
    call checknetcdfresult('Add nchymp dimension')
!
!   Global Attributes
!
    iretval = nf90_put_att(ncid,nf90_global,'Conventions','CF-1.6')
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global,'title',trim(schym(6)))
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global,'history', &
             'created on '//trim(schym(4))//' on '//trim(schym(5)))
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global,'institution',institution)
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global,'source',source)
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global,'comment', &
             'Input dataset is '//trim(schym(15)))
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global,'nslices',mchym(6))
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global,'is_restart',mchym(9))
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global,'landuse_code_for_sea', &
              mchym(10))
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global, &
             'landuse_code_for_internal_water',mchym(11))
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global, &
             'mm5_domain_for_museo_data',mchym(12))
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global,'dem_source_code', &
             mchym(13))
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global,'output_frequency',mchym(14))
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global, &
             'temperature_source_flag',mchym(15))
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global,'hourly_steps',mchym(16))
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global,'rainfall_calculation', &
             mchym(17))
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global,'river_selected', &
             mchym(18))
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global,'time_slices',mchym(22))
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global, &
             'integration_steps_per_hour',mchym(21))
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global, &
             'temperature_calculation',mchym(23))
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global,'modis_data_used', &
             mchym(24))
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global,'requested_variables', &
             schym(7))
    call checknetcdfresult('Add global attribute')
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global, &
             'radius_of_influence_for_sparse_data_interpolation', &
             rchym(10))
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global, &
             'approximate_resolution_in_meters',rchym(6))
    iretval = nf90_put_att(ncid,nf90_global, &
             'approximate_lat-lon_resolution_in_meters',rchym(5))
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global, &
             'Start_longitude',rchym(1))
    call checknetcdfresult('Add global attribute')
    iretval = nf90_put_att(ncid,nf90_global, &
             'Start_latitude',rchym(2))
    call checknetcdfresult('Add global attribute')
    iretval = nf90_def_var(ncid,'lon',nf90_float,id_dims(1:1), &
             lon_varid)
    call checknetcdfresult('Add lon var')
!#ifdef NETCDF4_HDF5
!    iretval = nf90_def_var_deflate(ncid, lon_varid, 1, 1, deflate_level)
!    call checknetcdfresult('Error setting deflate on lon')
!#endif
    iretval = nf90_put_att(ncid,lon_varid,'units','degrees_east')
    call checknetcdfresult('Add lon units')
    iretval = nf90_put_att(ncid,lon_varid,'long_name','Longitude')
    call checknetcdfresult('Add lon long_name')
    iretval = nf90_put_att(ncid,lon_varid,'standard_name','longitude')
    call checknetcdfresult('Add lon standard_name')
    iretval = nf90_def_var(ncid,'lat',nf90_float,id_dims(2:2), &
            lat_varid)
    call checknetcdfresult('Add lat var')
!#ifdef NETCDF4_HDF5
!    iretval = nf90_def_var_deflate(ncid, lat_varid, 1, 1, deflate_level)
!    call checknetcdfresult('Error setting deflate on lat')
!#endif
    iretval = nf90_put_att(ncid,lat_varid,'units','degrees_north')
    call checknetcdfresult('Add lat units')
    iretval = nf90_put_att(ncid,lat_varid,'long_name','Latitude')
    call checknetcdfresult('Add lat long_name')
    iretval = nf90_put_att(ncid,lat_varid,'standard_name','latitude')
    call checknetcdfresult('Add lat standard_name')
    iretval = nf90_def_var(ncid,'mchym',nf90_int,id_dims(4:4), &
           mchym_varid)
    call checknetcdfresult('Add mchym var')
    iretval = nf90_def_var(ncid,'rchym',nf90_int,id_dims(4:4), &
           rchym_varid)
    call checknetcdfresult('Add rchym var')
!
! Write time info if requested
!
    if (ncdata.ge.0) then
      write(amese,'(i10)')mchym(5)
      call gmafromindex(ncdata,hour0,day0,month0,year0)
      write(xyear,'(i4)')year0 ; write(xmonth,'(i2)')month0
      bisest = .false.
      if (MOD(year0,400)==0.or.(MOD(year0,4)==0.and.(.not. &
                 (MOD(year0,100)==0)))) then
        xmese = (/'013124','022924','033124','043024','053124','063024', &
                '073124','083124','093024','103124','113024','123124'/)
        bisest = .true.
      else
        xmese = (/'013124','022824','033124','043024','053124','063024', &
                '073124','083124','093024','103124','113024','123124'/)
      end if
      if ((xyear//xmonth=='199912').and.(mchym(5)<1000000) ) then
        fcentury = .false.
      end if
      if (xyear(3:4)//xmese(month0) > amese .and. fcentury) then
!        call index_to_date(mchym(5),year1,month1,day1,hour1)
         call gmafromindex(mchym(5),hour1,day1,month1,year1)
      else
!        call index_to_date(ncdata,year1,month1,day1,hour1)
         call gmafromindex(ncdata,hour1,day1,month1,year1)
        read(xmese(month0)(1:2),'(i2)')month1
        read(xmese(month0)(3:4),'(i2)')day1
        read(xmese(month0)(5:6),'(i2)')hour1
      end if

      write(isodate0,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a)') &
           year0,'-',month0,'-',day0,' ',hour0,':00:00'
      write(isodate1,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a)') &
           year1,'-',month1,'-',day1,' ',hour1,':00:00'
      write(timeunit,'(a,i0.4,a,i0.2,a,i0.2,a,i0.2,a)') &
           'seconds since ',year0,'-',month0,'-',day0,' ',hour0, &
          ':00:00 UTC'

      iretval = nf90_put_att(ncid,nf90_global,'date_start',isodate0)
      call checknetcdfresult('Add global attribute')
      iretval = nf90_put_att(ncid,nf90_global,'date_end',isodate1)
      call checknetcdfresult('Add global attribute')
      iretval = nf90_def_dim(ncid,'time',nf90_unlimited,id_dims(3))
      call checknetcdfresult('Add time dimension')
      iretval = nf90_def_var(ncid,'time',nf90_double,id_dims(3:3), &
             time_varid)
      call checknetcdfresult('Add time var')
      iretval = nf90_put_att(ncid,time_varid,'units', trim(timeunit))
      call checknetcdfresult('Add time units')
      iretval = nf90_put_att(ncid,time_varid,'long_name','Time')
      call checknetcdfresult('Add time long_name')
      iretval = nf90_put_att(ncid,time_varid,'standard_name','time')
      call checknetcdfresult('Add time standard_name')
      iretval = nf90_put_att(ncid,time_varid,'calendar','gregorian')
      call checknetcdfresult('Add time calendar')
      do ivar = 1 , ndynamic
        call add_variable(ivar,dynamic)
      end do
    else
      do ivar = 1 , nstatic
        call add_variable(ivar,static)
      end do
    endif
    iretval = nf90_enddef(ncid)
    call checknetcdfresult('Exit from netcdf define mode')
    allocate(lats(nlat),lons(nlon),vals(nlon,nlat))!,ivals(nlon,nlat))
    do j = 1 , nlat
      lats(j) = real(j-1)*dlat+lat0
    end do
    do i = 1 , nlon
      lons(i) = real(i-1)*dlon+lon0
    end do
    iretval = nf90_put_var(ncid,lon_varid,lons)
    call checknetcdfresult('Write longitudes')
    iretval = nf90_put_var(ncid,lat_varid,lats)
    call checknetcdfresult('Write latitudes')
    iretval = nf90_put_var(ncid,mchym_varid,mchym)
    call checknetcdfresult('Write mchym')
    iretval = nf90_put_var(ncid,rchym_varid,rchym)
    call checknetcdfresult('Write rchym')
    deallocate(lons,lats)
    istep = 0
  end subroutine createfile

  subroutine createfile_rst(restname,ncdata)
    implicit none
    character(len=*) , intent(in) :: restname
    integer, intent(in) :: ncdata
    logical :: fcentury
    integer :: i, j
    fcentury = .true.

    schymes(maxdyn+1:maxdyn+maxrst) = schym(56:60)

    nlon = mchym(2)
    nlat = mchym(3)

    lon0 = rchym(1)
    lat0 = rchym(2)
    dlon = rchym(8)
    dlat = rchym(9)

    nrestart = 5
    allocate(maprst_ids(nrestart))

    if (deflate_level>0) then
    iretval = nf90_create(restname,nf90_netcdf4,ncidrst)
    else if (deflate_level==0) then
    iretval = nf90_create(restname,nf90_clobber,ncidrst)
    else
      print*,"deflate_level not valid EXIT"
    end if

    call checknetcdfresult('create restart file')
    iretval = nf90_def_dim(ncidrst,'lon',nlon,rd_dims(1))
    call checknetcdfresult('Add lon dimension')
    iretval = nf90_def_dim(ncidrst,'lat',nlat,rd_dims(2))
    call checknetcdfresult('Add lat dimension')
    if (ncdata.ge.0) then
      write(amese,'(i10)')mchym(5)
      call gmafromindex(ncdata,hour0,day0,month0,year0)
      write(xyear,'(i4)')year0 ; write(xmonth,'(i2)')month0
      bisest = .false.
      if (MOD(year0,400)==0.or.(MOD(year0,4)==0.and.(.not. &
                 (MOD(year0,100)==0)))) then
        xmese = (/'013124','022924','033124','043024','053124','063024', &
                '073124','083124','093024','103124','113024','123124'/)
        bisest = .true.
      else
        xmese = (/'013124','022824','033124','043024','053124','063024', &
                '073124','083124','093024','103124','113024','123124'/)
      end if
      if ((xyear//xmonth=='199912').and.(mchym(5)<1000000) ) then
        fcentury = .false.
      end if
      if (xyear(3:4)//xmese(month0) > amese .and. fcentury) then
!        call index_to_date(mchym(5),year1,month1,day1,hour1)
         call gmafromindex(mchym(5),hour1,day1,month1,year1)
      else
!        call index_to_date(ncdata,year1,month1,day1,hour1)
         call gmafromindex(ncdata,hour1,day1,month1,year1)
        read(xmese(month0)(1:2),'(i2)')month1
        read(xmese(month0)(3:4),'(i2)')day1
        read(xmese(month0)(5:6),'(i2)')hour1
      end if

      write(isodate0,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a)') &
           year0,'-',month0,'-',day0,' ',hour0,':00:00'
      write(isodate1,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a)') &
           year1,'-',month1,'-',day1,' ',hour1,':00:00'
      write(timeunit,'(a,i0.4,a,i0.2,a,i0.2,a,i0.2,a)') &
           'seconds since ',year0,'-',month0,'-',day0,' ',hour0, &
          ':00:00 UTC'
      iretval = nf90_def_var(ncidrst,'lon',nf90_float,rd_dims(1:1), &
               lon_varid)
      call checknetcdfresult('Add lon var')
      iretval = nf90_put_att(ncidrst,lon_varid,'units','degrees_east')
      call checknetcdfresult('Add lon units')
      iretval = nf90_put_att(ncidrst,lon_varid,'long_name','Longitude')
      call checknetcdfresult('Add lon long_name')
      iretval = nf90_put_att(ncidrst,lon_varid,'standard_name','longitude')
      call checknetcdfresult('Add lon standard_name')
      iretval = nf90_def_var(ncidrst,'lat',nf90_float,rd_dims(2:2), &
              lat_varid)
      call checknetcdfresult('Add lat var')
      iretval = nf90_put_att(ncidrst,lat_varid,'units','degrees_north')
      call checknetcdfresult('Add lat units')
      iretval = nf90_put_att(ncidrst,lat_varid,'long_name','Latitude')
      call checknetcdfresult('Add lat long_name')
      iretval = nf90_put_att(ncidrst,lat_varid,'standard_name','latitude')
      call checknetcdfresult('Add lat standard_name')
      iretval = nf90_put_att(ncidrst,nf90_global,'Conventions','CF-1.6')
      call checknetcdfresult('Add global attribute')
      iretval = nf90_put_att(ncidrst,nf90_global,'date_start',isodate0)
      call checknetcdfresult('Add global attribute')
      iretval = nf90_put_att(ncidrst,nf90_global,'date_end',isodate1)
      call checknetcdfresult('Add global attribute')
      iretval = nf90_def_dim(ncidrst,'time',nf90_unlimited,rd_dims(3))
      call checknetcdfresult('Add time dimension')
      iretval = nf90_def_var(ncidrst,'time',nf90_double,rd_dims(3:3), &
             timerst_varid)
      call checknetcdfresult('Add time var')
      iretval = nf90_put_att(ncidrst,timerst_varid,'units', trim(timeunit))
      call checknetcdfresult('Add time units')
      iretval = nf90_put_att(ncidrst,timerst_varid,'long_name','Time')
      call checknetcdfresult('Add time long_name')
      iretval = nf90_put_att(ncidrst,timerst_varid,'standard_name','time')
      call checknetcdfresult('Add time standard_name')
      iretval = nf90_put_att(ncidrst,timerst_varid,'calendar','gregorian')
      call checknetcdfresult('Add time calendar')
      do ivar = 1 , nrestart
        call addrst_variable(ivar,dynamic)
      end do
    endif
    iretval = nf90_enddef(ncidrst)
    call checknetcdfresult('Exit from netcdf define mode')
    allocate(lats(nlat),lons(nlon))!,vals(nlon,nlat))!,ivals(nlon,nlat))
    do j = 1 , nlat
      lats(j) = real(j-1)*dlat+lat0
    end do
    do i = 1 , nlon
      lons(i) = real(i-1)*dlon+lon0
    end do
    iretval = nf90_put_var(ncidrst,lon_varid,lons)
    call checknetcdfresult('Write longitudes')
    iretval = nf90_put_var(ncidrst,lat_varid,lats)
    call checknetcdfresult('Write latitudes')
    deallocate(lons,lats)
    rstep = 0
  end subroutine createfile_rst

  subroutine write_statvar_real(vname,rval)
    implicit none
    character(len=*) , intent(in) :: vname
    real , dimension(:,:) , intent(in) :: rval
    integer :: indx
    if ( vname == 'lat' .or. vname == 'lon' ) return
    indx = search_indx(vname)
    iretval = nf90_put_var(ncid,var(indx)%varid,rval)
    call checknetcdfresult('Write '//var(indx)%vname)
  end subroutine write_statvar_real

  subroutine write_statvar_integer(vname,ival)
    implicit none
    character(len=*) , intent(in) :: vname
    integer , dimension(:,:) , intent(in) :: ival
    integer :: indx
    if ( vname == 'lat' .or. vname == 'lon' ) return
    indx = search_indx(vname)
    iretval = nf90_put_var(ncid,var(indx)%varid,ival)
    call checknetcdfresult('Write '//var(indx)%vname)
  end subroutine write_statvar_integer

  subroutine write_dynvar_real(vname,rval)
    implicit none
    character(len=*) , intent(in) :: vname
    real , dimension(:,:) , intent(in) :: rval
    integer :: indx
    istart(1) = 1
    istart(2) = 1
    istart(3) = istep
    icount(1) = nlon
    icount(2) = nlat
    icount(3) = 1
    indx = search_indx(vname)
    iretval = nf90_put_var(ncid,var(indx)%varid,rval,istart,icount)
    call checknetcdfresult('Write '//var(indx)%vname)
  end subroutine write_dynvar_real

  subroutine write_dynvar_integer2(vname,ival)
    implicit none
    character(len=*) , intent(in) :: vname
    integer(short) , dimension(:,:) , intent(in) :: ival
    integer :: indx
    istart(1) = 1
    istart(2) = 1
    istart(3) = istep
    icount(1) = nlon
    icount(2) = nlat
    icount(3) = 1
    indx = search_indx(vname)
    iretval = nf90_put_var(ncid,var(indx)%varid,ival,istart,icount)
    call checknetcdfresult('Write '//var(indx)%vname)
  end subroutine write_dynvar_integer2

  subroutine write_dynvar_integer4(vname,ival)
    implicit none
    character(len=*) , intent(in) :: vname
    integer(long) , dimension(:,:) , intent(in) :: ival
    integer :: indx
    istart(1) = 1
    istart(2) = 1
    istart(3) = istep
    icount(1) = nlon
    icount(2) = nlat
    icount(3) = 1
    indx = search_indx(vname)
    iretval = nf90_put_var(ncid,var(indx)%varid,ival,istart,icount)
    call checknetcdfresult('Write '//var(indx)%vname)
  end subroutine write_dynvar_integer4

  subroutine writerst_dynvar_real(vname,rval)
    implicit none
    character(len=*) , intent(in) :: vname
    real , dimension(:,:) , intent(in) :: rval
    integer :: indx
    istart(1) = 1
    istart(2) = 1
    istart(3) = rstep
    icount(1) = nlon
    icount(2) = nlat
    icount(3) = 1
    indx = search_indx(vname)
    iretval = nf90_put_var(ncidrst,var(indx)%varid,rval,istart,icount)
    call checknetcdfresult('Write '//var(indx)%vname)
  end subroutine writerst_dynvar_real

  subroutine add_timestep
    implicit none
    double precision , dimension(1) :: xtime
    istep = istep + 1
    xtime(1) = dble(istep-1)*3600.0D0*mchym(14)
    istart(1) = istep
    icount(1) = 1
    iretval = nf90_put_var(ncid,time_varid,xtime,istart(1:1), &
           icount(1:1))
    call checknetcdfresult('Write time')
  end subroutine add_timestep

  subroutine addrst_timestep
    implicit none
    double precision , dimension(1) :: xtime
    rstep = rstep + 1
    xtime(1) = dble(rstep)*3600.0D0*mchym(26)
    istart(1) = rstep
    icount(1) = 1
    iretval = nf90_put_var(ncidrst,timerst_varid,xtime,istart(1:1), &
           icount(1:1))
    call checknetcdfresult('Write time')
  end subroutine addrst_timestep

   subroutine add_variable(indx,ldyn)
     implicit none
     integer , intent(in) :: indx
     logical , intent(in) :: ldyn
     character(len=8) :: vname
     integer :: varid , int_indx
     if ( ldyn ) then
       call search_name(indx)
       int_indx = map_ids(indx)
       vname = var(int_indx)%vname
       if ( var(int_indx)%itype == 0 ) then
         iretval = nf90_def_var(ncid,vname,nf90_int,id_dims(1:3),varid)
         if (deflate_level>0) then
           iretval = nf90_def_var_deflate(ncid,varid,1,1,deflate_level)
           iretval = nf90_def_var_chunking(ncid, varid, 0, chunksizes)
         end if
!         iretval = NF90_DEF_VAR_FILL(ncid, varid, 0, NF90_FILL_INT)
!         print*,iretval
         iretval = nf90_put_att(ncid,varid,'_FillValue',NF90_FILL_INT)
         call checknetcdfresult('NF90_DEF_VAR_FILL fill_int')
       else if ( var(int_indx)%itype == 1 ) then
         iretval = nf90_def_var(ncid,vname,nf90_float,id_dims(1:3),varid)
         if (deflate_level>0) then
           iretval = nf90_def_var_deflate(ncid,varid,1,1,deflate_level)
           iretval = nf90_def_var_chunking(ncid, varid, 0, chunksizes)
         end if
!         iretval = NF90_DEF_VAR_FILL(ncid, varid, 0, NF90_FILL_FLOAT)
!         print*,iretval
         iretval = nf90_put_att(ncid,varid,'_FillValue',NF90_FILL_FLOAT)
         call checknetcdfresult('NF90_DEF_VAR_FILL fill_float')
       else if ( var(int_indx)%itype == 2 ) then
         iretval = nf90_def_var(ncid,vname,nf90_short,id_dims(1:3),varid)
         if (deflate_level>0) then
         iretval = nf90_def_var_deflate(ncid,varid,1,1,deflate_level)
         iretval = nf90_def_var_chunking(ncid, varid, 0, chunksizes)
         end if
!         iretval = NF90_DEF_VAR_FILL(ncid, varid, 0, NF90_FILL_INT)
!         print*,iretval
!         iretval = nf90_put_att(ncid,varid,'_FillValue',NF90_FILL_SHORT)
         iretval = nf90_put_att(ncid,varid,'_FillValue',-NF90_FILL_SHORT)
         call checknetcdfresult('NF90_DEF_VAR_FILL fill_int')
         iretval = NF90_PUT_ATT(ncid,varid,'scale_factor',scale_factor_por)
         call checknetcdfresult('scale_factor_por')
         iretval = NF90_PUT_ATT(ncid,varid,'add_offset',add_offset_por)
         call checknetcdfresult('add_offset_por')
       end if
     else
       int_indx = indx
       vname = var(indx)%vname
       if ( vname == 'mchym') return
       if ( vname == 'lat' .or. vname == 'lon' ) return
       if ( var(int_indx)%itype == 0 ) then
         iretval = nf90_def_var(ncid,vname,nf90_int,id_dims(1:2),varid)
         if (deflate_level>0) then
         iretval = nf90_def_var_deflate(ncid,varid,1,1,deflate_level)
         iretval = nf90_def_var_chunking(ncid, varid, 0, chunksizes)
         end if
!         iretval = NF90_DEF_VAR_FILL(ncid, varid, 0, NF90_FILL_INT)
!         print*,iretval
         iretval = nf90_put_att(ncid,varid,'_FillValue',NF90_FILL_INT)
         call checknetcdfresult('NF90_DEF_VAR_FILL fill_int')
       else if ( var(int_indx)%itype == 1 ) then
         iretval = nf90_def_var(ncid,vname,nf90_float,id_dims(1:2),varid)
         if (deflate_level>0) then
         iretval = nf90_def_var_deflate(ncid,varid,1,1,deflate_level)
         iretval = nf90_def_var_chunking(ncid, varid, 0, chunksizes)
         end if
!         iretval = NF90_DEF_VAR_FILL(ncid, varid, 0, NF90_FILL_FLOAT)
!         print*,iretval
         iretval = nf90_put_att(ncid,varid,'_FillValue',NF90_FILL_FLOAT)
         call checknetcdfresult('NF90_DEF_VAR_FILL fill_float')
       end if
     end if
     call checknetcdfresult('Add '//vname//' var')
     iretval = nf90_put_att(ncid,varid,'units',var(int_indx)%vunit)
     call checknetcdfresult('Add '//vname//' units')
     iretval = nf90_put_att(ncid,varid,'long_name',var(int_indx)%lname)
     call checknetcdfresult('Add '//vname//' long_name')
     iretval = nf90_put_att(ncid,varid,'standard_name',var(int_indx)%sname)
     call checknetcdfresult('Add '//vname//' standard_name')
     var(int_indx)%varid = varid
     if ( vname == 'ctr') then
        iretval = nf90_put_att(ncid,varid,"missing_value", real(-10))
        call checknetcdfresult('Put Attribute missing_value in ' // var(indx)%vname)
     end if
     if ( vname == 'bas') then
!        iretval = nf90_put_att(ncid,varid,"valid_range", (/-10,nsec/))
!        call checknetcdfresult('Put Attribute valid_range in ' // var(indx)%vname)
        iretval = nf90_put_att(ncid,varid,"missing_value", real(-10))
        call checknetcdfresult('Put Attribute missing_value in ' // var(indx)%vname)
     end if
     var(int_indx)%varid = varid
   end subroutine add_variable

   subroutine addrst_variable(indx,ldyn)
     implicit none
     integer , intent(in) :: indx
     logical , intent(in) :: ldyn
     character(len=8) :: vname
     integer :: varid , int_indx
     if ( ldyn ) then
       call searchrst_name(indx)
       int_indx = maprst_ids(indx)
       vname = var(int_indx)%vname
       if ( var(int_indx)%itype == 0 ) then
         iretval = nf90_def_var(ncidrst,vname,nf90_int,rd_dims(1:3),varid)
         if (deflate_level>0) then
         iretval = nf90_def_var_deflate(ncidrst,varid,1,1,deflate_level)
         iretval = nf90_def_var_chunking(ncidrst, varid, 0, chunksizes)
         end if
       else
         iretval = nf90_def_var(ncidrst,vname,nf90_float,rd_dims(1:3),varid)
         if (deflate_level>0) then
         iretval = nf90_def_var_deflate(ncidrst,varid,1,1,deflate_level)
         iretval = nf90_def_var_chunking(ncidrst, varid, 0, chunksizes)
         end if
       end if
     else
       int_indx = indx
       vname = var(indx)%vname
       if ( vname == 'mchym') return
       if ( vname == 'lat' .or. vname == 'lon' ) return
       if ( var(int_indx)%itype == 0 ) then
         iretval = nf90_def_var(ncidrst,vname,nf90_int,rd_dims(1:2),varid)
         if (deflate_level>0) then
         iretval = nf90_def_var_deflate(ncidrst,varid,1,1,deflate_level)
         iretval = nf90_def_var_chunking(ncidrst, varid, 0, chunksizes)
         end if
       else
         iretval = nf90_def_var(ncidrst,vname,nf90_float,rd_dims(1:2),varid)
         if (deflate_level>0) then
         iretval = nf90_def_var_deflate(ncidrst,varid,1,1,deflate_level)
         iretval = nf90_def_var_chunking(ncidrst, varid, 0, chunksizes)
         end if
       end if
     end if
     call checknetcdfresult('Add '//vname//' var')
     iretval = nf90_put_att(ncidrst,varid,'units',var(int_indx)%vunit)
     call checknetcdfresult('Add '//vname//' units')
     iretval = nf90_put_att(ncidrst,varid,'long_name',var(int_indx)%lname)
     call checknetcdfresult('Add '//vname//' long_name')
     iretval = nf90_put_att(ncidrst,varid,'standard_name',var(int_indx)%sname)
     call checknetcdfresult('Add '//vname//' standard_name')
     var(int_indx)%varid = varid
     if ( vname == 'ctr') then
        iretval = nf90_put_att(ncidrst,varid,"missing_value", real(-10))
        call checknetcdfresult('Put Attribute missing_value in ' //var(indx)%vname)
     end if
     if ( vname == 'bas') then
        iretval = nf90_put_att(ncidrst,varid,"valid_range", real((/-10,nsec/)))
        call checknetcdfresult('Put Attribute valid_range in ' //var(indx)%vname)
        iretval = nf90_put_att(ncidrst,varid,"missing_value", real(-10))
        call checknetcdfresult('Put Attribute missing_value in ' //var(indx)%vname)
     end if
     var(int_indx)%varid = varid
   end subroutine addrst_variable


   subroutine search_name(indx)
     implicit none
     integer, intent(in) :: indx
     integer :: ivar
     map_ids(indx) = -1
     do ivar = start_dyn , maxf
       if ( schymes(indx) == var(ivar)%vname ) then
         map_ids(indx) = ivar
         exit
       end if
     end do
     if ( map_ids(indx) == -1 ) then
       write(0,*) 'Variable ', schymes(indx), ' not found in defined vars'
       stop
     end if
   end subroutine search_name

   subroutine searchrst_name(indx)
     implicit none
     integer, intent(in) :: indx
     integer :: ivar
     maprst_ids(indx) = -1
     do ivar = start_rst , maxf
       if ( schymes(indx+maxdyn) == var(ivar)%vname ) then
         maprst_ids(indx) = ivar
         exit
       end if
     end do
     if ( maprst_ids(indx) == -1 ) then
       write(0,*) 'Restart variable ', schymes(indx), ' not found in defined vars'
       stop
     end if
   end subroutine searchrst_name

  subroutine read_restart_dyn_NC
    implicit none
    integer :: lun,idate,itime,icode,i,j,istatus
    integer, dimension(4) :: istart,icount
    idate=mchym(4)
    write(6,'(18x,a,a)') 'Opening Restart file ',schym(18)
    write(6,'(12x,a)') 'Now reading dynamic fields from '// &
        trim(schym(18))
    port=0.0 ; h2o=0.0 ; gh2o=0.0 ; snow=0.0 ; deepw=0.0
    call getlun(lun)
    istatus = nf90_open(trim(schym(18)),nf90_nowrite,lun)
    call chkncstatus(istatus,nf90_noerr,'open',schym(18),'')
    istatus = nf90_inq_dimid(lun,'time',icode)
    call chkncstatus(istatus,nf90_noerr,'time code',schym(18),'')
    istatus = nf90_inquire_dimension(lun,icode,len=itime)
    call chkncstatus(istatus,nf90_noerr,'time dim',schym(18),'')
    istart(1)=1
    istart(2)=1
    icount(1)=nlon
    icount(2)=nlat
    istart(3)=itime ; icount(3)=1
    istatus = nf90_inq_varid(lun,'rpor',icode)
    call chkncstatus(istatus,nf90_noerr,'rpor code',trim(schym(18)) &
       ,'')
    istatus = nf90_get_var(lun,icode,port,istart(1:3), &
                icount(1:3))
    call chkncstatus(istatus,nf90_noerr,'Reading rpor',trim(schym(18)), &
       '')
    call writetime('Initializing flow-discharge for',idate)
    istatus = nf90_inq_varid(lun,'rh2o',icode)
    call chkncstatus(istatus,nf90_noerr,'rh2o code',trim(schym(18)) &
       ,'')
    istatus = nf90_get_var(lun,icode,h2o,istart(1:3), &
                icount(1:3))
    call chkncstatus(istatus,nf90_noerr,'Reading rh2o',trim(schym(18)), &
       '')
    call writetime('Initializing Total water content for',idate)
!    bwet=0.0
!    do i=1,nlon ; do j=1,nlat
!       if (alfa(i,j)>0.0001) bwet(i,j)=port(i,j)/alfa(i,j)
!    enddo ; enddo
!    h2o=bwet*dx
    istatus = nf90_inq_varid(lun,'rgwt',icode)
    call chkncstatus(istatus,nf90_noerr,'rgwt code',trim(schym(18)) &
       ,'')
    istatus = nf90_get_var(lun,icode,gh2o,istart(1:3), &
                icount(1:3))
    call chkncstatus(istatus,nf90_noerr,'Reading rgwt',trim(schym(18)) &
       ,'')
    call writetime('Initializing ground water for',idate)
    istatus = nf90_inq_varid(lun,'rsno',icode)
    call chkncstatus(istatus,nf90_noerr,'rsno code',trim(schym(18)) &
       ,'')
    istatus = nf90_get_var(lun,icode,snow,istart(1:3), &
                icount(1:3))
    call chkncstatus(istatus,nf90_noerr,'Reading rsno',trim(schym(18)) &
       ,'')
    call writetime('Initializing snow accumulation for',idate)
    istatus = nf90_inq_varid(lun,'rdgw',icode)
    call chkncstatus(istatus,nf90_noerr,'rdgw code',trim(schym(18)) &
       ,'')
    istatus = nf90_get_var(lun,icode,deepw,istart(1:3), &
                icount(1:3))
    call chkncstatus(istatus,nf90_noerr,'Reading rdgw',trim(schym(18)) &
       ,'')
    call writetime('Initializing deep ground water for',idate)
    istatus = nf90_close(lun)
    call chkncstatus(istatus,nf90_noerr,'close',trim(schym(18)),'')
  end subroutine read_restart_dyn_NC

  subroutine closefile
    implicit none
    iretval = nf90_close(ncid)
    call checknetcdfresult('close output file')
    deallocate(vals)
    deallocate(map_ids)
  end subroutine closefile

  subroutine closerstfile
    implicit none
    iretval = nf90_close(ncidrst)
    call checknetcdfresult('close restart file')
    deallocate(maprst_ids)
  end subroutine closerstfile


  subroutine read_restart_stat_NC
    integer :: lun,nlon1,nlat1,nchymp1
    real :: slon,slat,dij1
    if (myid == 0) then
      call getlun(lun)
      istatus = nf90_open(trim(schym(12)),nf90_nowrite,lun)
      call chkncstatus(istatus,nf90_noerr,'open',trim(schym(12)),'')
!      istatus = nf90_inq_dimid(lun,'nchymp',icode)
!      call chkncstatus(istatus,nf90_noerr,'nchymp code',trim(schym(12)),'')
!      istatus = nf90_inquire_dimension(lun,icode,len=nchymp1)
!      call chkncstatus(istatus,nf90_noerr,'nchymp dim',trim(schym(12)),'')
      istatus = nf90_inq_dimid(lun,'lat',icode)
      call chkncstatus(istatus,nf90_noerr,'nlat code',trim(schym(12)),'')
      istatus = nf90_inquire_dimension(lun,icode,len=nlat1)
      call chkncstatus(istatus,nf90_noerr,'nlat dim',trim(schym(12)),'')
      istatus = nf90_inq_dimid(lun,'lon',icode)
      call chkncstatus(istatus,nf90_noerr,'nlon code',trim(schym(12)),'')
      istatus = nf90_inquire_dimension(lun,icode,len=nlon1)
      call chkncstatus(istatus,nf90_noerr,'nlon dim',trim(schym(12)),'')
      istatus = nf90_get_att(lun,nf90_global,'Start_longitude',slon)
      call chkncstatus(istatus,nf90_noerr,'Start_longitude',    &
        trim(schym(12)),'')
      istatus = nf90_get_att(lun,nf90_global,'Start_latitude',slat)
      call chkncstatus(istatus,nf90_noerr,'Start_latitude',trim(schym(12)),'')
      istatus = nf90_get_att(lun,nf90_global, &
            'approximate_lat-lon_resolution_in_meters',dij1)
      call chkncstatus(istatus,nf90_noerr,   &
        'approximate_lat-lon_resolution_in_meters',trim(schym(12)),'')
!      istatus = nf90_inq_varid(lun,'mchym',icode)
!      call chkncstatus(istatus,nf90_noerr,'mchym',trim(schym(12)),'')
!      istatus = nf90_get_var(lun,icode,mchym)
!      call chkncstatus(istatus,nf90_noerr,'mchym',trim(schym(12)),'')
!      istatus = nf90_inq_varid(lun,'rchym',icode)
!      call chkncstatus(istatus,nf90_noerr,'rchym',trim(schym(12)),'')
!      istatus = nf90_get_var(lun,icode,rchym)
!      call chkncstatus(istatus,nf90_noerr,'rchym',trim(schym(12)),'')
      if (nlon1.ne.mchym(2)) then
         call chymerror(31,0,0.0,'NLON')
      else if (nlat1.ne.mchym(3)) then
         call chymerror(31,0,0.0,'NLAT')
!      else if (nchymp1.ne.mchym(19)) then
!         call chymerror(31,0,0.0,'NCHYMP')
      else if (abs(slon-rchym(1))>0.0001) then
         call chymerror(31,0,0.0,'SLON')
      else if (abs(slat-rchym(2))>0.0001) then
         call chymerror(31,0,0.0,'SLAT')
      else if (abs(dij1-rchym(5))>0.00001) then
         call chymerror(31,0,0.0,'DIJ')
      endif
!      istatus = nf90_inq_varid(lun,'lat',icode)
!      call chkncstatus(istatus,nf90_noerr,'lat',trim(schym(12)),'')
!      istatus = nf90_get_var(lun,icode,lat)
!      call chkncstatus(istatus,nf90_noerr,'lat',trim(schym(12)),'')
!      istatus = nf90_inq_varid(lun,'lon',icode)
!      call chkncstatus(istatus,nf90_noerr,'lon',trim(schym(12)),'')
!      istatus = nf90_get_var(lun,icode,lon)
!      call chkncstatus(istatus,nf90_noerr,'lon',trim(schym(12)),'')
      istatus = nf90_inq_varid(lun,'dem',icode)
      call chkncstatus(istatus,nf90_noerr,'dem',trim(schym(12)),'')
      istatus = nf90_get_var(lun,icode,dem)
      call chkncstatus(istatus,nf90_noerr,'dem',trim(schym(12)),'')
      istatus = nf90_inq_varid(lun,'fdm',icode)
      call chkncstatus(istatus,nf90_noerr,'fdm',trim(schym(12)),'')
      istatus = nf90_get_var(lun,icode,fmap)
      call chkncstatus(istatus,nf90_noerr,'fdm',trim(schym(12)),'')
      istatus = nf90_inq_varid(lun,'acc',icode)
      call chkncstatus(istatus,nf90_noerr,'acc',trim(schym(12)),'')
      istatus = nf90_get_var(lun,icode,accl)
      call chkncstatus(istatus,nf90_noerr,'acc',trim(schym(12)),'')
      istatus = nf90_inq_varid(lun,'lus',icode)
      call chkncstatus(istatus,nf90_noerr,'lus',trim(schym(12)),'')
      istatus = nf90_get_var(lun,icode,luse)
      call chkncstatus(istatus,nf90_noerr,'lus',trim(schym(12)),'')
      istatus = nf90_inq_varid(lun,'aer',icode)
      call chkncstatus(istatus,nf90_noerr,'aer',trim(schym(12)),'')
      istatus = nf90_get_var(lun,icode,area)
      call chkncstatus(istatus,nf90_noerr,'aer',trim(schym(12)),'')
      istatus = nf90_inq_varid(lun,'dra',icode)
      call chkncstatus(istatus,nf90_noerr,'dra',trim(schym(12)),'')
      istatus = nf90_get_var(lun,icode,drai)
      call chkncstatus(istatus,nf90_noerr,'dra',trim(schym(12)),'')
      istatus = nf90_inq_varid(lun,'run',icode)
      call chkncstatus(istatus,nf90_noerr,'run',trim(schym(12)),'')
      istatus = nf90_get_var(lun,icode,runt)
      call chkncstatus(istatus,nf90_noerr,'run',trim(schym(12)),'')
      istatus = nf90_close(lun)
      call chkncstatus(istatus,nf90_noerr,'close',trim(schym(12)),'')
    end if
    call mpi_bcast(dem(1,1),nlon*nlat,MPI_REAL, 0,mycomm,mpierr)
    call mpi_bcast(fmap(1,1),nlon*nlat,MPI_REAL, 0,mycomm,mpierr)
    call mpi_bcast(accl(1,1),nlon*nlat,MPI_REAL, 0,mycomm,mpierr)
    call mpi_bcast(luse(1,1),nlon*nlat,MPI_REAL, 0,mycomm,mpierr)
    call mpi_bcast(area(1,1),nlon*nlat,MPI_REAL, 0,mycomm,mpierr)
    call mpi_bcast(drai(1,1),nlon*nlat,MPI_REAL, 0,mycomm,mpierr)
    call mpi_bcast(runt(1,1),nlon*nlat,MPI_REAL, 0,mycomm,mpierr)
    call mpi_barrier(mycomm,mpierr)
  end subroutine read_restart_stat_NC

  subroutine getasterdata(cfile, ndat, h, ifound)
    implicit none
    character(len=150) :: cfile
    integer :: ndat, lun, demlun, istatus, ifound
    integer(2) :: h(ndat,ndat)
    istatus = nf90_open(cfile,nf90_nowrite,lun)
    if (istatus /= nf90_noerr) then
       write(6,'(/,14x,a,a,a,a)')  'File '//cfile// &
        'not present in the Aster database'// &
        'continue running anyway'
        ifound = 1
        return
    end if
    call chkncstatus(istatus,nf90_noerr,'open',cfile, "")
    istatus = nf90_inq_varid(lun,'Band1',demlun)
    call chkncstatus(istatus,nf90_noerr,'dem inquire',cfile, "")
    istatus = nf90_get_var(lun,demlun,h)
    call chkncstatus(istatus,nf90_noerr,'dem data',cfile, "")
    istatus = nf90_close(lun)
    call chkncstatus(istatus,nf90_noerr,'close',cfile, "")
  end subroutine getasterdata

  subroutine gethydrodata(cfile, nlond, nlatd, h, lond, latd, ifound)
    implicit none
    character(len=150) :: cfile
    integer :: nlond, nlatd, lun, demlun, istatus , icode, ifound
!    integer(2) :: h(ndat,ndat)
    integer, allocatable, dimension(:,:) :: h
    real, allocatable, dimension(:) :: latd, lond
    istatus = nf90_open(cfile,nf90_nowrite,lun)
    if (istatus /= nf90_noerr) then
       write(6,'(/,14x,a,a,a,a)')  'File '//cfile// &
        'not present in the Hydroshed database'// &
        'continue running anyway'
        ifound = 1
        return
    end if
    call chkncstatus(istatus,nf90_noerr,'open',cfile, "")
    istatus = nf90_inq_dimid(lun,'lat',icode)
    call chkncstatus(istatus,nf90_noerr,'nlat code',cfile,'')
    istatus = nf90_inquire_dimension(lun,icode,len=nlatd)
    call chkncstatus(istatus,nf90_noerr,'nlat dim',cfile,'')
    istatus = nf90_inq_dimid(lun,'lon',icode)
    call chkncstatus(istatus,nf90_noerr,'nlon code',cfile,'')
    istatus = nf90_inquire_dimension(lun,icode,len=nlond)
    call chkncstatus(istatus,nf90_noerr,'nlon dim',cfile,'')
    if (allocated(h)) deallocate(h)
    if (allocated(lond)) deallocate(lond)
    if (allocated(latd)) deallocate(latd)
    allocate(h(nlond,nlatd))
    allocate(lond(nlond))
    allocate(latd(nlatd))
    istatus = nf90_inq_varid(lun,'Band1',demlun)
    call chkncstatus(istatus,nf90_noerr,'dem inquire',cfile, "")
    istatus = nf90_get_var(lun,demlun,h)
    call chkncstatus(istatus,nf90_noerr,'dem data',cfile, "")
    istatus = nf90_inq_varid(lun,'lon',demlun)
    call chkncstatus(istatus,nf90_noerr,'lon inquire',cfile, "")
    istatus = nf90_get_var(lun,demlun,lond)
    call chkncstatus(istatus,nf90_noerr,'lon data',cfile, "")
    istatus = nf90_inq_varid(lun,'lat',demlun)
    call chkncstatus(istatus,nf90_noerr,'lat inquire',cfile, "")
    istatus = nf90_get_var(lun,demlun,latd)
    call chkncstatus(istatus,nf90_noerr,'lat data',cfile, "")
    istatus = nf90_close(lun)
    call chkncstatus(istatus,nf90_noerr,'close',cfile, "")
  end subroutine gethydrodata

  subroutine gethydrodata_int(cfile, nlond, nlatd, h, lond, latd, ifound)
    implicit none
    character(len=150) :: cfile
    integer :: nlond, nlatd, lun, demlun, istatus , icode, ifound, ilon, ilat
!    integer(2) :: h(ndat,ndat)
    integer, allocatable, dimension(:,:) :: h
    double precision, allocatable, dimension(:) :: latd, lond
    istatus = nf90_open(cfile,nf90_nowrite,lun)
    if (istatus /= nf90_noerr) then
       write(6,'(/,14x,a,a,a,a)')  'File '//cfile// &
        'not present in the Hydroshed database'// &
        'continue running anyway'
        ifound = 1
        return
    end if
    call chkncstatus(istatus,nf90_noerr,'open',cfile, "")
    istatus = nf90_inq_dimid(lun,'lat',icode)
    call chkncstatus(istatus,nf90_noerr,'nlat code',cfile,'')
    istatus = nf90_inquire_dimension(lun,icode,len=nlatd)
    call chkncstatus(istatus,nf90_noerr,'nlat dim',cfile,'')
    istatus = nf90_inq_dimid(lun,'lon',icode)
    call chkncstatus(istatus,nf90_noerr,'nlon code',cfile,'')
    istatus = nf90_inquire_dimension(lun,icode,len=nlond)
    call chkncstatus(istatus,nf90_noerr,'nlon dim',cfile,'')
    if (allocated(h)) deallocate(h)
    if (allocated(lond)) deallocate(lond)
    if (allocated(latd)) deallocate(latd)
    allocate(h(mchym(2),mchym(3)))
    allocate(lond(nlond))
    allocate(latd(nlatd))
    istatus = nf90_inq_varid(lun,'lon',demlun)
    call chkncstatus(istatus,nf90_noerr,'lon inquire',cfile, "")
    istatus = nf90_get_var(lun,demlun,lond)
    call chkncstatus(istatus,nf90_noerr,'lon data',cfile, "")
    istatus = nf90_inq_varid(lun,'lat',demlun)
    call chkncstatus(istatus,nf90_noerr,'lat inquire',cfile, "")
    istatus = nf90_get_var(lun,demlun,latd)
    call chkncstatus(istatus,nf90_noerr,'lat data',cfile, "")
    if (trim(cfile) == "museo/HydroSHEDS/dir/011/dir_EuroCORDEX_small_final_011.nc" &
        .or. trim(cfile) == "museo/HydroSHEDS/dir/011/dir_Africa_011.nc" &
        .or. trim(cfile) == "output/Islanda_1kmHYDROSHED_scaleIslanda_remapcon_100_1000_21_011degree.nc" &
        .or. trim(cfile) == "output/EuroCORDEX_006degreeEUROPE_DEM_0.06degree_fdm.nc" &
        .or. trim(cfile) == "output/AfricaCORDEX_006degree_stkAfrica_DEM_0.06degree_fdm.nc" &
        .or. trim(cfile) == "output/NAMCORDEX_006degree_stkNAM_DEM_0.06degree_fdm.nc" &
        .or. trim(cfile) == "output/CAMCORDEX_006degree_stkCAM_DEM_0.06degree_fdm.nc" &
        .or. trim(cfile) == "output/SAMCORDEX_006degree_stkSAM_DEM_0.06degree_fdm.nc" &
        .or. trim(cfile) == "output/WASCORDEX_006degree_stkWAS_DEM_0.06degree_fdm.nc" &
        .or. trim(cfile) == "output/AUSCORDEX_006degree_stkAUS_DEM_0.06degree_fdm.nc" &
        .or. trim(cfile) == "output/SEACORDEX_006degree_stkSEA_DEM_0.06degree_fdm.nc" &
        .or. trim(cfile) == "output/EASCORDEX_006degree_stkEAS_DEM_0.06degree_fdm.nc" &
        .or. trim(cfile) == "output/NAMCORDEXWORLD_006degree_stkNA_CORDEXWORLD_DEM_0.06degree_fdm.nc" &
        .or. trim(cfile) == "output/SAMCORDEXWORLD_006degree_stkSA_CORDEXWORLD_DEM_0.06degree_fdm.nc" &
        .or. trim(cfile) == "output/EURCORDEXWORLD_006degree_stkEU_CORDEXWORLD_DEM_0.06degree_fdm.nc" &
        .or. trim(cfile) == "output/AFRCORDEXWORLD_006degree_stkAF_CORDEXWORLD_DEM_0.06degree_fdm.nc" &
        .or. trim(cfile) == "output/fdm_SAS.nc") then
      istatus = nf90_inq_varid(lun,'fdm',demlun)
    else
      istatus = nf90_inq_varid(lun,'Band1',demlun)
    end if
    call chkncstatus(istatus,nf90_noerr,'dem inquire',cfile, "")
!    istatus = nf90_get_var(lun,demlun,h)
    ilon=minloc(abs(lond-rchym(1)),dim=1)
    ilat=minloc(abs(latd-rchym(2)),dim=1)
    print*,ilon,ilat
    istatus = nf90_get_var(lun, demlun, h, start = (/ ilon, ilat /), count = (/mchym(2), mchym(3)/) )
    call chkncstatus(istatus,nf90_noerr,'dem data',cfile, "")
    istatus = nf90_close(lun)
    call chkncstatus(istatus,nf90_noerr,'close',cfile, "")
  end subroutine gethydrodata_int


  subroutine gethydrodata_real(cfile, nlond, nlatd, h, lond, latd, ifound)
    implicit none
    character(len=150) :: cfile
    integer :: nlond, nlatd, lun, demlun, istatus , icode, ifound
!    integer(2) :: h(ndat,ndat)
    real, allocatable, dimension(:,:) :: h
    double precision, allocatable, dimension(:) :: latd, lond
    istatus = nf90_open(cfile,nf90_nowrite,lun)
    if (istatus /= nf90_noerr) then
       write(6,'(/,14x,a,a,a,a)')  'File '//cfile// &
        'not present in the Hydroshed database'// &
        'continue running anyway'
        ifound = 1
        return
    end if
    call chkncstatus(istatus,nf90_noerr,'open',cfile, "")
    istatus = nf90_inq_dimid(lun,'lat',icode)
    call chkncstatus(istatus,nf90_noerr,'nlat code',cfile,'')
    istatus = nf90_inquire_dimension(lun,icode,len=nlatd)
    call chkncstatus(istatus,nf90_noerr,'nlat dim',cfile,'')
    istatus = nf90_inq_dimid(lun,'lon',icode)
    call chkncstatus(istatus,nf90_noerr,'nlon code',cfile,'')
    istatus = nf90_inquire_dimension(lun,icode,len=nlond)
    call chkncstatus(istatus,nf90_noerr,'nlon dim',cfile,'')
    if (allocated(h)) deallocate(h)
    if (allocated(lond)) deallocate(lond)
    if (allocated(latd)) deallocate(latd)
    allocate(h(nlond,nlatd))
    allocate(lond(nlond))
    allocate(latd(nlatd))
    istatus = nf90_inq_varid(lun,'Band1',demlun)
    call chkncstatus(istatus,nf90_noerr,'dem inquire',cfile, "")
    istatus = nf90_get_var(lun,demlun,h)
    call chkncstatus(istatus,nf90_noerr,'dem data',cfile, "")
    istatus = nf90_inq_varid(lun,'lon',demlun)
    call chkncstatus(istatus,nf90_noerr,'lon inquire',cfile, "")
    istatus = nf90_get_var(lun,demlun,lond)
    call chkncstatus(istatus,nf90_noerr,'lon data',cfile, "")
    istatus = nf90_inq_varid(lun,'lat',demlun)
    call chkncstatus(istatus,nf90_noerr,'lat inquire',cfile, "")
    istatus = nf90_get_var(lun,demlun,latd)
    call chkncstatus(istatus,nf90_noerr,'lat data',cfile, "")
    istatus = nf90_close(lun)
    call chkncstatus(istatus,nf90_noerr,'close',cfile, "")
  end subroutine gethydrodata_real

  subroutine ein75readnc(nlon,nlat,ora,giorno,mese,anno,ifl,pi, &
              la,lo,n)
    implicit none
    integer :: ora,giorno,mese,anno,ifl,ni,log,n,i,j,istatus,rec,itime
    integer :: lunt,len_trim,na,nrec
    real, dimension(nlon*nlat) :: pi, la, lo
    integer :: nlon,nlat,tcode,xcode
    logical :: bisest
    logical first
    data first /.true./
    save first
    integer :: xmese
    integer :: n2dim
    parameter (n2dim=120000)
    character(len=8) :: adata
    character(len=60) :: dirt
    character(len=132) :: cfl, oldcflt
    save :: tcode,lunt,oldcflt,dirt,xcode
    real, dimension(nlon,nlat) :: tmpm, temp
    real, dimension(nlon) :: xlon
    real, dimension(nlat) :: xlat
    real :: fill_value,tmp,offset,scalef
    integer, dimension(4) :: istart, icount
    xmese=0
    if (mese>1) xmese=julianday(monthlen(mese-1,anno), &
       mese-1,anno)
!    if (nlon*nlat>n2dim) then
!       write(6,'(/,14x,a,a)')  'Too big Era domain, increase '// &
!        'n2dim parameter inside ein75readnc'
!       call exit(1)
!    end if
    write(adata,'(i8)') anno*10000+mese*100+giorno
    if (first) then
      call mvgetiflags(70,log)
      if (log.le.0.or.log.ge.100) log=6
      if (mchym(15) == 6) then
        dirt="./museo/TEMP/EIN75/"
        write(cfl,'(a,i4,a,i4,a)') 'erai_tas_', &
           anno,'0101_',anno,'1231_00_0.75.nc'
        write(log,'(18x,a,a,i4,a,i4,a)')'Opening Era file ERAIN75', &
         'erai_tas_',anno,'0101_',anno,'1231_00_0.75.nc'
      else if (mchym(15) == 7 .or. mchym(15) == 8) then
        dirt="./museo/TEMP/ERA5/"
        write(cfl,'(a,a,a,a,a)') 'tas_',adata(1:4),'_',adata(5:6),'.nc'
        write(log,'(18x,a,a)')'Opening Era file ERA5', cfl
      end if
      call getlun(lunt)
      istatus = nf90_open(trim(dirt)//cfl,nf90_nowrite,lunt)
      call chkncstatus(istatus,nf90_noerr,'open',dirt,cfl)
      oldcflt=cfl
      istatus = nf90_inq_varid(lunt,'t2m',tcode)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_varid(lunt,'tas',tcode)
        call chkncstatus(istatus,nf90_noerr,'temp code',dirt,cfl)
      end if
      istatus = nf90_inq_varid(lunt,'time',itime)
      call chkncstatus(istatus,nf90_noerr,'time',dirt,cfl)
      first=.false.
    end if
    if (mchym(15) == 6) then
      write(cfl,'(a,i4,a,i4,a)') 'erai_tas_', &
       anno,'0101_',anno,'1231_00_0.75.nc'
    else if (mchym(15) == 7) then
      write(cfl,'(a,a,a,a,a)') 'tas_',adata(1:4),'_',adata(5:6),'.nc'
    end if
    if (cfl(1:len_trim(cfl)).ne.oldcflt(1:len_trim(oldcflt))) then
      istatus = nf90_close(lunt)
      call chkncstatus(istatus,nf90_noerr,'close',trim(dirt),oldcflt)
      call getlun(lunt)
      write(6,'(18x,a)') 'Opening ERA file '//cfl(1:len_trim(cfl))
      istatus = nf90_open(trim(dirt)//cfl,nf90_nowrite,lunt)
      call chkncstatus(istatus,nf90_noerr,'open',trim(dirt),cfl)
      oldcflt=cfl
    endif
    istart(1)=1 ; istart(2)=1 ; icount(1)=nlon ; icount(2)=nlat
    if (mchym(15) == 6) then
      nrec=(xmese*4+(giorno-1)*4+(ora/6))+1
    else if (mchym(15) == 7) then
      nrec=(giorno-1)*24+(ora)+1
    end if
    if (nrec.eq.0) nrec=1
    if (ifl.eq.1) then
      istart(3)=nrec
      icount(3)=1
      istatus = nf90_get_var(lunt,tcode,tmpm,istart(1:3), &
              icount(1:3))
      call chkncstatus(istatus,nf90_noerr,'temperature data', &
          dirt,cfl)
      istatus = nf90_get_att(lunt,tcode,'add_offset', &
             tmp)
      if ( istatus == nf90_noerr ) then
        offset = tmp
        istatus = nf90_get_att(lunt,tcode,'scale_factor',tmp)
        call chkncstatus(istatus,nf90_noerr,'scale factor',dirt,cfl)
        scalef = tmp
      else
        offset = 0.0
        scalef = 1.0
      end if
      temp = tmpm * scalef + offset - 273.15
      call erashiftlon(temp,nlon,nlat)
      call erashiftlat(temp,nlon,nlat)
      call regcmmat2vec(temp,nlon,nlat,1.,pi)
    else
      write (6,'(16x,a)') '  Bad flag passed to ein75readnc.Exiting...'
      call exit(1)
    endif
    n=nlon*nlat
    na = 0
    if (mchym(15) == 6) then
      do i=1,480
        do j=1,241
          na=na+1
          la(na)=eralatold(i,j) ; lo(na)=eralonold(i,j)
        enddo
      enddo
    else if (mchym(15) == 7) then
      do i=1,1440
        do j=1,721
          na=na+1
          la(na)=eralat5(i,j) ; lo(na)=eralon5(i,j)
        enddo
      enddo
    else if (mchym(15) == 8) then
      istatus = nf90_inq_varid(lunt,'lon',xcode)
      istatus = nf90_get_var(lunt,xcode,xlon)
      istatus = nf90_inq_varid(lunt,'lat',xcode)
      istatus = nf90_get_var(lunt,xcode,xlat)
      do i=1,nlon
        do j=1,nlat
          na=na+1
          la(na)=xlat(j) ; lo(na)=xlon(i)
        enddo
      enddo
    end if
    return
  end subroutine ein75readnc

  subroutine persiannreadncll(pnlon,pnlat,ora,giorno,mese,anno,ifl &
              ,pi,la,lo,n)
    implicit none
    integer :: ora,giorno,mese,anno,ifl,n,na
    integer :: pnlon,pnlat,i,j,itime
    real, dimension(pnlon*pnlat) :: pi, la, lo
    logical first
    data first /.true./
    save first
    character(len=26) :: dir
    character(len=132) :: cfl, oldcfl
    character(len=64) :: rainmis
    save rainmis,oldcfl
    integer :: len_trim,log,lun,code,pcode,ld,istatus,rec
    integer, dimension(4) :: istart,icount
    save log,lun,pcode,istart,icount,ld,dir
    real, dimension(pnlon,pnlat) :: tmpm, raini
    integer, parameter :: n2dim=700000
    real, dimension(n2dim) :: loclat,loclon
    real :: pfact
    integer :: xmese
    save loclat,loclon,pfact
    if (anno > 2010 .or. anno < 2000) then
      write(6,'(/,14x,a,a)') "Please insert a year between 2000", &
         " and 2010 for PERSIANN dataset"
      call exit(1)
    end if
    xmese=0
    if (mese.gt.1) xmese=julianday(monthlen(mese-1,anno), &
         mese-1,anno)
    if (anno == 2000) then
      xmese=julianday(monthlen(mese-1,anno), &
        mese-1,anno)-(31+29)
      if (mese <= 2) then
        write(6,'(/,14x,a)')'Persiann Dataset begin at '// &
          '2000/03/01 00:00'
        call exit(1)
      end if
    end if
    dir="./museo/PREC/PERSIANN/"

    if (pnlon*pnlat.gt.n2dim) then
      write(6,'(/,14x,a,a)') 'Too big PERSIANN domain, increase '// &
       'n2dim parameter inside persiannreadncll'
      call exit(1)
    endif
    if (first) then
      call mvgetiflags(70,log)
      if (log.le.0.or.log.ge.100) log=6
      write(cfl,'(a,i4,a)') 'PERSIANN-HOURLY_',anno,'.nc'
      call getlun(lun)
      write(log,'(18x,a,a,i4,a)') 'Opening Persiann file ', &
       'PERSIANN-HOURLY_',anno,'.nc'
      istatus = nf90_open(trim(dir)//cfl,nf90_nowrite,lun)
      call chkncstatus(istatus,nf90_noerr,'open',dir,cfl)
      oldcfl=cfl
      tmpm = 0
      istatus = nf90_inq_varid(lun,'pre',pcode)
      if (istatus.ne.nf90_noerr) istatus = nf90_inq_varid(lun,'pre', &
         pcode)
      call chkncstatus(istatus,nf90_noerr,'rain code',dir,cfl)
      rainmis=' '
      istatus = nf90_get_att(lun,pcode,'units',rainmis)
      call chkncstatus(istatus,nf90_noerr,'rain units',dir,cfl)
      pfact=1
      istatus = nf90_inq_varid(lun,'time',itime)
      call chkncstatus(istatus,nf90_noerr,'time',dir,cfl)
      istart(1)=1 ; istart(2)=1 ; icount(1)=pnlon ; icount(2)=pnlat
      first=.false.
    endif
    write(cfl,'(a,i4,a)') 'PERSIANN-HOURLY_',anno,'.nc'
    if (cfl(1:len_trim(cfl)).ne.oldcfl(1:len_trim(oldcfl))) then
      istatus = nf90_close(lun)
      call chkncstatus(istatus,nf90_noerr,'close',trim(dir),oldcfl)
      call getlun(lun)
      write(log,'(18x,a)') 'Opening PERSIANN file '// &
                          cfl(1:len_trim(cfl))
      istatus = nf90_open(trim(dir)//cfl,nf90_nowrite,lun)
      call chkncstatus(istatus,nf90_noerr,'open',trim(dir),cfl)
      oldcfl=cfl
    endif
    if (anno == 2000) then
      rec=(xmese*24+(giorno-1)*24+ora)+1
      if (rec.eq.0) rec=1
    else
      rec=(xmese*24+(giorno-1)*24+ora)+1 ; if (rec.eq.0) rec=1
    end if
    if (ifl.eq.1) then
      istart(3)=rec ; icount(3)=1
      istatus = nf90_get_var(lun,pcode,raini,istart(1:3), &
              icount(1:3))
      call chkncstatus(istatus,nf90_noerr,'rain data',dir,cfl)
      call erashiftlon(raini,pnlon,pnlat)
      call regcmmat2vec(raini,pnlon,pnlat,pfact,pi)
    else
      write (6,'(16x,a)') '  Bad flag passed to persiannreadnc.' &
              //'Exiting...'
      call exit(1)
    endif
    n=pnlon*pnlat
    na = 0
    do i=1,1440 ; do j=1,480
      na=na+1
      la(na)=persiannlat0_25(i,j)
      lo(na)=persiannlon0_25(i,j)
    enddo ; enddo
    call mvsetflags('Calendar',0.0)
  end subroutine persiannreadncll

  subroutine trimmreadncll(tnlon,tnlat,ora,giorno,mese,anno,ifl,pi &
             ,la,lo,n)
    implicit none
    integer :: ora,giorno,mese,anno,ifl,n,na,nrec,iora,igiorno,imese
    integer :: ianno
    integer :: tnlon,tnlat,i,j,itime,rdata
    real, dimension(tnlon*tnlat) :: pi, la, lo
    logical :: first
    data first /.true./
    save first
    character(len=60) :: dir
    character(len=132) :: cfl, oldcfl
    character(len=64) :: rainmis
    save rainmis,oldcfl
    integer :: len_trim,log,lun,code,pcode,scode,ld,istatus,rec
    integer, dimension(4) :: istart, icount
    save log,lun,pcode,istart,icount,scode,ld,dir
    real, dimension(tnlon,tnlat) :: tmpm, raini
    integer, parameter :: n2dim=576000
    real, dimension(n2dim) :: loclat, loclon
    real :: pfact
    integer :: xmese
    save loclat,loclon,pfact
    xmese = monthlen(mese,anno)
    imese = mese
    ianno = anno
    igiorno = giorno
    if (ora < 3 ) then
       iora = 3
    else if (ora > 2 .and. ora < 6 ) then
       iora = 6
    else if (ora > 5 .and. ora < 9 ) then
       iora = 9
    else if (ora > 8 .and. ora < 12 ) then
       iora = 12
    else if (ora > 11 .and. ora < 15 ) then
       iora = 15
    else if (ora > 14 .and. ora < 18 ) then
       iora = 18
    else if (ora > 17 .and. ora < 21 ) then
       iora = 21
    else if (ora > 20  ) then
       iora = 00
       if (giorno >= xmese) then
         igiorno = 1
           if (imese == 12) then
             imese = 1
             ianno = ianno + 1
           else
             imese = imese + 1
           end if
       else
         igiorno = igiorno + 1
       end if
    end if
    if (tnlon*tnlat.gt.n2dim) then
       write(6,'(/,14x,a,a)')  'Too big Trimm domain, increase '// &
         'n2dim parameter inside trimmreadncll'
       call exit(1)
    end if
    call mvgetiflags(70,log) ; if (log.le.0.or.log.ge.100) log=6
    write(dir,'(a,i4,a)')'./museo/PREC/TRMM/',ianno,'/'
    if (iora < 10) then
       if (igiorno < 10) then
          if (imese < 10) then
            write(cfl,'(a,i4,a,i1,a,i1,a,i1,a)') 'trmm_3B42_', &
              ianno,'0',imese,'0',igiorno,'.0',iora,'.7.nc'
          else
            write(cfl,'(a,i4,i2,a,i1,a,i1,a)') 'trmm_3B42_', &
              ianno,imese,'0',igiorno,'.0',iora,'.7.nc'
          end if
       else
          if (imese < 10) then
            write(cfl,'(a,i4,a,i1,i2,a,i1,a)') 'trmm_3B42_', &
              ianno,'0',imese,igiorno,'.0',iora,'.7.nc'
          else
            write(cfl,'(a,i4,i2,i2,a,i1,a)') 'trmm_3B42_', &
              ianno,imese,igiorno,'.0',iora,'.7.nc'
          end if
       end if
    else
       if (igiorno < 10) then
          if (imese < 10) then
            write(cfl,'(a,i4,a,i1,a,i1,a,i2,a)') 'trmm_3B42_', &
              ianno,'0',imese,'0',igiorno,'.',iora,'.7.nc'
          else
            write(cfl,'(a,i4,i2,a,i1,a,i2,a)') 'trmm_3B42_', &
              ianno,imese,'0',igiorno,'.',iora,'.7.nc'
          end if
       else
          if (imese < 10) then
            write(cfl,'(a,i4,a,i1,i2,a,i2,a)') 'trmm_3B42_', &
              ianno,'0',imese,igiorno,'.',iora,'.7.nc'
          else
            write(cfl,'(a,i4,i2,i2,a,i2,a)') 'trmm_3B42_', &
              ianno,imese,igiorno,'.',iora,'.7.nc'
          end if
       end if
    end if
    call getlun(lun)
    pfact=1.0
    do i=1,5
      istatus = nf90_open(trim(dir)//cfl,nf90_nowrite,lun)
      if (istatus == nf90_noerr) exit
    end do
    call chkncstatus(istatus,nf90_noerr,'open',dir,cfl)
    if (first) then
      oldcfl=cfl
      istatus = nf90_inq_varid(lun,'precip',pcode)
      call chkncstatus(istatus,nf90_noerr,'rain code',dir,cfl)
      istatus = nf90_inq_varid(lun,'time',itime)
      call chkncstatus(istatus,nf90_noerr,'time',dir,cfl)
      first = .false.
    end if
    if (cfl(1:len_trim(cfl)).ne.oldcfl(1:len_trim(oldcfl))) then
      istatus = nf90_close(lun)
      call chkncstatus(istatus,nf90_noerr,'close',trim(dir),oldcfl)
      call getlun(lun)
      write(log,'(18x,a)') 'Opening Trimm file '// &
         cfl(1:len_trim(cfl))
      do i=1,5
        istatus = nf90_open(trim(dir)//cfl,nf90_nowrite,lun)
        if (istatus == nf90_noerr) exit
      end do
      call chkncstatus(istatus,nf90_noerr,'open',trim(dir),cfl)
      istatus = nf90_inq_varid(lun,'precip',pcode)
      call chkncstatus(istatus,nf90_noerr,'rain code',dir,cfl)
      istatus = nf90_inq_varid(lun,'time',itime)
      call chkncstatus(istatus,nf90_noerr,'time',dir,cfl)
      oldcfl=cfl
    endif
    istart(1)=1 ; istart(2)=1 ; istart(3)=1
    icount(1)=tnlon ; icount(2)=tnlat ; icount(3)=1
    istatus = nf90_get_var(lun,pcode,tmpm,istart(1:3), &
      icount(1:3))
    call chkncstatus(istatus,nf90_noerr,'rain data',dir,cfl)
    !Added a close statement which was missing (AF 27jan15)
    istatus = nf90_close(lun)
    call chkncstatus(istatus,nf90_noerr,'close',trim(dir),cfl)
    !End edit
!    call erashiftlat(tmpm,tnlon,tnlat)
    call regcmmat2vec(tmpm,tnlon,tnlat,pfact,pi)
    n=tnlon*tnlat
    na = 0
    do i=1,1440
      do j=1,400
        na=na+1
        la(na)=trimmlat(i,j)
        lo(na)=trimmlon(i,j)
      enddo
    enddo
    return
  end subroutine trimmreadncll

  subroutine erareadncll(nlon,nlat,ora,giorno,mese,anno,ifl, &
             pi,la,lo,n,rdata)
    implicit none
    integer :: ora,giorno,mese,anno,ifl,n,na,nrec
    integer :: nlon,nlat,i,j,itime,rdata,ifin,jfin
    logical first
    data first /.true./
    save first
    character(len=60) :: dir
    character(len=132) :: cfl,oldcfl
    character(len=64) :: rainmis
    save rainmis,oldcfl
    integer :: len_trim,log,lun,code,pcode,ld,istatus,rec,ivarid
    integer,dimension(4) :: istart,icount
    save log,lun,pcode,istart,icount,ld,dir
!oldwindows    real, dimension(nlon,nlat) :: tmpm, raini !Fabio-nc
    real, allocatable, dimension(:,:) :: tmpm,raini !Fabio-nc
    real, dimension(nlon) :: glon !Fabio-nc
    real, dimension(nlat) :: glat !Fabio-nc
    integer, parameter :: n2dim=120000
    real, dimension(n2dim) :: loclat,loclon
    real :: pfact
    integer :: xmese, iband
    logical :: bisest
    real :: mancanti,fill_value,tmp,offset,scalef
    save loclat,loclon,pfact,fill_value,mancanti,gdomain
    type(global_domain) :: gdomain
    integer :: itile , iti , itf
    real, dimension(:) :: pi,la,lo
    xmese=0
!    if (ora == 0) then
!      ora = 23
!      if (giorno == 1) then
!        if (mese == 1) then
!          anno = anno - 1
!          mese = 12
!          giorno = 31
!        else
!          mese = mese - 1
!          giorno = monthlen(mese,anno)
!        end if
!      else
!        giorno = giorno - 1
!      end if
!    end if
    if (mese.gt.1) xmese=julianday(monthlen(mese-1,anno), &
       mese-1,anno)
    if (nlon*nlat.gt.n2dim) then
       write(6,'(/,14x,a,a)')  'Too big Era domain, increase '// &
        'n2dim parameter inside erareadncll'
       call exit(1)
    end if
    if (first) then
      call mvgetiflags(70,log) ; if (log.le.0.or.log.ge.100) log=6
      if (rdata == 8) then
        dir="./museo/PREC/EIN15/"
        write(cfl,'(a,i4,a,i4,a)') 'erai_precip_', &
          anno,'0101_',anno,'1231_00_1.50.nc'
        call getlun(lun)
        write(log,'(18x,a,a,i4,a,i4,a)')'Opening Era file ERAIN15', &
          'erai_precip_',anno,'0101_',anno,'1231_00_1.50.nc'
        pfact=1000./3.
      else if (rdata == 10) then
        dir="./museo/PREC/EIN75/"
        write(cfl,'(a,i4,a,i4,a)') 'erai_precip_', &
          anno,'0101_',anno,'1231_00_0.75.nc'
        write(log,'(18x,a,a,i4,a,i4,a)')'Opening Era file ERAIN75', &
          'erai_precip_',anno,'0101_',anno,'1231_00_0.75.nc'
        pfact=1000./3.
      end if
      istatus = nf90_open(trim(dir)//cfl,nf90_nowrite,lun)
      call chkncstatus(istatus,nf90_noerr,'open',dir,cfl)
      oldcfl=cfl
      istatus = nf90_inq_varid(lun,'latitude',ivarid)
      call chkncstatus(istatus,nf90_noerr,'Missing latitude variable in file ',dir,cfl)
      istatus = nf90_get_var(lun,ivarid,glat)
      call chkncstatus(istatus,nf90_noerr,'Error reading lat var in file ',dir,cfl)
      istatus = nf90_inq_varid(lun,'longitude',ivarid)
      call chkncstatus(istatus,nf90_noerr,'Missing longitude variable in file ',dir,cfl)
      istatus = nf90_get_var(lun,ivarid,glon)
      call chkncstatus(istatus,nf90_noerr,'Error reading lon var in file ',dir,cfl)
      !
      ! Find window to read
      !
      call get_window(glat,glon,lat,lon,iband,gdomain)

      allocate (tmpm(sum(gdomain%ni),gdomain%nj), raini(sum(gdomain%ni),gdomain%nj), &
                rainiold15(sum(gdomain%ni),gdomain%nj),rainiold75(sum(gdomain%ni),gdomain%nj))
      tmpm = 0
      raini = 0

      istatus = nf90_inq_varid(lun,'tp',pcode)
      if (istatus.ne.nf90_noerr) then
          istatus = nf90_inq_varid(lun,'TP',pcode)
      end if
      call chkncstatus(istatus,nf90_noerr,'rain code',dir,cfl)
      istatus = nf90_inq_varid(lun,'time',itime)
      call chkncstatus(istatus,nf90_noerr,'time',dir,cfl)

!      iti = 1
!      do itile = 1 , gdomain%ntiles
!        istart(1) = gdomain%igstart(itile)
!        icount(1) = gdomain%ni(itile)
!        ! Latitudes are reversed in original file
!        istart(2) = gdomain%jgstart
!        icount(2) = gdomain%nj
!        itf = iti + gdomain%ni(itile) - 1
!!        istatus = nf90_get_var(lun,pcode,tmpm(iti:itf,:,:),istart,icount)
!!        call chkncstatus(istatus,nf90_noerr,'Error read var tp era',dir,cfl)
!        iti = iti + gdomain%ni(itile)
!      print*,'gdomain%ntiles',gdomain%ntiles
!      print*,'gdomain%igstart(itile)',gdomain%igstart(itile)
!      print*,'gdomain%ni(itile)',gdomain%ni(itile)
!      print*,'gdomain%jgstart',gdomain%jgstart
!      print*,'gdomain%nj',gdomain%nj
!      print*,'itf',itf
!      print*,'iti',iti
!      end do
!
! Reg the old value of precipitation if hour is greater than 2
!
!oldwindows      istart(1)=1 ; istart(2)=1 ; icount(1)=nlon ; icount(2)=nlat
      if (ora <= 3) then
        if (rdata == 10) then
          rainiold75 = 0.
        else
          rainiold15 = 0.
        end if
      else
        if (.not.allocated(tmpm)) allocate (tmpm(sum(gdomain%ni),gdomain%nj))
        tmpm = 0
        rec=(xmese*8+(giorno-1)*8+((ora)/3))+1
        if (rec.eq.0) rec=1
        istart(3) = rec ; icount(3) = 1
!oldwindows        istatus = nf90_get_var(lun,pcode,tmpm,istart(1:3), &
!oldwindows            icount(1:3))
!oldwindows        call chkncstatus(istatus,nf90_noerr,'rain data',dir,cfl)
        iti = 1
        do itile = 1 , gdomain%ntiles
          istart(1) = gdomain%igstart(itile)
          icount(1) = gdomain%ni(itile)
          ! Latitudes are reversed in original file
          istart(2) = gdomain%jgstart
          icount(2) = gdomain%nj
          itf = iti + gdomain%ni(itile) - 1
!        print*,'itf',itf
!        print*,'iti',iti
!        print*,'istart(1)',istart(1)
!        print*,'icount(1)',icount(1)
!        print*,'istart(2)',istart(2)
!        print*,'icount(2)',icount(2)
          istatus = nf90_get_var(lun,pcode,tmpm(iti:itf,:),istart,icount)
          call chkncstatus(istatus,nf90_noerr,'Error read var tp era',dir,cfl)
          iti = iti + gdomain%ni(itile)
!        print*,'gdomain%ntiles',gdomain%ntiles
!        print*,'gdomain%igstart(itile)',gdomain%igstart(itile)
!        print*,'gdomain%ni(itile)',gdomain%ni(itile)
!        print*,'gdomain%jgstart',gdomain%jgstart
!        print*,'gdomain%nj',gdomain%nj
!        print*,'itf',itf
!        print*,'iti',iti
        end do
!       print*,"schiatti qui5?"

        istatus = nf90_get_att(lun,pcode,'add_offset', tmp)
        if ( istatus == nf90_noerr ) then
          offset = tmp
          istatus = nf90_get_att(lun,pcode,'scale_factor',tmp)
          call chkncstatus(istatus,nf90_noerr,'scale factor',dir,cfl)
          scalef = tmp
        else
          offset = 0.0
          scalef = 1.0
        end if
        tmpm = tmpm * scalef + offset
!       print*,"schiatti qui4?"
!        print*,maxval(tmpm)
        call erashiftlon(tmpm,sum(gdomain%ni),gdomain%nj)
        call erashiftlat(tmpm,sum(gdomain%ni),gdomain%nj)
        if (rdata == 8) then
           rainiold15=tmpm
        else if (rdata == 10) then
           rainiold75=tmpm
        end if
      end if
      first=.false.
    end if
    if (rdata == 8) then
      write(cfl,'(a,i4,a,i4,a)') 'erai_precip_', &
          anno,'0101_',anno,'1231_00_1.50.nc'
    end if
    if (rdata == 10) then
      write(cfl,'(a,i4,a,i4,a)') 'erai_precip_', &
          anno,'0101_',anno,'1231_00_0.75.nc'
    end if
    if (cfl(1:len_trim(cfl)).ne.oldcfl(1:len_trim(oldcfl))) then
       istatus = nf90_close(lun)
       call chkncstatus(istatus,nf90_noerr,'close',trim(dir),oldcfl)
       call getlun(lun)
       write(log,'(18x,a)') 'Opening ERA file '//cfl(1:len_trim(cfl))
       istatus = nf90_open(trim(dir)//cfl,nf90_nowrite,lun)
       call chkncstatus(istatus,nf90_noerr,'open',trim(dir),cfl)
       oldcfl=cfl
    endif
    nrec=(xmese*8+(giorno-1)*8+((ora)/3))+1
    if (nrec.eq.0) nrec=1
    if (ifl.eq.1) then
       if (.not.allocated(raini)) allocate (raini(sum(gdomain%ni),gdomain%nj))
       raini = 0
!OldFF         nrec = nrec-976 !FF
       istart(3)=nrec ; icount(3)=1
!oldwindows       istatus = nf90_get_var(lun,pcode,raini,istart(1:3), &
!oldwindows               icount(1:3))
!oldwindows       call chkncstatus(istatus,nf90_noerr,'rain data',dir,cfl)
       iti = 1
       do itile = 1 , gdomain%ntiles
         istart(1) = gdomain%igstart(itile)
         icount(1) = gdomain%ni(itile)
         ! Latitudes are reversed in original file
         istart(2) = gdomain%jgstart
         icount(2) = gdomain%nj
         itf = iti + gdomain%ni(itile) - 1
!        print*,'itf',itf
!        print*,'iti',iti
!        print*,'istart(1)',istart(1)
!        print*,'icount(1)',icount(1)
!        print*,'istart(2)',istart(2)
!        print*,'icount(2)',icount(2)

         istatus = nf90_get_var(lun,pcode,raini(iti:itf,:),istart,icount)
         call chkncstatus(istatus,nf90_noerr,'Error read var tp era',dir,cfl)
         iti = iti + gdomain%ni(itile)
!       print*,'gdomain%ntiles',gdomain%ntiles
!       print*,'gdomain%igstart(itile)',gdomain%igstart(itile)
!       print*,'gdomain%ni(itile)',gdomain%ni(itile)
!       print*,'gdomain%jgstart',gdomain%jgstart
!       print*,'gdomain%nj',gdomain%nj
!       print*,'itf',itf
!       print*,'iti',iti
       end do
!       print*,"schiatti qui6?"
!       print*,'pcode',pcode

       istatus = nf90_get_att(lun,pcode,'add_offset',tmp)
       if ( istatus == nf90_noerr ) then
         offset = tmp
         istatus = nf90_get_att(lun,pcode,'scale_factor',tmp)
         call chkncstatus(istatus,nf90_noerr,'scale factor',dir,cfl)
         scalef = tmp
       else
         offset = 0.0
         scalef = 1.0
       end if
!       print*,'offsetscalef',offset,scalef
       raini = raini * scalef + offset
!       print*,'maxval',maxval(raini)
!         print*,'ma davero schiatti qui?????'
!       call erashiftlon(raini,sum(gdomain%ni),gdomain%nj)
!         print*,'ma davero schiatti qui?????'
!       call erashiftlat(raini,sum(gdomain%ni),gdomain%nj)
!       print*,"schiatti qui2?"
       do i=1,sum(gdomain%ni) ; do j=1,gdomain%nj
         if ( raini(i,j) < 0 ) raini(i,j) = 0      !Soluzione temporanea a valori negativi di precipitazione dovuti probabilmente ad arrotondamenti
       end do ; end do
       if (ora <=3) then
         if (rdata == 10) then
           rainiold75 = 0.
         else if (rdata == 8) then
           rainiold15 = 0.
         end if
       else
!         if (.not.allocated(tmpm)) allocate (tmpm(sum(gdomain%ni),gdomain%nj))
!         tmpm = 0
!         rec=(xmese*8+(giorno-1)*8+((ora)/3))
!         if (rec.eq.0) rec=1
!!         print*,rec
!!OLDFF         rec = rec-976 !FF
!         istart(3) = rec ; icount(3) = 1
!!oldwindows         istatus = nf90_get_var(lun,pcode,tmpm,istart(1:3), &
!!oldwindows            icount(1:3))
!!oldwindows         call chkncstatus(istatus,nf90_noerr,'rain data',dir,cfl)
!         iti = 1
!         do itile = 1 , gdomain%ntiles
!           istart(1) = gdomain%igstart(itile)
!           icount(1) = gdomain%ni(itile)
!           ! Latitudes are reversed in original file
!           istart(2) = gdomain%jgstart
!           icount(2) = gdomain%nj
!           itf = iti + gdomain%ni(itile) - 1
!           istatus = nf90_get_var(lun,pcode,tmpm(iti:itf,:),istart,icount)
!           call chkncstatus(istatus,nf90_noerr,'Error read var tp era',dir,cfl)
!           iti = iti + gdomain%ni(itile)
!!         print*,'gdomain%ntiles',gdomain%ntiles
!!         print*,'gdomain%igstart(itile)',gdomain%igstart(itile)
!!         print*,'gdomain%ni(itile)',gdomain%ni(itile)
!!         print*,'gdomain%jgstart',gdomain%jgstart
!!         print*,'gdomain%nj',gdomain%nj
!!         print*,'itf',itf
!!         print*,'iti',iti
!         end do
!
!         istatus = nf90_get_att(lun,pcode,'add_offset',tmp)
!         call chkncstatus(istatus,nf90_noerr,'add_offset',dir,cfl)
!         offset = tmp
!         istatus = nf90_get_att(lun,pcode,'scale_factor',tmp)
!         call chkncstatus(istatus,nf90_noerr,'scale factor',dir,cfl)
!         scalef = tmp
!         tmpm = tmpm * scalef + offset
!!         print*,'maxval tmpm',maxval(tmpm)
!         call erashiftlon(tmpm,sum(gdomain%ni),gdomain%nj)
!         call erashiftlat(tmpm,sum(gdomain%ni),gdomain%nj)
!         raini = tmpm
!       print*,"schiatti qui1?"
         if (rdata == 8) then
           raini = raini - rainiold15
           rainiold15=tmpm
         else if (rdata == 10) then
           raini = raini - rainiold75
           rainiold75=tmpm
         end if
       end if
!       print*,"schiatti qui?"
!       if (allocated(pi)) deallocate(pi)
!       allocate(pi(sum(gdomain%ni)*gdomain%nj))
       call regcmmat2vec(raini,sum(gdomain%ni),gdomain%nj,pfact,pi)
       do i=1,sum(gdomain%ni) ; do j=1,gdomain%nj
         if ( raini(i,j) < 0 ) raini(i,j) = 0      !Soluzione temporanea a valori negativi di precipitazione dovuti probabilmente ad arrotondamenti
       end do ; end do
    else
       write (6,'(16x,a)') '  Bad flag passed to erareadnc.  &
               Exiting...'
       call exit(1)
    endif
    n=sum(gdomain%ni)*gdomain%nj
    na = 0
    if (rdata == 8) then
      do i=1,240 ; do j=1,121
        na=na+1
        la(na)=eralat1_5(i,j) ; lo(na)=eralon1_5(i,j)
      enddo ; enddo
    end if
!OLDFF    do i = 1 , nlon !FF
!OLDFF      do j = 1 , nlat !FF
!OLDFF        na = na+1 !FF
!OLDFF        lo(na) =slon + dij*i !FF
!OLDFF        la(na) = slat + dij*j !!FFF
!OLDFF      end do   !!FF
!OLDFF    end do   !!FF
!       if (allocated(la)) deallocate(la)
!       allocate(la(sum(gdomain%ni)*gdomain%nj))
!       if (allocated(lo)) deallocate(lo)
!       allocate(lo(sum(gdomain%ni)*gdomain%nj))
    ifin = gdomain%igstart(1)+gdomain%ni(1)
    jfin = gdomain%jgstart+gdomain%nj
!    print*,'ifin',ifin
!    print*,'jfin',jfin
!    print*,'gdomain%igstart(1)',gdomain%igstart(1)
!    print*,'gdomain%ni(1)',gdomain%ni(1)
!    print*,'gdomain%jgstart',gdomain%jgstart
!    print*,'gdomain%nj',gdomain%nj
    if (rdata == 10) then
      do i=gdomain%igstart(1),ifin
         do j=gdomain%jgstart,jfin
        na=na+1
        la(na)=eralat(i,j) ; lo(na)=eralon(i,j)
      enddo ; enddo
    end if
!    print*,'maxval la',maxval(la)
!    print*,'maxval lo',maxval(lo)
!    print*,'MAXVAL',maxval(raini)
    if (allocated(raini)) deallocate(raini)
    if (allocated(tmpm)) deallocate(tmpm)
    return
  end subroutine erareadncll

  subroutine erareadncll2(nlon,nlat,ora,giorno,mese,anno,ifl, &
             pi,la,lo,n,rdata)
    implicit none
    integer :: ora,giorno,mese,anno,ifl,n,na,nrec
    integer :: nlon,nlat,i,j,itime,rdata,ifin,jfin
    logical first
    data first /.true./
    save first
    character(len=60) :: dir
    character(len=132) :: cfl,oldcfl
    character(len=64) :: rainmis
    save rainmis,oldcfl
    integer :: len_trim,log,lun,code,pcode,ld,istatus,rec,ivarid
    integer,dimension(4) :: istart,icount
    save log,lun,pcode,istart,icount,ld,dir
    real, allocatable, dimension(:,:) :: tmpm,raini !Fabio-nc
    real, dimension(nlon) :: glon !Fabio-nc
    real, dimension(nlat) :: glat !Fabio-nc
    integer, parameter :: n2dim=120000
    real, dimension(n2dim) :: loclat,loclon
    real :: pfact
    integer :: xmese, iband
    logical :: bisest
    real :: mancanti,fill_value,tmp,offset,scalef
    save loclat,loclon,pfact,fill_value,mancanti,gdomain
    type(global_domain) :: gdomain
    integer :: itile , iti , itf
    real, dimension(:) :: pi,la,lo
    xmese=0
    if (mese.gt.1) xmese=julianday(monthlen(mese-1,anno), &
       mese-1,anno)
    if (nlon*nlat.gt.n2dim) then
       write(6,'(/,14x,a,a)')  'Too big Era domain, increase '// &
        'n2dim parameter inside erareadncll'
       call exit(1)
    end if
    if (first) then
      call mvgetiflags(70,log) ; if (log.le.0.or.log.ge.100) log=6
      if (rdata == 8) then
        dir="./museo/PREC/EIN15/"
        write(cfl,'(a,i4,a,i4,a)') 'erai_precip_', &
          anno,'0101_',anno,'1231_00_1.50.nc'
        call getlun(lun)
        write(log,'(18x,a,a,i4,a,i4,a)')'Opening Era file ERAIN15', &
          'erai_precip_',anno,'0101_',anno,'1231_00_1.50.nc'
        pfact=1000./3.
      else if (rdata == 10) then
        dir="./museo/PREC/EIN75/"
!        write(cfl,'(a,i4,a)') 'erainterim_075_',anno,'.nc'
        write(cfl,'(a,i4,a,i4,a)') 'erai_precip_',anno,'0101_',anno,'1231_00_0.75.nc'
        write(log,'(18x,a,a,i4,a)')'Opening Era file ERAIN75', &
        'erainterim_075_',anno,'.nc'
        pfact=1.
      end if
      istatus = nf90_open(trim(dir)//cfl,nf90_nowrite,lun)
      call chkncstatus(istatus,nf90_noerr,'open',dir,cfl)
      oldcfl=cfl
      istatus = nf90_inq_varid(lun,'lat',ivarid)
      call chkncstatus(istatus,nf90_noerr,'Missing latitude variable in file ',dir,cfl)
      istatus = nf90_get_var(lun,ivarid,glat)
      call chkncstatus(istatus,nf90_noerr,'Error reading lat var in file ',dir,cfl)
      istatus = nf90_inq_varid(lun,'lon',ivarid)
      call chkncstatus(istatus,nf90_noerr,'Missing longitude variable in file ',dir,cfl)
      istatus = nf90_get_var(lun,ivarid,glon)
      call chkncstatus(istatus,nf90_noerr,'Error reading lon var in file ',dir,cfl)
      !
      ! Find window to read
      !
      call get_window(glat,glon,lat,lon,iband,gdomain)

      allocate (tmpm(sum(gdomain%ni),gdomain%nj), raini(sum(gdomain%ni),gdomain%nj), &
                rainiold15(sum(gdomain%ni),gdomain%nj),rainiold75(sum(gdomain%ni),gdomain%nj))
      raini = 0

      istatus = nf90_inq_varid(lun,'rain',pcode)
      if (istatus.ne.nf90_noerr) then
          istatus = nf90_inq_varid(lun,'TP',pcode)
      end if
      call chkncstatus(istatus,nf90_noerr,'rain code',dir,cfl)
      istatus = nf90_inq_varid(lun,'time',itime)
      call chkncstatus(istatus,nf90_noerr,'time',dir,cfl)
      if (.not.allocated(raini)) allocate (raini(sum(gdomain%ni),gdomain%nj))
      raini = 0
      rec=(xmese*8+(giorno-1)*24+(ora))+1
      if (rec.eq.0) rec=1
      istart(3) = rec ; icount(3) = 1
      iti = 1
      do itile = 1 , gdomain%ntiles
        istart(1) = gdomain%igstart(itile)
        icount(1) = gdomain%ni(itile)
        ! Latitudes are reversed in original file
        istart(2) = gdomain%jgstart
        icount(2) = gdomain%nj
        itf = iti + gdomain%ni(itile) - 1
        istatus = nf90_get_var(lun,pcode,raini(iti:itf,:),istart,icount)
        call chkncstatus(istatus,nf90_noerr,'Error read var tp era',dir,cfl)
        iti = iti + gdomain%ni(itile)
      end do
!      call erashiftlon(raini,sum(gdomain%ni),gdomain%nj)
!      call erashiftlat(raini,sum(gdomain%ni),gdomain%nj)
      first=.false.
    end if
    if (rdata == 8) then
      write(cfl,'(a,i4,a,i4,a)') 'erai_precip_', &
          anno,'0101_',anno,'1231_00_1.50.nc'
    end if
    if (rdata == 10) then
      write(cfl,'(a,i4,a)') 'erainterim_075_',anno,'.nc'
    end if
    if (cfl(1:len_trim(cfl)).ne.oldcfl(1:len_trim(oldcfl))) then
       istatus = nf90_close(lun)
       call chkncstatus(istatus,nf90_noerr,'close',trim(dir),oldcfl)
       call getlun(lun)
       write(log,'(18x,a)') 'Opening ERA file '//cfl(1:len_trim(cfl))
       istatus = nf90_open(trim(dir)//cfl,nf90_nowrite,lun)
       call chkncstatus(istatus,nf90_noerr,'open',trim(dir),cfl)
       oldcfl=cfl
    endif
    nrec=(xmese*8+(giorno-1)*24+(ora))+1
    if (nrec.eq.0) nrec=1
    if (ifl.eq.1) then
       if (.not.allocated(raini)) allocate (raini(sum(gdomain%ni),gdomain%nj))
       raini = 0
       istart(3)=nrec ; icount(3)=1
       iti = 1
       do itile = 1 , gdomain%ntiles
         istart(1) = gdomain%igstart(itile)
         icount(1) = gdomain%ni(itile)
         istart(2) = gdomain%jgstart
         icount(2) = gdomain%nj
         itf = iti + gdomain%ni(itile) - 1

         istatus = nf90_get_var(lun,pcode,raini(iti:itf,:),istart,icount)
         call chkncstatus(istatus,nf90_noerr,'Error read var tp era',dir,cfl)
         iti = iti + gdomain%ni(itile)
       end do
!       call erashiftlat(raini,sum(gdomain%ni),gdomain%nj)

       do i=1,sum(gdomain%ni) ; do j=1,gdomain%nj
         if ( raini(i,j) < 0 ) raini(i,j) = 0      !Soluzione temporanea a valori negativi di precipitazione dovuti probabilmente ad arrotondamenti
       end do ; end do
       call regcmmat2vec(raini,sum(gdomain%ni),gdomain%nj,pfact,pi)
    else
       write (6,'(16x,a)') '  Bad flag passed to erareadnc.  &
               Exiting...'
       call exit(1)
    endif
    n=sum(gdomain%ni)*gdomain%nj
    na = 0
    if (rdata == 8) then
      do i=1,240 ; do j=1,121
        na=na+1
        la(na)=eralat1_5(i,j) ; lo(na)=eralon1_5(i,j)
      enddo ; enddo
    end if
    ifin = gdomain%igstart(1)+gdomain%ni(1)
    jfin = gdomain%jgstart+gdomain%nj
    if (rdata == 10) then
      do i=gdomain%igstart(1),ifin
         do j=gdomain%jgstart,jfin
        na=na+1
        la(na)=eralat(i,j) ; lo(na)=eralon(i,j)
      enddo ; enddo
    end if
    if (allocated(raini)) deallocate(raini)
    if (allocated(tmpm)) deallocate(tmpm)
    return
  end subroutine erareadncll2

  subroutine erareadncll3(nlon,nlat,ora,giorno,mese,anno,ifl, &
             pi,la,lo,n,rdata)
    implicit none
    integer :: ora,giorno,mese,anno,ifl,n,na,nrec
    integer :: nlon,nlat,i,j,itime,rdata,ifin,jfin
    logical first
    data first /.true./
    save first
    character(len=60) :: dir
    character(len=132) :: cfl,oldcfl
    character(len=64) :: rainmis
    character(len=8) :: adata
    save rainmis,oldcfl
    integer :: len_trim,log,lun,code,pcode,ld,istatus,rec,ivarid
    integer,dimension(4) :: istart,icount
    save log,lun,pcode,istart,icount,ld,dir
    real, allocatable, dimension(:,:) :: tmpm,raini !Fabio-nc
    real, dimension(nlon) :: glon !Fabio-nc
    real, dimension(nlat) :: glat !Fabio-nc
    integer, parameter :: n2dim=120000
    real, dimension(n2dim) :: loclat,loclon
    real :: pfact
    integer :: xmese, iband
    logical :: bisest
    real :: mancanti,fill_value,tmp,offset,scalef
    save loclat,loclon,pfact,fill_value,mancanti,gdomain
    type(global_domain) :: gdomain
    integer :: itile , iti , itf
    real, dimension(:) :: pi,la,lo
    xmese=0
    if (mese.gt.1) xmese=julianday(monthlen(mese-1,anno), &
       mese-1,anno)
!    if (nlon*nlat.gt.n2dim) then
!       write(6,'(/,14x,a,a)')  'Too big Era domain, increase '// &
!        'n2dim parameter inside erareadncll'
!       call exit(1)
!    end if
    write(adata,'(i8)') anno*10000+mese*100+giorno
    if (first) then
      call mvgetiflags(70,log) ; if (log.le.0.or.log.ge.100) log=6
      if (rdata == 8) then
        dir="./museo/PREC/EIN15/"
        write(cfl,'(a,i4,a,i4,a)') 'erai_precip_', &
          anno,'0101_',anno,'1231_00_1.50.nc'
        call getlun(lun)
        write(log,'(18x,a,a,i4,a,i4,a)')'Opening Era file ERAIN15', &
          'erai_precip_',anno,'0101_',anno,'1231_00_1.50.nc'
        pfact=1000./3.
      else if (rdata == 10) then
        dir="./museo/PREC/EIN75/"
!        write(cfl,'(a,i4,a)') 'erainterim_075_',anno,'.nc'
        write(cfl,'(a,i4,a,i4,a)') 'erai_precip_',anno,'0101_',anno,'1231_00_0.75.nc'
        write(log,'(18x,a,a,i4,a)')'Opening Era file ERAIN75', &
        'erainterim_075_',anno,'.nc'
        pfact=1000./3.
      else if (rdata == 14 .or. rdata == 16) then
        dir="./museo/PREC/ERA5/"
        write(cfl,'(a,a,a,a,a)') 'pr_',adata(1:4),'_',adata(5:6),'.nc'
        write(log,'(18x,a,a)')'Opening Era file ERA5', &
        cfl
        pfact=1000.
      end if
      istatus = nf90_open(trim(dir)//cfl,nf90_nowrite,lun)
      call chkncstatus(istatus,nf90_noerr,'open',dir,cfl)
      oldcfl=cfl
      istatus = nf90_inq_varid(lun,'latitude',ivarid)
      call chkncstatus(istatus,nf90_noerr,'Missing latitude variable in file ',dir,cfl)
      istatus = nf90_get_var(lun,ivarid,glat)
      call chkncstatus(istatus,nf90_noerr,'Error reading lat var in file ',dir,cfl)
      istatus = nf90_inq_varid(lun,'longitude',ivarid)
      call chkncstatus(istatus,nf90_noerr,'Missing longitude variable in file ',dir,cfl)
      istatus = nf90_get_var(lun,ivarid,glon)
      call chkncstatus(istatus,nf90_noerr,'Error reading lon var in file ',dir,cfl)
      !
      ! Find window to read
      !
!      call get_window(glat,glon,lat,lon,iband,gdomain)

      istatus = nf90_inq_varid(lun,'rain',pcode)
      if (istatus.ne.nf90_noerr) then
          istatus = nf90_inq_varid(lun,'TP',pcode)
      end if
      if (istatus.ne.nf90_noerr) then
          istatus = nf90_inq_varid(lun,'tp',pcode)
      end if
      call chkncstatus(istatus,nf90_noerr,'rain code',dir,cfl)
      istatus = nf90_inq_varid(lun,'time',itime)
      call chkncstatus(istatus,nf90_noerr,'time',dir,cfl)
      if (.not.allocated(raini)) allocate (raini(nlon,nlat))
      raini = 0
      if (.not.allocated(tmpm) .and. rdata /= 14) then
        allocate (tmpm(nlon,nlat))
        tmpm = 0
      end if
      if (rdata == 14) then
        rec=(giorno-1)*24+(ora)+1
      else
        rec=(xmese*8+(giorno-1)*8+(ora/3))+1
      end if
      if (rec.eq.0) rec=1
      istart(3) = rec ; icount(3) = 1
      istart(1) = 1
      icount(1) = nlon
      ! Latitudes are reversed in original file
      istart(2) = 1
      icount(2) = nlat
      istatus = nf90_get_var(lun,pcode,raini(:,:),istart,icount)
      call chkncstatus(istatus,nf90_noerr,'Error read var tp era',dir,cfl)
      istatus = nf90_get_att(lun,pcode, &
            'scale_factor',scalef)
      if ( istatus == nf90_noerr ) then
        istatus = nf90_get_att(lun,pcode,  &
             'add_offset',offset)
        call chkncstatus(istatus,nf90_noerr,'Error find att add_offset',dir,cfl)
      else
        offset = 0.0
        scalef = 1.0
      end if
      if (ora>=3 .and. rdata /= 14) then
        istart(3) = rec-1 ; icount(3) = 1
        istatus = nf90_get_var(lun,pcode,tmpm(:,:),istart,icount)
        call chkncstatus(istatus,nf90_noerr,'Error read var tp era',dir,cfl)
        tmpm = (tmpm*scalef)+offset
        raini = ((raini*scalef)+offset)-tmpm
      else
        raini = ((raini*scalef)+offset)
      end if
      call erashiftlon(raini,nlon,nlat)
      call erashiftlat(raini,nlon,nlat)
      first=.false.
    end if
    if (rdata == 8) then
      write(cfl,'(a,i4,a,i4,a)') 'erai_precip_', &
          anno,'0101_',anno,'1231_00_1.50.nc'
    end if
    if (rdata == 10) then
!      write(cfl,'(a,i4,a)') 'erainterim_075_',anno,'.nc'
      write(cfl,'(a,i4,a,i4,a)') 'erai_precip_',anno,'0101_',anno,'1231_00_0.75.nc'
    end if
    if (rdata == 14) then
      write(cfl,'(a,a,a,a,a)') 'pr_',adata(1:4),'_',adata(5:6),'.nc'
    end if
    if (cfl(1:len_trim(cfl)).ne.oldcfl(1:len_trim(oldcfl))) then
       istatus = nf90_close(lun)
       call chkncstatus(istatus,nf90_noerr,'close',trim(dir),oldcfl)
       call getlun(lun)
       write(log,'(18x,a)') 'Opening ERA file '//cfl(1:len_trim(cfl))
       istatus = nf90_open(trim(dir)//cfl,nf90_nowrite,lun)
       call chkncstatus(istatus,nf90_noerr,'open',trim(dir),cfl)
       oldcfl=cfl
    endif
    if (rdata == 14) then
      nrec=(giorno-1)*24+(ora)+1
    else
      nrec=(xmese*8+(giorno-1)*8+(ora/3))+1
    end if
    if (nrec.eq.0) nrec=1
    if (ifl.eq.1) then
       if (.not.allocated(raini)) allocate (raini(nlon,nlat))
       raini = 0
       if (.not.allocated(tmpm) .and. rdata /= 14) then
         allocate (tmpm(nlon,nlat))
         tmpm = 0
       end if
       istart(3)=nrec ; icount(3)=1
       istart(1) = 1
       icount(1) = nlon
       istart(2) = 1
       icount(2) = nlat

       istatus = nf90_get_var(lun,pcode,raini(:,:),istart,icount)
       call chkncstatus(istatus,nf90_noerr,'Error read var tp era',dir,cfl)
       istatus = nf90_get_att(lun,pcode, &
             'scale_factor',scalef)
       if ( istatus == nf90_noerr ) then
         istatus = nf90_get_att(lun,pcode,  &
              'add_offset',offset)
         call chkncstatus(istatus,nf90_noerr,'Error find att add_offset',dir,cfl)
       else
         offset = 0.0
         scalef = 1.0
       end if
       if (ora>=3 .and. rdata /= 14 ) then
         istart(3) = nrec-1 ; icount(3) = 1
         istatus = nf90_get_var(lun,pcode,tmpm(:,:),istart,icount)
         call chkncstatus(istatus,nf90_noerr,'Error read var tp era',dir,cfl)
         tmpm = (tmpm*scalef)+offset
         raini = ((raini*scalef)+offset)-tmpm
       else
         raini = ((raini*scalef)+offset)
       end if
       call erashiftlon(raini,nlon,nlat)
       call erashiftlat(raini,nlon,nlat)

       do i=1,nlon ; do j=1,nlat
         if ( raini(i,j) < 0 ) raini(i,j) = 0      !Soluzione temporanea a valori negativi di precipitazione dovuti probabilmente ad arrotondamenti
       end do ; end do
       call regcmmat2vec(raini,nlon,nlat,pfact,pi)
    else
       write (6,'(16x,a)') '  Bad flag passed to erareadnc.  &
               Exiting...'
       call exit(1)
    endif
    n=nlon*nlat
    na = 0
    if (rdata == 8) then
      do i=1,240 ; do j=1,121
        na=na+1
        la(na)=eralat1_5(i,j) ; lo(na)=eralon1_5(i,j)
      enddo ; enddo
    end if
    if (rdata == 10) then
      do i=1,nlon
         do j=1,nlat
        na=na+1
        la(na)=eralatold(i,j) ; lo(na)=eralonold(i,j)
      enddo ; enddo
    end if
    if (rdata == 14) then
      do i=1,nlon
         do j=1,nlat
        na=na+1
        la(na)=eralat5(i,j) ; lo(na)=eralon5(i,j)
      enddo ; enddo
    end if
    if (rdata == 16) then
      do i=1,nlon
         do j=1,nlat
        na=na+1
        la(na)=glat(j) ; lo(na)=glon(i)
      enddo ; enddo
    end if
    if (allocated(raini)) deallocate(raini)
    if (allocated(tmpm)) deallocate(tmpm)
    return
  end subroutine erareadncll3

  subroutine hourlyrainnc(ora,giorno,mese,anno,pi,la,lo,n)
    implicit none
    integer :: ora,giorno,mese,anno
    real, dimension(:) :: pi,la,lo
    character(len=80) :: dir
    character(len=132) :: cfl,oldcfl
    logical first
    data first /.true./
    save oldcfl,first
    save log,lun,pcode,istart,icount,dir
    integer :: len_trim,log,lun,ivarid,istatus,rec,i,j,itime,n
    integer :: pcode,nrec
    integer , dimension(2) :: istart , icount
    real, dimension(3712) :: glon !Fabio-nc
    real, dimension(3712) :: glat !Fabio-nc
    real, dimension(3712) :: gpi !Fabio-nc
    call mvgetiflags(70,log) ; if (log.le.0.or.log.ge.100) log=6
    dir = "/home/clima-archive/afantini/flood/museodat2nc/clean_dewrain/cleaned/"
    if (mese .le. 9) then
      write(cfl,'(i4,a,i1,a)') anno,'0',mese,'.nc'
    else
      write(cfl,'(i4,i2,a)') anno,mese,'.nc'
    end if
    if (first) then
      call getlun(lun)
      write(log,'(18x,a,a)') 'Opening Hourly Precipitation file',cfl
      istatus = nf90_open(trim(dir)//cfl,nf90_nowrite,lun)
      call chkncstatus(istatus,nf90_noerr,'open',dir,cfl)
      oldcfl=cfl
      istatus = nf90_inq_varid(lun,'lat',ivarid)
      call chkncstatus(istatus,nf90_noerr,'Missing latitude variable in file ',dir,cfl)
      istatus = nf90_get_var(lun,ivarid,glat)
      call chkncstatus(istatus,nf90_noerr,'Error reading lat var in file ',dir,cfl)
      istatus = nf90_inq_varid(lun,'lon',ivarid)
      call chkncstatus(istatus,nf90_noerr,'Missing longitude variable in file ',dir,cfl)
      istatus = nf90_get_var(lun,ivarid,glon)
      call chkncstatus(istatus,nf90_noerr,'Error reading lon var in file ',dir,cfl)
      istatus = nf90_inq_varid(lun,'pr',pcode)
      call chkncstatus(istatus,nf90_noerr,'rain code',dir,cfl)
      istatus = nf90_inq_varid(lun,'time',itime)
      call chkncstatus(istatus,nf90_noerr,'time',dir,cfl)
      first=.false.
    end if
    if (cfl(1:len_trim(cfl)).ne.oldcfl(1:len_trim(oldcfl))) then
      istatus = nf90_close(lun)
      call chkncstatus(istatus,nf90_noerr,'close',trim(dir),oldcfl)
      call getlun(lun)
      write(log,'(18x,a,a)') 'Opening Hourly Precipitation file',cfl
      istatus = nf90_open(trim(dir)//cfl,nf90_nowrite,lun)
      call chkncstatus(istatus,nf90_noerr,'open',trim(dir),cfl)
      oldcfl=cfl
      istatus = nf90_inq_varid(lun,'lat',ivarid)
      call chkncstatus(istatus,nf90_noerr,'Missing latitude variable in file ',dir,cfl)
      istatus = nf90_get_var(lun,ivarid,glat)
      call chkncstatus(istatus,nf90_noerr,'Error reading lat var in file ',dir,cfl)
      istatus = nf90_inq_varid(lun,'lon',ivarid)
      call chkncstatus(istatus,nf90_noerr,'Missing longitude variable in file ',dir,cfl)
      istatus = nf90_get_var(lun,ivarid,glon)
      call chkncstatus(istatus,nf90_noerr,'Error reading lon var in file ',dir,cfl)
      istatus = nf90_inq_varid(lun,'pr',pcode)
      call chkncstatus(istatus,nf90_noerr,'rain code',dir,cfl)
      istatus = nf90_inq_varid(lun,'time',itime)
      call chkncstatus(istatus,nf90_noerr,'time',dir,cfl)
    endif
    rec=((giorno-1)*24+(ora))*4+1
    if (nrec.eq.0) nrec=1
    if (rec.eq.0) rec=1
    istart(2) = rec ; icount(2) = 1
    istart(1) = 1
    icount(1) = 3712
    ! Latitudes are reversed in original file
    istatus = nf90_get_var(lun,pcode,gpi,istart,icount)
    call chkncstatus(istatus,nf90_noerr,'Error read var pr hourlyrainnc',dir,cfl)
    do i=1,3712
      la(i)=glat(i) ; lo(i)=glon(i)
      pi(i)=gpi(i)
    enddo
    n = 3712

  end subroutine hourlyrainnc

!inmuseo  subroutine chkncstatus(istatus,nf90_noerr,field,dir,cfl)
!inmuseo    implicit none
!inmuseo    integer :: istatus,nf90_noerr
!inmuseo    character(len=*) :: field,dir,cfl
!inmuseo    character(len=80) :: nf90_strerror
!inmuseo    if (istatus.eq.nf90_noerr) return
!inmuseo    if (len_trim(field).eq.0.or.field(1:len_trim(field)).eq.'open') then
!inmuseo       write(6,'(/,10x,a)') 'NetCDF Library error while opening file'
!inmuseo       write(6,'(10x,a)') dir(1:len_trim(dir))//cfl(1:len_trim(cfl))
!inmuseo    else if (field(1:len_trim(field)).eq.'close') then
!inmuseo       write(6,'(/,10x,a)') 'NetCDF Library error while closing file'
!inmuseo       write(6,'(10x,a)') dir(1:len_trim(dir))//cfl(1:len_trim(cfl))
!inmuseo    else if (field(1:14).eq.'reading record') then
!inmuseo       write(6,'(/,10x,a)') 'NetCDF Library error while '// &
!inmuseo           field(1:len_trim(field))
!inmuseo       write(6,'(10x,a)') 'from'//dir(1:len_trim(dir))//cfl(1:len_trim(cfl))
!inmuseo    else
!inmuseo       write(6,'(/,10x,a)') 'NetCDF Library error while reading '// &
!inmuseo           field(1:len_trim(field))//' field from'
!inmuseo       write(6,'(10x,a)') dir(1:len_trim(dir))//cfl(1:len_trim(cfl))
!inmuseo    endif
!inmuseo    write(6,'(10x,a)') nf90_strerror(istatus)
!inmuseo    write(6,'(10x,a)') 'Exiting with status = 1'
!inmuseo    call exit(1)
!inmuseo    return
!inmuseo  end subroutine chkncstatus

  subroutine chymreadirec(lun,field,vsource,mat,idate)
    implicit none
    integer :: lun,idate,att,isize
    integer, dimension(nlon,nlat) :: mat
    character(len=*) :: field,vsource
    character(len=2) :: ctype
    character(len=3) :: chkfld
    integer i,j
    if (chymcrec.lt.0) call chymerror(16,-9999,-9999.0,'chymreadrec')
    att=0 ; call setmvlibinti(1,-9999)
    do while (att.le.1)
       do while (att.le.1.or.att.ge.1)
          read (lun,end=100) chkfld,ctype,isize,idate,vsource
          chymcrec=chymcrec+1
          if (chkfld(1:3).eq.field(1:3).and.ctype(1:1).eq.'i') then
             read(lun) mat
             return
          else
             read(lun)
          endif
       enddo
100    continue
       rewind(lun) ; read (lun) ; att=att+1 ; chymcrec=0 ; hourstep=0
       call setmvlibinti(1,att)
    enddo
    call chymerror(10,lun,0.0,field)
  end subroutine chymreadirec

  subroutine chymreadrrec(lun,field,vsource,mat,idate)
    implicit none
    integer :: lun,idate,att,isize
    real, dimension(nlon,nlat) :: mat
    character(len=*) :: field, vsource
    character(len=2) :: ctype
    character(len=3) :: chkfld
    integer i,j
    if (chymcrec.lt.0) call chymerror(16,-9999,-9999.0,'chymreadrec')
    att=0 ; call setmvlibinti(1,-9999)
    do while (att.le.1)
       do while (att.le.1.or.att.ge.1)
          read (lun,end=100,err=100) chkfld,ctype,isize,idate,vsource
          chymcrec=chymcrec+1
          if (chkfld(1:3).eq.field(1:3).and.ctype(1:1).eq.'r') then
             read(lun,end=100,err=100) mat
             return
          else
             read(lun,end=100,err=100)
          endif
       enddo
100    continue
       rewind(lun) ; read (lun) ; att=att+1 ; chymcrec=0
       call setmvlibinti(1,att)
    enddo
    call chymerror(10,lun,0.0,field)
  end subroutine chymreadrrec

  integer function search_indx(vname) result(ivar)
    character(len=*) , intent(in) :: vname
    integer :: iv
    ivar = -1
    do iv = 1 , maxf
      if ( var(iv)%vname == vname ) then
        ivar = iv
        exit
      end if
    end do
    if ( ivar < 0 ) then
      write(6,*) 'Variable ', trim(vname), &
      ' not found in defined vars'
      stop
    end if
  end function search_indx

  real function persiannlon0_25(i,j)
    implicit none
    integer i,j
    persiannlon0_25=-180.125+0.25*i
    return
  end function persiannlon0_25

  real function persiannlat0_25(i,j)
    implicit none
    integer i,j
    persiannlat0_25=60.125-0.25*j
    return
  end function persiannlat0_25

  real function trimmlon(i,j)
    implicit none
    integer i,j
    trimmlon=-180.125+0.25*i
!    trimmlon=0.125+0.25*i
    return
  end function trimmlon

  real function trimmlat(i,j)
    implicit none
    integer i,j
!    trimmlat=49.325-0.25*j
    trimmlat=-49.875+0.25*j
    return
  end function trimmlat

  subroutine checknetcdfresult(token)
    implicit none
    character(len=*) , intent(in) :: token
    if ( iretval /= nf90_noerr ) then
      write(6,*) 'Error in netcdf file in ',trim(token)
      write(6,*) 'NetCDF Error :',nf90_strerror(iretval)
      stop
    end if
  end subroutine checknetcdfresult

  ! Assumptions :
  !  GLAT is regular ordered array of latitudes in a GLOBAL grid
  !  GLON is regular ordered array of longitudes in a GLOBAL grid
  !  XLAT is local grid latitudes going South to North
  !  XLON is local grid longitudes  going West to East
  ! ALL LONGITUDES ARE in the range -180.0 <-> 180.0
  subroutine get_window(glat,glon,xlat,xlon,i_band,domain)
    implicit none
    real , dimension(:) , intent(in) :: glat
    real , dimension(:) , intent(in) :: glon
    real , dimension(:,:) , intent(in) :: xlat
    real , dimension(:,:) , intent(in) :: xlon
    integer , intent(in) :: i_band
    type(global_domain) , intent(out) :: domain

    real :: dlat , dlon
    real , allocatable , dimension(:,:) :: xlon360
    real :: maxlat
    real :: minlat
    real :: maxlon
    real :: minlon
    integer :: gi , gj , xi , xj , l1 , l2 , i , j , itmp

    xi = size(xlon,1)
    xj = size(xlat,2)
    maxlat = maxval(xlat)
    minlat = minval(xlat)
    maxlon = maxval(xlon)
    minlon = minval(xlon)

    gi = size(glon)
    gj = size(glat)
    dlat = abs(glat(2) - glat(1))
    dlon = abs(glon(2) - glon(1))
    domain%global_ni = gi
    domain%global_nj = gj
    if ( i_band == 1 ) then
      domain%ntiles = 1
      domain%igstart(1) = 1
      domain%igstop(1) = gi
      domain%igstart(2) = 0
      domain%igstop(2) = 0
    else if ( glon(gi) > 350.0 ) then
      ! Input data is     0 : 360 , xlon is -180 : 180
      allocate(xlon360(xi,xj))
      xlon360 = xlon
      where ( xlon < 0.0 )
        xlon360 = 360.0 + xlon
      end where
      if ( minval(xlon360(1,:)) < maxval(xlon360(xi,:)) ) then
        domain%ntiles = 1
        minlon = minval(xlon360)
        maxlon = maxval(xlon360)
        domain%igstart(1) = int((minlon-glon(1))/dlon) - 2
        domain%igstop(1) = int((maxlon-glon(1))/dlon) + 3
        domain%igstart(2) = 0
        domain%igstop(2) = 0
      else
        ! Cross Greenwich line
        minlon = minval(xlon360(1,:))
        maxlon = maxval(xlon(xi,:))
        domain%ntiles = 2
        domain%igstart(1) = int((minlon-glon(1))/dlon) - 2
        domain%igstop(1) = gi
        domain%igstart(2) = 1
        domain%igstop(2) = int((maxlon-glon(1))/dlon) + 3
      end if
      deallocate(xlon360)
    else
      ! Input Data is -180 : 180 , xlon is -180 : 180
      if ( xlon(1,xj)   <= xlon(xi,xj)   .and. &
           xlon(1,xj/2) <= xlon(xi,xj/2) .and. &
           xlon(1,1)    <= xlon(xi,1) ) then
        ! it is not crossing timeline
        domain%ntiles = 1
        domain%igstart(1) = int((minlon-glon(1))/dlon) - 2
        domain%igstop(1) = int((maxlon-glon(1))/dlon) + 3
        domain%igstart(2) = 0
        domain%igstop(2) = 0
      else
        domain%ntiles = 2
        minlon = 180.0
        do j = 1 , xj
          if ( xlon(1,j) > 0.0 ) minlon = min(minlon,xlon(1,j))
        end do
        maxlon = -180.0
        do j = 1 , xj
          if ( xlon(xi,j) < 0.0 ) maxlon = max(maxlon,xlon(xi,j))
        end do
        domain%igstart(1) = int((minlon-glon(1))/dlon) - 2
        domain%igstop(1) = gi
        domain%igstart(2) = 1
        domain%igstop(2) = int((maxlon-glon(1))/dlon) + 3
      end if
    end if
    if ( has_north_pole(xlat,xi/2) ) then
      ! North pole inside
      if ( glat(1) < glat(gj) ) then
        l1 = int((minval(xlat(:,1))-glat(1))/dlat) - 2
        l2 = int((minval(xlat(:,xj))-glat(1))/dlat) - 2
        domain%jgstart = min(l1,l2)
        domain%jgstop = gj
      else
        l1 = int((maxval(glat(1)-xlat(:,1)))/dlat) + 3
        l2 = int((maxval(glat(1)-xlat(:,xj)))/dlat) + 3
        domain%jgstart = 1
        domain%jgstop = max(l1,l2)
      end if
      domain%ntiles = 1
      domain%igstart(1) = 1
      domain%igstop(1) = gi
      domain%igstart(2) = 0
      domain%igstop(2) = 0
    else if ( has_south_pole(xlat,xi/2) ) then
      ! South Pole inside
      if ( glat(1) < glat(gj) ) then
        l1 = int((maxval(xlat(:,1))-glat(1))/dlat) + 3
        l2 = int((maxval(xlat(:,xj))-glat(1))/dlat) + 3
        domain%jgstart = 1
        domain%jgstop = max(l1,l2)
      else
        l1 = int((maxval(glat(1)-xlat(:,1)))/dlat) - 2
        l2 = int((maxval(glat(1)-xlat(:,xj)))/dlat) - 2
        domain%jgstart = min(l1,l2)
        domain%jgstop = gj
      end if
      domain%ntiles = 1
      domain%igstart(1) = 1
      domain%igstop(1) = gi
      domain%igstart(2) = 0
      domain%igstop(2) = 0
    else
      if ( glat(1) < glat(gj) ) then
        domain%jgstart = int((minlat-glat(1))/dlat) - 2
        domain%jgstop  = int((maxlat-glat(1))/dlat) + 3
      else
        domain%jgstart = int((minlat-glat(gj))/dlat) - 2
        domain%jgstop  = int((maxlat-glat(gj))/dlat) + 3
        domain%nj =  domain%jgstop - domain%jgstart + 1
        domain%jgstart = (gj+1)-(domain%jgstart+domain%nj-1)
        domain%jgstop = domain%jgstart + domain%nj - 1
      end if
    end if

    if ( domain%igstart(1) < 1 .and. domain%ntiles == 1 ) then
      domain%ntiles = 2
      itmp = domain%igstop(1)
      domain%igstart(1) = gi + domain%igstart(1) - 1
      domain%igstop(1) = gi
      domain%igstart(2) = 1
      domain%igstop(2) = itmp
    end if
    if ( domain%igstop(1) > gi .and. domain%ntiles == 1 ) then
      domain%ntiles = 2
      itmp = domain%igstop(1) - gi
      domain%igstart(1) = domain%igstart(1)
      domain%igstop(1) = gi
      domain%igstart(2) = 1
      domain%igstop(2) = itmp
    end if
    domain%ni = 0
    do i = 1 , domain%ntiles
      domain%ni(i) =  domain%igstop(i) - domain%igstart(i) + 1
    end do
    domain%jgstart = min(max(1,domain%jgstart),gj)
    domain%jgstop = max(min(domain%jgstop,gj),1)
    domain%nj =  domain%jgstop - domain%jgstart + 1

    contains

      logical function has_north_pole(l,i)
        real , intent(in) , dimension(:,:) :: l
        integer , intent(in) :: i
        integer :: j
        has_north_pole = .false.
        if ( all(l(i,:) < 0.0) ) return
        do j = 2 , size(l,2)
          if ( l(i,j) < l(i,j-1) ) then
            has_north_pole = .true.
            exit
          end if
        end do
      end function has_north_pole

      logical function has_south_pole(l,i)
        real , intent(in) , dimension(:,:) :: l
        integer , intent(in) :: i
        integer :: j
        has_south_pole = .false.
        if ( all(l(i,:) > 0.0) ) return
        do j = 2 , size(l,2)
          if ( l(i,j) < l(i,j-1) ) then
            has_south_pole = .true.
            exit
          end if
        end do
      end function has_south_pole

  end subroutine get_window

end module mod_ncio
