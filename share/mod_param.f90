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
module mod_param

     use mod_statparams
     use mod_runparams
     use mod_internal
     use mod_phys
     use mod_strman
     use mod_time
     use mod_mssg
     use mod_libmv
     use mod_mpimess
     use mod_varandtypes

     contains

   subroutine setparam
     implicit none
     integer :: i1,j1
     if (myid==0) then
       call read_namelist
!       namelist /chymconfig/ nlon , nlat , slon , slat , dij , chym_radius , &
!               nsli , nsave , demf , angiocycle , chym_sdate , chym_steps ,  &
!               chym_tempfl , chym_river , chym_savet , chym_savep ,          &
!               chym_restart , chym_grads , chym_modis , chym_tplot ,         &
!               chym_iplot , chym_zoom , chym_verbose , chym_mfile ,          &
!               chym_ofile , chym_rfile , chym_pfile , chym_tfile ,           &
!               chym_dsource , chym_ifile1 , chym_ifile2 , chem_symtype ,     &
!               chym_savefld , chym_netcdf , angionp , ncyc1 , ncyc2 , ncyc3, &
!               integrflag , cpar1 , cpar2 , cpar3 , cpar4 , cpar5 , cpar6 ,  &
!               cpar7 , cpar8 , cpar9 , cpar10 , infiltr , infi_lago ,        &
!               infi_fiume , infi_ice , chym_manning
!
!       call getarg(1, namelistfile)
!       call getlun(lun)
!       open(lun,file=namelistfile,status='old',action='read',iostat=iretval)
!       if ( iretval /= 0 ) then
!         write(0,*) 'Error opening namelist file'//trim(namelistfile)
!         stop
!       end if
!       write(6,'(7x,a)') 'Reading model namelist'
!       rewind(lun)
!       read(lun,nml=chymconfig,iostat=iretval)
!       if ( iretval /= 0 ) then
!         write(0,*) 'Error reading namelist'
!         stop
!       end if
!       if (demf == 10) dij = 0.0083333
       n2df = nlon*nlat
       n4df = 2*n2df
       isdate = chym_sdate
       ienddate = chym_edate
     end if
     call mpi_bcast(nlon,1,MPI_integer, 0,mycomm,mpierr)
     call mpi_bcast(nlat,1,MPI_integer, 0,mycomm,mpierr)
     call mpi_bcast(n2df,1,MPI_integer, 0,mycomm,mpierr)
     call mpi_bcast(n4df,1,MPI_integer, 0,mycomm,mpierr)
!     call mpi_bcast(slon,1,MPI_integer, 0,comm,ierr)
!     call mpi_bcast(slat,1,MPI_integer, 0,comm,ierr)
!     call mpi_bcast(dij,1,MPI_real, 0,comm,ierr)
!     call mpi_bcast(chym_radius,1,MPI_integer, 0,comm,ierr)
!     call mpi_bcast(nsli,1,MPI_integer, 0,comm,ierr)
!     call mpi_bcast(nsave,1,MPI_integer, 0,comm,ierr)
!     call mpi_bcast(demf,1,MPI_integer, 0,comm,ierr)
!     call mpi_bcast(angiocycle,1,MPI_integer, 0,comm,ierr)
!     call mpi_bcast(chym_ofile,len_trim(chym_ofile),MPI_char, 0,comm,ierr)
!     print*,chym_ofile
!     call mpi_barrier(comm,ierr)
!     call exit(0)

     ! Get work space from system

     allocate(accl(nlon,nlat) , alfa(nlon,nlat) , area(nlon,nlat) ,       &
              bwet(nlon,nlat) , ddeepw(nlon,nlat) , deepw(nlon,nlat) ,    &
              dem(nlon,nlat) , dh2o(nlon,nlat) , drai(nlon,nlat) ,        &
              dx(nlon,nlat) , evap(nlon,nlat) , gh2o(nlon,nlat) ,         &
              h2o(nlon,nlat) , lat(nlon,nlat) , lon(nlon,nlat) ,          &
              modis(nlon,nlat) , port(nlon,nlat) , rain(nlon,nlat) ,      &
              rsrm(nlon,nlat) , runt(nlon,nlat) , snow(nlon,nlat) ,       &
              temp(nlon,nlat) , wk(nlon,nlat) , wrk1(nlon,nlat) ,         &
              wrk2(nlon,nlat) , wrk4(n2df) , ca(nlon,nlat) ,              &
              fmap(nlon,nlat) , luse(nlon,nlat) , noflow(nlon,nlat) ,     &
              iwk(nlon,nlat) , caplot(nlon,nlat) , lca(nlon,nlat) ,       &
              canc(nlon,nlat) , lavonc(nlon,nlat) , allrivb(nlon,nlat) ,  &
              avgrain(nlon,nlat),avgh2o(nlon,nlat),                       &
              avgrsrm(nlon,nlat), avgport(nlon,nlat),avgbwet(nlon,nlat),  &
              avggh2o(nlon,nlat), avgevap(nlon,nlat), avgsnow(nlon,nlat) ,&
              avgarai(nlon,nlat), avgtemp(nlon,nlat), avgdeepw(nlon,nlat),&
              avgddeepw(nlon,nlat),                                       & !rainiold15(nlon,nlat),rainiold75(nlon,nlat),
              seqx(nlon+nlat), seqy(nlon+nlat) ,   &
              iriv(n4df) , jriv(n4df) , seqi(n4df) , seqj(n4df))
     i1 = nlon
     j1 = (nlat/nproc)+4

!     if (allocated(portsub)) deallocate(portsub)
!     allocate(portsub(jde1gb:jde2gb,ide1gb:ide2gb))
!     if (allocated(wkm1sub)) deallocate(wkm1sub)
!     allocate(wkm1sub(jde1gb:jde2gb,ide1gb:ide2gb))
!     if (allocated(h2osub)) deallocate(h2osub)
!     allocate(h2osub(jde1gb:jde2gb,ide1gb:ide2gb))
!     if (allocated(bwetsub)) deallocate(bwetsub)
!     allocate(bwetsub(jde1gb:jde2gb,ide1gb:ide2gb))
!     portsub = 0.
!     wkm1sub = 0.
!     h2osub = 0.
!     bwetsub = 0.
     call MPI_BARRIER(mycomm,mpierr)

   end subroutine setparam

   subroutine acquirescriptpar
    implicit none
    logical :: chktemfl
    integer :: i,ih,ii,j,n,nf,nsts,ndyn,k
    integer , dimension(8) :: tvals
    integer :: rainflag = 0
    integer :: iunit = -1
    character(len=80) :: sstr1,nameluse
    character(len=8) :: chcentur
    integer :: test1
    if (myid == 0) then
    write (6,'(/12x,a)') 'Acquiring script parameters.'
    do i = 1 , nlon
      do j = 1 , nlat
        lon(i,j) = slon + (i-1)*dij
        lat(i,j) = slat + (j-1)*dij
      end do
    end do
    elon = lon(nlon,nlat)
    elat = lat(nlon,nlat)
    avr = 0.5*(distance(lat(1,1),lon(1,1),lat(2,1),lon(2,1)) + &
          distance(lat(1,1),lon(1,1),lat(1,2),lon(1,2)))
    write (sstr1,'(i10,a,f6.2,a,f6.2)') &
      nlon , ' longitudes in the range: ' , slon , '-' , elon
    call noinspace(sstr1)
    write (6,'(15x,a)') trim(sstr1)
    write (sstr1,'(i10,a,f6.2,a,f6.2)') &
       nlat , ' latitudes in the range: ' , slat , '-' , elat
    call noinspace(sstr1)
    write (6,'(15x,a)') trim(sstr1)
    write (6,'(15x,a,f6.1,a,f6.1,a)') 'Geographical domain is ' ,        &
                                      0.001*distance(lat(1,1),lon(1,1),  &
                                      lat(nlon,1),lon(nlon,1)) , ' x' ,  &
                                      0.001*distance(lat(1,1),lon(1,1),  &
                                      lat(1,nlat),lon(1,nlat)) , ' Km.'
    write (6,'(15x,a,f6.1,a)') 'Approximate Resolution: ' , avr , ' meters.'
    mm5domain = 1
    if ( mm5d.lt.1.or.mm5d.lt.3 ) then
      if ( slon>=12.51742 .and. slon<=14.95717 .and. elon>=12.51742 .and.  &
           elon<=14.95717 .and. slat>=41.18784 .and. slat<=42.91781 .and.  &
           elat>=41.18784 .and. elat<=42.91781 ) then
        mm5domain = 3
      else if ( slon>=10.21659 .and. slon<=17.82529 .and. elon>=10.21659 .and. &
                elon<=17.82529 .and. slat>=39.65918 .and. slat<=44.06796 .and. &
                elat>=39.65918 .and. elat<=44.06796 ) then
        mm5domain = 2
      else
        mm5domain = 1
      end if
    end if
!    call mvsetflags('Data style',2.0)
!    mm5file = chym_ifile1
!    if ( ilat2<1 ) ilat2 = 1
!    if ( ilat2>nlat ) ilat2 = nlat
!    izoom = chym_zoom
!    tplot = chym_tplot
!    iplot = chym_iplot
    if ( izoom<=0 .or. ilon2<=ilon1 .or. ilat2<=ilat1 ) then
      plon1 = slon
      plon2 = elon
      plat1 = slat
      plat2 = elat
      izoom = 0
    else
      plon1 = lon(ilon1,ilat1)
      plon2 = lon(ilon2,ilat2)
      plat1 = lat(ilon1,ilat1)
      plat2 = lat(ilon2,ilat2)
    end if
!    rainflag = chym_savep
    schym(1)=chym_dsource ; schym(2)=chym_ifile1     ! needed by definerainsources
    call definerainsources
    mchym(14) = 1
    call date_and_time(values=tvals)
    ora = tvals(5)
    giorno = tvals(3)
    mese = tvals(2)
    anno = tvals(1)
    isdate = chym_sdate
    ienddate = chym_edate
!    ii = increasemm5index(isdate)
    ii = increasetime(isdate)
    if ( ii==-9999 ) then
      isdate = mm5index(ora,giorno,mese,anno)
      do i = 1 , 72 + ora
!        isdate = increasemm5index(isdate)
        isdate = increasetime(isdate)
      end do
      write (6,'(15x,a,i10)') 'Assuming operational run starting from ' , isdate
    end if
    mchym(15) = 1
    chktemfl = .true.
    if ( mchym(15)>=1 .and. mchym(15)<=5 ) chktemfl = .false.
    if ( mchym(15)>=21 .and. mchym(15)<=29 ) chktemfl = .false.     ! ACQWA codes
    if ( chktemfl ) call chymerror(22,mchym(15),rchym(1),schym(1))
    nsts = 0
    ii = isdate
    write (6,'(15x,a,i10)') 'Assuming operational run starting from ' ,isdate
!    iend = chym_sdate
    nsli = 0
    do while (ii <= ienddate)
!!      ii = increasemm5index(ii)
      ii = increasetime(ii)
!      if ( mod(ih,1)==0 ) then
!        iend = ii
        nsli = nsli + 1
!        nsts = nsts + 1
!      end if
    end do
    write (6,'(15x,a,i10)') 'Assuming operational run ending at ' , ienddate

    avr = 0.5*(distance(lat(1,1),lon(1,1),lat(2,1),lon(2,1))                     &
          +distance(lat(1,1),lon(1,1),lat(1,2),lon(1,2)))
    chunksizes(1) = nlon/num_chunky
    chunksizes(2) = nlat/num_chunkx
    chunksizes(3) = 1
    call caparameters(avr)

    mchym = 0
    rchym = 0.0
    schym = ' '
    mchym(01) = 1              ! who write the output
    mchym(02) = nlon           ! number of longitudes
    mchym(03) = nlat           ! number of latitudes
    mchym(04) = isdate         ! code of start date
    mchym(05) = ienddate           ! Integer code for end date
    mchym(06) = nsts           ! # of time slices
    mchym(07) = 9              ! # of static fields
    mchym(08) = 0              ! # of dynamical fields - set below
    mchym(09) = chym_restart   ! Restart flag
    mchym(10) = mare           ! landuse code for sea
    mchym(11) = lago           ! landuse code for internal water
    mchym(12) = mm5domain      ! mm5 domain used for museo data
    mchym(13) = demf           ! DEM source code
    mchym(14) = nsave          ! Every how many steps data are saved
    mchym(15) = chym_tempfl    ! Temperature source flag
    mchym(16) = nsli           ! # of hourly step for this run
    mchym(17) = rainflag       ! Calcluate(0), write(1) or read(2) rainfall
    mchym(18) = chym_river     ! Selected river
    mchym(19) = 50             ! # of mchym component in the current version
    mchym(20) = chym_grads     ! If > 0 Out for grads utility is also created
    mchym(21) = nsts           ! # of time slices (same as 06 if restart=0,1,2)
    mchym(22) = chym_steps     ! Number of integration steps per hour
    mchym(23) = chym_savet     ! Calcluate(0), write(1) or read(2) temperature
    mchym(24) = chym_modis     ! Flags to establish use of ModIS data
    mchym(25) = chym_netcdf    ! If > 0 Out in netCDF format is also created
    mchym(26) = rsave          ! Every how many steps data are saved
    mchym(29) = packnetcdf     ! If > 0 scale_factor and add_offset added to netCDF

    if (mchym(29) > 0) call compute_scale_and_offset(0., 60000., 16)

    rchym(1) = slon
    rchym(2) = slat
    rchym(3) = lon(nlon,nlat) ! last longitude
    rchym(4) = lat(nlon,nlat) ! last latitude
    rchym(5) = dij            ! lat-lon resolution
    rchym(6) = avr
    rchym(7) = chym_version   ! CHyM version
    rchym(8) = dij            ! lon resolution
    rchym(9) = dij            ! lat resolution
    rchym(10) = chym_radius   ! Radius of influence for sparse data interpol.
    rchym(11) = chym_regcm_rad! Radius of influence for regcm data interpol.

!    if (mchym(13) == 10) rchym(5) = 0.0083333

    schym(1) = sdata         ! Rain data sources
    schym(2) = chym_ifile1   ! Chym or MM5 file
    write(schym(3),'(i0.10)') isdate
    schym(4) = get_nowdate() ! date of production of file
    call hostnm(schym(5))    ! hostname which produced file
    schym(6) = chem_symtype  ! title of experiment
    schym(7) = chym_savefld  ! Fields to be saved
!    call datafrommm5index(mchym(4),schym(8))
!    call datafrommm5index(mchym(5),schym(9))
    call datafromidx(mchym(4),schym(8))
    call datafromidx(mchym(5),schym(9))
    schym(10) = trim(chym_mfile)         ! Museo db for discharge comparison
    schym(11) = trim(chym_ofile)         ! Output file
    schym(12) = trim(chym_sfile)         ! Static fields file
    schym(13) = trim(chym_pfile)         ! Precipitation file
    schym(14) = trim(chym_tfile)         ! Temperature file
    schym(15) = trim(chym_dsource)       ! Rain data sources
    schym(18) = trim(chym_rfile)         ! Restart file
    schym(16) = chym_ifile2
!    scst = schym(2)
!    if ( len_trim(scst)>50 ) then
!      schym(16) = scst(1:50)    ! RegCM directory files
!      schym(17) = scst(51:)     ! RegCM directory files
!    else
!      schym(16) = scst          ! RegCM directory files
!      schym(17) = ' '           ! RegCM directory files
!    end if
    schym(25) = 'dem'                 ! Digital Elevation Model
    schym(26) = 'lat'                 ! Latitudes Matrix
    schym(27) = 'lon'                 ! Longitudes Matrix
    schym(28) = 'fdm'                 ! Flow direction Matrix (0-8)
    schym(29) = 'acc'                 ! slope (tangent angle for flow direc.)
    schym(30) = 'lus'                 ! Land use Code
    schym(31) = 'aer'                 ! area (Km2) of each cell
    schym(32) = 'dra'                 ! Total Drained Area of each cell (Km2)
    schym(33) = 'run'                 ! Runoff Time
    schym(56) = 'rgwt'                ! Restart ground water field
    schym(57) = 'rh2o'                ! Restart total water content
    schym(58) = 'rsno'                ! Restart accumulated snow field
    schym(59) = 'rpor'                ! Restart discharge
    schym(60) = 'rdgw'                ! Restart deep water content
    ndyn=0
    do i=1,(len_trim(dflcodes)+1)/4
       call field2save(mchym,schym,dflcodes((i-1)*4+1:(i-1)*4+3),j)
       if (j.eq.1) then
          ndyn=ndyn+1
          schym(38+ndyn)=dflcodes((i-1)*4+1:(i-1)*4+3)
          write(6,'(20x,a)') schym(38+ndyn) &
            (1:len_trim(schym(38+ndyn)))//' field will be saved'
       endif
    enddo
    mchym(8)=ndyn
!    call mvsetflags('Century',float(savecen))
!    write (chcentur,'(i10)') mchym(04) !Fabio
!    print*,'CENTURYYYYYYYY:', chcentur
!    call chkcentury(chcentur) !Fabio
!    call chymstdomain(0)            ! Need to comunicate with mvilbdata
    cpar( 1)=cpar1
    cpar( 2)=cpar2
    cpar( 3)=cpar3
    cpar( 4)=cpar4
    cpar( 5)=cpar5
    cpar( 6)=cpar6
    cpar( 7)=cpar7
    cpar( 8)=cpar8
    cpar( 9)=cpar9
    cpar(10)=cpar10

    infi = infiltr
    infi(lago)=infi_lago
    infi(fiume)=infi_fiume
    infi(12)=infi_ice
!
!  Reading manning coefficienties
!
    if ( iunit<0 ) then
      call getlun(iunit)
    else
      close (iunit)
    end if
    open (iunit,file=chym_manning,status='old',err=101,action='read')
    read (iunit,*,end=102)
    read (iunit,*,end=102)
    read (iunit,*,end=102)
    do k=1,101
       read (iunit,'(2x,i3,2x,f5.3,1x,a)',end=103) test1,manning(k),nameluse
!       write(6,'(2x,i3,2x,f5.3,1x,a)'),test1,manning(k),nameluse
    enddo
    end if
    call mpi_bcast(lon(1,1),nlon*nlat,MPI_REAL, 0,mycomm,mpierr)
    call mpi_bcast(lat(1,1),nlon*nlat,MPI_REAL, 0,mycomm,mpierr)
    call mpi_barrier(mycomm,mpierr)
    close(iunit)
    write (6,'(12x,a)') 'Done.'
    return
101 write (6,'(12x,a)') 'Error opening manning coefficienties.'
102 write (6,'(12x,a)') 'Error reading manning coefficienties header.'
103 write (6,'(12x,a)') 'Error reading manning coefficienties.'
    call exit(0)
  end subroutine acquirescriptpar

  subroutine definerainsources
    implicit none
    integer :: i , istat , nn, n
    real :: r
    call takefields(schym(1),',',sources,n)
    do i = 1 , n
      call cv2lower(sources(i),sources(i))
      call nospace(sources(i))
      if ( trim(sources(i))=='intdb' ) then
        srcflag(1) = .true.
      else if ( trim(sources(i))=='radar' ) then
        srcflag(2) = .true.
      else if ( trim(sources(i))=='era' ) then
        srcflag(3) = .true.
      else if ( trim(sources(i))=='wrf1' ) then
        srcflag(4) = .true.
      else if ( trim(sources(i))=='wrf2' ) then
        srcflag(5) = .true.
      else if ( trim(sources(i))=='regcm' ) then
        srcflag(7) = .true.
      else if ( trim(sources(i))=='persiann' ) then
        srcflag(8) = .true.
      else if ( trim(sources(i))=='trmm' ) then
        srcflag(9)=.true.
      else if ( trim(sources(i))=='ein75' ) then
        srcflag(10)=.true.
      else if ( trim(sources(i))=='ein15' ) then
        srcflag(11)=.true.
      else if ( trim(sources(i))=='era5' ) then
        srcflag(14)=.true.
      else if ( trim(sources(i))=='txt' ) then
        srcflag(15)=.true.
      else if ( trim(sources(i))=='mm5file' ) then
        call getlun(lu(4))
        open(lu(4),status='old',form='unformatted',file=schym(2), &
    iostat=istat)
	    if (istat.eq.0) then
	      srcflag(12)=.true.
	      write (6,'(16x,a)') schym(2)(1:len_trim(schym(2)))// &
        ' MM5 output file will be used.'
	    else
	      close(lu(4))
	      write (6,'(16x,a)') 'Cannot open '// &
        schym(2)(1:len_trim(schym(2)))//', MM5 output file will NOT be used.'
        endif
        srcflag(12) = .true.
      else if ( trim(sources(i))=='museo' ) then
        srcflag(13) = .true.
      else if ( trim(sources(i))=='acqwa01' ) then
        srcflag(21) = .true.
      else if ( trim(sources(i))=='acqwa02' ) then
        srcflag(22) = .true.
      else if ( trim(sources(i))=='acqwa03' ) then
        srcflag(23) = .true.
      else if ( trim(sources(i))=='acqwa04' ) then
        srcflag(24) = .true.
      else if ( trim(sources(i))=='acqwa5a' ) then
        srcflag(25) = .true.
      else if ( trim(sources(i))=='acqwa5b' ) then
        srcflag(26) = .true.
      else if ( trim(sources(i))=='acqwa6a' ) then
        srcflag(27) = .true.
      else if ( trim(sources(i))=='acqwa6b' ) then
        srcflag(28) = .true.
      else if ( trim(sources(i))=='acqwac1' ) then
        srcflag(29) = .true.
      else if ( trim(sources(i))=='intdb-fil' ) then
        srcflag(30) = .true.
      else if ( trim(sources(i))=='intdb-fil-nc' ) then
        srcflag(31) = .true.
      else if ( sources(i)(1:1) /= ' ' .and. len_trim(sources(i)) > 0 ) then
        call chymerror(2,i,r,trim(sources(i)))
      end if
    end do
  end subroutine definerainsources
  subroutine caparameters(res)
    implicit none
    real :: res,d
    integer :: i,j
    write (6,'(15x,a)') 'Defining Weights for Cellular routines.'
    i=nlon/2
    j=nlat/2
    d=distance(lat(i,j),lon(i,j),lat(i+ir(1),j+jr(1)),lon(i+ir(1),j+jr(1)))
    call mvsetflags('Cellular Weights NO',1./d)
    d=distance(lat(i,j),lon(i,j),lat(i+ir(2),j+jr(2)),lon(i+ir(2),j+jr(2)))
    call mvsetflags('Cellular Weights N',1./d)
    d=distance(lat(i,j),lon(i,j),lat(i+ir(3),j+jr(3)),lon(i+ir(3),j+jr(3)))
    call mvsetflags('Cellular Weights NE',1./d)
    d=distance(lat(i,j),lon(i,j),lat(i+ir(4),j+jr(4)),lon(i+ir(4),j+jr(4)))
    call mvsetflags('Cellular Weights E',1./d)
    d=distance(lat(i,j),lon(i,j),lat(i+ir(5),j+jr(5)),lon(i+ir(5),j+jr(5)))
    call mvsetflags('Cellular Weights SE',1./d)
    d=distance(lat(i,j),lon(i,j),lat(i+ir(6),j+jr(6)),lon(i+ir(6),j+jr(6)))
    call mvsetflags('Cellular Weights S',1./d)
    d=distance(lat(i,j),lon(i,j),lat(i+ir(7),j+jr(7)),lon(i+ir(7),j+jr(7)))
    call mvsetflags('Cellular Weights SO',1./d)
    d=distance(lat(i,j),lon(i,j),lat(i+ir(8),j+jr(8)),lon(i+ir(8),j+jr(8)))
    call mvsetflags('Cellular Weights O',1./d)
    return
  end subroutine caparameters

  subroutine field2save(mchym,schym,field,flag)
    implicit none
    character field*3,schym(nchymp)*150 ; integer flag,mchym(nchymp)
    character fields(12)*3 ; integer nfields,i
    flag=0
    if (mchym(8).le.0) then
       call takefields(schym(7),',',fields,nfields)
       do i=1,nfields
          if (fields(i).eq.field) then
             flag=1 ; return
          endif
       enddo
    else
       do i=1,mchym(8)
          if (schym(38+i).eq.field) then
             flag=1 ; return
          endif
       enddo
    endif
    return
  end subroutine field2save

  subroutine calibration
    implicit none

!    cpar( 1)=4e-07! Return flow factor (4.8e-07)
!!   cpar( 1)=4.8e-07! Return flow factor (4.8e-07)
!    cpar( 2)=0.0015 ! Alpha coefficients for hydraulic radius (0.0015)
!    cpar( 3)=0.050  ! Beta coefficients for hydraulic radius (0.050)
!    cpar( 4)=0.050  ! Melting temperature factor (0.050)
!    cpar( 5)=0.0094 ! Melting shortwave rad. factor (0.0094)
!    cpar( 6)=500.0  ! River/land threshold (Km2) (500.0)
!    cpar( 7)=90.0   ! Number of days to consider for return flow (90)
!    cpar( 8)=4.5    ! Reduction of land/channel manning coefficient
!    cpar( 9)=200.0  ! River/land threshold (Km2) for returnflow
!    cpar(10)=0.0    ! Not yet used
!
!    infi= 40.0
!    infi(lago)=0.0
!    infi(fiume)=0.0
!    infi(12)=0.0    ! Ice

    perc=infi*0.01
    call basincalibration
    return
  end subroutine calibration

  subroutine basincalibration
    implicit none
    integer len_trim
    if (schym(10)(1:len_trim(schym(10))).eq.'acqwapo') then
       write(6,'(12x,a)') 'Adopting calibration parameters for Po basin '    &
           //'- ACQWA Project'
    else if (schym(10)(1:len_trim(schym(10))).eq.'acqwarhone') then
       write(6,'(12x,a)') 'Adopting calibration parameters for Rhone basin ' &
          //'- ACQWA Project'
       cpar(1)=2.0e-06      ! Return flow factor (4.8e-07)
       cpar(4)=0.03         ! Melting temperature factor (0.050)
       cpar(5)=0.00540      ! Melting shortwave rad. factor (0.0094)
       cpar(9)=200.0        ! River/land threshold (Km2) for returnflow
    endif
    return
  end subroutine basincalibration

  subroutine read_namelist
       namelist /chymconfig/ nlon , nlat , slon , slat , dij , chym_radius ,       &
               nsave , demf , angiocycle , chym_sdate , chym_edate , chym_steps ,  &
               chym_tempfl , chym_river , chym_savet , chym_savep ,                &
               chym_restart , chym_grads , chym_modis , chym_tplot ,               &
               chym_iplot , chym_zoom , chym_verbose , chym_mfile ,                &
               chym_ofile , chym_sfile , chym_rfile , rsave ,                      &
               chym_pfile , chym_tfile , chym_dsource , chym_ifile1 ,              &
               chym_ifile2 , chym_regcm_rad , chem_symtype ,                       &
               chym_savefld , chym_netcdf , angionp , threshdr , numrivdr ,        &
               numrunave , uphill , ncyc1 , ncyc2 , ncyc3,                         &
               integrflag , cpar1 , cpar2 , cpar3 , cpar4 , cpar5 , cpar6 ,        &
               cpar7 , cpar8 , cpar9 , cpar10 , infiltr , infi_lago ,              &
               infi_fiume , infi_ice , chym_manning , deflate_level ,              &
               num_chunky , num_chunkx , packnetcdf

       call getarg(1, namelistfile)
       call getlun(lun)
       open(lun,file=namelistfile,status='old',action='read',iostat=iretval)
       if ( iretval /= 0 ) then
         write(6,*) 'Error opening namelist file'//trim(namelistfile)
         stop
       end if
       write(6,'(7x,a)') 'Reading model namelist'
       rewind(lun)
       read(lun,nml=chymconfig,iostat=iretval)
       if ( iretval /= 0 ) then
         write(6,*) 'Error reading namelist'
         stop
       end if
  end subroutine read_namelist

  subroutine compute_scale_and_offset(minv, maxv, nbit)
    implicit none
    real, intent(in) :: minv,maxv
    integer, intent(in) :: nbit
    scale_factor_por = (maxv - minv) / (2 ** nbit - 1)
    ! translate the range to be symmetric about zero
    add_offset_por = minv + 2 ** (nbit - 1) * scale_factor_por
    print*,"The scale_factor_por and add_offset_por are: ",scale_factor_por,add_offset_por
  end subroutine compute_scale_and_offset

end module mod_param
