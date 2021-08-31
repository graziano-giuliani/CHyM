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

module mod_crtstatic


  use mod_libmv
  use mod_internal
  use mod_museo
  use mod_statparams
  use mod_vector
  use mod_runparams
  use mod_ncio
  use mod_strman
  use mod_cellaut
  use mod_param
  use mod_chymlib
  use mod_phys

  contains

  subroutine builddem
    implicit none
    dem = -1.0
    write (6,'(/12x,a)') 'Builting Digital Elevation model from MuSEO db.'
    if ( mchym(13)==2 .or. mchym(13)==3 .or. mchym(13)>=21 ) call worlddem
    if ( mchym(13)==1 .or. mchym(13)==3 ) call italydem
    if ( mchym(13)==4 ) call nasadem
    if ( mchym(13)==7 ) call asterdem
    if ( mchym(13)==9 ) call hydrodem
    if ( mchym(13)==10 ) then 
       call worlddem
       call italydem
!       call hydrodem_gen
!       call hydrodem_1k
!       call chkdemhole(dem)
!       call areamatrix
!       call buildlandusemap
!       call buildacclivitymap             ! Fill accl matrix
!!       call smoothHydroDEM
!       call demsmoothing
!       call chkdemhole(dem)
!       call buildacclivitymap             ! Fill accl matrix
!       call plotriverbasin
    end if
!    if ( mchym(13)==21 ) call rhonebasin
!    if ( mchym(13)/=10 ) then 
       call chkdemhole(dem)
!    end if
    write (6,'(12x,a)') 'Done.'
  end subroutine builddem

  subroutine italydem
    implicit none
    integer , parameter :: nmax = 145161
    real , dimension(nmax) :: xla , xlo , h
    real , dimension(8) :: hc
    integer :: i , j , k , lfile , lun , n , np , nz, len_trim
    character(len=240) :: cfile
    character(len=3) :: str
    write (6,'(15x,a)') 'Italian DEM will be used.'
    call getlun(lun)
    call openmuseofiles(lun,'demita.bin',0)
    np=0
    cfile='Used data of zones: '
    lfile=20
    do nz=1,300
       read(lun,end=100) n,(xla(i),i=1,n),(xlo(i),i=1,n),(h(i),i=1,n)
       if (n.gt.0) then
          if (.not.(minval(xlo,n).gt.rchym(3).or.maxval(xlo,n).le.rchym(1)  &
            .or.   minval(xla,n).gt.rchym(4).or.maxval(xla,n).le.rchym(2)   &
                   )) then
             do k=1,n
                call locateij(xla(k),xlo(k),rchym(2),rchym(1),rchym(5),     &
                    rchym(5),mchym(3),mchym(2),i,j)
                if (i.gt.0) then
                   dem(i,j)=h(k)
                   np=np+1
                endif
             end do
             write(str,'(i3)') nz
             cfile=cfile(1:lfile)//str//' '
             lfile=lfile+4
             if (lfile.gt.60) then
                write (6,'(15x,a)') cfile(1:len_trim(cfile))
                cfile='        '
                lfile=8
             endif
          endif
       endif
    end do
100 close(lun)
    write (6,'(15x,a)') cfile(1:len_trim(cfile)-1)
    return
  end subroutine italydem

  subroutine worlddem
    implicit none
    character(len=7) :: cfile
    integer :: ilx , ily
    real :: la1 , la2 , lo1 , lo2
    character(len=1) :: ns , we
    write (6,'(15x,a)') 'World USGS DEM will be used.'
    do ilx = -180 , 140 , 40
      if ( ilx>0 ) then
        we = 'E'
      else
        we = 'W'
      end if
      do ily = 90 , -10 , -50
        if ( ily>0 ) then
          ns = 'N'
        else
          ns = 'S'
        end if
        if ( abs(ilx)>=100 ) then
          write (cfile,'(a1,i3,a1,i2)') we , iabs(ilx) , ns , iabs(ily)
        else
          write (cfile,'(a2,i2,a1,i2)') we//'0' , iabs(ilx) , ns , iabs(ily)
        end if
        lo1 = ilx
        lo2 = ilx + 40
        la1 = ily - 50
        la2 = ily
        if ( lo1<=rchym(3) .and. lo2>rchym(1) .and. la1<=rchym(4) .and. la2>rchym(2) ) then
          write (6,'(18x,a,4f9.1)') 'Using tile '//cfile , lo1 , lo2 , la1 , la2
          call readtiles(cfile)
        end if
      end do
    end do
    la1 = -90
    la2 = -60
    do ilx = -180 , 120 , 60
      if ( ilx>0 ) then
        we = 'E'
      else
        we = 'W'
      end if
      if ( abs(ilx)>=100 ) then
        write (cfile,'(a1,i3,a3)') we , iabs(ilx) , 'S80'
      else if ( abs(ilx)>=10 ) then
        write (cfile,'(a2,i2,a3)') we//'0' , iabs(ilx) , 'S80'
      else
        cfile = 'W000S80'
      end if
      lo1 = ilx
      lo2 = ilx + 60
      if ( lo1<=rchym(3) .and. lo2>rchym(1) .and. la1<=rchym(4) .and. la2>rchym(2) ) then
        write (6,'(18x,a,4f9.1)') 'Using tile '//cfile , lo1 , lo2 , la1 , la2
        call readtiles(cfile)
      end if
    end do
  end subroutine worlddem

  subroutine nasadem
    implicit none
    integer , parameter :: ndat = 1201
    real :: dxy , xla , xlo
    character(len=60) :: cfile
    integer(2) , dimension(ndat,ndat) :: h , nm
    integer :: i , ii , ila , ilo , j , jj , lat1 , lat2 , lon1 , lon2, lun, nadat
    write (6,'(15x,a)') 'NASA high resolution DEM will be used.'
    call getlun(lun)
    dxy = 1./ndat
    if (rchym(2) >= 0) then
      lat1=int(rchym(2))
    else
      lat1=int(rchym(2))-1
    end if
    if (rchym(4) >= 0) then
      lat2=int(rchym(4))
    else
      lat2=int(rchym(4))-1
    end if
    if (rchym(1) >= 0) then
      lon1=int(rchym(1))
    else
      lon1=int(rchym(1))-1
    end if
    if (rchym(3) >= 0) then
      lon2=int(rchym(3))
    else
      lon2=int(rchym(3))-1
    end if
    nm = 0
    do ila = lat1 , lat2
      do ilo = lon1 , lon2
        call readnasadem(ila,ilo,h,ndat,nadat)
        do ii = 1 , nadat
          do jj = 1 , nadat
            xla = ila + (jj-1)*dxy
            xlo = ilo + (ii-1)*dxy
            call locateij(xla,xlo,rchym(2),rchym(1),rchym(5),rchym(5), &
               mchym(3),mchym(2),i,j)
            if ( i>0 .and. j>0 .and. i<=nlon .and. j<=nlat .and. h(ii,jj)>=0 ) &
                 then
              if ( nm(i,j)==0 ) then
                nm(i,j) = 1
                dem(i,j) = h(ii,jj)
              else
                nm(i,j) = nm(i,j) + 1
                dem(i,j) = dem(i,j) + h(ii,jj)
              end if
            end if
          end do
        end do
      end do
    end do
    do i = 1 , nlon
      do j = 1 , nlat
        if ( nm(i,j)>1 ) dem(i,j) = dem(i,j)/nm(i,j)
      end do
    end do
  end subroutine nasadem

  subroutine asterdem
    implicit none
    integer , parameter :: ndat = 3601
    real :: dxy , xla , xlo
    character(len=60) :: cfile
    integer(2) , dimension(ndat,ndat) :: h
    integer , allocatable , dimension (:,:) :: nm
    integer :: i , ii , ila , ilo , j , jj , lat1 , lat2 , lon1 , lon2, lun,nadat,ifound
    allocate(nm(mchym(2),mchym(3)))
    write (6,'(15x,a)') 'ASTER 1-arcosec (30m) high resolution DEM will be used.'
    call getlun(lun)
    dxy = 1./ndat
    if (rchym(2) >= 0) then
      lat1=int(rchym(2))
    else
      lat1=int(rchym(2))-1
    end if
    if (rchym(4) >= 0) then
      lat2=int(rchym(4))
    else
      lat2=int(rchym(4))-1
    end if
    if (rchym(1) >= 0) then
      lon1=int(rchym(1))
    else
      lon1=int(rchym(1))-1
    end if
    if (rchym(3) >= 0) then
      lon2=int(rchym(3))
    else
      lon2=int(rchym(3))-1
    end if
    nm = 0
    do ila = lat1 , lat2
      do ilo = lon1 , lon2
        ifound = 0
        call readasterdem(ila,ilo,h,ndat,ifound)
        if (ifound == 1) cycle
        do ii = 1 , ndat
          do jj = 1 , ndat
            xla = ila + (jj-1)*dxy
            xlo = ilo + (ii-1)*dxy
            call locateij(xla,xlo,rchym(2),rchym(1),rchym(5),rchym(5), &
               mchym(3),mchym(2),i,j)
            if ( i>0 .and. j>0 .and. i<=nlon .and. j<=nlat .and. h(ii,jj)>=0 ) &
                 then
              if ( nm(i,j)==0 ) then
                nm(i,j) = 1
                dem(i,j) = h(ii,jj)
              else
                nm(i,j) = nm(i,j) + 1
                dem(i,j) = dem(i,j) + h(ii,jj)
              end if
            end if
          end do
        end do
      end do
    end do
    do i = 1 , nlon
      do j = 1 , nlat
        if ( nm(i,j)>1 ) dem(i,j) = dem(i,j)/nm(i,j)
      end do
    end do
    deallocate(nm)
  end subroutine asterdem

  subroutine readnasadem(ila,ilo,h,ndat,nadat)
    implicit none
    integer :: ndat,nadat,lun,len_trim,ila,ilo,i,j
    integer(2) , dimension(ndat,ndat) :: h
    character(len=80) :: cfile
    call getlun(lun)
    call demnasa(float(ila),float(ilo),-lun,cfile)
    write(6,'(21x,a)') 'Reading '//cfile(1:len_trim(cfile))
    nadat=0
    open(lun,file=cfile,status='old',form='unformatted',access='stream', &
           err=100)
    read (lun) ((h(i,j),i=1,ndat),j=ndat,1,-1)
    nadat=ndat
100 close(lun)
    return
  end subroutine readnasadem

  subroutine readasterdem(ila,ilo,h,ndat,ifound)
    implicit none
    integer :: ndat,lun,len_trim,ila,ilo,i,j,ifound
    integer(2), dimension(ndat,ndat) ::h
    character(len=150) :: cfile
    call getlun(lun)
    call getastername(float(ila),float(ilo),cfile)
    write(6,'(21x,a)') 'Reading '//cfile(1:len_trim(cfile))
    call getasterdata(cfile, ndat, h, ifound)
    return
  end subroutine readasterdem

  subroutine getastername(xlat,xlon,lfile)
    implicit none
    character(len=150) :: lfile
    character(len=100) :: dir
    character(len=4) :: sn, ew, lats, lons
    real :: xlat,xlon 
    integer :: lat,lon
    lat=nint(abs(xlat)) ; lon=nint(abs(xlon))
    if (lat.gt.90.or.lon.gt.180) then
       write(6,'(10x,a)') 'Invalid Lat/Lon Range passed to demASTER routine'
       call exit(0)
    endif
    !call museodir(dir)
!    dir="/home/clima-archive/afantini/flood/DEM/ASTER_DEM/data/"
    dir="museo/DEM/ASTER/"
    call integer2string(lat,2,lats) ; call integer2string(lon,3,lons)
    sn='N' ; if (xlat.lt.0.0) sn='S'
    ew='E' ; if (xlon.lt.0.0) ew='W'
    lfile=trim(dir)//"ASTGTM2_"//sn(1:1)//lats(1:2)// &
         ew(1:1)//lons(1:3)//'_dem.nc'
    return
  end subroutine getastername

  subroutine readtiles(tile)
    implicit none
    character(len=7) :: tile
    intent (in) tile
    integer , parameter :: nxm = 7200
    integer(2) , dimension(nxm) :: idem
    real , dimension(nlon,nlat) :: w
    real :: delx , xslat , xslon , demv
    integer :: ny , nx , i , j
    character(len=60) :: cfile
    integer :: k , l , lun , np
    real :: xla , xlo
    call getlun(lun)
    call openmuseofiles(lun,'DEM/USGS_1km/'//tile//'.dem',0)
    read (lun) nx , ny , delx , xslon , xslat
    np = 0
    w = -1
    do j = 1 , ny
      read (lun) (idem(i),i = 1,nx)
      do i = 1 , nx
        if ( idem(i)>0 ) then
          demv = idem(i)
        else if ( idem(i)>-1000 ) then
          demv = 0.1
        else
          demv = -0.1
        end if
        xla = xslat - (j-1)*delx
        xlo = xslon + (i-1)*delx
        k = nint((xlo-rchym(1))/rchym(5)) + 1
        l = nint((xla-rchym(2))/rchym(5)) + 1
        if ( k>0 .and. l>0 .and. k<=nlon .and. l<=nlat ) then
          if ( w(k,l)<0.0 ) then
            dem(k,l) = demv
            w(k,l) = distance(xla,xlo,lat(k,l),lon(k,l))
            np = np + 1
          else if ( w(k,l)>distance(xla,xlo,lat(k,l),lon(k,l)) ) then
            dem(k,l) = demv
            w(k,l) = distance(xla,xlo,lat(k,l),lon(k,l))
          end if
!--------------------------  Extra lines -------------------------------------
!       if (k.eq.293.and.l.eq.176) then
!       write (6,'(7f11.5))') dem(k,l),xla,xlo,lat(k,l),lon(k,l),
!       2                distance(xla,xlo,lat(k,l),lon(k,l)),w(k,l)
!       endif
!-----------------------------------------------------------------------------
        end if
      end do
    end do
    close (lun)
    write (cfile,'(i10,a)') np , ' points used to fill CHYM matrix'
    call noinspace(cfile)
    write (6,'(21x,a)') trim(cfile)
  end subroutine readtiles 

  subroutine hydrodem_gen
    implicit none
    character(len=150) :: cfile
    real , allocatable , dimension (:,:) :: h
!    integer , allocatable , dimension (:,:) :: dir
    double precision , allocatable , dimension (:) :: lond ,latd
    integer :: ii, jj, i, j, nlond , nlatd, ifound
    if (.not.allocated(fmap1)) allocate(fmap1(nlon,nlat))
    print*,"we are going to read the direction matrix from file"
    cfile = "output/SAMCORDEXWORLD_006degree_stkSA_CORDEXWORLD_DEM_0.06degree_fdm.nc"
!    call gethydrodata_int(cfile, nlond, nlatd,dir,lond,latd,ifound)
    
!    print*,"dir matrix read",maxval(dir),minval(dir)
    ii = 1 ; jj = 1
!    print*,"int((rchym(1)-rchym(1))/rchym(5))+1", int((rchym(1)-rchym(1))/rchym(5))+1
!    print*,"int((rchym(3)-rchym(1))/rchym(5))+1", int((rchym(3)-rchym(1))/rchym(5))+1
!    print*,"int((rchym(2)-rchym(2))/rchym(5))+1", int((rchym(2)-rchym(2))/rchym(5))+1
!    print*,"int((rchym(4)-rchym(2))/rchym(5))+1", int((rchym(4)-rchym(2))/rchym(5))+1
!      do i = int((rchym(1)-rchym(1))/rchym(5))+1, int((rchym(3)-rchym(1))/rchym(5))+1
!        do j =int((rchym(2)-rchym(2))/rchym(5))+1, int((rchym(4)-rchym(2))/rchym(5))+1
!  !        fmap1(ii,jj) = int(dir(i,j))
!  !        fmap1(ii,jj) = int(demd(i,j))
!          fmap1(ii,jj) = int(fmap(i,j))
!          if (fmap1(ii,jj) == 255) print*,ii,jj,i,j
!          jj = jj + 1
!        end do
!        jj = 1
!        ii = ii + 1
!      end do
    fmap1 = fmap
    do i = 1 , nlon
      fmap1(i,1) = -1
      fmap1(i,nlat) = -1
      noflow(i,1) = 0
      noflow(i,nlat) = 0
    end do
    do j = 1 , nlat
      fmap1(1,j) = -1
      fmap1(nlon,j) = -1
      noflow(1,j) = 0
      noflow(nlon,j) = 0
    end do
    fmap = fmap1
  end subroutine hydrodem_gen

  subroutine hydrodem_1k
    implicit none
    character(len=150) :: cfile,nfile
    real , allocatable , dimension (:,:) :: h
    real , allocatable , dimension (:,:) :: dir
    double precision , allocatable , dimension (:) :: lond ,latd
    integer :: ii, jj, i, j, nlond , nlatd, ifound
    real :: minlon,maxlon,minlat,maxlat
    if (rchym(1) < -32.00417 .and. rchym(1) > -144.9958  .and.                 &
        rchym(3) < -32.00417 .and. rchym(3) > -144.9958  .and.                 & 
        rchym(2) <  59.99583 .and. rchym(2) >  -55.99583 .and.                 &
        rchym(4) <  59.99583 .and. rchym(4) >  -55.99583 ) then
      write (6,'(15x,a)') 'America HydroSHEDS (1KM) DEM and DIR will be used.'
      nfile = "ca_na_sa.nc"
      minlon=-144.9958 ;maxlon=-32.00417 ;minlat=-55.99583 ;maxlat=59.99583
    elseif ( rchym(1) < 179.9958  .and. rchym(1) > 112.0042  .and.             &
             rchym(3) < 179.9958  .and. rchym(3) > 112.0042  .and.             & 
             rchym(2) < -10.00417 .and. rchym(2) > -55.99583 .and.             &
             rchym(4) < -10.00417 .and. rchym(4) > -55.99583 ) then
      write (6,'(15x,a)') 'Australia HydroSHEDS (1KM) DEM and DIR will be used.'
      nfile = "au.nc"
      minlon=112.0042 ;maxlon=179.9958 ;minlat=-55.99583 ;maxlat=-10.00417
    elseif ( rchym(1) < 179.9958  .and. rchym(1) > -18.99583  .and.             &
             rchym(3) < 179.9958  .and. rchym(3) > -18.99583  .and.             & 
             rchym(2) <  59.99583 .and. rchym(2) > -34.99583  .and.             &
             rchym(4) <  59.99583 .and. rchym(4) > -34.99583 ) then
      write (6,'(15x,a)') 'EurAfrAsia HydroSHEDS (1KM) DEM and DIR will be used.'
      nfile = "af_as_eu.nc"
      minlon=-18.99583 ;maxlon=179.9958 ;minlat=-34.99583 ;maxlat=59.99583
    else
      write (6,'(15x,a)') 'The domain limits are outside the HydroSHEDS domain '
      write (6,'(15x,a)') 'availability, please use a different DEM or keep the ' 
      write (6,'(15x,a)') 'domain boundaries among the following limits: '
      write (6,'(15x,a)') 'Eurasia and Africa (minlon=-18.99583 ;maxlon=179.9958 ;minlat=-34.99583 ;maxlat=59.99583)'
      write (6,'(15x,a)') 'America (minlon=-144.9958 ;maxlon=-32.00417 ;minlat=-55.99583 ;maxlat=59.99583)'
      write (6,'(15x,a)') 'Australia and NewZeland (minlon=112.0042 ;maxlon=179.9958 ;minlat=-55.99583 ;maxlat=-10.00417)'
      call exit(1)
    end if 
    cfile = "museo/DEM/HydroSHEDS_30s/"//nfile
    ii = 1
    jj = 1
    if (.not.allocated(fmap1)) allocate(fmap1(nlon,nlat))
    call gethydrodata_real(cfile, nlond, nlatd,h,lond,latd,ifound)
!    write (6,'(15x,a,a)') cfile,' reading completed'
    cfile = "museo/FLOWDIR/HydroSHEDS_30s/"//nfile
    call gethydrodata_real(cfile, nlond, nlatd,dir,lond,latd,ifound)
!    write (6,'(15x,a,a)') cfile,' reading completed'
      do i = int((rchym(1)-minlon)/0.0083333) , int((rchym(3)-minlon)/0.0083333)
        do j =int((rchym(2)-minlat)/0.0083333),int((rchym(4)-minlat)/0.0083333)
          dem(ii,jj) = h(i,j)
          jj = jj + 1
        end do
        jj = 1
        ii = ii + 1
      end do
      if (nfile == "ca_na_sa.nc") minlon=-137.9958
      ii = 1
      jj = 1
      do i = int((rchym(1)-minlon)/0.0083333) , int((rchym(3)-minlon)/0.0083333)
        do j =int((rchym(2)-minlat)/0.0083333),int((rchym(4)-minlat)/0.0083333)
          fmap1(ii,jj) = int(dir(i,j))
          jj = jj + 1
        end do
        jj = 1
        ii = ii + 1
      end do
      write (6,'(15x,a,a)') 'HydroSHEDS reading completed'
      where(dem .le. 0) dem = -0.1
      where(fmap1 .eq. 8) fmap1 = 7
      where(fmap1 .eq. 4) fmap1 = 6
      where(fmap1 .eq. 1) fmap1 = 4
      where(fmap1 .eq. 2) fmap1 = 5
      where(fmap1 .eq. 16) fmap1 = 8
      where(fmap1 .eq. 32) fmap1 = 1
      where(fmap1 .eq. 64) fmap1 = 2
      where(fmap1 .eq. 128) fmap1 = 3
      where(fmap1 .eq. 247) fmap1 = 0
      where(fmap1 .eq. 255) fmap1 = 0
      do i = 1 , nlon
        fmap1(i,1) = -1
        fmap1(i,nlat) = -1
        noflow(i,1) = 0
        noflow(i,nlat) = 0
      end do
      do j = 1 , nlat
        fmap1(1,j) = -1
        fmap1(nlon,j) = -1
        noflow(1,j) = 0
        noflow(nlon,j) = 0
      end do
      fmap = fmap1
  end subroutine hydrodem_1k

  subroutine hydrodem
    implicit none
    integer , parameter :: ndat = 6000
    real :: dxy , xla , xlo
    character(len=60) :: cfile
!    integer(2) , dimension(ndat,ndat) :: h
    integer , allocatable , dimension (:,:) :: nm , h
    real , allocatable , dimension (:) :: lond ,latd
    integer :: i , ii , ila , ilo , j , jj , lat1 , lat2 , lon1 , lon2, lun,nadat
    integer :: ila1, ilo1 , nlond , nlatd, ifound
    allocate(nm(mchym(2),mchym(3)))
    write (6,'(15x,a)') 'HydroSHEDS (90m) high resolution DEM will be used.'
    call getlun(lun)
    dxy = 5./ndat
    if (rchym(2) >= 0) then
      lat1=int(rchym(2))
    else
      lat1=int(rchym(2))-1
    end if
    if (rchym(4) >= 0) then
      lat2=int(rchym(4))
    else
      lat2=int(rchym(4))-1
    end if
    if (rchym(1) >= 0) then
      lon1=int(rchym(1))
    else
      lon1=int(rchym(1))-1
    end if
    if (rchym(3) >= 0) then
      lon2=int(rchym(3))
    else
      lon2=int(rchym(3))-1
    end if
    nm = 0
    do ila = lat1 , lat2
      do ilo = lon1 , lon2
        ifound = 0
        if (mod(ila,5)>=1 .and. mod(ila,5)<=4 .or. mod(ila,5)>=-4 .and. mod(ila,5)<=-1) then
           if (ila >= 0) then
             ila1 = ila - mod(ila,5)
           else
             ila1 = ila - mod(ila,5)-5
           end if
        else
           ila1 = ila
        end if
        if (mod(ilo,5)>=1 .and. mod(ilo,5)<=4 .or. mod(ilo,5)>=-4 .and. mod(ilo,5)<=-1) then
           if (ilo >= 0) then
             ilo1 = ilo - mod(ilo,5)
           else
             ilo1 = ilo - mod(ilo,5)-5
           end if
        else
           ilo1 = ilo
        end if
        call readhydrodem(ila1,ilo1,h,nlond,nlatd,lond,latd,ifound)
        if ( ifound == 1 ) cycle
        do ii = 1 , nlond
          do jj = 1 , nlatd
!            xla = ila1 + (jj-1)*dxy
!            xlo = ilo1 + (ii-1)*dxy
!            call locateij(xla,xlo,rchym(2),rchym(1),rchym(5),rchym(5), &
            call locateij(latd(jj),lond(ii),rchym(2),rchym(1),rchym(5),rchym(5), &
               mchym(3),mchym(2),i,j)
            if ( i>0 .and. j>0 .and. i<=nlon .and. j<=nlat .and. h(ii,jj)>=0 ) &
                 then
              if ( nm(i,j)==0 ) then
                nm(i,j) = 1
                dem(i,j) = h(ii,jj)
              else
                nm(i,j) = nm(i,j) + 1
                dem(i,j) = dem(i,j) + h(ii,jj)
              end if
            end if
          end do
        end do
      end do
    end do
    do i = 1 , nlon
      do j = 1 , nlat
        if ( nm(i,j)>1 ) dem(i,j) = dem(i,j)/nm(i,j)
      end do
    end do
    deallocate(nm)
  end subroutine hydrodem

  subroutine readhydrodem(ila,ilo,h,nlond,nlatd,lond,latd,ifound)
    implicit none
    integer :: nlond,nlatd,lun,len_trim,ila,ilo,i,j,ifound
!    integer(2), dimension(ndat,ndat) ::h
    integer, allocatable, dimension(:,:) ::h
    real, allocatable, dimension(:) ::lond,latd
    character(len=150) :: cfile
    call getlun(lun)
    call gethydroname(float(ila),float(ilo),cfile)
    write(6,'(21x,a)') 'Reading '//cfile(1:len_trim(cfile))
    call gethydrodata(cfile, nlond, nlatd, h,lond,latd,ifound)
    return
  end subroutine readhydrodem

  subroutine gethydroname(xlat,xlon,lfile)
    implicit none
    character(len=150) :: lfile
    character(len=100) :: dir
    character(len=4) :: sn, ew, lats, lons
    real :: xlat,xlon
    integer :: lat,lon
    lat=nint(abs(xlat)) ; lon=nint(abs(xlon))
    if (lat.gt.90.or.lon.gt.180) then
       write(6,'(10x,a)') 'Invalid Lat/Lon Range passed to dem HydroSHEDS routine'
       call exit(1)
    endif
    !call museodir(dir)
!    dir="/home/clima-archive/afantini/flood/DEM/HydroSHEDS_DEM/netcdf/"
    dir="museo/DEM/HydroSHEDS_CON_3s/"
    call integer2string(lat,2,lats) ; call integer2string(lon,3,lons)
    sn='n' ; if (xlat.lt.0.0) sn='s'
    ew='e' ; if (xlon.lt.0.0) ew='w'
    lfile=trim(dir)//sn(1:1)//lats(1:2)// &
         ew(1:1)//lons(1:3)//'.nc'
!    write(6,'(10x,a)') lfile
    return
  end subroutine gethydroname

  subroutine chkdemhole(dem)
    implicit none
    real , dimension(nlon,nlat) :: dem
    real , dimension(nlon,nlat) :: work , wk1
    integer , dimension(nlon,nlat) :: icl , luse
    real :: err , errmax , errmean
    character(len=60) :: cfile , title
    integer :: i , ic , j , k , nc , nhole , np
    luse = 5
    write (6,'(15x,a)') 'Checking DEM in order to avoid undefined cells.'
    nhole = 0
    do j = 1 , nlat
      do i = 1 , nlon
        icl(i,j) = 0
        if ( dem(i,j)<0.0 ) then
          nc = 0
          if ( i==1 .or. j==1 .or. i==nlon .or. j==nlat ) then
            np = 0
            do k = 1 , 8
              if ( i+ir(k)>=1 .and. i+ir(k)<=nlon .and. j+jr(k)>=1 .and. j+jr(k)   &
                   <=nlat ) then
                if ( dem(i+ir(k),j+jr(k))<0.0 ) nc = nc + 1
                np = np + 1
              end if
            end do
          else
            np = 8
            do k = 1 , 8
              if ( dem(i+ir(k),j+jr(k))<0.0 ) nc = nc + 1
            end do
          end if
          if ( nc<np ) then
            icl(i,j) = 8
            nhole = nhole + 1
          end if
        end if
      end do
    end do
    if ( nhole>5*nlat ) then
      write (title,'(i10,a)') nhole , ' undefined cells found.'
      call no2space(title)
      call noinspace(title)
      write (6,'(15x,a)') title
      write (6,'(15x,a)') 'Calling CellAut algorithm.'
      hmax = maxval(dem)
      ic = 0
      errmean = 1000
      errmax = 1000
      do while ( ic<1000 .and. errmean>0.01 .and. errmax>0.5 )
        wk1 = dem
        call d2cellcycle(dem,icl,work,nlon,nlat,0.1)
        errmean = 0
        errmax = 0
        np = 0
        do i = 1 , nlon
          do j = 1 , nlat
            if ( icl(i,j)>0 ) then
              np = np + 1
              err = (wk1(i,j)-dem(i,j))**2
              if ( err>errmax ) errmax = err
              errmean = errmean + err
            end if
          end do
        end do
        errmean = errmean/np
!     write(6,'(i10,2f15.8)') ic,errmean,errmax
        ic = ic + 1
      end do
      idemcor = 1
    else
      idemcor = 0
    end if
    write (6,'(15x,a)') 'Done.'
  end subroutine chkdemhole

  subroutine buildlandusemap
    implicit none
    integer , parameter :: ncols = 43200
    integer , parameter :: nrows = 21600
    integer(1) , dimension(ncols) :: ivals
    integer , dimension(nlon,nlat) :: iflg
    real :: cost , ylat , ylon
    character(len=70) :: cfile
    integer :: i , irec , irec1 , irec2 , j , lun , n , ncic , nzero
    write (6,'(/12x,a)') 'Building Land Use Map using USGS data.'
    do j = 1 , nlat
      do i = 1 , nlon
        luse(i,j) = 0
      end do
    end do
    cost = 1./120.
    irec1 = nint(1.0+(90.0-rchym(4))/cost) - 4
    irec2 = nint(1.0+(90.0-rchym(2))/cost) + 4
    call getlun(lun)
  ! call openmuseofiles(lun,'bats2_0.marco',0)
    call openmuseofiles(lun,'oge2_0.marco',0)
    do irec = 1 , irec1 - 1
      read (lun)
    end do
    do irec = irec1 , irec2
      read (lun) ivals
      ylat = 90.0 - (irec-1)*cost - 3*cost                     ! correzione
      j = nint((ylat-rchym(2))/rchym(5)) + 1
      if ( j>=1 .and. j<=nlat ) then
        call correctusgs(ivals,irec)
        do n = 1 , ncols
          ylon = -180.0 + (n-1)*cost
          i = nint((ylon-rchym(1))/rchym(5)) + 1
          if ( i>=1 .and. i<=nlon ) then
            luse(i,j) = ivals(n)
          end if
        end do
      end if
    end do
    close (lun)
    ncic = 0
    nzero = 1
    do while ( nzero>0 .and. ncic<10 )
      nzero = 0
      do j = 1 , nlat
        do i = 1 , nlon
          if ( luse(i,j)==0 ) then
            nzero = nzero + 1
            iflg(i,j) = 8
          else
            iflg(i,j) = 0
          end if
        end do
      end do
      if ( nzero>0 ) then
        write (cfile,'(a,i2,a,i7,a)')  &
          'After ',ncic,' d2icellcycle cicles ',nzero, ' undefined points.'
        call noinspace(cfile)
        call no2space(cfile)
        write (6,'(15x,a)') trim(cfile)
        call d2icellcycle(luse,iflg,iwk,nlon,nlat)
        ncic = ncic + 1
      end if
    end do
    write (cfile,'(a,i2,a,i7,a)') &
      'After ',ncic,' d2icellcycle cicles ',nzero,' undefined points.'
    call noinspace(cfile)
    call no2space(cfile)
    write (6,'(15x,a)') trim(cfile)
    write (6,'(15x,a)') 'Correcting LandUse map using DEM model.'
    if (mchym(13) /= 9) then
      do j = 1 , nlat
        do i = 1 , nlon
          if ( dem(i,j)<=0.00001 ) luse(i,j) = mare
        end do
      end do
    end if
    write (6,'(12x,a)') 'Done.'
  end subroutine buildlandusemap

  subroutine correctusgs(ivals,irec)
    implicit none
    integer , parameter :: ncols = 43200
    integer :: irec
    integer(1) , dimension(ncols) :: ivals
    intent (in) irec
    intent (out) ivals
    integer :: ii
    if ( irec<0 .or. irec>ncols ) then
      return
    else if ( irec==5671 ) then                   ! vicino adriatico Prov. di TE
      ivals(23279) = 31
    else if ( irec==5672 ) then
      ivals(23277) = 31
      ivals(23278) = 31
      ivals(23280) = 31
    else if ( irec==5673 ) then
      ivals(23279) = 31
    else if ( irec==5713 ) then                   ! vicino L'Aquila
      ivals(23203) = 30
    else if ( irec==5714 ) then
      ivals(23204) = 30
      ivals(23205) = 30
    else if ( irec==5752 ) then                   ! tra le provincie di CB e CH
      ivals(23378) = 2
    else if ( irec==5753 ) then
      ivals(23378) = 58
      ivals(23379) = 2
      ivals(23381) = 2
    else if ( irec==5754 ) then
      ivals(23380) = 2
    else if ( irec==5755 ) then
      ivals(23377) = 31
    else if ( irec==5756 ) then
      ivals(23379) = 31
    else if ( irec==5757 ) then
      ivals(23374) = 31
      ivals(23375) = 93
    else if ( irec==5758 ) then
      ivals(23376) = 93
      ivals(23377) = 58
    else if ( irec==5759 ) then
      ivals(23369) = 30
    else if ( irec==5760 ) then
      ivals(23371) = 58
    else if ( irec==5761 ) then
      ivals(23368) = 58
    else if ( irec==5762 ) then
      ivals(23367) = 31
      ivals(23369) = 31
      ivals(23370) = 31
    else if ( irec==5763 ) then
      ivals(23365) = 30
      ivals(23366) = 30
      ivals(23368) = 31
      ivals(23369) = 31
    else if ( irec==5764 ) then
      ivals(23364) = 30
      ivals(23367) = 31
    else if ( irec==5765 ) then
      ivals(23364) = 30
      ivals(23365) = 31
      ivals(23367) = 58
    else if ( irec==5766 ) then
      ivals(23363) = 30
      ivals(23366) = 30
    else if ( irec==5767 ) then
      ivals(23364) = 30
      ivals(23365) = 31
    else if ( irec==5769 ) then
      ivals(23362) = 56
      ivals(23363) = 30
    else if ( irec==5770 ) then
      ivals(23361) = 30
      ivals(23363) = 30
      ivals(23364) = 30
    else if ( irec==5771 ) then
      ivals(23362) = 30
      ivals(23363) = 30
    else if ( irec==5772 ) then
      ivals(23357) = 46
    else if ( irec==5773 ) then
      ivals(23358) = 31
      ivals(23359) = 31
    else if ( irec==5774 ) then
      ivals(23355) = 30
      ivals(23356) = 31
    else if ( irec==5775 ) then
      ivals(23354) = 26
      ivals(23355) = 30
      ivals(23357) = 30
    else if ( irec==5776 ) then
      ivals(23356) = 30
    else if ( irec==5777 ) then
      ivals(23353) = 26
    else if ( irec==5778 ) then
      ivals(23353) = 26
      ivals(23354) = 26
    else if ( irec==5779 ) then
      ivals(23352) = 30
      ivals(23355) = 31
    else if ( irec==5780 ) then
      ivals(23353) = 30
      ivals(23354) = 31
    end if

    if ( irec==5345 ) then                   ! cambiato da me!! (Erika Coppola)
      ivals(22551) = 31
    else if ( irec==5346 ) then
      ivals(22551) = 31
      ivals(22552) = 31
      ivals(22553) = 31
    else if ( irec==5347 ) then
      ivals(22551) = 31
      ivals(22552) = 31
      ivals(22553) = 31
    else if ( irec==5348 ) then
      ivals(22553) = 31
      ivals(22554) = 31
    else if ( irec==5349 ) then
      ivals(22565) = 31
    else if ( irec==5395 ) then
      ivals(22931) = 31
      ivals(22932) = 31
      ivals(22933) = 31
      ivals(22934) = 31
      ivals(22882) = 56
      ivals(22883) = 56
      ivals(22884) = 56
      ivals(22885) = 56
      ivals(22892) = 56
      ivals(22893) = 56
      ivals(22894) = 56
      ivals(22904) = 56
      ivals(22943) = 31
      ivals(22953) = 31
    else if ( irec==5396 ) then
      ivals(22885) = 31
    else if ( irec==5397 ) then
      ivals(22885) = 93
      ivals(22956) = 31
      ivals(22957) = 31
      ivals(22958) = 31
    else if ( irec==5398 ) then
      ivals(22885) = 56
      ivals(22886) = 56
      ivals(22887) = 56
      ivals(22958) = 31
      ivals(22959) = 31
      ivals(22960) = 31
    else if ( irec==5399 ) then
      ivals(22885) = 31
      ivals(22886) = 31
      ivals(22887) = 31
      ivals(22961) = 31
      ivals(22962) = 31
    else if ( irec==5400 ) then
      ivals(22884) = 31
      ivals(22885) = 31
      ivals(22887) = 31
      ivals(22888) = 31
      ivals(22963) = 31
    else if ( irec==5401 ) then
      ivals(22883) = 31
      ivals(22884) = 31
      ivals(22885) = 31
      ivals(22886) = 31
      ivals(22965) = 31
      ivals(22966) = 31
      ivals(22967) = 31
      ivals(23020) = 31
      ivals(23035) = 31
      ivals(23036) = 31
      ivals(23042) = 31
      ivals(23043) = 31
      ivals(23044) = 31
      ivals(23045) = 31
      ivals(23046) = 31
      ivals(23047) = 31
    else if ( irec==5402 ) then
      ivals(22965) = 31
      ivals(22966) = 31
      ivals(22967) = 31
      ivals(22883) = 31
      ivals(22884) = 31
      ivals(22885) = 31
      ivals(22886) = 31
      ivals(23020) = 31
      ivals(23035) = 31
      ivals(23036) = 31
      ivals(23037) = 31
      ivals(23042) = 31
      ivals(23043) = 31
      ivals(23044) = 31
      ivals(23045) = 31
      ivals(23046) = 31
      ivals(23047) = 31
    else if ( irec==5403 ) then
      ivals(22968) = 31
      ivals(22969) = 31
      ivals(22882) = 31
      ivals(22883) = 31
      ivals(22884) = 31
      ivals(22885) = 31
      ivals(22886) = 31
      ivals(22887) = 31
      ivals(22970) = 31
      ivals(22971) = 31
      ivals(23017) = 31
      ivals(23018) = 31
      ivals(23022) = 31
      ivals(23023) = 31
      ivals(23024) = 31
      ivals(23036) = 31
      ivals(23037) = 31
      ivals(23038) = 31
      ivals(23044) = 31
      ivals(23045) = 31
      ivals(23046) = 31
      ivals(23047) = 31
      ivals(23048) = 31
    else if ( irec==5404 ) then
      ivals(22880) = 31
      ivals(22881) = 31
      ivals(22882) = 31
      ivals(22883) = 31
      ivals(22884) = 31
      ivals(22885) = 31
      ivals(22971) = 31
      ivals(22972) = 31
      ivals(22973) = 31
      ivals(23019) = 61
      ivals(23025) = 31
    else if ( irec==5405 ) then
      ivals(22881) = 31
      ivals(22882) = 31
      ivals(22883) = 31
      ivals(22884) = 31
      ivals(23012) = 94
    else if ( irec==5406 ) then
      ivals(22879) = 56
      ivals(22880) = 56
      ivals(22882) = 56
      ivals(22883) = 56
      ivals(22973) = 31
      ivals(22974) = 31
      ivals(23013) = 31
      ivals(23014) = 31
    else if ( irec==5407 ) then
      ivals(22879) = 56
      ivals(22881) = 56
      ivals(22882) = 56
      ivals(22878) = 56
      ivals(22973) = 31
      ivals(22974) = 31
      ivals(22975) = 31
      ivals(22987) = 31
    else if ( irec==5408 ) then
      ivals(22874) = 31
      ivals(22875) = 31
      ivals(22876) = 31
      ivals(22877) = 31
      ivals(22880) = 31
      ivals(22975) = 31
      ivals(22988) = 31
      ivals(22989) = 31
      ivals(22990) = 31
      ivals(23008) = 31
    else if ( irec==5409 ) then
      ivals(22865) = 31
      ivals(22866) = 31
      ivals(22867) = 31
      ivals(22868) = 31
      ivals(22869) = 31
      ivals(22870) = 31
      ivals(22871) = 31
      ivals(22872) = 31
      ivals(22873) = 31
      ivals(22874) = 31
      ivals(22875) = 31
      ivals(22876) = 31
      ivals(22877) = 31
      ivals(22878) = 31
      ivals(22879) = 31
      ivals(22992) = 31
      ivals(23003) = 31
      ivals(23004) = 31
      ivals(23009) = 31
      ivals(23010) = 31
    else if ( irec==5410 ) then
      ivals(22867) = 31
      ivals(22868) = 31
      ivals(22869) = 31
      ivals(22870) = 31
      ivals(22871) = 31
      ivals(22872) = 31
      ivals(22873) = 31
      ivals(22874) = 31
      ivals(22875) = 31
      ivals(22876) = 31
      ivals(22993) = 31
      ivals(22994) = 31
      ivals(23005) = 57
      ivals(23006) = 57
    else if ( irec==5411 ) then
      ivals(22999) = 31
      ivals(23000) = 31
    else if ( irec==5412 ) then
      ivals(23001) = 31
      ivals(23002) = 31
    else if ( irec==5413 ) then
      ivals(23001) = 31
      ivals(23002) = 31
    else if ( irec==5394 ) then
      ivals(22902) = 56
      ivals(22903) = 56
      ivals(22905) = 56
      ivals(22906) = 56
      ivals(22907) = 56
      ivals(22927) = 31
      ivals(22928) = 31
      ivals(22929) = 31
      ivals(22930) = 31
      ivals(22931) = 31
      ivals(22932) = 31
      ivals(22936) = 31
      ivals(22937) = 31
      ivals(22910) = 56
      ivals(22911) = 56
      ivals(22912) = 56
      ivals(22913) = 56
      ivals(22880) = 31
      ivals(22881) = 31
      ivals(22882) = 31
      ivals(22883) = 31
      ivals(22886) = 31
      ivals(22887) = 31
      ivals(22888) = 31
      ivals(22889) = 31
      ivals(22890) = 31
      ivals(22891) = 31
      ivals(22892) = 31
      ivals(22941) = 31
      ivals(22941) = 31
      ivals(22945) = 31
      ivals(22946) = 31
      ivals(22947) = 31
      ivals(22948) = 31
      ivals(22949) = 31
      ivals(22950) = 31
      ivals(22951) = 31
      ivals(22952) = 31
      ivals(22956) = 94
    else if ( irec==5393 ) then
      ivals(22926) = 31
      ivals(22935) = 31
      ivals(22909) = 56
      ivals(22909) = 56
      ivals(22910) = 56
      ivals(22911) = 56
      ivals(22912) = 56
      ivals(22913) = 56
      ivals(22884) = 31
      ivals(22885) = 31
      ivals(22886) = 31
      ivals(22887) = 31
      ivals(22888) = 31
      ivals(22889) = 31
      ivals(22903) = 56
      ivals(22904) = 56
      ivals(22905) = 56
      ivals(22908) = 56
      ivals(22938) = 31
      ivals(22939) = 31
      ivals(22940) = 31
      ivals(22947) = 31
      ivals(22948) = 31
      ivals(22949) = 31
      ivals(22944) = 31
      ivals(22945) = 31
      ivals(22946) = 31
      ivals(22954) = 31
      ivals(22955) = 31
    else if ( irec==5392 ) then
      ivals(22936) = 31
      ivals(22911) = 56
      ivals(22914) = 56
      ivals(22917) = 56
      ivals(22918) = 56
      ivals(22919) = 56
      ivals(22920) = 56
      ivals(22921) = 56
      ivals(22922) = 31
      ivals(22923) = 31
      ivals(22924) = 31
      ivals(22925) = 31
      ivals(22906) = 56
      ivals(22907) = 56
    else if ( irec==5391 ) then
      ivals(22936) = 31
      ivals(22910) = 56
      ivals(22911) = 56
      ivals(22913) = 56
      ivals(22914) = 56
      ivals(22916) = 56
      ivals(22917) = 56
      ivals(22918) = 56
      ivals(22919) = 56
      ivals(22920) = 56
      ivals(22921) = 56
      ivals(22922) = 31
      ivals(22923) = 31
      ivals(22924) = 31
      ivals(22906) = 56
    else if ( irec==5390 ) then
      ivals(22912) = 56
      ivals(22915) = 56
      ivals(22916) = 56
      ivals(22917) = 56
      ivals(22918) = 56
      ivals(22919) = 56
      ivals(22920) = 56
      ivals(22921) = 56
    else if ( irec==5353 ) then
      ivals(22931) = 56
    else if ( irec==5354 ) then
      ivals(22933) = 31
      ivals(22940) = 56
      ivals(22941) = 56
      ivals(22942) = 56
      ivals(22943) = 56
    else if ( irec==5355 ) then
      ivals(22942) = 56
      ivals(22943) = 56
      ivals(22944) = 56
      ivals(22945) = 56
    else if ( irec==5357 ) then
      ivals(22953) = 56
    else if ( irec==5358 ) then
      ivals(22954) = 31
      ivals(22955) = 31
    else if ( irec==5364 ) then
      ivals(22954) = 31
      ivals(22955) = 31
    else if ( irec==5365 ) then
      ivals(22954) = 31
      ivals(22955) = 31
      ivals(22956) = 31
    else if ( irec==5366 ) then
      ivals(22956) = 31
    else if ( irec==5370 ) then
      ivals(22955) = 31
      ivals(22956) = 31
    else if ( irec==5371 ) then
      ivals(22957) = 31
    else if ( irec==5372 ) then
      ivals(22957) = 31
      ivals(22958) = 31
    else if ( irec==5378 ) then
      ivals(22963) = 31
    else if ( irec==5379 ) then
      ivals(22965) = 31
    end if

! correzioni di Marco Verdecchia - Nord Ovest - Bacini del Po e Ticino
    if ( irec==5295 ) then
      ivals(22611) = 60
    else if ( irec==5296 ) then
      do ii = 22611 , 22613
        ivals(ii) = 60
      end do
    else if ( irec==5297 ) then
      do ii = 22610 , 22613
        ivals(ii) = 30
      end do
    else if ( irec==5298 ) then
      do ii = 22610 , 22612
        ivals(ii) = 60
      end do
    else if ( irec==5299 ) then
      do ii = 22610 , 22612
        ivals(ii) = 60
      end do
    else if ( irec==5300 ) then
      do ii = 22610 , 22612
        ivals(ii) = 60
      end do
    else if ( irec==5301 ) then
      do ii = 22610 , 22612
        ivals(ii) = 60
      end do
    else if ( irec==5302 ) then
      do ii = 22611 , 22613
        ivals(ii) = 60
      end do
    else if ( irec==5303 ) then
      do ii = 22609 , 22614
        ivals(ii) = 25
      end do
    else if ( irec==5304 ) then
      do ii = 22610 , 22613
        ivals(ii) = 60
      end do
    else if ( irec==5305 ) then
      do ii = 22611 , 22615
        ivals(ii) = 31
      end do
    else if ( irec==5306 ) then
      do ii = 22613 , 22615
        ivals(ii) = 31
      end do
    else if ( irec==5307 ) then
      ivals(22615) = 30
    else if ( irec==5319 ) then
      do ii = 22643 , 22644
        ivals(ii) = 25
      end do
    else if ( irec==5320 ) then
      do ii = 22645 , 22646
        ivals(ii) = 56
      end do
    else if ( irec==5321 ) then
      ivals(22644) = 25
    else if ( irec==5322 ) then
      do ii = 22644 , 22645
        ivals(ii) = 60
      end do
    else if ( irec==5323 ) then
      do ii = 22641 , 22646
        ivals(ii) = 60
      end do
    else if ( irec==5324 ) then
      do ii = 22643 , 22647
        ivals(ii) = 93
      end do
    else if ( irec==5325 ) then
      do ii = 22643 , 22644
        ivals(ii) = 93
      end do
    else if ( irec==5326 ) then
      do ii = 22645 , 22646
        ivals(ii) = 93
      end do
    else if ( irec==5327 ) then
      ivals(22647) = 93
    else if ( irec==5328 ) then
      do ii = 22648 , 22649
        ivals(ii) = 30
      end do
    else if ( irec==5329 ) then
      ivals(22646) = 93
    else if ( irec==5330 ) then
      do ii = 22646 , 22648
        ivals(ii) = 56
      end do
    else if ( irec==5331 ) then
      ivals(22648) = 93
    else if ( irec==5332 ) then
      ivals(22648) = 93
    else if ( irec==5333 ) then
      do ii = 22647 , 22649
        ivals(ii) = 31
      end do
    else if ( irec==5334 ) then
      do ii = 22648 , 22650
        ivals(ii) = 31
      end do
    else if ( irec==5335 ) then
      do ii = 22649 , 22651
        ivals(ii) = 31
      end do
    else if ( irec==5336 ) then
      do ii = 22647 , 22651
        ivals(ii) = 31
      end do
    else if ( irec==5338 ) then
      do ii = 22649 , 22651
        ivals(ii) = 31
      end do
    else if ( irec==5339 ) then
      do ii = 22651 , 22655
        ivals(ii) = 93
      end do
    else if ( irec==5340 ) then
      do ii = 22654 , 22656
        ivals(ii) = 93
      end do
    else if ( irec==5337 ) then
      do ii = 22648 , 22651
        ivals(ii) = 31
      end do
    else if ( irec==5344 ) then
      do ii = 22658 , 22659
        ivals(ii) = 31
      end do
    else if ( irec==5345 ) then
      do ii = 22659 , 22661
        ivals(ii) = 31
      end do
    else if ( irec==5346 ) then
      do ii = 22660 , 22661
        ivals(ii) = 93
      end do
    else if ( irec==5347 ) then
      do ii = 22660 , 22662
        ivals(ii) = 93
      end do
    else if ( irec==5348 ) then
      do ii = 22661 , 22664
        ivals(ii) = 93
      end do
    else if ( irec==5349 ) then
      do ii = 22661 , 22664
        ivals(ii) = 31
      end do
    else if ( irec==5350 ) then
      do ii = 22663 , 22665
        ivals(ii) = 31
      end do
    else if ( irec==5351 ) then
      do ii = 22663 , 22665
        ivals(ii) = 31
      end do
    else if ( irec==5352 ) then
      do ii = 22663 , 22665
        ivals(ii) = 31
      end do
    else if ( irec==5353 ) then
      do ii = 22663 , 22667
        ivals(ii) = 31
      end do
    else if ( irec==5354 ) then
      do ii = 22664 , 22666
        ivals(ii) = 31
      end do
    else if ( irec==5355 ) then
      do ii = 22664 , 22666
        ivals(ii) = 31
      end do
    else if ( irec==5356 ) then
      ivals(22666) = 31
    else if ( irec==5357 ) then
      do ii = 22666 , 22668
        ivals(ii) = 31
      end do
    else if ( irec==5358 ) then
      do ii = 22668 , 22670
        ivals(ii) = 31
      end do
    else if ( irec==5359 ) then
      do ii = 22669 , 22671
        ivals(ii) = 5
      end do
    else if ( irec==5360 ) then
      do ii = 22671 , 22674
        ivals(ii) = 31
      end do
    else if ( irec==5361 ) then
      do ii = 22674 , 22676
        ivals(ii) = 31
      end do
    else if ( irec==5362 ) then
      do ii = 22674 , 22676
        ivals(ii) = 31
      end do
    else if ( irec==5363 ) then
      do ii = 22676 , 22678
        ivals(ii) = 31
      end do
    else if ( irec==5364 ) then
      do ii = 22676 , 22678
        ivals(ii) = 31
      end do
      ivals(22753) = 30
    else if ( irec==5365 ) then
      do ii = 22677 , 22679
        ivals(ii) = 31
      end do
      do ii = 22754 , 22755
        ivals(ii) = 56
      end do
    else if ( irec==5366 ) then
      do ii = 22678 , 22680
        ivals(ii) = 31
      end do
    else if ( irec==5367 ) then
      do ii = 22680 , 22682
        ivals(ii) = 56
      end do
    else if ( irec==5368 ) then
      do ii = 22680 , 22682
        ivals(ii) = 5
      end do
    else if ( irec==5369 ) then
      do ii = 22681 , 22683
        ivals(ii) = 5
      end do
    else if ( irec==5370 ) then
      do ii = 22682 , 22686
        ivals(ii) = 31
      end do
    else if ( irec==5371 ) then
      do ii = 22684 , 22688
        ivals(ii) = 31
      end do
    else if ( irec==5373 ) then
      ivals(22689) = 31
    else if ( irec==5375 ) then
      do ii = 22692 , 22694
        ivals(ii) = 31
      end do
    else if ( irec==5376 ) then
      do ii = 22694 , 22696
        ivals(ii) = 31
      end do
    else if ( irec==5377 ) then
      do ii = 22587 , 22588
        ivals(ii) = 31
      end do
      do ii = 22700 , 22701
        ivals(ii) = 31
      end do
    else if ( irec==5378 ) then
      do ii = 22702 , 22703
        ivals(ii) = 93
      end do
      do ii = 22568 , 22591
        ivals(ii) = 31
      end do
    else if ( irec==5379 ) then
      do ii = 22705 , 22706
        ivals(ii) = 31
      end do
      do ii = 22570 , 22592
        ivals(ii) = 31
      end do
    else if ( irec==5380 ) then
      do ii = 22707 , 22709
        ivals(ii) = 31
      end do
      do ii = 22622 , 22627
        ivals(ii) = 31
      end do
    else if ( irec==5381 ) then
      do ii = 22707 , 22711
        ivals(ii) = 31
      end do
      do ii = 22703 , 22705
        ivals(ii) = 94
      end do
      do ii = 22624 , 22632
        ivals(ii) = 31
      end do
      do ii = 22610 , 22617
        ivals(ii) = 31
      end do
    else if ( irec==5382 ) then
      do ii = 22707 , 22711
        ivals(ii) = 31
      end do
      do ii = 22703 , 22705
        ivals(ii) = 57
      end do
      do ii = 22624 , 22632
        ivals(ii) = 31
      end do
      do ii = 22611 , 22617
        ivals(ii) = 31
      end do
    else if ( irec==5383 ) then
      do ii = 22706 , 22707
        ivals(ii) = 31
      end do
      do ii = 22709 , 22713
        ivals(ii) = 31
      end do
      do ii = 22702 , 22705
        ivals(ii) = 31
      end do
      do ii = 22630 , 22634
        ivals(ii) = 31
      end do
      do ii = 22612 , 22619
        ivals(ii) = 21
      end do
      do ii = 22790 , 22795
        ivals(ii) = 31
      end do
    else if ( irec==5384 ) then
      do ii = 22714 , 22720
        ivals(ii) = 31
      end do
      do ii = 22700 , 22704
        ivals(ii) = 57
      end do
      do ii = 22630 , 22635
        ivals(ii) = 31
      end do
      do ii = 22757 , 22759
        ivals(ii) = 31
      end do
      do ii = 22788 , 22799
        ivals(ii) = 31
      end do
    else if ( irec==5385 ) then
      do ii = 22716 , 22722
        ivals(ii) = 31
      end do
      do ii = 22701 , 22703
        ivals(ii) = 57
      end do
      do ii = 22632 , 22636
        ivals(ii) = 31
      end do
      do ii = 22756 , 22759
        ivals(ii) = 31
      end do
      do ii = 22789 , 22794
        ivals(ii) = 31
      end do
      ivals(22796) = 31
    else if ( irec==5386 ) then
      do ii = 22723 , 22727
        ivals(ii) = 31
      end do
      do ii = 22690 , 22697
        ivals(ii) = 31
      end do
      do ii = 22702 , 22703
        ivals(ii) = 57
      end do
      do ii = 22634 , 22637
        ivals(ii) = 31
      end do
      do ii = 22741 , 22759
        ivals(ii) = 31
      end do
      do ii = 22772 , 22795
        ivals(ii) = 31
      end do
      ivals(22798) = 37
    else if ( irec==5387 ) then
      do ii = 22725 , 22733
        ivals(ii) = 31
      end do
      do ii = 22689 , 22694
        ivals(ii) = 19
      end do
      do ii = 22697 , 22699
        ivals(ii) = 31
      end do
      do ii = 22636 , 22638
        ivals(ii) = 31
      end do
      do ii = 22734 , 22745
        ivals(ii) = 31
      end do
      do ii = 22749 , 22757
        ivals(ii) = 31
      end do
      do ii = 22759 , 22761
        ivals(ii) = 94
      end do
      do ii = 22774 , 22797
        ivals(ii) = 31
      end do
    else if ( irec==5388 ) then
      do ii = 22673 , 22679
        ivals(ii) = 31
      end do
      do ii = 22636 , 22639
        ivals(ii) = 31
      end do
      do ii = 22680 , 22691
        ivals(ii) = 57
      end do
    else if ( irec==5389 ) then
      do ii = 22675 , 22678
        ivals(ii) = 57
      end do
      do ii = 22681 , 22692
        ivals(ii) = 57
      end do
      do ii = 22638 , 22640
        ivals(ii) = 31
      end do
      do ii = 22763 , 22765
        ivals(ii) = 31
      end do
      do ii = 22746 , 22750
        ivals(ii) = 31
      end do
      do ii = 22772 , 22793
        ivals(ii) = 31
      end do
    else if ( irec==5390 ) then
      do ii = 22675 , 22679
        ivals(ii) = 57
      end do
      do ii = 22680 , 22680
        ivals(ii) = 31
      end do
      do ii = 22685 , 22686
        ivals(ii) = 31
      end do
      ivals(22638) = 31
      do ii = 22764 , 22766
        ivals(ii) = 31
      end do
      do ii = 22770 , 22771
        ivals(ii) = 31
      end do
      ivals(22748) = 31
      do ii = 22773 , 22790
        ivals(ii) = 31
      end do
    else if ( irec==5391 ) then
      do ii = 22669 , 22673
        ivals(ii) = 31
      end do
      do ii = 22638 , 22640
        ivals(ii) = 31
      end do
      ivals(22677) = 93
      ivals(22768) = 31
      ivals(22772) = 31
      ivals(22773) = 31
      do ii = 22782 , 22783
        ivals(ii) = 31
      end do
    else if ( irec==5392 ) then
      do ii = 22670 , 22673
        ivals(ii) = 31
      end do
      do ii = 22638 , 22640
        ivals(ii) = 31
      end do
      ivals(22677) = 93
    else if ( irec==5393 ) then
      do ii = 22667 , 22675
        ivals(ii) = 31
      end do
      do ii = 22639 , 22640
        ivals(ii) = 31
      end do
    else if ( irec==5394 ) then
      do ii = 22663 , 22665
        ivals(ii) = 31
      end do
      do ii = 22639 , 22641
        ivals(ii) = 31
      end do
      do ii = 22669 , 22670
        ivals(ii) = 93
      end do
    else if ( irec==5395 ) then
      do ii = 22640 , 22644
        ivals(ii) = 31
      end do
      do ii = 22660 , 22666
        ivals(ii) = 57
      end do
      ivals(22648) = 31
      ivals(22639) = 31
    else if ( irec==5396 ) then
      do ii = 22641 , 22664
        ivals(ii) = 31
      end do
    else if ( irec==5397 ) then
      do ii = 22642 , 22662
        ivals(ii) = 31
      end do
    else if ( irec==5398 ) then
      do ii = 22649 , 22660
        ivals(ii) = 31
      end do
    else if ( irec==5399 ) then
      do ii = 22653 , 22658
        ivals(ii) = 57
      end do
    else if ( irec==5400 ) then
      do ii = 22651 , 22656
        ivals(ii) = 57
      end do
    else if ( irec==5401 ) then
      do ii = 22649 , 22653
        ivals(ii) = 93
      end do
    else if ( irec==5402 ) then
      do ii = 22650 , 22653
        ivals(ii) = 93
      end do
    else if ( irec==5403 ) then
      do ii = 22648 , 22652
        ivals(ii) = 57
      end do
    end if

! correzioni di Marco Verdecchia - Nord Ovest - Bacino del Po
    if ( irec==5383 ) then
      do ii = 22800 , 22802
        ivals(ii) = 31
      end do
    else if ( irec==5384 ) then
      do ii = 22802 , 22804
        ivals(ii) = 31
      end do
    else if ( irec==5385 ) then
      do ii = 22800 , 22804
        ivals(ii) = 93
      end do
    else if ( irec==5386 ) then
      do ii = 22804 , 22806
        ivals(ii) = 31
      end do
    else if ( irec==5387 ) then
      do ii = 22805 , 22807
        ivals(ii) = 31
      end do
    else if ( irec==5388 ) then
      do ii = 22806 , 22808
        ivals(ii) = 31
      end do
    else if ( irec==5389 ) then
      do ii = 22807 , 22808
        ivals(ii) = 31
      end do
    else if ( irec==5390 ) then
      do ii = 22808 , 22810
        ivals(ii) = 31
      end do
    else if ( irec==5391 ) then
      do ii = 22808 , 22810
        ivals(ii) = 31
      end do
    else if ( irec==5392 ) then
      do ii = 22808 , 22810
        ivals(ii) = 31
      end do
    else if ( irec==5393 ) then
      do ii = 22808 , 22810
        ivals(ii) = 31
      end do
    else if ( irec==5394 ) then
      do ii = 22822 , 22825
        ivals(ii) = 31
      end do
      ivals(22810) = 31
    else if ( irec==5395 ) then
      do ii = 22810 , 22814
        ivals(ii) = 56
      end do
      do ii = 22819 , 22820
        ivals(ii) = 56
      end do
      do ii = 22824 , 22830
        ivals(ii) = 31
      end do
    else if ( irec==5396 ) then
      do ii = 22812 , 22817
        ivals(ii) = 56
      end do
      ivals(22821) = 93
      do ii = 22828 , 22833
        ivals(ii) = 31
      end do
    else if ( irec==5397 ) then
      do ii = 22814 , 22818
        ivals(ii) = 56
      end do
      do ii = 22831 , 22835
        ivals(ii) = 31
      end do
    else if ( irec==5398 ) then
      do ii = 22835 , 22837
        ivals(ii) = 31
      end do
    else if ( irec==5399 ) then
      do ii = 22835 , 22838
        ivals(ii) = 31
      end do
    else if ( irec==5400 ) then
      do ii = 22850 , 22852
        ivals(ii) = 31
      end do
      do ii = 22838 , 22841
        ivals(ii) = 31
      end do
    else if ( irec==5401 ) then
      do ii = 22839 , 22852
        ivals(ii) = 31
      end do
    else if ( irec==5402 ) then
      do ii = 22841 , 22854
        ivals(ii) = 31
      end do
    else if ( irec==5403 ) then
      do ii = 22843 , 22855
        ivals(ii) = 56
      end do
    else if ( irec==5403 ) then
      do ii = 22846 , 22856
        ivals(ii) = 30
      end do
    else if ( irec==5404 ) then
      do ii = 22846 , 22856
        ivals(ii) = 30
      end do
    else if ( irec==5405 ) then
      do ii = 22854 , 22859
        ivals(ii) = 30
      end do
    else if ( irec==5407 ) then
      do ii = 22856 , 22861
        ivals(ii) = 30
      end do
    else if ( irec==5408 ) then
      do ii = 22859 , 22862
        ivals(ii) = 31
      end do
    else if ( irec==5409 ) then
      ivals(22864) = 31
    end if
  end subroutine correctusgs

  subroutine d2icellcycle(luse,fl,work,ix,jx)
    implicit none
    integer :: ix , jx
    integer , dimension(ix,jx) :: fl , luse , work
    intent (in) fl
    intent (inout) luse , work
    integer , dimension(9) :: di , dj
    integer :: i , imax , ival , j , k , l , npoint
    integer , dimension(8) :: wk1 , wk2 , wk3
    data di/0 , 1 , 0 , -1 , 1 , 1 , -1 , -1 , 0/
    data dj/1 , 0 , -1 , 0 , 1 , -1 , -1 , 1 , 0/
    do i = 1 , ix
      do j = 1 , jx
        work(i,j) = luse(i,j)
        if ( fl(i,j)>0 ) then
          if ( inside(i,2,ix-2) .and. inside(j,2,jx-1) ) then
            npoint = 8
            do l = 1 , npoint
              wk1(l) = luse(i+di(l),j+dj(l))
              wk3(l) = fl(i+di(l),j+dj(l))
            end do
          else
            npoint = 0
            do l = 1 , 8
              if ( inside(i+di(l),1,ix) .and. inside(j+dj(l),1,jx) ) then
                npoint = npoint + 1
                wk1(npoint) = luse(i+di(l),j+dj(l))
                wk3(npoint) = fl(i+di(l),j+dj(l))
              end if
            end do
          end if
          imax = 0
          do k = 1 , npoint
            if ( wk3(k)==0 ) then
              do l = 1 , npoint
                if ( wk1(k)==wk1(l) .and. wk3(l)==0 ) wk2(k) = wk2(k) + 1
              end do
              if ( imax<wk2(k) ) then
                imax = wk2(k)
                ival = wk1(k)
              end if
            end if
          end do
          if ( imax>0 ) work(i,j) = ival
        end if
      end do
    end do
    do i = 1 , ix
      do j = 1 , jx
        luse(i,j) = work(i,j)
      end do
    end do
  end subroutine d2icellcycle

  logical function inside(n,n1,n2)
    implicit none
    integer :: n , n1 , n2
    intent (in) n , n1 , n2
    inside = .true.
    if ( n<n1 .or. n>n2 ) inside = .false.
  end function inside

  subroutine buildhyddirmap
    character(len=256) :: fdemd,outf
!    real(8), dimension(nlond,nlatd) :: demdd
    double precision , allocatable , dimension (:) :: lond ,latd
    integer :: ii, jj, i, j, nlond , nlatd, ifound
    fdemd = 'museo/FLOWDIR/HYDROSHEDS_001/w001001_001degree_mulc-1_missval.nc'
    outf = 'EU_CORDEXWORLD_DEM_0.06degree_fdm.nc'
    call gethydrodata_int(fdemd, nlond, nlatd, demd, lond, latd, ifound)
    luse = 1
    do jj = 1 , nlat
      do ii = 1 , nlon
        if ( demd(ii,jj)>=50000 ) then
          luse(ii,jj) = mare
          demd(ii,jj) = -500000000
        end if
        if ( demd(ii,jj)<=-150000000 ) then
          luse(ii,jj) = mare
          demd(ii,jj) = -150000000
        end if
      end do
    end do
    call buildflowdirmap(.False.)
    call areamatrix
!    call exit(0)
  end subroutine buildhyddirmap

  subroutine buildflowdirmap(flagd)
    implicit none
    character(len=2) , dimension(8) :: col
    real :: deltah , diff , hmin , xmax , xmin
    integer :: i , ii , j , nbin , ncyc , nzero , nzero2
    integer , dimension(8) :: icol
    real , dimension(0:8) :: rk
    real , dimension(nlon,nlat) :: plot , work , s
    integer , dimension(nlon,nlat) :: icl
    real , dimension(100) :: vv
    character(len=80) :: stitle
    character(len=80) :: title
    character(len=60) :: cfile
    logical :: plotgeop , plotlake, flagd
    real , dimension(nlon,nlat) :: savedem
    data icol/5 , 7 , 9 , 11 , 3 , 17 , 4 , 2/
    data col/'NW' , 'N' , 'NE' , 'E' , 'SE' , 'S' , 'SW' , 'W'/
    data plotlake/.false./
    data plotgeop/.true./
    write (6,'(/12x,a)') 'Building Flow Direction Map.'
    rk(0) = 0.0
    nzero = 0
    do i = 1 , nlon
      fmap(i,1) = -1
      fmap(i,nlat) = -1
      noflow(i,1) = 0
      noflow(i,nlat) = 0
    end do
    do j = 1 , nlat
      fmap(1,j) = -1
      fmap(nlon,j) = -1
      noflow(1,j) = 0
      noflow(nlon,j) = 0
    end do
!---------------------------------------------------------- Add feb 2005 ------
    do i = 2 , nlon - 1
      do j = 2 , nlat - 1
        if ( luse(i,j)/=mare ) then
          noflow(i,j) = 8
        else
          noflow(i,j) = 0
        end if
      end do
    end do
    if (flagd) then
    savedem = dem
    else
    savedem = demd
    end if

    if (flagd) then
      ncyc = ncyc1
      call writeint(6,15,'DEM Smoothing by CA algorithm:',ncyc,'cycles')
      do ii = 1 , ncyc
        call d2cellcycle(dem,noflow,plot,nlon,nlat,0.005)
      end do
      ncyc = ncyc2
      call writeint(6,15,'DEM Smoothing by filling algorithm:',ncyc,'cycles')
      do ii = 1 , ncyc
        call demholefilling
      end do
    end if

    call estabilishfowdir(flagd)
    nzero = noflowcounter(noflow)
    call writeint(6,15,'No-flow points are',nzero,' ')
    ii = 0
    if (flagd) then
    do while ( ii<ncyc3 .and. nzero>0 )
      ii = ii + 1
      call d2cellcycle(dem,noflow,plot,nlon,nlat,0.035)             !Scommentami
      call estabilishfowdir(flagd)
      nzero = noflowcounter(noflow)
      if ( mod(ii,5)==0 ) then
        write (title,'(i10,a)') nzero , ' No-flow points after'
        call writeint(6,15,title,ii,'CA Algorithm Correction')
      end if
    end do
!------------------------------------------------------------------------------
    if ( nzero>0 ) then
      write (6,'(15x,a)') 'Correcting flow-dir with direction-algorithm.'
      call dircorrflow(icl)
      call dircorrflow(icl)
      nzero2 = noflowcounter(icl)
      call writeint(6,15,'No-flow points are now ',nzero2,' ')
    end if
    if ( demf==1 ) then
      stitle = 'DEM Smooting Algorithm 1 (DSA1)'
    else
      dem = savedem
      call demsmoothing
      stitle = 'DEM Smooting Algorithm 2 (DSA2)'
    end if
    plot=1.0
    call flowcheck(plot,work)
    if ( mchym(13)==10 ) fmap = fmap1
    end if
    write (6,'(12x,a)') 'Done.'
  end subroutine buildflowdirmap


  subroutine areamatrix
    implicit none
    real :: d1 , d2
    integer :: i , j
    write (6,'(/12x,a)') 'Building Area and Total Drained Area fields.'
    do i = 1 , nlon - 1
      do j = 1 , nlat - 1
        d1 = 0.001*distance(lat(i,j),lon(i,j),lat(i+1,j),lon(i+1,j))
        d2 = 0.001*distance(lat(i,j),lon(i,j),lat(i,j+1),lon(i,j+1))
        area(i,j) = d1*d2
      end do
    end do
    do i = 1 , nlon
      area(i,nlat) = area(i,nlat-1)
    end do
    do j = 1 , nlat
      area(nlon,j) = area(nlon-1,j)
    end do
    wrk1 = area
    call flowcheck(wrk1,drai)
    write (6,'(12x,a)') 'Done.'
  end subroutine areamatrix

  subroutine buildacclivitymap
    implicit none
    real :: cfact , delta , thresh , vmax
    logical :: first
    integer :: i , i1 , i2 , j , j1 , j2 , k , n , ncor , nnofl
    character(len=80) :: title
    real , dimension(nlon,nlat) :: wk1
    integer , dimension(100) :: icol
    character(len=10) , dimension(100) :: col
    data first/.true./
    save first
    if ( first ) then
      write (6,'(/12x,a)') 'Building Incline Map.'
    else
      write (6,'(/12x,a)') 'Correcting Incline Map.'
    end if
    accl(:,:) = 0.0
    wk1(:,:) = 0.0
    ncor = 0
    nnofl = 0
    do j = 2 , nlat - 1
      do i = 2 , nlon - 1
        k = fmap(i,j)
        if ( luse(i,j)/=mare .and. k>0 .and. k<=8 ) then
          accl(i,j) = (dem(i,j)-dem(i+ir(k),j+jr(k)))  &
                      /distance(lat(i,j),lon(i,j),lat(i+ir(k),j+jr(k)), &
                      lon(i+ir(k),j+jr(k)))
          if ( accl(i,j)<=1.0E-05 ) then
            if ( mchym(13)/=10 ) then
              noflow(i,j) = 8
              ncor = ncor + 1
              accl(i,j) = 0.0
            else
              accl(i,j) = 1.2E-05
            end if
          else
            noflow(i,j) = 0
          end if
        else if ( luse(i,j)/=mare .and. luse(i,j)/=lago .and. k==0 ) then
          nnofl = nnofl + 1
          noflow(i,j) = 8
          accl(i,j) = 0.0
        else if ( luse(i,j)==fiume .or. luse(i,j)==lago ) then
          accl(i,j) = 0.0
          noflow(i,j) = 8
        else if ( luse(i,j)==mare ) then
          accl(i,j) = 0.0
          noflow(i,j) = 0
        else
          write (6,'(10x,a,4i8)') &
            ' Flux error inside buildacclivitymap.' , i , j , k , luse(i,j)
          nnofl = nnofl + 1
          noflow(i,j) = 8
          accl(i,j) = 0.0
        end if
        wk1(i,j) = fmap(i,j)
      end do
    end do
    do j = 1 , nlat
      noflow(1,j) = -1
      noflow(nlon,j) = -1
      accl(1,j) = 0.0
      accl(nlon,j) = 0.0
    end do
    do i = 1 , nlon
      noflow(i,1) = -1
      noflow(i,nlat) = -1
      accl(i,1) = 0.0
      accl(i,nlat) = 0.0
    end do
    write (6,'(13x,i5,a)') ncor , ' cells have acclivity lower than 1E-05.'
    write (6,'(13x,i5,a)') nnofl , ' points have no flow direction.'
    write (6,'(15x,a)') 'Correcting Incline field with CA algorithm.'

! The following CA cycles affect only the grid points for which
! acclivity has not been defined (noflow=8)
!    if ( mchym(13)/=10 ) then
      do i = 1 , 100
        call d2cellcycle(accl,noflow,wk,nlon,nlat,0.1)
      end do
!    end if

! This control has been eliminated, in few situation at least, it seems
! that increase the number of steps of integration is enough to
! solve numerica singularities along  the river path.
    ncor = 0
    thresh = 0.00030
    do i = 1 , nlon
      do j = 1 , nlat
        if ( accl(i,j)>thresh .and. drai(i,j)>200.0 ) then 
          ncor = ncor + 1
!          accl(i,j)=thresh
        end if
      end do
    end do
    write (6,'(13x,i5,a,f8.4)') &
      ncor,' cells had acclivity greater than ',thresh
    first = .false.
    write (6,'(12x,a)') 'Done.'
  end subroutine buildacclivitymap

  subroutine reconnectdem
    implicit none
    integer :: i , iter , j
    logical :: modify
    real , dimension(nlon,nlat) :: plot , work
    character(len=60) :: str1 , str2
    write (6,'(/12x,a)') 'Reconnecting severe DEM singularities.'
    modify = .false.
    do iter = 1 , angiocycle
      write (6,'(15x,a,i1)') 'Iteration number ' , iter
      call findthemouths(drai,100.,0)
      if ( selriver>nsec .or. selriver<=0 ) selriver = 1
      call basinpaint(dem,fmap,luse,work,nlon,nlat, &
                      isec(selriver),jsec(selriver),plot,0)
      do i = 2 , nlon - 1
        do j = 2 , nlat - 1
          if ( work(i,j)>10.0 .and. nint(plot(i,j))==-5 ) then
            write (str1,'(i3,a,i3)') i , '-' , j
            call nospace(str1)
            write (str2,'(a,f10.1,a)') 'CHyM Angioplasty to cell: '//   &
                     trim(str1)//'(' , work(i,j) ,')'
            write (6,'(15x,a)') str2
            call angioplasty(dem,fmap,luse,plot,wk,nlon,nlat,i,j)
            modify = .true.
          end if
        end do
      end do
    end do
    write (6,'(12x,a)') 'Done.'
    if ( modify ) then
      call areamatrix
      call buildacclivitymap
    end if
  end subroutine reconnectdem

  subroutine angioplasty(dem,fmap,luse,wk1,wk3,nlon,nlat,i,j)
    implicit none
    integer :: i , j , nlat , nlon
    real , dimension(nlon,nlat) :: dem , wk1 , wk3
    integer , dimension(nlon,nlat) :: fmap , luse
    intent (in) i , j , luse
    intent (inout) dem , wk1
    real :: h1 , h2 , slope , xmax , xmin
    integer :: i1 , i2 , ibet , idir , idist , ii , iskip , istep , j1 , j2 , &
               jbet , jdir , jj , jskip , jstep , mindist , ncal , np , nstep
    logical :: plot
    character(len=60) :: str1 , str2
    data plot/.false./
    data ncal/0/

    ncal = ncal + 1
    if ( ncal==10 .and. plot ) then
      plot = .false.
      write (6,'(12x,a)') 'Too many calls to angioplasty. No plot produced.'
    end if
    np = angionp    ! Seeing several "Cannot solve singularity" increase this number
    i1 = i - np
    if ( i1<1 ) i1 = 1
    j1 = j - np
    if ( j1<1 ) j1 = 1
    i2 = i + np
    if ( i2>nlon ) i2 = nlon
    j2 = j + np
    if ( j2>nlat ) j2 = nlat
    ibet = 0
    jbet = 0
    mindist = 1000
    do ii = i1 , i2
      do jj = j1 , j2
        idist = iabs(i-ii) + iabs(j-jj)
        if ( dem(ii,jj)<dem(i,j) .and. &
             nint(wk1(ii,jj))/=-5 .and. idist<mindist ) then
          ibet = ii
          jbet = jj
          mindist = iabs(i-ii) + iabs(j-jj)
        end if
      end do
    end do
    if ( ibet==0 ) then
      write (6,'(15x,a)') 'Now Applying less restrictive algorithm'
      do ii = i1 , i2
        do jj = j1 , j2
          idist = iabs(i-ii) + iabs(j-jj)
          if ( nint(wk1(ii,jj))/=-5 .and. idist<mindist ) then
            ibet = ii
            jbet = jj
            mindist = iabs(i-ii) + iabs(j-jj)
          end if
        end do
      end do
      if ( ibet/=0 ) dem(i,j) = dem(ibet,jbet) + mindist
    end if
    if ( ibet==0 ) then
      write (6,'(12x,a)') 'Cannot solve singularity.'
      return
    end if
    if ( ibet<i ) then
      istep = -1
      idir = 8
      iskip = 1
    else
      istep = 1
      idir = 4
      iskip = -1
    end if
    if ( jbet<j ) then
      jstep = -1
      jdir = 6
      jskip = 1
    else
      jstep = 1
      jdir = 2
      jskip = -1
    end if
    nstep = 0
    h1 = dem(i,j)
    h2 = dem(ibet,jbet)
    slope = (dem(i,j)-dem(ibet,jbet))/mindist
    if ( ibet<i ) then
      fmap(i,j) = 8
    else
      fmap(i,j) = 4
    end if
    wk1(i,j) = dem(i,j)
    do ii = i + istep , ibet + iskip , istep
      nstep = nstep + 1
      fmap(ii,j) = idir
      dem(ii,j) = h1 - nstep*slope
      wk1(ii,j) = dem(ii,j)
    end do
    do jj = j , jbet + jskip , jstep
      nstep = nstep + 1
      fmap(ibet,jj) = jdir
      dem(ibet,jj) = h1 - nstep*slope
      wk1(ibet,jj) = dem(ibet,jj)
    end do
  end subroutine angioplasty

  subroutine findthemouths(w,cut,ilog)
    implicit none
    real :: cut
    integer :: ilog
    real , dimension(nlon,nlat) :: w
    intent (in) cut , ilog , w
    real :: dist , xxlat , xxlon
    integer :: i , idir , j , k , lun
    character(len=25) :: name , test
    if ( ilog==1 ) then
      write (6,'(21x,a)') 'Rivers mouths finding...'
      write (6,'(8x,68(''-''))')
      write (6,'(10x,a)') '#    I   J     Lat      Lon   Area(Km2)   Name'
      write (6,'(8x,68(''-''))')
    end if
    nsec = 0
    do j = nlat - 1 , 2 , -1
      do i = 2 , nlon - 1
        idir = fmap(i,j)
        if ( idir>=1 .and. idir<=8 ) then
          if ( luse(i+ir(idir),j+jr(idir))==mare .and. w(i,j)>cut ) then
            call chymrivername(lat(i,j),lon(i,j),name)
            nsec = nsec + 1
            if ( nsec>mxsc ) call chymerror(1,0,0.0,'findthemouths')
            isec(nsec) = i
            jsec(nsec) = j
            if ( ilog==1 ) &
              write (6,'(9x,i2,2x,2i4,2f9.4,f11.1,3x,a25)') nsec ,  &
                i , j , lat(i,j) , lon(i,j) , w(i,j) , name
          end if
        end if
      end do
    end do
  end subroutine findthemouths

  subroutine demholefilling
    implicit none
    real , dimension(8) :: hc
    integer :: i , ihole , j , k
    do i = 2 , nlon - 1
      do j = 2 , nlat - 1
        if ( luse(i,j)/=mare ) then
          ihole = 0
          do k = 1 , 8
            hc(k) = dem(i+ir(k),j+jr(k))
            if ( dem(i,j)<dem(i+ir(k),j+jr(k)) ) ihole = ihole + 1
          end do
          if ( ihole==8 ) dem(i,j) = minval(hc(1:8)) + 1
        end if
      end do
    end do
  end subroutine demholefilling

  subroutine estabilishfowdir(flagd)
    implicit none
    integer :: i , j , k
    integer , dimension(1) :: iloc
    real , dimension(0:8) :: rk
    logical :: flagd
    rk(0) = 0.0
    do j = 2 , nlat - 1
      do i = 2 , nlon - 1
        if ( luse(i,j)/=mare ) then  ! land use chk
          do k = 1 , 8
            if (flagd) then
            if ( dem(i+ir(k),j+jr(k))<dem(i,j) ) then
              rk(k) = (dem(i,j)-dem(i+ir(k),j+jr(k))) / &
                  distance(lat(i,j),lon(i,j), &
                           lat(i+ir(k),j+jr(k)),lon(i+ir(k),j+jr(k)))
            else
              rk(k) = -1.0
            end if
            else
            if ( demd(i+ir(k),j+jr(k))<demd(i,j) ) then
              rk(k) = (demd(i,j)-demd(i+ir(k),j+jr(k))) / &
                  distance(lat(i,j),lon(i,j), &
                           lat(i+ir(k),j+jr(k)),lon(i+ir(k),j+jr(k)))
            else
              rk(k) = -1.0
            end if
            end if
          end do
          iloc = maxloc(rk(0:8))
          fmap(i,j) = iloc(1) - 1
        end if
      end do
    end do
  end subroutine estabilishfowdir

  function noflowcounter(cam)
    implicit none
    integer , dimension(nlon,nlat) :: cam
    integer :: noflowcounter
    intent (out) cam
    integer :: i , j
    noflowcounter = 0
    do i = 2 , nlon - 1
      do j = 2 , nlat - 1
        if ( luse(i,j)==mare .or. luse(i,j)==lago ) then
          cam(i,j) = 0
        else if ( fmap(i,j)==0 ) then
          noflowcounter = noflowcounter + 1
          cam(i,j) = 8
        else
          cam(i,j) = 0
        end if
      end do
    end do
  end function noflowcounter

  subroutine dircorrflow(icl)
    implicit none
    integer , dimension(nlon,nlat) :: icl
    intent (inout) icl
    integer :: i , j , k
    do i = 2 , nlon - 1
      do j = 2 , nlat - 1
        if ( icl(i,j)==8 ) then
          do k = 1 , 4
            if ( fmap(i+ir(k),j+jr(k))==fmap(i+ir(k+4),j+jr(k+4)) ) then
              fmap(i,j) = fmap(i+ir(k),j+jr(k))
              icl(i,j) = 0
            end if
          end do
        end if
        if ( icl(i,j)==8 ) then
          do k = 1 , 8
            if ( fmap(i+ir(k),j+jr(k))==k ) then
              fmap(i,j) = k
              icl(i,j) = 0
            end if
          end do
        end if
      end do
    end do
  end subroutine dircorrflow

  subroutine writeint(unit,spaces,before,n,after)
    implicit none
    character(len=*) :: after , before
    integer :: n , spaces , unit
    intent (in) after , before , n , spaces , unit
    character(len=20) :: form
    character(len=132) :: str
    write (str,'(a,1x,i8,1x,a)') before , n , after
    call no2space(str)
    call noinspace(str)
    write (form,'(a,i3,a)') '(' , spaces , 'x,a)'
    write (unit,form) trim(str)
  end subroutine writeint

  subroutine demsmoothing
    implicit none
    integer :: icycle , i , idir , j , maxcycle , nsin
    nsin = 10
    icycle = 0
    maxcycle = 300
    do while ( icycle<=maxcycle .and. nsin>0 )
      nsin = 0
      icycle = icycle + 1
      do i = 2 , nlon - 1
        do j = 2 , nlat - 1
          if ( fmap(i,j)>=1 .and. fmap(i,j)<=8 ) then
            idir = fmap(i,j)
            if ( dem(i,j)<=dem(i+ir(idir),j+jr(idir)) ) then
              dem(i,j) = dem(i+ir(idir),j+jr(idir)) + 0.1
              nsin = nsin + 1
            end if
          end if
        end do
      end do
      if ( icycle==1 .or. mod(icycle,10)==0 ) &
        write (6,'(i20,a,i3)') &
          nsin , ' DEM cells modified by demsmoothing during cycle ' , icycle
    end do
    if ( mod(icycle,10)/=0 ) &
      write (6,'(i20,a,i3)') &
        nsin , ' DEM cells modified by demsmoothing during cycle ', icycle
  end subroutine demsmoothing

  subroutine flowcheck(water,wk1)
    implicit none
    real , dimension(nlon,nlat) :: water , wk1
    intent (inout) water , wk1
    integer :: i , idir , ifact , ii , isa , j , jsa , niter
    real :: xmax
    if (.not.allocated(wk)) allocate(wk(nlon,nlat))
    wk1(:,:) = 0.0
!    wk1 = (rchym(6)*0.001)**2 !Fabio
!    wk1 = 0.0
    write (6,'(18x,a)') 'Flow Check Module Called.'
    niter = (nlon+nlat)*2                        ! Deve dipendere da Ris.
    do ii = 1 , niter
      wk(:,:) = 0.0
      do j = 2 , nlat - 1
        do i = 2 , nlon - 1
          
          idir = fmap(i,j)
          if ( idir>0 ) then
            wk(i+ir(idir),j+jr(idir)) = wk(i+ir(idir),j+jr(idir)) + water(i,j)
            water(i,j) = 0
          end if
        end do
      end do
!   do i=2,nlon-1
!   do j=2,nlat-1
      do j = 1 , nlat
        do i = 1 , nlon
          water(i,j) = water(i,j) + wk(i,j)
          wk1(i,j) = wk1(i,j) + wk(i,j)
        end do
      end do
    end do
    wk(1:nlon,1:nlat) = wk1(1:nlon,1:nlat)
    xmax = 0.0
    write (6,'(18x,a)') 'Flow Check Module First part ended.'
    do j = 1 , nlat
      do i = 1 , nlon
        if ( wk1(i,j)>xmax ) then
          isa = i
          jsa = j
          xmax = wk1(i,j)
        end if
        if ( (rchym(3)-rchym(1))>5.0 .or. (rchym(4)-rchym(2))>5.0 ) then
          ifact = 600                                  ! Spiegare
        else
          ifact = 300
        end if
        if ( luse(i,j)/=mare .and. luse(i,j)/=lago ) then
          if ( wk1(i,j)>float((nlon+nlat)*ifact) ) then
            write (6,'(20x,a,2i4,f9.1)') 'Discarded: ' , i , j , wk1(i,j)
            noflow(i,j) = 8
          end if
        end if
      end do
    end do
    write (6,'(21x,a,f9.1)') 'Maximum value of "Rolling Stones" is: ' , xmax
!    write (6,'(21x,a,i3,a,i3,a,f9.4,a,f9.4,a)') 'on the point ' , isa , '-' , &
!           jsa , ' (' , lat(isa,jsa) , 'N-' , lon(isa,jsa) , 'E)'
  end subroutine flowcheck

  subroutine chymrivername(lat,lon,title)
    implicit none
    real :: lat,lon,xxlat,xxlon
    integer :: lun,i
    character(len=*) :: title
    character(len=50) :: river
    write (title,'(a)') 'Unknown'
    call getlun(lun)
    call openmuseofiles(lun,'river.basin',1)
    do i=1,10000
       read (lun,'(a50,2f9.5)',end=100) river,xxlat,xxlon
       if (distance(lat,lon,xxlat,xxlon).le.1000.0) then
          title=river
          exit
       endif
    enddo
100 close(lun)
    return
  end subroutine chymrivername

  subroutine riveronlanduse
    implicit none
    integer :: i , j
    do i = 1 , nlon
      do j = 1 , nlat
        if ( drai(i,j)>cpar(6) ) then 
           if (luse(i,j)/=mare .and. mchym(13)/=10) luse(i,j) = fiume
           if ( mchym(13)==10 ) luse(i,j) = fiume
        end if
        if ( luse(i,j)<=0 .or. luse(i,j)>lntypes ) then
          write (6,'(9x,a,3i3)') 'Undefined land use type: ' , i , j , luse(i,j)
          write (6,'(9x,a)') 'Exiting.'
          call exit(0)
        end if
      end do
    end do
  end subroutine riveronlanduse

  subroutine runoffspeed
    implicit none
    real :: alfamin , delta , xgamma , mann , tresh , vmax
    real :: enne , hrad
    integer :: i , idir , intstep , j , land
    !!!
    wrk2 = 0
    alfa = 0.0
    xgamma = 0.33
    delta = cpar(8)           ! Param. for land/channel flow
    tresh = cpar(6)
    alfamin = 0.1           ! Minimum value of surface runoff speed
    do i = 2 , nlon - 1
      do j = 2 , nlat - 1
        idir = fmap(i,j)
        land = luse(i,j)
        mann = manning(land)
        if ( idir>=1 .and. idir<=8 .and. land/=mare .and. land>0 ) then
          if ( land>100 .or. land<=0 ) then
            write (6,'(10x,a,i5)') 'Wrong value for landuse code: ' , land
            stop 'flux error inside runoffspeed'
          end if
          dx(i,j) = distance(lat(i,j),lon(i,j), &
                             lat(i+ir(idir),j+jr(idir)), &
                             lon(i+ir(idir),j+jr(idir)))
          if ( drai(i,j)>tresh ) then
            enne = mann/delta
          else
            enne = mann/(1+(delta-1)*(1+(drai(i,j)-tresh)/tresh))
          end if
          hrad = cpar(2) + cpar(3)*((drai(i,j)*1.E00)**xgamma)
          alfa(i,j) = ((hrad**0.6666*accl(i,j)**0.5)/(enne))
          if ( alfa(i,j)<alfamin ) alfa(i,j) = alfamin
        end if
      end do
    end do
 
    vmax = 0.0
    do i = 1 , nlon
      do j = 1 , nlat
        if ( alfa(i,j)>vmax .and. drai(i,j)>100.0 ) vmax = alfa(i,j)
      end do
    end do
    intstep = nint(vmax*6.0)
    if ( intstep<10 ) intstep = 10
! Following line is commented because algorithm does not work
! write(6,'(12x,a,i4)') 'Suggested number of steps/hour=',intstep
    alfamin = 0.1
    do i = 1 , nlon
      do j = 1 , nlat
        if ( luse(i,j)==mare ) then
          wrk2(i,j) = -10.0
        else if ( alfa(i,j)<alfamin ) then
          wrk2(i,j) = dx(i,j)/alfamin
        else
          wrk2(i,j) = dx(i,j)/alfa(i,j)
        end if
      end do
    end do
    call runofftime(wrk2,fmap,runt,nlon,nlat)
    wk = 1
    call rollingstones2(wk,fmap,wrk2,nlon,nlat)
    call rollingstones2(runt,fmap,wk,nlon,nlat)
    do i = 1 , nlon
      do j = 1 , nlat
        runt(i,j) = wk(i,j)/wrk2(i,j) - runt(i,j)
      end do
    end do
  end subroutine runoffspeed

  subroutine write_stat_NC
    implicit none
    integer savecen
    if ( mchym(13)/=10 )call plotriverbasin
    call createfile(trim(schym(11))//'.static_fields.nc', -1) 
!    call createfile(trim(schym(11))//'.static_fields.nc', &
!      mchym,rchym,schym,-1)
    call mvgetiflags(57,savecen)
    if (mchym(9).eq.2) then
      call runoffspeed
      call read_restart_dyn_NC
      call write_stat_NC0
    else
       call write_stat_NC0
       bwet=0.0 ; port=0.0 ; gh2o=0.0 ; snow=0.0 ; deepw=0.0
    endif
    call closefile
    call mvsetflags('Century',float(savecen))
    return
  end subroutine write_stat_NC

  subroutine write_stat_NC0
!    use chymdata , only : mchym,rchym,schym,nlon,nlat,chymwriterec,lavonc
!    use chymdata , only : dem,lat,lon,fmap,accl,luse,area,drai,runt,logun
!    use chymdata , only : allrivb
!    use mo_ncio
    implicit none
    integer lenfn
    lenfn=len_trim(schym(11))
    write(6,'(/,15x,a)') 'Initializing output file: '//schym(11)(1:lenfn)
    call write_statvar('dem',dem)
    call write_statvar('fdm',fmap)
    call write_statvar('acc',accl)
    call write_statvar('lus',luse)
    call write_statvar('aer',area)
    call write_statvar('dra',drai)
    call write_statvar('run',runt)
    call write_statvar('alf',alfa)
    call write_statvar('ctr',lavonc)
    call write_statvar('bas',allrivb)
    write(6,'(15x,a)') 'Done.'
  end subroutine write_stat_NC0

end module mod_crtstatic
