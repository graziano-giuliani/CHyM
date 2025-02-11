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
module mod_museo

  use mod_time
  use mod_strman
  use mod_internal
  use mod_vector
  use mod_statparams

  integer , parameter :: npt = 198912

  contains

  subroutine snowfrommodis(rec,year,modis,nlon,nlat,n1,n2,r1,r2,r3,r4)
    implicit none
    integer :: n1 , n2 , nlat , nlon , rec , year
    real :: r1 , r2 , r3 , r4
    real , dimension(nlon,nlat) :: modis
    intent (in) nlat , nlon
    intent (out) modis , r1 , r2 , r3 , r4
    intent (inout) n1 , n2
    !
    ! Local variables
    !
    integer :: i , j , lun , oldrec , oldyear
    !
    !*** End of declarations rewritten by SPAG
    !
    data oldyear , lun/ - 1 , -1/
    save lun,oldrec,oldyear
    if ( year/=oldyear ) then
     if ( lun>0 ) close (lun)
     call openmuseodb(lun,'modis',year)
     oldrec = 0
    end if
    call skipandreadgame(lun,rec,oldrec,0)
    modis = 255
    read (lun) n1 , n2 , r1 , r2 , r3 , r4 , ((modis(i,j),i=1,n1),j=1,n2)
  end subroutine snowfrommodis

  subroutine wrfrainv(ora,giorno,mese,anno,dom,pi,la,lo,n)
    implicit none
    integer :: anno , dom , giorno , mese , n , ora
    real , dimension(:) :: la , lo , pi
    intent (in) dom
    intent (out) la , lo , pi
    intent (inout) n
    integer :: i , last , lun , nlat , nlon , nrec , odom , orec , xflag ,  &
               year
    real , dimension(npt) , save :: lat , lon , rain
    data lun , year , odom/3* - 1/   ! Corresponding to 444x448 grid
    !     integer nlo,nla
    !     parameter (nlo=262,nla=184)     ! Dominio 1
    !     integer nlo,nla
    !     parameter (nlo=444,nla=448)     ! Dominio 2
    !     real xmat(nlo,nla),xlat(nlo,nla),xlon(nlo,nla)
    if ( year/=anno .or. dom/=odom ) then
      call openmuseodb(lun,'wrfrai',anno*10+dom)
      read (lun)
      last = 0
      do i = 1 , index1h(23,31,12,anno)
        read (lun) nlon , nlat , xflag
        if ( xflag==1 ) last = i
      end do
      rewind (lun)
      read (lun) year , nlon , nlat , (lat(i),i=1,nlon*nlat) ,               &
                (lon(i),i=1,nlon*nlat)
    ! read (lun) year,nlon,nlat,((xlat(i,j),i=1,nlon),j=1,nlat), &
    !           ((xlon(i,j),i=1,nlon),j=1,nlat)
    ! call plottastamatrice(xlon,xlon,xlat,nlo,nla,'WRF Latitudes')
      orec = 1
      year = anno
      odom = dom
    end if
    nrec = index1h(ora,giorno,mese,anno) + 1
    if ( nrec>last ) then
      n = -1
      return
    end if
    call skipandreadgame(lun,nrec,orec,0)
    read (lun) nlon , nlat , xflag , (rain(i),i=1,nlon*nlat)
    if ( xflag==1 ) then
      n = nlon*nlat
      do i = 1 , n
        pi(i) = rain(i)
        la(i) = lat(i)
        lo(i) = lon(i)
      end do
    else
      n = 0
    end if
  end subroutine wrfrainv

  subroutine demnasa(xlat,xlon,iunit,efile)
    implicit none
    character(len=*) :: efile
    integer :: iunit
    real :: xlat , xlon
    intent (in) iunit , xlat , xlon
    intent (out) efile
    character(60) :: dir
    character(4) :: ew , lats , lons , sn
    integer :: lat , lon
    character(80) :: lfile
    integer , intrinsic :: nint
    lat = nint(abs(xlat))
    lon = nint(abs(xlon))
    if ( lat>90 .or. lon>180 ) then
      write (6,'(10x,a)') 'Invalid Lat/Lon Range passed to demnasa routine'
      call exit(0)
    end if
    call museodir(dir)
    call integer2string(lat,2,lats)
    call integer2string(lon,3,lons)
    sn = 'N'
    if ( xlat<0.0 ) sn = 'S'
    ew = 'E'
    if ( xlon<0.0 ) ew = 'W'
    lfile = trim(dir)//'DEM/NASA/'//sn(1:1)//lats(1:2)        &
            //ew(1:1)//lons(1:3)//'.hgt'
    if ( iunit<=0 ) then
      efile = lfile
    else
      open (iunit,file=lfile,status='old',form='unformatted',access='stream')
    end if
  end subroutine demnasa
     
  integer function numberoflines(lun)
    implicit none
    integer :: lun
    intent (in) lun
    character(80) :: a
    integer :: nl
    nl = 0
    rewind (lun)
    do while ( nl<=0 .or. nl>0 )
      read (lun,'(a)',end=100) a
      nl = nl + 1
    end do
  100  rewind (lun)
    numberoflines = nl
  end function numberoflines
     
  subroutine skipandreadgame(lun,curr,last,flag)
    implicit none
    integer :: curr , flag , last , lun
    intent (in) curr , flag
    intent (inout) last
    integer :: i , i1 , i2 , ilog
    character(80) :: msg
    call mvgetiflags(1,ilog)
    if ( curr==last+1 ) then
      last = curr
      if ( ilog>0 ) write (6,'(12x,a,i3)') &
        'No records skipped from skipandreadgame on unit', lun
      return
    else if ( curr<=last .or. last<0 ) then
      rewind (lun)
      i1 = 1
      i2 = curr - 1
    else if ( curr>last ) then
      i1 = last + 1
      i2 = curr - 1
    end if
    if ( ilog>0 ) then
      write (msg,'(a,i6,a,i6,a,i3)') 'skipandreadgame skipped record from ' ,&
                                    i1 , ' to ' , i2 , ' from unit ' , lun
      call no2space(msg)
      write (6,'(12x,a)') trim(msg)
    end if
    if ( flag==0 ) then
      do i = i1 , i2
        read (lun,err=100,end=100)
      end do
    else
      do i = i1 , i2
        read (lun,*,err=100,end=100)
      end do
    end if
    last = curr
    return
  100  call mvliberror('skipandreadgame','Error while skipping',lun,     &
                          -9999.)
  end subroutine skipandreadgame
     
  subroutine openmuseodb(iunit,museodb,anno)
    implicit none
    integer :: anno , iunit
    character(*) :: museodb
    intent (in) anno
    character(80) :: file
    integer :: flag
    if ( museodb(1:3)=='mm5' .or. museodb(1:6)=='wrfrai' ) then
      write (file,'(a,i5,a)') trim(museodb) , anno , '.dat'
      flag = 0
    else
      write (file,'(a,i4,a)') trim(museodb) , anno , '.dat'
      flag = 0
    end if
    call openmuseofiles(iunit,file,flag)
  end subroutine openmuseodb
     
  subroutine openmuseoarchive(iunit,cfile,anno,dom)
    implicit none
    integer :: anno , dom , iunit
    character(len=*) :: cfile
    intent (in) anno , dom , iunit
    character(5) :: ext
    call museopath
    write (filename,'(i10,a4)') dom , '.dat'
    if ( cfile(1:3)=='mm5' ) then
      write (ext,'(i1,a4)',err=200) dom , '.dat'
    else
      ext = '.dat'
    end if
    write (filename,'(a,i4,a,i4,a)') trim(pathm) , anno ,        &
                              '/'//trim(cfile) , anno , ext
    open (iunit,file=filename,status='old',form='unformatted',err=100)
    return
  100  write (filename,'(a,i4,a)') trim(pathm)                &
                              //trim(cfile) , anno , ext
    open (iunit,file=filename,status='old',form='unformatted',err=200)
    return
  200  write (6,'(10x,a)') 'MVLib: error opening '//trim(filename)
    stop ' '
  end subroutine openmuseoarchive
     
  subroutine museoformat(iflag,form)
    implicit none
    character(len=*) :: form
    integer :: iflag
    intent (in) iflag
    intent (out) form
    if ( iflag==1 ) then
      form = '(a35,a2,i9,2f10.5,i5)'      ! comuni.tutti,
      if ( iflg(1)>=10 ) then
        write (6,'(a)') ' '
        write (6,'(a)') '          use comuni.tutti format as follow:'
        write (6,'(a)') '     character com*35,prov*2'
        write (6,'(a)') '     integer quota,flag ; real lat,lon'
        write (6,'(a)') '     read (unit,form) com,prov,quota,lat,lon,flag'
        write (6,'(a)') ' '
      end if
    else if ( iflag==2 ) then
      form = '(i3,1x,a45,a20,a2)'         ! comuni.idrografico riga 1
    else if ( iflag==3 ) then
      form = '(i10,2x,a30,a10,2f8.4)'     ! comuni.idrografico riga 2
    else if ( iflag==4 ) then
      form = '(i4,1x,a35,a25,a2,i7,1x,a25,a10,2f8.4,1x,3i1)'  ! comuni.iira
    else if ( iflag==5 ) then
      form = '(i12,1x,a50,1x,a2,1x,a50,2f8.3)'           ! dewetra
    else if ( iflag==6 ) then
      form = '(i3,2x,a40,a8,2f12.4)'                     ! Simone Fatichi
    end if
  end subroutine museoformat
     
  subroutine museodir(ofile)
    implicit none
    character(len=*) :: ofile
    intent (out) ofile
    call museopath
    ofile = pathm
  end subroutine museodir
     
  subroutine openmuseofiles(iunit,cfile,flag)
    implicit none
    character(len=*) :: cfile
    character(len=64) :: safe
    integer :: flag , iunit
    integer :: ii
    character(132) :: lfile
    call resetmvlibint
    lfile = cfile
    safe = cfile
    if ( safe(1:13)=='northitalyrai' .or. safe(1:13)=='northitalytem' ) then
      lfile = lfile(14:17)//'/'//lfile
    else if ( safe(1:3)=='nox' ) then
      lfile = lfile(4:7)//'/'//lfile
    else if ( safe(1:24)=='iiraportate1930-2002.dat' ) then
      lfile = '0000/'//lfile
    else if ( safe(1:4)=='iira' ) then
      lfile = lfile(5:8)//'/'//lfile
    else if ( safe(1:11)=='idrografico' ) then
      lfile = lfile(12:15)//'/'//lfile
    else if ( safe(1:8)=='idro-all' ) then
      lfile = lfile(9:12)//'/'//lfile
    else if ( safe(1:2)=='Po' ) then
      lfile = lfile(3:6)//'/'//lfile
    else if ( safe(1:5)=='arssa' ) then
      lfile = lfile(6:9)//'/'//lfile
    else if ( safe(1:5)=='lazio' ) then
      lfile = lfile(6:9)//'/'//lfile
    else if ( safe(1:5)=='radar' .and. safe(1:7)/='radar.l' ) then
      lfile = lfile(6:9)//'/'//lfile
    else if ( safe(1:13)=='pioggiaoraria' ) then
      if ( safe(1:17)=='pioggiaoraria-fil' ) then
        lfile = lfile(18:21)//'/'//lfile
      else
        lfile = lfile(14:17)//'/'//lfile
      end if
    else if ( safe(1:10)=='temporaria' ) then
      lfile = lfile(11:14)//'/'//lfile
    else if ( safe(1:6)=='scerni' ) then
      lfile = lfile(7:10)//'/'//lfile
    else if ( safe(1:15)=='MonteMidia.grid' ) then
      lfile = 'radar/'//lfile
    else if ( safe(1:10)=='MonteMidia' ) then
      lfile = lfile(11:14)//'/'//lfile
    else if ( safe(1:5)=='Micra' ) then
      lfile = lfile(6:9)//'/'//lfile
    else if ( safe(1:6)=='wrfrai' ) then
      lfile = lfile(7:10)//'/'//lfile
    else if ( safe(1:3)=='wrf' ) then                ! No longer used
      lfile = lfile(4:7)//'/'//lfile
    else if ( safe(1:11)=='dewrain-fil' ) then
      lfile = lfile(12:15)//'/'//lfile
    else if ( safe(1:7)=='dewrain' ) then
      lfile = lfile(8:11)//'/'//lfile
    else if ( safe(1:7)=='dewtemp' ) then
      lfile = lfile(8:11)//'/'//lfile
    else if ( safe(1:7)=='dewidro' ) then
      lfile = lfile(8:11)//'/'//lfile
    else if ( safe(1:7)=='dewdisc' ) then
      lfile = lfile(8:11)//'/'//lfile
    else if ( safe(1:9)=='worldtemp' .and. safe(1:13)/='worldtemp.dat' ) then
      lfile = lfile(10:13)//'/'//lfile
    else if ( safe(1:9)=='worldrain' ) then
      lfile = lfile(10:13)//'/'//lfile
    else if ( safe(1:7)=='comuni.' ) then
      lfile = 'comuni/'//lfile
    else if ( safe(5:19)=='0100_eraint.dat' ) then
      lfile = 'era-interim/'//lfile
    else if ( safe(1:10)=='acqwarhone' ) then
      lfile = lfile(11:14)//'/'//lfile
    else if ( safe(1:7)=='acqwapo' ) then
      lfile = lfile(8:11)//'/'//lfile
    else if ( safe(1:8)=='acqwachu' ) then
      lfile = lfile(9:12)//'/'//lfile
    else if ( safe(1:5)=='modis' ) then
      lfile = lfile(6:9)//'/'//lfile
    else if ( safe(1:14)=='mm5terrain.dat' .or. safe(1:7)=='mm5.dat' ) then
      lfile = 'mm5/'//lfile
    else if ( (safe(1:3)=='mm5' .or. safe(1:4)=='/mm5') ) then
      ii = index(safe,'.dat')
      lfile = safe(ii-5:ii-2)//'/'//safe
    else if ( safe(1:9)=='BAtot.dat' ) then
      lfile = 'EVAP/'//lfile
    else if ( safe(1:12)=='oge2_0.marco' ) then
      lfile = 'LUSE/'//lfile
    else if ( safe(1:13)=='worldtemp.dat' ) then
      lfile = 'TEMP/'//lfile
    else if ( safe(1:10)=='demita.bin' ) then
      lfile = 'DEM/'//lfile
    else if ( safe(1:6)=='river.' ) then
      lfile = 'RIVER/'//lfile
    else if ( safe(1:16)=='alarmabruzzo.cnt' .or. &
              safe(1:12)=='demlakes.bin' ) then
      lfile = 'dem/'//lfile
    end if
    call museofile(lfile,lfile)
    if ( iunit<0 ) then
      call getlun(iunit)
    else
      close (iunit)
    end if
    if ( flag==0 ) then
      open (iunit,file=lfile,status='old',form='unformatted',err=100,        &
            action='read')
    else if ( flag==1 ) then
      open (iunit,file=lfile,status='old',err=100,action='read')
    else if ( flag==2 ) then
      open (iunit,file=lfile,status='unknown',form='unformatted',err=100)
    else if ( flag==3 ) then
      open (iunit,file=lfile,status='unknown',err=100)
    else
      call mvliberror('openmuseofiles','Bad flag passed while opening '//     &
                      trim(safe),flag,-9999.0)
    end if
    if ( iflg(1)>0 ) write (6,'(2x,a,i2,a)')                       &
                           'openmuseofiles successfully open '//   &
                           trim(lfile)//' (' , iunit , ')'
    return
  100  write (6,'(10x,a)') 'MVLib: error opening '//trim(lfile)
    if ( iflg(62)==1 ) then
      close (iunit)
      call setmvlibinti(1,1)
      iflg(63) = 1
      return
    else
      call exit(0)
    end if
  end subroutine openmuseofiles
     
  subroutine museofile(ifile,ofile)
    implicit none
    character(len=*) :: ifile , ofile
    intent (out) ofile
    call museopath
    if ( trim(ifile)=='river.contour' ) then
      filename = pathm(1:npathm)//'dem/'//ifile
    else
      filename = pathm(1:npathm)//ifile
    end if
    ofile = filename
  end subroutine museofile
  
  ! questa routine potrebbe essere chiamata da museomm52d
  subroutine museomm52drec(u,ih,v,lat,lon,ixm,jxm,lix,ljx)
    implicit none
    integer :: ih , ixm , jxm , lix , ljx , u
    real , dimension(ixm,jxm) :: lat , lon
    real , dimension(ixm,jxm,25) :: v
    intent (in) ih , ixm , jxm , u
    intent (out) lat , lon , v
    intent (inout) lix , ljx
    integer :: ii , jj , lkx , lns , n
    if ( ih>0 ) then
      read (u) lix , ljx , lkx , lns , ((lat(ii,jj),ii=1,lix),jj=1,ljx) ,    &
               ((lon(ii,jj),ii=1,lix),jj=1,ljx) ,                            &
               (((v(ii,jj,1),ii=1,lix),jj=1,ljx),n=1,ih)
    else
      read (u) lix , ljx , lkx , lns , ((lat(ii,jj),ii=1,lix),jj=1,ljx) ,    &
               ((lon(ii,jj),ii=1,lix),jj=1,ljx) ,                            &
               (((v(ii,jj,n),ii=1,lix),jj=1,ljx),n=1,lns)
    end if
  end subroutine museomm52drec
     
  subroutine museomm52dw(iy,im,ig,ih,id,u,v,lat,lon,ixm,jxm,lix,ljx)
    implicit none
    integer :: id , ig , ih , im , ixm , iy , jxm , lix , ljx
    real , dimension(ixm,jxm) :: lat , lon
    real , dimension(ixm,jxm,25) :: u , v
    intent (in) id
    intent (inout) lix , ljx
    character(60) :: cfile
    integer :: iday , irec , irec0 , iyold , lix1 , ljx1 , lun91 , lun92
    data irec0 , iyold/0 , -1/
    data lun91 , lun92/ - 1 , -1/
    if ( iy/=iyold ) then
      write (cfile,'(a,i4,i1,a)') 'mm5ugr' , iy , id , '.dat'
      call openmuseofiles(lun91,cfile,0)
      write (cfile,'(a,i4,i1,a)') 'mm5vgr' , iy , id , '.dat'
      call openmuseofiles(lun92,cfile,0)
      irec0 = 0
      iyold = iy
    end if
    irec = indexofyear(ig,im,iy)
    if ( irec>366 .or. irec<1 .or. ih>24 .or. ih<0 ) then
      write (6,'(7x,a)') 'MVLib museomm52dw: bad parameters passed. Return'
      return
    end if
    if ( irec<=irec0 ) then
      irec0 = 0
      rewind (lun91)
      rewind (lun92)
    end if
    do iday = irec0 + 1 , irec - 1
      read (lun91)
      read (lun92)
    end do
    call museomm52drec(lun91,ih,u,lat,lon,ixm,jxm,lix,ljx)
    call museomm52drec(lun92,ih,v,lat,lon,ixm,jxm,lix1,ljx1)
    if ( lix1/=lix .or. ljx1/=ljx ) then
      lix = 0
      ljx = 0
    end if
    irec0 = irec
  end subroutine museomm52dw
     
  subroutine museomm52d(fld,iy,im,ig,ih,id,v,lat,lon,ixm,jxm,lix,ljx)
    implicit none
    character(len=*) :: fld
    integer :: id , ig , ih , im , ixm , iy , jxm , lix , ljx
    real , dimension(ixm,jxm) :: lat , lon
    real , dimension(ixm,jxm,24) :: v
    intent (in) id , ih , ixm , jxm
    intent (out) lat , lon , v
    intent (inout) lix , ljx
    character(60) :: file1 , file2
    integer :: iday , ii , irec , irec0 , irec1 , isl , jj , lkx , lns ,    &
               lun , n
    data file1/'pippo'/
    data irec0/0/
    data lun/ - 1/
    write (file2,'(a,i4,i1,a)') 'mm5'//trim(fld) , iy , id , '.dat'
    if ( file1/=file2 ) then
      if ( lun>0 ) then
        close (lun)
      else
        call getlun(lun)
      end if
      call openmuseofiles(lun,file2,0)
      file1 = file2
      irec0 = 0
      if ( iflg(1)>=100 ) write (6,'(7x,a)') 'MVLib museomm52d reopen'
    end if
    irec = indexofyear(ig,im,iy)
    if ( irec>366 .or. irec<1 .or. abs(ih)>25 ) then
      write (6,'(7x,a)') 'MVLib museomm52d: bad parameters passed. Return'
      write (6,'(20x,a,2i10)') 'irec and ih were ' , irec , ih
      return
    end if
    if ( irec<=irec0 ) then
      irec0 = 0
      rewind (lun)
      if ( iflg(1)>=100 ) write (6,'(7x,a)') 'MVLib museomm52d rewind'
    end if
    irec1 = irec0 + 1
    if ( iflg(1)>=100 ) &
      write (6,'(7x,a,i5)') 'museomm52d reading rec: ', irec
    do iday = irec1 , irec - 1
      read (lun)
    end do
    if ( ih>0 ) then
      read (lun) lix , ljx , lkx , lns , ((lat(ii,jj),ii=1,lix),jj=1,ljx) ,  &
                 ((lon(ii,jj),ii=1,lix),jj=1,ljx) ,                          &
                 (((v(ii,jj,1),ii=1,lix),jj=1,ljx),n=1,ih)
    else if ( ih==0 ) then
      read (lun) lix , ljx , lkx , lns , ((lat(ii,jj),ii=1,lix),jj=1,ljx) ,  &
                 ((lon(ii,jj),ii=1,lix),jj=1,ljx) ,                          &
                 (((v(ii,jj,n),ii=1,lix),jj=1,ljx),n=1,24)
    else
      isl = -ih
      read (lun) lix , ljx , lkx , lns , ((lat(ii,jj),ii=1,lix),jj=1,ljx) ,  &
                 ((lon(ii,jj),ii=1,lix),jj=1,ljx) ,                          &
                 (((v(ii,jj,n),ii=1,lix),jj=1,ljx),n=1,isl)
    end if
    irec0 = irec
    call mvsetflags('Dominio',float(id))
  end subroutine museomm52d
     
  subroutine museosynop(iop,year,synop)
    implicit none
    integer :: iop , year
    real , dimension(366,24,200) :: synop
    intent (in) iop , year
    intent (out) synop
    character(60) :: cfile
    integer :: i , lun
    if ( iflg(1)>=10 ) &
      write (6,'(7x,a,i2,a,i4,a)')  &
         'MVLib: Reading record ' , iop , ' of year ' , year , &
         ' from MuSEO files of SYNOP.'
    write (cfile,'(a,i4,a)') 'synop' , year , '.dat'
    call getlun(lun)
    call openmuseofiles(lun,cfile,0)
    do i = 1 , iop - 1
      read (lun)
    end do
    read (lun) synop
    close (lun)
  end subroutine museosynop
     
  subroutine synopgrid(zlat,zlon,maxs,nz)
    implicit none
    integer :: maxs , nz
    real , dimension(maxs) :: zlat , zlon
    intent (in) maxs
    intent (out) nz , zlat , zlon
    integer :: i , lun
    call getlun(lun)
    call openmuseofiles(lun,'comuni.synop',1)
    do i = 1 , maxs
      read (lun,'(30x,f7.4,4x,f7.4)',end=100) zlat(i) , zlon(i)
      nz = i
    end do
  100  close (lun)
  end subroutine synopgrid
     
  subroutine museoarssa(record,anno,vettore)
    implicit none
    integer :: anno , record
    real , dimension(366,24,100) :: vettore
    intent (in) anno , record
    intent (out) vettore
    character(60) :: cfile
    integer :: i , lun
    write (cfile,'(a,i4,a)') 'scerni' , anno , '.dat'
    call getlun(lun)
    call openmuseofiles(lun,cfile,0)
    do i = 1 , record - 1
      read (lun)
    end do
    read (lun) vettore
    close (lun)
  end subroutine museoarssa
     
  subroutine arssagrid(zlat,zlon,maxs,nz)
    implicit none
    integer :: maxs , nz
    real , dimension(maxs) :: zlat , zlon
    intent (in) maxs
    intent (out) nz , zlat , zlon
    integer :: i , lun
    call getlun(lun)
    call openmuseofiles(lun,'comuni.arssa',1)
    do i = 1 , maxs
      read (lun,'(33x,2f12.4)',end=100) zlat(i) , zlon(i)
      nz = i
    end do
   100  close (lun)
  end subroutine arssagrid
     
  subroutine arssacode(lat,lon)
    implicit none
    integer , parameter :: maxcom = 100
    real , dimension(maxcom) :: lat , lon
    intent (out) lat , lon
    integer :: code , i , lun
    do i = 1 , maxcom
      lat(i) = 0.0
      lon(i) = 0.0
    end do
    call getlun(lun)
    call openmuseofiles(lun,'comuni.arssa',1)
    do i = 1 , maxcom
      read (lun,'(i2,36x,f7.4,5x,f7.4)',end=100) code , lat(code) , lon(code)
    end do
   100  close (lun)
  end subroutine arssacode
     
  subroutine worldtemp(ora,giorno,mese,anno,p,lat,lon,n)
    implicit none
    integer :: anno , giorno , mese , n , ora
    real , dimension(:) :: lat , lon , p
    intent (out) lat , lon , p
    intent (inout) n
    logical :: fbady , first
    integer :: i , ii , isource , last , lastyear , lun , nn , rec
    data first/.true./
    data fbady/.false./
    save lun,lastyear,last,first,fbady
    if ( first ) then
                     ! Find bad year?
      lastyear = anno
      first = .false.
      call getlun(lun)
      call mvsaveflags
      call mvsetflags('openmuseofiles behaviour',1.0)
      call openmuseodb(lun,'worldtemp',anno)
      call mvgetiflags(63,ii)
      call mvrestflags
      if ( ii==1 ) then
        fbady = .true.
        return
      end if
      last = 0
    else if ( anno/=lastyear ) then
      lastyear = anno
      if ( fbady ) then
        call getlun(lun)
        fbady = .false.
      else
        close (lun)
      end if
      call mvsaveflags
      call mvsetflags('openmuseofiles behaviour',1.0)
      call openmuseodb(lun,'worldtemp',anno)
      call mvgetiflags(63,ii)
      call mvrestflags
      if ( ii==1 ) then
        fbady = .true.
        return
      end if
      last = 0
    else if ( fbady ) then
      return
    end if
    rec = index1h(ora,giorno,mese,anno)
    if ( rec==last+1 ) then
      read (lun,err=100,end=100) nn , &
           (lat(i+n),lon(i+n),isource,p(i+n),i=1,nn)
    else
      rewind (lun)
      do ii = 1 , rec - 1
        read (lun,err=100,end=100)
      end do
      read (lun,err=100,end=100) nn , &
           (lat(i+n),lon(i+n),isource,p(i+n),i=1,nn)
    end if
    last = rec
    n = n + nn
    return
   100  close (lun)
    write (6,'(10x,a)') 'Flux error inside worldtemp. Exiting...'
    call exit(0)
  end subroutine worldtemp
     
  subroutine worldrain(ora,giorno,mese,anno,p,lat,lon,n)
    implicit none
    integer :: anno , giorno , mese , n , ora
    real , dimension(:) :: lat , lon , p
    intent (out) lat , lon , p
    intent (inout) n
    logical :: fbady , first
    integer :: i , ii , isource , last , lastyear , lun , nn , rec
    data first/.true./
    data fbady/.false./
    save lun,lastyear,last,first,fbady
    if ( first ) then
                     ! Find bad year?
      lastyear = anno
      first = .false.
      call getlun(lun)
      call mvsaveflags
      call mvsetflags('openmuseofiles behaviour',1.0)
      call openmuseodb(lun,'worldrain',anno)
      call mvgetiflags(63,ii)
      call mvrestflags
      if ( ii==1 ) then
        fbady = .true.
        return
      end if
      last = 0
    else if ( anno/=lastyear ) then
      lastyear = anno
      if ( fbady ) then
        call getlun(lun)
        fbady = .false.
      else
        close (lun)
      end if
      call mvsaveflags
      call mvsetflags('openmuseofiles behaviour',1.0)
      call openmuseodb(lun,'worldrain',anno)
      call mvgetiflags(63,ii)
      call mvrestflags
      if ( ii==1 ) then
        fbady = .true.
        return
      end if
      last = 0
    else if ( fbady ) then
      return
    end if
    rec = index1h(ora,giorno,mese,anno)
    if ( rec==last+1 ) then
      read (lun,err=100,end=100) nn , &
           (lat(i+n),lon(i+n),isource,p(i+n),i=1,nn)
    else
      rewind (lun)
      do ii = 1 , rec - 1
        read (lun,err=100,end=100)
      end do
      read (lun,err=100,end=100) nn , &
          (lat(i+n),lon(i+n),isource,p(i+n),i=1,nn)
    end if
    last = rec
    n = n + nn
    return
  100  close (lun)
    write (6,'(10x,a)') 'Flux error inside worldrain. Exiting...'
    call exit(0)
  end subroutine worldrain
     
  subroutine hourlyrain(ora,giorno,mese,eanno,p,lat,lon,n)
    implicit none
    integer :: eanno , giorno , mese , n , ora
    real , dimension(:) :: lat , lon , p
    intent (out) lat , lon , p
    intent (inout) n
    integer :: anno , i , ii , irec , isource , last , lastyear , lun
    character(40) :: cfile
    logical :: first
    data first/.true./
    save lun,first,lastyear,last
    anno = i4digityear(eanno)
    if ( anno<mvlibmagicnum(2) .or. anno>mvlibmagicnum(3) ) then
      n = 0
      return
    end if
    if ( first ) then
      call getlun(lun)
      if (srcflag(30)) then 
        write (cfile,'(a,i4,a)') 'pioggiaoraria-fil' , anno , '.dat'
      else
        write (cfile,'(a,i4,a)') 'pioggiaoraria' , anno , '.dat'
      end if
      call openmuseofiles(lun,cfile,0)
      last = 0
      lastyear = anno
      first = .false.
    else if ( anno/=lastyear ) then
      close (lun)
      if (srcflag(30)) then 
        write (cfile,'(a,i4,a)') 'pioggiaoraria-fil' , anno , '.dat'
      else
        write (cfile,'(a,i4,a)') 'pioggiaoraria' , anno , '.dat'
      end if
      call openmuseofiles(lun,cfile,0)
      last = 0
      lastyear = anno
    end if
    irec = index1h(ora,giorno,mese,anno)
    if ( irec>=1 .and. irec<=366*24 ) then
      if ( irec<=last ) then
        rewind (lun)
        last = 0
      end if
      do ii = last + 1 , irec - 1
        read (lun,err=100,end=100) n
      end do
      n = 0
      read (lun,err=100,end=100) n , (lat(i),lon(i),isource,p(i),i=1,n)
      last = irec
      return
    end if
  100  close (lun)
    write (6,'(10x,a)') 'Flux error inside hourlyrain. Exiting...'
    write (6,'(10x,a,4i8)') 'hh/dd/mm/yy are:' , ora , giorno , mese , anno
    call exit(0)
  end subroutine hourlyrain
     
  subroutine hourlytemp(ora,giorno,mese,eanno,p,lat,lon,n)
    implicit none
    integer :: eanno , giorno , mese , n , ora
    real , dimension(:) :: lat , lon , p
    intent (out) lat , lon , p
    intent (inout) n
    integer :: anno , i , ii , irec , isource , last , lastyear , lun
    logical :: first
    data first/.true./
    save lun,lastyear,last,first
    anno = i4digityear(eanno)
    if ( anno<mvlibmagicnum(2) .or. anno>mvlibmagicnum(3) ) then
      n = 0
      return
    end if
    if ( first ) then
      call getlun(lun)
      call openmuseodb(lun,'temporaria',anno)
      last = 0
      lastyear = anno
      first = .false.
    else if ( anno/=lastyear ) then
      close (lun)
      call openmuseodb(lun,'temporaria',anno)
      last = 0
      lastyear = anno
    end if
    irec = index1h(ora,giorno,mese,anno)
    if ( irec>=1 .and. irec<=366*24 ) then
      if ( irec<=last ) then
        rewind (lun)
        last = 0
      end if
      do ii = last + 1 , irec - 1
        read (lun,err=100,end=100) n
      end do
      n = 0
      read (lun,err=100,end=100) n , (lat(i),lon(i),isource,p(i),i=1,n)
      last = irec
      return
    end if
   100  close (lun)
    write (6,'(10x,a)') 'Flux error inside hourltemp. Exiting...'
    write (6,'(10x,a,4i8)') 'hh/dd/mm/yy are:' , ora , giorno , mese , anno
    call exit(0)
  end subroutine hourlytemp
     
  subroutine idrograficosensori(iu,flag,rsens,com,lat,lon,irec,nc)
    implicit none
    integer :: flag , iu , nc
    character(len=*) :: rsens
    character(len=*) , dimension(:) :: com
    integer , dimension(:) :: irec
    real , dimension(:) :: lat , lon
    intent (in) flag
    intent (out) irec
    intent (inout) nc
    integer :: i , n , nch
    character(60) :: sens , sensore , st
    if ( flag==1 ) then
      call openmuseofiles(iu,'comuni.idrografico',1)
    else if ( flag==2 ) then
      call openmuseofiles(iu,'comuni.idrografico-all',1)
    else if ( flag==3 ) then
      call openmuseofiles(iu,'comuni.iira',1)
    else
      stop '  Stop inside idrograficosensori - Bad value of flag parameter'
    end if
    nc = 0
    n = 1
    call cv2lower(rsens,sensore)
    nch = len_trim(sensore)
    do while ( n/=0 )
      call idrograficosensore(iu,n,i,st,com(nc+1),sens,lat(nc+1),lon(nc+1))
      if ( n>0 .and. (sens(1:nch)==sensore(1:nch) .or. nch==1) ) then
        nc = nc + 1
        irec(nc) = n
      end if
    end do
    close (iu)
  end subroutine idrograficosensori
 
  subroutine idrograficosensore(iunit,n,icd,loc1,loc2,sensore,lat,lon)
    implicit none
    integer :: icd , iunit , n
    real :: lat , lon
    character(len=*) :: loc1 , loc2 , sensore
    intent (in) iunit
    intent (out) icd , lat , loc1 , loc2 , lon , n
    character(20) :: com
    character(10) :: mis
    character(2) :: prov
    character(30) :: sens
    read (iunit,'(i3,1x,a45,a20,a2)',end=100) n , loc1 , com , prov
    loc2 = trim(com)//' ('//prov//')'
    read (iunit,'(4x,i6,2x,a30,a10,2f8.4)') icd , sens , mis , lat , lon
    sensore = trim(sens)//' ('//trim(mis)//')'
    call cv2lower(sensore,sensore)
    return
  100  n = 0
  end subroutine idrograficosensore
     
  subroutine idrograficomisure(u,sensore,ih,id,im,iy,comuni,lat,lon,va,ng)
    implicit none
    integer :: id , ih , im , iy , ng , u
    character(len=*) :: sensore
    character(len=*) , dimension(:) :: comuni
    real , dimension(:) :: lat , lon , va
    intent (inout) comuni , lat , lon , ng , va
    integer :: i , ip , nc
    call idrograficosensori(u,2,sensore,comuni,lat,lon,i1dum,nc)
    call idrograficomisura(u,x1dum,0,ih,id,im,iy,va(i),ip)
    do i = 1 , nc
      call idrograficomisura(u,x1dum,i1dum(i),ih,id,im,iy,va(i),ip)
    end do
    ng = nc
    i = 1
    do while ( ng>i )
      if ( nint(va(i))==-9999 ) then
        comuni(i) = comuni(ng)
        lat(i) = lat(ng)
        lon(i) = lon(ng)
        va(i) = va(ng)
        ng = ng - 1
      end if
      if ( nint(va(i))/=-9999 ) i = i + 1
    end do
    if ( nint(va(ng))==-9999 ) ng = ng - 1
    write (6,'(i20,a)') nc , ' sensors found.'
    write (6,'(i20,a)') ng , ' good measurements found.'
  end subroutine idrograficomisure
     
  subroutine idrograficomisura(u,v,irec,ora,giorno,mese,anno,va,ip)
    implicit none
    integer :: anno , giorno , ip , irec , mese , ora , u
    real :: va
    real , dimension(:) :: v
    intent (in) irec , u
    intent (inout) ip , v , va
    integer :: i , ir , lasr , lasy , nn , np , nrec
    integer , dimension(1000) :: iflag
    character(11) :: sensore
    data lasy/1800/
    save lasy,lasr,nrec
    if ( lasy/=anno .or. irec<=0 ) then
      call museopath
      write (filename,'(a,i4,a)') pathm(1:npathm)//'comuni.idrografico-all'
      if ( lasy==1800 ) then
        open (u,file=filename,status='old')
        do i = 1 , 1000
          read (u,'(a11)',end=20)
          read (u,'(12x,a11)') sensore
          if ( sensore(1:11)=='pluviometro' ) then
            iflag(i) = 1
          else
            iflag(i) = 0
          end if
          nrec = i
        end do
   20   close (u)
      end if
      write (filename,'(a,i4,a)') pathm(1:npathm)//'idro-all' , anno , '.dat'
      open (u,file=filename,status='old',form='unformatted')
      if ( irec<=0 ) return
      lasr = 0
      lasy = anno
    else if ( lasr>=irec ) then
      rewind (u)
      lasr = 0
    end if
    if ( irec>nrec ) then
      write (6,'(12x,a,i10)') &
        'idrograficomisura: Access to unexisting record: ', irec
      stop ' '
    end if
    do ir = lasr + 1 , irec
      read (u) nn , (v(i),i=1,nn)
    end do
    lasr = irec
    ip = index15ofyear(00,ora,giorno,mese,anno)
    if ( ip-1>nn ) then
      va = -9999
      do i = nn + 1 , ip
        v(i) = -9999
      end do
      return
    end if
    va = 0.0
    np = 0
    do i = ip - 1 , ip + 2
      if ( nint(v(i))/=-9999 ) then
        va = va + v(i)
        np = np + 1
      end if
    end do
    if ( np==0 ) then
      va = -9999
    else if ( iflag(irec)==0 ) then
      va = va/np
    end if
  end subroutine idrograficomisura
     
  subroutine montemidiavec(irec,year,v1,v2,v3,n)
    implicit none
    integer , parameter :: nlon = 480 , nlat = 480 , maxdata = 180938
    integer :: n , irec , year
    real , dimension(maxdata) :: v1 , v2 , v3
    intent (out) v1 , v2 , v3
    intent (inout) n
    integer :: i , j , map , n1
    real , dimension(nlon,nlat) :: rain
    real :: ylat , ylon
    call montemidiamap(irec,year,rain,n1)
    if ( n1>=0 ) then
      n = 0
      do i = 1 , nlon
        do j = 1 , nlat
          call montemidiagrid(i,j,ylon,ylat,map)
          if ( map==1 ) then
            n = n + 1
            v1(n) = rain(i,j)
            v2(n) = ylat
            v3(n) = ylon
          end if
        end do
      end do
    else
      n = 0
    end if
  end subroutine montemidiavec

  subroutine montemidiamap(irec,year,rain,n)
    implicit none
    integer , parameter :: nlon = 480 , nlat = 480
    integer :: n , irec , year
    real , dimension(nlon,nlat) :: rain
    call radarlikemap(irec,year,rain,n,1)
  end subroutine montemidiamap
     
  subroutine wrfmap(irec,year,rain,n)                ! No longer used
    implicit none
    integer , parameter :: nlon = 480 , nlat = 480
    integer :: n , irec , year
    real , dimension(nlon,nlat) :: rain
    call radarlikemap(irec,year,rain,n,2)
  end subroutine wrfmap
     
  subroutine micramap(irec,year,rain,n)                ! No longer used
    implicit none
    integer , parameter :: nlon = 480 , nlat = 480
    integer :: n , irec , year
    real , dimension(nlon,nlat) :: rain
    call radarlikemap(irec,year,rain,n,3)
  end subroutine micramap
     
  subroutine radarlikemap(irec,year,rain,n,db)
    implicit none
    integer , parameter :: nlon = 480 , nlat = 480 , maxdata = 180938
    integer :: db , n , irec , year
    real , dimension(nlon,nlat) :: rain
    intent (in) db , irec
    intent (out) rain
    intent (inout) n
    logical :: first
    integer :: i , j , k , lun , olddb , oldrec , oldyear
    integer , intrinsic :: nint
    real , dimension(maxdata) :: v1 , v2 , v3
    data first/.true./
    data olddb/ - 1/
    save oldyear,oldrec,lun,first,olddb
    if ( first ) then
      call getlun(lun)
      if ( db==1 ) then
        call openmuseodb(lun,'MonteMidia',year)
      else if ( db==2 ) then
        call openmuseodb(lun,'wrf',year)                   ! No longer used
      else if ( db==3 ) then
        call openmuseodb(lun,'Micra',year)              ! No longer used
      else
        stop ' Flux error inside radarlikemap'
      end if
      olddb = db
      oldyear = year
      oldrec = 0
      first = .false.
    else if ( oldyear/=year .or. db/=olddb ) then
      close (lun)
      if ( db==1 ) then
        call openmuseodb(lun,'MonteMidia',year)
      else if ( db==2 ) then
        call openmuseodb(lun,'wrf',year)                   ! No longer used
      else if ( db==3 ) then
        call openmuseodb(lun,'Micra',year)              ! No longer used
      end if
      oldyear = year
      oldrec = 0
    end if
    if ( irec/=oldrec+1 ) then
      rewind (lun)
      do i = 1 , irec - 1
        read (lun)
      end do
    end if
    read (lun) n , (v1(i),i=1,n) , (v2(i),i=1,n) , (v3(i),i=1,n)
    oldrec = irec
    rain = 0.0
    do k = 1 , n
      i = nint(v2(k))
      j = nint(v3(k))
      rain(i,j) = v1(k)
    end do
  end subroutine radarlikemap
     
  subroutine montemidiagrid(i,k,ylon,ylat,map)
    implicit none
    integer :: i , k , map
    real :: ylat , ylon
    intent (in) i , k
    intent (out) map
    intent (inout) ylat
    real , save :: d , pi , re , rlat , rlon , x , y
    logical :: first
    integer :: j , nlat , nlon
    data first/.true./
    save nlon,nlat,first
    if ( first ) then
      d = 0.5                ! Passo della griglia (km)
      rlat = 42.057                ! Latitudine della posizionde del radar
      rlon = 13.177                ! Longitudine della posizionde del radar
      re = 6371.0                ! Raggio terrestre
      pi = (4.*atan(1.0))/180.     ! Pi greco / 180
      nlon = 480                ! Numero di Longitudini
      nlat = 480                ! Numero di Latitudini
      first = .false.
    end if
    !  La riga seguente ci sta perchè nelle convenzioni dei radaristi
    !  la prima latitudine è quella in alto, ciò che chiaramente non è
    !  comodo per i grafici. A ciò corrisponde il fatto che in radar2museo
    !  i files di Picciotti vengono letti scambiando le latitudini.
    !  do j=nlat,1,-1 ; read(lun,*,end=10,err=10) (p(i,j),i=1,nlon) ; enddo
    !  Vedi la routine takeonehour.
    !  Per la riga seguente, questa routine fornisce la latitudine
    !  secondo la nostra (di museo) convenzione.
    j = 1 + nlat - k
    x = (i-nlon/2-0.5)*d   ! Distanza x dal radar (km) da -119.75 a 119.75
    y = (nlat/2-j+0.5)*d   ! Distanza y dal radar (km) da -119.75 a 119.75
    ylat = rlat + y/((pi)*re)
    ylon = rlon + x/((pi)*re*cos(ylat*pi))
    map = 0
    if ( distance(ylat,ylon,rlat,rlon)<120000. ) map = 1
  end subroutine montemidiagrid

  integer function nrectimeseries(db)
    implicit none
    character(len=*) :: db
    if ( trim(db)=='iira' ) then
      nrectimeseries = 464
    else if ( trim(db)=='lazio' ) then
      nrectimeseries = 137
    else if ( trim(db)=='arssa' ) then
      nrectimeseries = 300
    else if ( trim(db)=='portateiira' ) then
      nrectimeseries = 54
    else if ( trim(db)=='Po' ) then
      nrectimeseries = 7
    else if ( trim(db)=='dewetra' .or. &
              trim(db)=='dewrain' .or. &
              trim(db)=='dewrain-fil' .or. &
              trim(db)=='dewtemp' ) then
      nrectimeseries = mvlibmagicnum(1)
    else if ( trim(db)=='dewidro' ) then
      nrectimeseries = mvlibmagicnum(7)
    else if ( trim(db)=='acqwarhone' ) then
      nrectimeseries = 38
    else if ( trim(db)=='acqwapo' ) then
      nrectimeseries = 501
    else if ( trim(db)=='acqwachu' ) then
      nrectimeseries = 2
    else
      nrectimeseries = 0
    end if
  end function nrectimeseries
     
  subroutine museoanagrafica(edb,irec,com,sens,xlat,xlon)
    implicit none
    character(len=*) :: com , edb , sens
    integer :: irec
    real :: xlat , xlon
    intent (in) edb , irec
    intent (out) xlat , xlon
    intent (inout) com , sens
    character(80) :: com2
    character(40) :: db , cfile , oldfile
    logical :: first
    character(64) :: cform
    integer :: i , icd , ir , jj , ldb , lun , nn , nrec , nst , oldrec
    character(45) :: loc
    character(10) :: mis
    character(2) :: prov
    data first/.true./
    save :: oldrec,lun,first
    db = edb
    ldb = len_trim(db)
    if ( db(1:ldb)=='iira' .or. db(1:ldb)=='lazio' .or. &
         db(1:ldb)=='arssa' .or. db(1:ldb)=='portateiira' .or. &
         db(1:ldb)=='Po' .or. db(1:ldb)=='acqwarhone' .or. &
         db(1:ldb)=='acqwapo' .or. db(1:ldb)=='acqwachu' ) then
      nrec = nrectimeseries(db)
      write (cfile,'(a)') 'comuni.'//db(1:ldb)
    else if ( db(1:ldb)=='dewetra' .or. &
              db(1:ldb)=='dewrain' .or. &
              db(1:ldb)=='dewrain-fil' .or. &
              db(1:ldb)=='dewtemp' ) then
      nrec = mvlibmagicnum(1)
      cfile = 'comuni.dewetra'
      sens = ' '
      if ( db(1:ldb)=='dewrain' .or. db(1:ldb)=='dewrain-fil' ) sens = 'Precipitation (mm)'
      if ( db(1:ldb)=='dewtemp' ) sens = 'Temperature (C)'
      db = 'dewetra'
    else if ( db(1:ldb)=='dewidro' ) then
      nrec = mvlibmagicnum(7)
      cfile = 'comuni.dewidro'
      sens = 'Hydrometric level (m)'
      db = 'dewetra'
    else if ( db(1:ldb)=='dewdisc' ) then
      nrec = mvlibmagicnum(1)
      cfile = 'comuni.dewdisc'
      sens = 'Discharge (m3/s)'
      db = 'dewetra'
    else
      write (6,'(10x,a)') 'museoanagrafica error, unknown db: '//db
      call exit(0)
    end if
    ldb = len_trim(db)
    if ( irec<1 .or. irec>nrec ) then
      write (6,'(10x,a)') 'museoanagrafica error, invalid record for '//     &
                        trim(db)//' database'
      write (6,'(10x,a,i4)') 'It must be in the range 1-' , nrec
      call exit(0)
    end if
    if ( first ) then
      lun = -1
      call openmuseofiles(lun,cfile,1)
      oldfile = cfile
      oldrec = 0
      first = .false.
    else if ( oldfile/=cfile ) then
      close (lun)
      call openmuseofiles(lun,cfile,1)
      oldfile = cfile
      oldrec = 0
    end if
    if ( irec<=oldrec ) then
      if ( iflg(1)>0 ) write (6,'(10x,a,i3)') 'museoanagrafica rewinds '//   &
             trim(cfile)//' and skips until record ' , irec
      rewind (lun)
      if ( db(1:ldb)=='lazio' .or. &
           db(1:ldb)=='portateiira' .or. &
           db(1:2)=='Po' .or. db(1:ldb)=='dewetra' .or. &
           db(1:10)=='acqwarhone' .or. db(1:7)=='acqwapo' .or. &
           db(1:8)=='acqwachu' ) then
        do i = 1 , irec - 1
          read (lun,'(a)')
        end do
      else
        do i = 1 , irec - 1
          read (lun,'(a)')
          read (lun,'(a)')
        end do
      end if
    else if ( irec/=oldrec+1 ) then
      if ( iflg(1)>0 ) write (6,'(10x,a,i3,a,i3)') 'museoanagrafica skips '//&
        trim(cfile)//' from record ' ,oldrec + 1 , ' to ' , irec
      if ( db(1:ldb)=='lazio' .or. &
           db(1:ldb)=='portateiira' .or. &
           db(1:ldb)=='Po' .or. db(1:ldb)=='dewetra' .or. &
           db(1:10)=='acqwarhone' .or. db(1:7)=='acqwapo' .or. &
           db(1:8)=='acqwachu' ) then
        do i = oldrec + 1 , irec - 1
          read (lun,'(a)')
        end do
      else
        do i = oldrec + 1 , irec - 1
          read (lun,'(a)')
          read (lun,'(a)')
        end do
      end if
    end if
    if ( db(1:ldb)=='lazio' ) then
      read (lun,'(a35,2x,i9,2f10.5,1x,a2)') com , icd , xlat , xlon , prov
      com = trim(com)//' ('//prov//')'
      sens = 'Pluviometro (mm)'
    else if ( db(1:ldb)=='Po' ) then
      read (lun,'(a35,a2,i9,2f10.5)') com , prov , jj , xlat , xlon
      com = trim(com)//' ('//prov//')'
      sens = 'Discharge (m3/sec)'
    else if ( db(1:ldb)=='portateiira' ) then
      read (lun,'(13x,a40,5x,i2,2f10.6)') com , nn , xlat , xlon
      write (sens,'(i2,a)') nn , ' Discharge (m3/sec)'
    else if ( db(1:ldb)=='dewetra' ) then
      read (lun,'(64x,a2,1x,a50,2f8.3)') prov , com2 , xlat , xlon
      com = trim(com2)//' ('//prov//')'
    else if ( db(1:ldb)=='acqwarhone' .or. db(1:ldb)=='acqwapo' .or.        &
              db(1:ldb)=='acqwachu' ) then
      call museoformat(6,cform)
      read (lun,cform) nst , com , sens , xlat , xlon
    else
      read (lun,'(i3,1x,a45,a20,a2)') ir , loc , com , prov
      read (lun,'(4x,i6,2x,a30,a10,2f8.4)') icd , sens , mis , xlat , xlon
      call nameadjust(com)
      com = trim(com)//' ('//prov//')'
      sens = trim(sens)//' ('//trim(mis)//')'
    end if
    oldrec = irec
  end subroutine museoanagrafica
     
  subroutine museogoodseries(x1,y1,ifirst,ndat,x2,y2,n)
    implicit none
    integer :: ifirst , n , ndat
    real , dimension(:) :: x1 , x2 , y1 , y2
    intent (in) ifirst , ndat , x1 , y1
    intent (out) x2 , y2
    intent (inout) n
    integer :: i
    n = 0
    do i = ifirst , ifirst + ndat - 1
      if ( nint(y1(i))/=-9999 ) then
        n = n + 1
        x2(n) = x1(i)
        y2(n) = y1(i)
      end if
    end do
  end subroutine museogoodseries
     
  subroutine museotimeseries(db,irec,anno,x,y,n)
    implicit none
    integer :: anno , n , irec
    character(*) :: db
    real , dimension(:) :: x , y
    intent (in) irec
    intent (out) x , y
    intent (inout) n
    real :: dt
    character(40) :: cfile , oldfile
    logical :: first
    integer :: i , lun , n1 , ncomp , nrec , oldrec , ldb
    character(80) :: str1 , str2
    data first/.true./
    save first,oldrec,lun,oldfile
    ldb = len_trim(db)
    if ( db(1:ldb)=='iira' ) then
      call chktimeseriesrange(db,anno,2003,2011)
      nrec = 464
      dt = 1.0/(24*4)
      ncomp = 366*24*4
    else if ( db(1:ldb)=='lazio' ) then
      call chktimeseriesrange(db,anno,2007,2011)
      nrec = 137
      dt = 1.0/(24*4)
      ncomp = 366*24*4
    else if ( db(1:ldb)=='arssa' ) then
      call chktimeseriesrange(db,anno,1997,2007)
      nrec = 300
      dt = 1.0/24
      ncomp = 366*24
    else if ( db(1:ldb)=='Po' ) then
      call chktimeseriesrange(db,anno,1995,2005)
      nrec = 7
      dt = 1.0/24
      ncomp = 365*24
      if ( mod(anno,4)==0 ) ncomp = 366*24
    else if ( db(1:ldb)=='dewrain' .or. db(1:ldb)=='dewtemp' &
            .or. db(1:ldb)=='dewrain-fil' ) then
      call chktimeseriesrange(db,anno,2000,-1)
      nrec = mvlibmagicnum(1)
      dt = 1.0/(24*4)
      ncomp = 366*24*4
    else if ( db(1:ldb)=='dewidro' ) then
      call chktimeseriesrange(db,anno,2000,-1)
      nrec = mvlibmagicnum(7)
      dt = 1.0/(24*4)
      ncomp = 366*24*4
    else if ( db(1:ldb)=='dewdisc' ) then
      call chktimeseriesrange(db,anno,2000,-1)
      nrec = mvlibmagicnum(1)
      dt = 1.0/(24*4)
      ncomp = 366*24*4
    else if ( db(1:ldb)=='acqwarhone' ) then
      call chktimeseriesrange(db,anno,1990,2008)
      nrec = 38
      dt = 1.0/24
      ncomp = 365*24
      if ( mod(anno,4)==0 ) ncomp = 366*24
    else if ( db(1:ldb)=='acqwapo' ) then
      call chktimeseriesrange(db,anno,1999,2010)
      nrec = 501
      dt = 1.0/24
      ncomp = 365*24
      if ( mod(anno,4)==0 ) ncomp = 366*24
    else
      write (6,'(10x,a)') 'museotimeseries error, unknown db: '//trim(db)
      call exit(0)
    end if
    if ( irec<0 .or. irec>nrec ) then
      write (6,'(10x,a)') 'museotimeseries error, invalid record for '//     &
                  db(1:ldb)//' database'
      write (6,'(10x,a,i4)') 'It must be in the range 1-' , nrec
      call exit(0)
    end if
    write (cfile,'(a,i4,a)') db(1:ldb) , anno , '.dat'
    if ( first ) then
      lun = -1
      call openmuseofiles(lun,cfile,0)
      oldfile = cfile
      oldrec = 0
      first = .false.
    else if ( oldfile/=cfile ) then
      close (lun)
      call openmuseofiles(lun,cfile,0)
      oldfile = cfile
      oldrec = 0
    end if
    if ( irec==0 ) then
      rewind (lun)
      write (6,'(2x,7('' Rec compn.''))')
      str1 = ' '
      do i = 1 , nrec
        read (lun) n1
        write (str2,'(a,i5,i6)') trim(str1) , i , n1
        if ( mod(i,70)==0 ) write (6,'(/2x,7('' Rec compn.''))')
        if ( mod(i,7)==0 ) then
          write (6,'(1x,a)') trim(str2)
          str1 = ' '
        else
          str1 = str2
        end if
      end do
      write (6,'(1x,a)') trim(str1)
      rewind (lun)
      oldrec = 0
      return
    end if
    do i = 1 , ncomp
      x(i) = (i-1)*dt
    end do
    if ( irec<=oldrec ) then
      if ( iflg(1)>0 ) write (6,'(10x,a,i3)') 'Rewinding '// &
         trim(cfile)//' and skipping until record ' , irec
      rewind (lun)
      do i = 1 , irec - 1
        read (lun)
      end do
    else if ( irec/=oldrec+1 ) then
      if ( iflg(1)>0 ) write (6,'(10x,a,i3,a,i3)') 'Skipping '//             &
          trim(cfile)//' from record ' ,oldrec + 1 , ' to ' , irec
      do i = oldrec + 1 , irec - 1
        read (lun)
      end do
    end if
    if ( db(1:ldb)=='Po' ) then
      n = ncomp
      read (lun) (y(i),i=1,n)
    else
      do i = 1 , ncomp
        y(i) = -9999
      end do
      read (lun) n , (y(i),i=1,n)
    end if
    oldrec = irec
  end subroutine museotimeseries
     
  subroutine chktimeseriesrange(db,anno,year1,eyear2)
    implicit none
    integer :: anno , eyear2 , year1
    character(len=*) :: db
    intent (in) anno , eyear2 , year1
    integer :: day , hour , minu , month , year2
    if ( eyear2>0 ) then
      year2 = eyear2
    else
      call whattimeisit(minu,hour,day,month,year2)
    end if
    if ( anno<year1 .or. anno>year2 ) then
      write (6,'(10x,a)') &
        'museotimeseries error reading '//trim(db)//' database'
      write (6,'(10x,a,i4,a,i4)') &
        'data are available since year ' , year1 , ' to ' , year2
      call exit(0)
    end if
  end subroutine chktimeseriesrange

  subroutine eradata(lun,ora,giorno,mese,anno,rain)
    implicit none
    integer , parameter :: nlon = 480 , nlat = 241
    integer :: anno , giorno , lun , mese , ora
    real , dimension(nlon,nlat) :: rain
    integer :: annol , i , j , mesel , nskip , slice
    real , dimension(nlat) :: work
    data slice/0/
    data mesel , annol/0 , 0/
    save :: slice,mesel,annol
    if ( lun<1 ) call getlun(lun)
    if ( mesel/=mese .or. annol/=anno ) then
      close (lun)
      call eraheader(lun,ora,giorno,mese,anno,rain,rain)
      annol = anno
      mesel = mese
      slice = 0
    end if
    nskip = (giorno-1)*24 + ora
    if ( nskip+1<=slice ) then
      rewind (lun)
      read (lun)
      read (lun)
      read (lun)
      slice = nskip + 1
    else
      nskip = nskip - slice
      slice = (giorno-1)*24 + ora + 1
    end if
    do i = 1 , nskip
      read (lun)
    end do
    read (lun) ((rain(i,j),i=1,nlon),j=nlat,1,-1)
    call erashiftlon(rain,nlon,nlat)
  end subroutine eradata
     
  real function eralon(i,j)
    implicit none
    integer :: i , j
    intent (in) i
    if (i <= 239) then
      eralon = 0.75 + 0.75*i
    else
      eralon = -180. + 0.75 * (i - 239) 
    end if
!    eralon = -180.75 + 0.75*i
  end function eralon
     
  real function eralat(i,j)
    implicit none
    integer :: i , j
    intent (in) j
    if (j <= 121) then
      eralat = 90.75 - 0.75*j
    else
      eralat = 0. - 0.75* (j - 121)
    end if
!    eralat = -90.75 + 0.75*j
  end function eralat
     
  real function eralon5(i,j)
    implicit none
    integer :: i , j
    intent (in) i
    eralon5 = -180.25 + 0.25*i
  end function eralon5
     
  real function eralat5(i,j)
    implicit none
    integer :: i , j
    intent (in) j
    eralat5 = -90.25 + 0.25*j
  end function eralat5

  real function eralonold(i,j)
    implicit none
    integer :: i , j
    intent (in) i
    eralonold = -180.75 + 0.75*i
  end function eralonold
     
  real function eralatold(i,j)
    implicit none
    integer :: i , j
    intent (in) j
    eralatold = -90.75 + 0.75*j
  end function eralatold

  real function eralon1_5(i,j)
    implicit none
    integer i,j
    eralon1_5=-181.5+1.5*i
    return
  end function eralon1_5

  real function eralat1_5(i,j)
    implicit none
    integer i,j
    eralat1_5=-91.5+1.5*j
    return
  end function eralat1_5

     
  subroutine eraheader(lun,ora,giorno,mese,anno,lat,lon)
    implicit none
    integer , parameter :: nlon = 480 , nlat = 241
    integer :: anno , giorno , lun , mese , ora
    real , dimension(nlon,nlat) :: lat , lon
    intent (out) lat
    intent (inout) lon
    character(60) :: cfile
    integer :: i , ii , isdate , j , n , ngridp , ntimes
    real :: res
    real , dimension(nlat) :: work
    ii = mm5index(ora,giorno,mese,anno)
    call integer2string(ii,8,cfile)
    if ( giorno<=monthlen(mese,anno) ) then
      if ( ora<=23 .and. ora>=0 ) then
        if ( anno<1989 .or. anno>2005 ) then
          write (6,'(a)')                                                    &
               ' MVLib: ERA data not available for the requested data.'
          write (6,'(a)') ' Requeasted data was for '//cfile(1:8)
          call exit(0)
        else
          ii = mm5index(ora,giorno,mese,anno)
          if ( lun<1 ) call getlun(lun)
          call openmuseofiles(lun,cfile(1:4)//'0100_eraint.dat',0)
          read (lun) n , n , ngridp , ntimes , isdate , res
          read (lun) ((lat(i,j),i=1,nlon),j=nlat,1,-1) ,&
                     ((lon(i,j),i=1,nlon),j=nlat,1,-1)
          call erashiftlon(lon,nlon,nlat)
          do i = 1 , nlon/2
            do j = 1 , nlat
              lon(i,j) = lon(i,j) - 360.0
            end do
          end do
    !     call worldplot(lon,nlon,nlat,'Longitudes
    !     matrix','$TITLE','dlc6p10') call
    !     worldplot(lat,nlon,nlat,'Latitudes matrix','$TITLE','dlc6p10')
          return
        end if
      end if
    end if
    write (6,'(a)')                                                         &
      ' MVLib: ERA Error reading parameter in YYMMDDHH format.'
    write (6,'(a)') ' Requested data was for '//cfile(1:8)
    call exit(0)
  end subroutine eraheader
     
  subroutine erashiftlon(mat,nlon,nlat)
    implicit none
    integer :: nlat , nlon
    real , dimension(nlon,nlat) :: mat
    real , dimension(nlat) :: work
    intent (in) nlat , nlon
    intent (inout) mat
    integer :: i , nloop
    nloop = nlon/2
    do i = 1 , nloop
      work(1:nlat) = mat(nloop+i,1:nlat)
      mat(nloop+i,1:nlat) = mat(i,1:nlat)
      mat(i,1:nlat) = work(1:nlat)
    end do
  end subroutine erashiftlon

  subroutine erashiftlat(mat,nlon,nlat)
    implicit none
    integer :: nlat , nlon
    real , dimension(nlon,nlat) :: mat
    real , dimension(nlon) :: work
    intent (in) nlat , nlon
    intent (inout) mat
    integer :: i , nloop
    nloop = nlat/2
    do i = 1 , nloop
      work(1:nlon)=mat(1:nlon,i)
      mat(1:nlon,i)=mat(1:nlon,nlat-i+1)
      mat(1:nlon,nlat-i+1)=work(1:nlon)
    end do
  end subroutine erashiftlat

  subroutine iiraflowdis(enome,ianno,numanni,port)
    implicit none
    integer , parameter :: nstaz = 54
    character(len=*) :: enome
    integer :: ianno , numanni
    real , dimension(366) :: port
    intent (out) port
    intent (inout) enome , ianno , numanni
    integer , dimension(80) :: anni
    integer :: anno , i , ian , inn , irec , lnm , lun1 , nrec
    character(40) :: nome
    character(5) :: sens
    character(44) , dimension(nstaz) :: sensore
    character(50) :: stazione
    real :: xlat , xlon
    nome = enome
    call no2space(nome)
    call noinspace(nome)
    call cv2lower(nome,nome)
    call extractint(nome,irec)
    lnm = len_trim(nome)
    inn = 0
    do i = 1 , nstaz
      call museoanagrafica('portateiira',i,stazione,sens,xlat,xlon)
      sensore(i) = stazione
      call no2space(stazione)
      call noinspace(stazione)
      call cv2lower(stazione,stazione)
      read (sens(1:3),'(i2)') numanni
      if ( nome(1:lnm)==stazione(1:lnm) ) go to 100
      if ( i==irec ) then
        enome = stazione
        call nameadjust(enome)
        call strsub(' A ',' a ',enome)
        call strsub(' Ad ',' ad ',enome)
        go to 100
      end if
      inn = inn + numanni
    end do
    write (6,'(10x,a)') 'MVLib severe error. iiraflowdis called with '//    &
                        'invalid station name'
    write (6,'(10x,a)') 'The complete list follow.'
    do i = 1 , nstaz
      write (6,'(i12,2x,a)') i , trim(sensore(i))
    end do
    call exit
     
  100  call getlun(lun1)
    call openmuseofiles(lun1,'iiraportate1930-2002.dat',0)
    do i = 1 , inn
      read (lun1)
    end do
    if ( ianno>0 ) then
      do i = 1 , numanni
        read (lun1) anno , port
        anni(i) = anno
        if ( anno==ianno ) go to 200
      end do
      write (6,'(10x,a,i5)') 'Invalid year passed to iiraflowdis: ' , ian
      write (6,'(10x,a)') 'Valid list for "'//trim(enome)//'" follows.'
      write (6,'(13i6)') (anni(i),i=1,numanni)
      call exit
    else
      nrec = abs(ianno)
      if ( nrec>=1 .and. nrec<=numanni ) then
        do i = 1 , nrec
          read (lun1) anno , port
        end do
        ianno = anno
      else
        write (6,'(10x,a,i5)') 'Invalid record passed to iiraflowdis: ' ,    &
                               abs(ianno)
        write (6,'(10x,a,i2)') &
                'Valid range for "'//trim(enome)//'" is 1-' , numanni
        call exit
      end if
    end if
  200  close (lun1)
  end subroutine iiraflowdis

  subroutine acqwascen5b(ora,giorno,mese,anno,flag,pi,la,lo,n)
    implicit none
    integer , parameter :: lonm = 630 , latm = 384
    integer :: anno , flag , giorno , mese , n , ora
    real , dimension(:) :: la , lo , pi
    intent (out) la , lo
    intent (inout) n , pi
    integer :: crec , i , j , ilog , lun , nlat , nlon , plun , prec ,       &
               pyear , tlun , trec , tyear
    character(40) , save :: dir
    character(48) :: cfile
    logical :: first
    real , dimension(lonm,latm) , save :: lat , lon , rain , temp
    !
    !*** End of declarations rewritten by SPAG
    !
    data first/.true./
    data trec , prec/ - 1 , -1/
    data tyear , pyear/ - 1 , -1/
    data tlun , plun , lun/ - 1 , -1 , -1/
    save nlon,nlat,first,trec,tyear,pyear,tlun,plun
    if ( first ) then
      call mvgetiflags(70,ilog)
      if ( ilog<=0 .or. ilog>=100 ) ilog = 6
      call mvgetiflags(71,i)
      if ( i==22 ) then
        nlon = lonm
        nlat = latm
        dir = 'acqwascen5b/Rh'
        cfile = 'acqwascen5b/latlon5brhone.dat'
      else
        nlon = 216
        nlat = 286
        dir = 'acqwascen5b/Po'
        cfile = 'acqwascen5b/latlon5bpo.dat'
      end if
      lun = -1
      call openmuseofiles(lun,cfile,0)
      read (lun) ((lat(i,j),i=1,nlon),j=1,nlat) ,  &
                 ((lon(i,j),i=1,nlon),j=1,nlat)
      close(lun)
      first = .false.
    end if
    if ( flag==1 ) then
      if ( pyear/=anno ) then
        write (ilog,'(18x,a,i4,a)') 'Opening '//trim(dir)//'Rain' ,  &
                                  anno , '.dat'
        call openmuseodb(plun,trim(dir)//'Rain',anno)
        pyear = anno
        prec = 0
      end if
      crec = index1d(giorno,mese,anno)
      if ( crec/=prec ) then
        call skipandreadgame(plun,crec,prec,0)
        read (plun) ((rain(i,j),i=1,nlon),j=1,nlat)
      end if
      n = 0
      do i = 1 , nlon
        do j = 1 , nlat
          if ( nint(rain(i,j))/=-9999 ) then
            n = n + 1
            pi(n) = rain(i,j)/24.
            if ( pi(n)<0.2 ) pi(n) = 0.0
            la(n) = lat(i,j)
            lo(n) = lon(i,j)
          end if
        end do
      end do
    else if ( flag==2 ) then
      if ( tyear/=anno ) then
        write (ilog,'(18x,a,i4,a)') 'Opening '//trim(dir)//'Temp' ,  &
                                anno , '.dat'
        call openmuseodb(tlun,trim(dir)//'Temp',anno)
        tyear = anno
        trec = 0
      end if
      n = 0
      crec = index1d(giorno,mese,anno)
      if ( crec==trec ) return
      call skipandreadgame(tlun,crec,trec,0)
      read (tlun) ((temp(i,j),i=1,nlon),j=1,nlat)
      do i = 1 , nlon
        do j = 1 , nlat
          if ( nint(temp(i,j))/=-9999 ) then
            n = n + 1
            pi(n) = temp(i,j)
            la(n) = lat(i,j)
            lo(n) = lon(i,j)
          end if
        end do
      end do
    else
     call mvliberror('acqwascen5b','Bad flag passed',flag,-9999.0)
    end if
  end subroutine acqwascen5b
     
  subroutine acqwascen6b(ora,giorno,mese,anno,flag,pi,la,lo,n)
    implicit none
    integer , parameter :: lonm = 630 , latm = 384
    integer :: anno , flag , giorno , mese , n , ora
    real , dimension(:) :: la , lo , pi
    intent (out) la , lo
    intent (inout) n , pi
    integer :: crec , i , j , ilog , lun , nlat , nlon , plun , prec ,       &
               pyear , tlun , trec , tyear
    character(40) , save :: dir
    character(48) :: cfile
    logical :: first
    real , dimension(lonm,latm) , save :: lat , lon , rain , temp
    data first/.true./
    data trec , prec/ - 1 , -1/
    data tyear , pyear/ - 1 , -1/
    data tlun , plun , lun/ - 1 , -1 , -1/
    if ( first ) then
      call mvgetiflags(70,ilog)
      if ( ilog<=0 .or. ilog>=100 ) ilog = 6
      call mvgetiflags(71,i)
      if ( i==22 ) then
        nlon = lonm
        nlat = latm
        dir = 'acqwascen6b/Rh'
        cfile = 'acqwascen6b/latlon6brhone.dat'
      else
        nlon = 216
        nlat = 286
        dir = 'acqwascen6b/Po'
        cfile = 'acqwascen6b/latlon6bpo.dat'
      end if
      lun = -1
      call openmuseofiles(lun,cfile,0)
      read (lun) ((lat(i,j),i=1,nlon),j=1,nlat) , &
                 ((lon(i,j),i=1,nlon),j=1,nlat)
      close(lun)
      first = .false.
    end if
    if ( flag==1 ) then
      if ( pyear/=anno ) then
        write (ilog,'(18x,a,i4,a)') 'Opening '//trim(dir)//'Rain' ,  &
                                  anno , '.dat'
        call openmuseodb(plun,trim(dir)//'Rain',anno)
        pyear = anno
        prec = 0
      end if
      crec = index1d(giorno,mese,anno)
      if ( crec/=prec ) then
        call skipandreadgame(plun,crec,prec,0)
        read (plun) ((rain(i,j),i=1,nlon),j=1,nlat)
      end if
      n = 0
      do i = 1 , nlon
        do j = 1 , nlat
          if ( nint(rain(i,j))/=-9999 ) then
            n = n + 1
            pi(n) = rain(i,j)/24.
            if ( pi(n)<0.2 ) pi(n) = 0.0
            la(n) = lat(i,j)
            lo(n) = lon(i,j)
          end if
        end do
      end do
    else if ( flag==2 ) then
      if ( tyear/=anno ) then
        write (ilog,'(18x,a,i4,a)') 'Opening '//trim(dir)//'Temp' ,  &
                                  anno , '.dat'
        call openmuseodb(tlun,trim(dir)//'Temp',anno)
        tyear = anno
        trec = 0
      end if
      n = 0
      crec = index1d(giorno,mese,anno)
      if ( crec==trec ) return
      call skipandreadgame(tlun,crec,trec,0)
      read (tlun) ((temp(i,j),i=1,nlon),j=1,nlat)
      do i = 1 , nlon
        do j = 1 , nlat
          if ( nint(temp(i,j))/=-9999 ) then
            n = n + 1
            pi(n) = temp(i,j)
            la(n) = lat(i,j)
            lo(n) = lon(i,j)
          end if
        end do
      end do
    else
      call mvliberror('acqwascen6b','Bad flag passed',flag,-9999.0)
    end if
  end subroutine acqwascen6b
     
  subroutine acqwascen03(ora,giorno,mese,anno,flag,pi,la,lo,n)
    implicit none
    integer , parameter :: nlon = 30 , nlat = 36
    integer :: anno , flag , giorno , mese , n , ora
    real , dimension(:) :: la , lo , pi
    intent (out) la , lo
    intent (inout) n , pi
    integer :: crec , i , j , ilog , lun , plun , prec , pyear , tlun ,      &
               trec , tyear
    character(48) :: cfile
    logical :: first
    real , dimension(nlon,nlat) , save :: lat , lon , mat
    data first/.true./
    data trec , prec/ - 1 , -1/
    data tyear , pyear/ - 1 , -1/
    data tlun , plun , lun/ - 1 , -1 , -1/
    if ( first ) then
      call mvgetiflags(70,ilog)
      if ( ilog<=0 .or. ilog>=100 ) ilog = 6
      call getlun(lun)
      call openmuseofiles(lun,'acqwascen03/latlon.dat',0)
      read (lun) i , j
      if ( i/=nlon .or. j/=nlat ) &
        call mvliberror('acqwascen03','Bad nlon/nlat parameters',nlon,   &
        float(nlat))
      read (lun) lat , lon
      close (lun)
      first = .false.
    end if
    if ( flag==1 ) then
      if ( pyear/=anno ) then
        write (cfile,'(a,i4,a)') 'museo/acqwascen03/rain' , anno , '.dat'
        write (ilog,'(18x,a)') 'Opening '//trim(cfile)
        call openmuseodb(plun,'acqwascen03/rain',anno)
        pyear = anno
        prec = 0
      end if
      crec = index1h(ora,giorno,mese,anno)
      call skipandreadgame(plun,crec,prec,0)
      read (plun) mat
      n = 0
      do i = 1 , nlon
        do j = 1 , nlat
          n = n + 1
          pi(n) = mat(i,j)
          if ( pi(n)<1.0 ) pi(n) = 0.0
          la(n) = lat(i,j)
          lo(n) = lon(i,j)
        end do
      end do
    else if ( flag==2 ) then
      if ( tyear/=anno ) then
        write (cfile,'(a,i4,a)') 'museo/acqwascen03/temp' , anno , '.dat'
        write (ilog,'(18x,a)') 'Opening '//trim(cfile)
        call openmuseodb(tlun,'acqwascen03/temp',anno)
        tyear = anno
        trec = 0
      end if
      crec = index1h(ora,giorno,mese,anno)
      call skipandreadgame(tlun,crec,trec,0)
      read (tlun) mat
      n = 0
      do i = 1 , nlon
        do j = 1 , nlat
          n = n + 1
          pi(n) = mat(i,j)
          la(n) = lat(i,j)
          lo(n) = lon(i,j)
        end do
      end do
    else
      call mvliberror('acqwascen03','Bad flag passed',flag,-9999.0)
    end if
  end subroutine acqwascen03
     
  subroutine acqwascen5a(ora,giorno,mese,anno,flag,pi,la,lo,n)
    implicit none
    integer , parameter :: nlon = 70 , nlat = 40
    integer :: anno , flag , giorno , mese , n , ora
    real , dimension(:) :: la , lo , pi
    intent (out) la , lo
    intent (inout) n , pi
    integer :: crec , i , j , ilog , lun , plun , prec , pyear , tlun ,      &
               trec , tyear
    logical :: first
    integer , dimension(nlon,nlat) , save :: mask
    character(22) :: ms
    real , dimension(nlon,nlat) , save :: rain , temp
    data first/.true./
    data trec , prec/ - 1 , -1/
    data tyear , pyear/ - 1 , -1/
    data tlun , plun , lun/ - 1 , -1 , -1/
    data ms/'acqwa Scen 5a database'/
    if ( first ) then
      call mvgetiflags(70,ilog)
      if ( ilog<=0 .or. ilog>=100 ) ilog = 6
      call getlun(lun)
      call openmuseofiles(lun,'acqwascen5a/mask.dat',0)
      read (lun) mask
      close (lun)
      first = .false.
    end if
    if ( flag==1 ) then
      if ( pyear/=anno ) then
        call openmuseodb(plun,'acqwascen5a/rain',anno)
        pyear = anno
        prec = 0
      end if
      crec = index3h(ora,giorno,mese,anno)
      if ( crec/=prec ) then
        call skipandreadgame(plun,crec,prec,0)
        read (plun) rain
      end if
      n = 0
      do i = 1 , nlon
        do j = 1 , nlat
          if ( mask(i,j)==1 ) then
            n = n + 1
            pi(n) = rain(i,j)
            if ( pi(n)<1.0 ) pi(n) = 0.0
            la(n) = 41.5 + (j-1)*0.25
            lo(n) = 3.0 + (i-1)*0.25
          end if
        end do
      end do
    else if ( flag==2 ) then
      if ( tyear/=anno ) then
        call openmuseodb(tlun,'acqwascen5a/temp',anno)
        tyear = anno
        trec = 0
      end if
      crec = index3h(ora,giorno,mese,anno)
      n = 0
      if ( crec/=trec ) then
        call skipandreadgame(tlun,crec,trec,0)
        read (tlun) temp
      else
        return
      end if
      do i = 1 , nlon
        do j = 1 , nlat
          if ( mask(i,j)==1 ) then
            n = n + 1
            pi(n) = temp(i,j)
            la(n) = 41.5 + (j-1)*0.25
            lo(n) = 3.0 + (i-1)*0.25
          end if
        end do
      end do
    else
      call mvliberror('acqwascen5a','Bad flag passed',flag,-9999.0)
    end if
  end subroutine acqwascen5a
     
  subroutine acqwascen6a(ora,giorno,mese,anno,flag,pi,la,lo,ndat)
    implicit none
    integer :: anno , flag , giorno , mese , ndat , ora
    real , dimension(:) :: la , lo , pi
    if ( flag==1 ) then
      call acqwascen6ap(ora,giorno,mese,anno,pi,la,lo,ndat)
    else if ( flag==2 ) then
      call acqwascen6at(ora,giorno,mese,anno,pi,la,lo,ndat)
    else
      call mvliberror('acqwascen6a','Bad flag passed',flag,-9999.0)
    end if
  end subroutine acqwascen6a
     
  subroutine acqwascen6at(ora,giorno,mese,anno,pi,la,lo,ndat)
    implicit none
    integer , parameter :: nlon = 70 , nlat = 40
    integer :: anno , giorno , mese , ndat , ora
    real , dimension(:) :: la , lo , pi
    intent (out) la , lo , pi
    intent (inout) ndat
    character(64) :: cfile , dire , ofile
    character(4) :: field
    logical :: first
    integer :: i , j , lun , oldrec , irec
    real , dimension(nlon,nlat) , save :: lat , lon , map
    data field/'temp'/
    data first/.true./
    data ofile/'nofile'/
    data dire/'museo/acqwascen6a/'/
    data oldrec/ - 1/
    if ( anno>2050 .or. anno<1951 ) then
      ndat = 0
      return
    end if
    if ( first ) then
      do i = 1 , nlon
        do j = 1 , nlat
          lat(i,j) = (j-1)*0.25 + 41.5
          lon(i,j) = (i-1)*0.25 + 3.0
        end do
      end do
      call getlun(lun)
      first = .false.
    end if
    write (cfile,'(a,i4,a)') trim(dire)//field(1:4) , anno , '.dat'
    if ( trim(cfile)/=trim(ofile) ) then
      close (lun)
      open (lun,file=cfile,status='old',form='unformatted')
      ofile = cfile
      oldrec = 0
    end if
    irec = index3h(ora,giorno,mese,anno)
    if ( irec==oldrec ) then
      ndat = 0
      return
    else if ( irec==oldrec+1 ) then
      read (lun) map
    else
      rewind (lun)
      do i = 1 , irec - 1
        read (lun)
      end do
      read (lun) map
    end if
    oldrec = irec
    ndat = 0
    do i = 1 , nlon
      do j = 1 , nlat
        if ( nint(map(i,j))>-9000 ) then
          ndat = ndat + 1
          pi(ndat) = map(i,j)
          la(ndat) = lat(i,j)
          lo(ndat) = lon(i,j)
        end if
      end do
    end do
  end subroutine acqwascen6at

  subroutine acqwascen6ap(ora,giorno,mese,anno,pi,la,lo,ndat)
    implicit none
    integer , parameter :: nlon = 70 , nlat = 40
    integer :: anno , giorno , mese , ndat , ora
    real , dimension(:) :: la , lo , pi
    intent (out) la , lo
    intent (inout) ndat , pi
    character(64) :: cfile , dire , ofile
    character(4) :: field
    logical :: first
    integer :: i , j , lun , oldrec , irec
    real , dimension(nlon,nlat) , save :: lat , lon , map
    data field/'rain'/
    data first/.true./
    data ofile/'nofile'/
    data dire/'museo/acqwascen6a/'/
    data oldrec/ - 1/
    if ( anno>2050 .or. anno<1951 ) then
      ndat = 0
      return
    end if
    if ( first ) then
      do i = 1 , nlon
        do j = 1 , nlat
          lat(i,j) = (j-1)*0.25 + 41.5
          lon(i,j) = (i-1)*0.25 + 3.0
        end do
      end do
      call getlun(lun)
      first = .false.
    end if
    write (cfile,'(a,i4,a)') trim(dire)//field(1:4) , anno ,'.dat'
    if ( trim(cfile)/=trim(ofile) ) then
      close (lun)
      open (lun,file=cfile,status='old',form='unformatted')
      ofile = cfile
      oldrec = 0
    end if
    irec = index3h(ora,giorno,mese,anno)
    if ( irec==oldrec ) then
    else if ( irec==oldrec+1 ) then
      read (lun) map
    else
      rewind (lun)
      do i = 1 , irec - 1
        read (lun)
      end do
      read (lun) map
    end if
    oldrec = irec
    ndat = 0
    do i = 1 , nlon
      do j = 1 , nlat
        if ( nint(map(i,j))>-9000 ) then
          ndat = ndat + 1
          pi(ndat) = map(i,j)/3.0
          if ( pi(ndat)<1.0 ) pi(ndat) = 0.0
          la(ndat) = lat(i,j)
          lo(ndat) = lon(i,j)
        end if
      end do
    end do
  end subroutine acqwascen6ap
     
  subroutine regcmgetdims(ora,giorno,mese,anno,nlon,nlat,simname)
    !
    ! Legge le dimensioni nlat,nlon (iy,jx) da un file regcm
    !
    implicit none
    include 'netcdf.inc'
    integer :: anno , giorno , mese , nlat , nlon , ora
    integer :: anno2 , mese2
    intent (in) anno , mese
    character(300) :: cfl
    character(132) :: dir
    character(:), allocatable :: cname
    character(:), allocatable :: cvar
    character(30) :: cdat
    character(*) :: simname
    integer :: istatus , iyd , jxd , ld , ilog , lun
    data lun/ - 1/
!    ilog = 6
    if ( trim(schym(2)) == 'AFR44-ERAI') then
      if (anno >= 1979 .and. anno <= 1980) then
        write (cfl,'(a,a)') trim(simname)//'1979010112-1980123112.nc'
      else if (anno >= 1981 .and. anno <=1985) then
        write (cfl,'(a,a)') trim(simname)//'1981010112-1985123112.nc'
      else if (anno >= 1986 .and. anno <=1990) then
        write (cfl,'(a,a)') trim(simname)//'1986010112-1990123112.nc'
      else if (anno >= 1991 .and. anno <=1995) then
        write (cfl,'(a,a)') trim(simname)//'1991010112-1995123112.nc'
      else if (anno >= 1996 .and. anno <=2000) then
        write (cfl,'(a,a)') trim(simname)//'1996010112-2000123112.nc'
      else if (anno >= 2001 .and. anno <=2005) then
        write (cfl,'(a,a)') trim(simname)//'2001010112-2005123112.nc'
      else if (anno >= 2006 .and. anno <=2009) then
        write (cfl,'(a,a)') trim(simname)//'2006010112-2009123112.nc'
      end if
    else if ( trim(schym(2)) == 'SAM44-ERAI') then
      if (anno >= 1979 .and. anno <= 1980) then
        write (cfl,'(a,a)') trim(simname)//'1979010112-1980123112.nc'
      else if (anno >= 1981 .and. anno <=1985) then
        write (cfl,'(a,a)') trim(simname)//'1981010112-1985123112.nc'
      else if (anno >= 1986 .and. anno <=1990) then
        write (cfl,'(a,a)') trim(simname)//'1986010112-1990123112.nc'
      else if (anno >= 1991 .and. anno <=1995) then
        write (cfl,'(a,a)') trim(simname)//'1991010112-1995123112.nc'
      else if (anno >= 1996 .and. anno <=2000) then
        write (cfl,'(a,a)') trim(simname)//'1996010112-2000123112.nc'
      else if (anno >= 2001 .and. anno <=2005) then
        write (cfl,'(a,a)') trim(simname)//'2001010112-2005123112.nc'
      else if (anno >= 2006 .and. anno <=2008) then
        write (cfl,'(a,a)') trim(simname)//'2006010112-2008123112.nc'
      end if
    else if ( trim(schym(2)) == 'MED44-ERAIv4') then
      if (anno >= 1979 .and. anno <= 1980) then
        write (cfl,'(a,a)') trim(simname)//'1979010112-1980123112.nc'
      else if (anno >= 1981 .and. anno <=1985) then
        write (cfl,'(a,a)') trim(simname)//'1981010112-1985123112.nc'
      else if (anno >= 1986 .and. anno <=1990) then
        write (cfl,'(a,a)') trim(simname)//'1986010112-1990123112.nc'
      else if (anno >= 1991 .and. anno <=1995) then
        write (cfl,'(a,a)') trim(simname)//'1991010112-1995123112.nc'
      else if (anno >= 1996 .and. anno <=2000) then
        write (cfl,'(a,a)') trim(simname)//'1996010112-2000123112.nc'
      else if (anno >= 2001 .and. anno <=2005) then
        write (cfl,'(a,a)') trim(simname)//'2001010112-2005123112.nc'
      else if (anno >= 2006 .and. anno <=2010) then
        write (cfl,'(a,a)') trim(simname)//'2006010112-2010123112.nc'
      else if (anno >= 2011 .and. anno <=2012) then
        write (cfl,'(a,a)') trim(simname)//'2011010112-2012123112.nc'
      end if
    else if ( trim(schym(2)) == 'MED11-ERAI') then
      if (anno >= 1979 .and. anno <= 1980) then
        write (cfl,'(a,a)') trim(simname)//'1979010112-1980123112.nc'
      else if (anno >= 1981 .and. anno <=1985) then
        write (cfl,'(a,a)') trim(simname)//'1981010112-1985123112.nc'
      else if (anno >= 1986 .and. anno <=1990) then
        write (cfl,'(a,a)') trim(simname)//'1986010112-1990123112.nc'
      else if (anno >= 1991 .and. anno <=1995) then
        write (cfl,'(a,a)') trim(simname)//'1991010112-1995123112.nc'
      else if (anno >= 1996 .and. anno <=2000) then
        write (cfl,'(a,a)') trim(simname)//'1996010112-2000123112.nc'
      else if (anno >= 2001 .and. anno <=2005) then
        write (cfl,'(a,a)') trim(simname)//'2001010112-2005123112.nc'
      else if (anno >= 2006 .and. anno <=2010) then
        write (cfl,'(a,a)') trim(simname)//'2006010112-2010123112.nc'
      end if
    else if ( trim(schym(2)) == 'EUR44-ERAI') then
      if (anno >= 1979 .and. anno <= 1980) then
        write (cfl,'(a,a)') trim(simname)//'1979010112-1980123112.nc'
      else if (anno >= 1981 .and. anno <=1985) then
        write (cfl,'(a,a)') trim(simname)//'1981010112-1985123112.nc'
      else if (anno >= 1986 .and. anno <=1990) then
        write (cfl,'(a,a)') trim(simname)//'1986010112-1990123112.nc'
      else if (anno >= 1991 .and. anno <=1995) then
        write (cfl,'(a,a)') trim(simname)//'1991010112-1995123112.nc'
      else if (anno >= 1996 .and. anno <=2000) then
        write (cfl,'(a,a)') trim(simname)//'1996010112-2000123112.nc'
      else if (anno >= 2001 .and. anno <=2005) then
        write (cfl,'(a,a)') trim(simname)//'2001010112-2005123112.nc'
      else if (anno >= 2006 .and. anno <=2010) then
        write (cfl,'(a,a)') trim(simname)//'2006010112-2010123112.nc'
      else if (anno >= 2011 .and. anno <=2012) then
        write (cfl,'(a,a)') trim(simname)//'2011010112-2012123112.nc'
      end if
    else if ( trim(schym(2)) == 'EUR44-MOHC-his') then
      if (anno >= 1950 .and. anno <=1950) then
        write (cfl,'(a,a)') trim(simname)//'1950010112-1950123012.nc'
      else if (anno >= 1951 .and. anno <=1955) then
        write (cfl,'(a,a)') trim(simname)//'1951010112-1955123012.nc'
      else if (anno >= 1956 .and. anno <=1960) then
        write (cfl,'(a,a)') trim(simname)//'1956010112-1960123012.nc'
      else if (anno >= 1961 .and. anno <=1965) then
        write (cfl,'(a,a)') trim(simname)//'1961010112-1965123012.nc'
      else if (anno >= 1966 .and. anno <=1970) then
        write (cfl,'(a,a)') trim(simname)//'1966010112-1970123012.nc'
      else if (anno >= 1971 .and. anno <=1975) then
        write (cfl,'(a,a)') trim(simname)//'1971010112-1975123012.nc'
      else if (anno >= 1976 .and. anno <=1980) then
        write (cfl,'(a,a)') trim(simname)//'1976010112-1980123012.nc'
      else if (anno >= 1981 .and. anno <=1985) then
        write (cfl,'(a,a)') trim(simname)//'1981010112-1985123012.nc'
      else if (anno >= 1986 .and. anno <=1990) then
        write (cfl,'(a,a)') trim(simname)//'1986010112-1990123012.nc'
      else if (anno >= 1991 .and. anno <=1995) then
        write (cfl,'(a,a)') trim(simname)//'1991010112-1995123012.nc'
      else if (anno >= 1996 .and. anno <=2000) then
        write (cfl,'(a,a)') trim(simname)//'1996010112-2000123012.nc'
      else if (anno >= 2001 .and. anno <=2005) then
        write (cfl,'(a,a)') trim(simname)//'2001010112-2005123012.nc'
      end if
    else if ( trim(schym(2)) == 'CAM44-MPI-his') then
      if (anno >= 1970 .and. anno <=1970) then
        write (cfl,'(a,a)') trim(simname)//'1970010112-1970123112.nc'
      else if (anno >= 1971 .and. anno <=1975) then
        write (cfl,'(a,a)') trim(simname)//'1971010112-1975123112.nc'
      else if (anno >= 1976 .and. anno <=1980) then
        write (cfl,'(a,a)') trim(simname)//'1976010112-1980123112.nc'
      else if (anno >= 1981 .and. anno <=1985) then
        write (cfl,'(a,a)') trim(simname)//'1981010112-1985123112.nc'
      else if (anno >= 1986 .and. anno <=1990) then
        write (cfl,'(a,a)') trim(simname)//'1986010112-1990123112.nc'
      else if (anno >= 1991 .and. anno <=1995) then
        write (cfl,'(a,a)') trim(simname)//'1991010112-1995123112.nc'
      else if (anno >= 1996 .and. anno <=2000) then
        write (cfl,'(a,a)') trim(simname)//'1996010112-2000123112.nc'
      else if (anno >= 2001 .and. anno <=2005) then
        write (cfl,'(a,a)') trim(simname)//'2001010112-2005123112.nc'
      end if
    else if ( trim(schym(2)) == 'CAM44-MPI-rcp85') then
      if (anno >= 2006 .and. anno <=2010) then
        write (cfl,'(a,a)') trim(simname)//'2006010112-2010123112.nc'
      else if (anno >= 2011 .and. anno <=2015) then
        write (cfl,'(a,a)') trim(simname)//'2011010112-2015123112.nc'
      else if (anno >= 2016 .and. anno <=2020) then
        write (cfl,'(a,a)') trim(simname)//'2016010112-2020123112.nc'
      else if (anno >= 2021 .and. anno <=2025) then
        write (cfl,'(a,a)') trim(simname)//'2021010112-2025123112.nc'
      else if (anno >= 2026 .and. anno <=2030) then
        write (cfl,'(a,a)') trim(simname)//'2026010112-2030123112.nc'
      else if (anno >= 2031 .and. anno <=2035) then
        write (cfl,'(a,a)') trim(simname)//'2031010112-2035123112.nc'
      else if (anno >= 2036 .and. anno <=2040) then
        write (cfl,'(a,a)') trim(simname)//'2036010112-2040123112.nc'
      else if (anno >= 2041 .and. anno <=2045) then
        write (cfl,'(a,a)') trim(simname)//'2041010112-2045123112.nc'
      else if (anno >= 2046 .and. anno <=2050) then
        write (cfl,'(a,a)') trim(simname)//'2046010112-2050123112.nc'
      else if (anno >= 2051 .and. anno <=2055) then
        write (cfl,'(a,a)') trim(simname)//'2051010112-2055123112.nc'
      else if (anno >= 2056 .and. anno <=2060) then
        write (cfl,'(a,a)') trim(simname)//'2056010112-2060123112.nc'
      else if (anno >= 2061 .and. anno <=2065) then
        write (cfl,'(a,a)') trim(simname)//'2061010112-2065123112.nc'
      else if (anno >= 2066 .and. anno <=2070) then
        write (cfl,'(a,a)') trim(simname)//'2066010112-2070123112.nc'
      else if (anno >= 2071 .and. anno <=2075) then
        write (cfl,'(a,a)') trim(simname)//'2071010112-2075123112.nc'
      else if (anno >= 2076 .and. anno <=2080) then
        write (cfl,'(a,a)') trim(simname)//'2076010112-2080123112.nc'
      else if (anno >= 2081 .and. anno <=2085) then
        write (cfl,'(a,a)') trim(simname)//'2081010112-2085123112.nc'
      else if (anno >= 2086 .and. anno <=2090) then
        write (cfl,'(a,a)') trim(simname)//'2086010112-2090123112.nc'
      else if (anno >= 2091 .and. anno <=2095) then
        write (cfl,'(a,a)') trim(simname)//'2091010112-2095123112.nc'
      else if (anno >= 2096 .and. anno <=2099) then
        write (cfl,'(a,a)') trim(simname)//'2096010112-2099123012.nc'
      end if
    else if ( trim(schym(2)) == 'AFR44-MPI-his' .or. trim(schym(2)) == 'AFR44-MPIMR-his') then
      if (anno >= 1970 .and. anno <=1970) then
        write (cfl,'(a,a)') trim(simname)//'19700101-19701231.nc'
      else if (anno >= 1971 .and. anno <=1975) then
        write (cfl,'(a,a)') trim(simname)//'19710101-19751231.nc'
      else if (anno >= 1976 .and. anno <=1980) then
        write (cfl,'(a,a)') trim(simname)//'19760101-19801231.nc'
      else if (anno >= 1981 .and. anno <=1985) then
        write (cfl,'(a,a)') trim(simname)//'19810101-19851231.nc'
      else if (anno >= 1986 .and. anno <=1990) then
        write (cfl,'(a,a)') trim(simname)//'19860101-19901231.nc'
      else if (anno >= 1991 .and. anno <=1995) then
        write (cfl,'(a,a)') trim(simname)//'19910101-19951231.nc'
      else if (anno >= 1996 .and. anno <=2000) then
        write (cfl,'(a,a)') trim(simname)//'19960101-20001231.nc'
      else if (anno >= 2001 .and. anno <=2005) then
        write (cfl,'(a,a)') trim(simname)//'20010101-20051130.nc'
      end if
    else if ( trim(schym(2)) == 'AFR44-MPI-rcp85' .or. trim(schym(2)) == 'AFR44-MPIMR-rcp85') then
      if (anno >= 2007 .and. anno <=2010) then
        write (cfl,'(a,a)') trim(simname)//'20070101-20101231.nc'
      else if (anno >= 2011 .and. anno <=2015) then
        write (cfl,'(a,a)') trim(simname)//'20110101-20151231.nc'
      else if (anno >= 2016 .and. anno <=2020) then
        write (cfl,'(a,a)') trim(simname)//'20160101-20201231.nc'
      else if (anno >= 2021 .and. anno <=2025) then
        write (cfl,'(a,a)') trim(simname)//'20210101-20251231.nc'
      else if (anno >= 2026 .and. anno <=2030) then
        write (cfl,'(a,a)') trim(simname)//'20260101-20301231.nc'
      else if (anno >= 2031 .and. anno <=2035) then
        write (cfl,'(a,a)') trim(simname)//'20310101-20351231.nc'
      else if (anno >= 2036 .and. anno <=2040) then
        write (cfl,'(a,a)') trim(simname)//'20360101-20401231.nc'
      else if (anno >= 2041 .and. anno <=2045) then
        write (cfl,'(a,a)') trim(simname)//'20410101-20451231.nc'
      else if (anno >= 2046 .and. anno <=2050) then
        write (cfl,'(a,a)') trim(simname)//'20460101-20501231.nc'
      else if (anno >= 2051 .and. anno <=2055) then
        write (cfl,'(a,a)') trim(simname)//'20510101-20551231.nc'
      else if (anno >= 2056 .and. anno <=2060) then
        write (cfl,'(a,a)') trim(simname)//'20560101-20601231.nc'
      else if (anno >= 2061 .and. anno <=2065) then
        write (cfl,'(a,a)') trim(simname)//'20610101-20651231.nc'
      else if (anno >= 2066 .and. anno <=2070) then
        write (cfl,'(a,a)') trim(simname)//'20660101-20701231.nc'
      else if (anno >= 2071 .and. anno <=2075) then
        write (cfl,'(a,a)') trim(simname)//'20710101-20751231.nc'
      else if (anno >= 2076 .and. anno <=2080) then
        write (cfl,'(a,a)') trim(simname)//'20760101-20801231.nc'
      else if (anno >= 2081 .and. anno <=2085) then
        write (cfl,'(a,a)') trim(simname)//'20810101-20851231.nc'
      else if (anno >= 2086 .and. anno <=2090) then
        write (cfl,'(a,a)') trim(simname)//'20860101-20901231.nc'
      else if (anno >= 2091 .and. anno <=2095) then
        write (cfl,'(a,a)') trim(simname)//'20910101-20951231.nc'
      else if (anno >= 2096 .and. anno <=2099) then
        write (cfl,'(a,a)') trim(simname)//'20960101-20991231.nc'
      end if
    else if ( trim(schym(2)(1:6)) == 'CORDEX') then
    !else if ( trim(schym(2)(1:6)) == 'CORDEX_EUR-11_HadGEM_hist_MARCONI') then
!      print*,"SIAMO ENTRATI QUI"
!      anno2 = anno ; mese2 = mese + 1
!      if ( mese == 12 ) then
!         anno2=anno+1 ; mese2 = 01
!      end if
!      write (cdat,'(i6,a,i6,a)') anno*100 + mese ,'010300-',anno2*100 + mese2, &
!             '010000.nc'
!      print*,"CDAT",cdat
!      cvar = "tas"
!      call cordexnamefile(trim(schym(16)),cdat,cvar,cname)
!      cfl = cname 
      cfl = simname
    else if ( trim(schym(2)) == 'standard' ) then
     call chymregcmoutput(dir)
     write (cfl,'(a,i6,a)') trim(schym(16))//'SRF.' , anno*100 + mese , '0100.nc'
    end if
!    print*,trim(cfl)
!    print*,trim(cfl)
!    ld = len_trim(dir)
!    call mvgetiflags(70,ilog)
!    if ( ilog<=0 .or. ilog>=100 ) ilog = 6
!    write (ilog,'(10x,a)') 'Acquiring RegCM parameter from the file: '
!    write (ilog,'(10x,a)') dir(1:ld)//trim(cfl)
!    call openncregcmfile(lun,dir(1:ld)//cfl,ilog)
    call openncregcmfile(lun,trim(cfl),ilog)
!    istatus = nf_inq_dimid(lun,'iy',iyd)
!    if ( istatus/=nf_noerr ) then
!      write (ilog,'(16x,a)') &
!        '  ERROR Reading IY dimension from '//dir(1:ld)//cfl
!      stop
!    end if
    istatus = nf_inq_dimid(lun,'y',iyd)
    if ( istatus/=nf_noerr ) then
      istatus = nf_inq_dimid(lun,'iy',iyd)
      if ( istatus/=nf_noerr ) then
        write (ilog,'(16x,a)') &
          '  ERROR Reading IY dimension from '//dir(1:ld)//cfl
        stop
      end if
    end if
    istatus = nf_inq_dimlen(lun,iyd,nlat)
    if ( istatus/=nf_noerr ) then
      write (ilog,'(16x,a)') &
        '  ERROR Reading IY dimension from '//dir(1:ld)//cfl
      stop
    end if
!    istatus = nf_inq_dimid(lun,'jx',jxd)
!    if ( istatus/=nf_noerr ) then
!      write (ilog,'(16x,a)') &
!        '  ERROR Reading JX dimension from '//dir(1:ld)//cfl
!      stop
!    end if
!    istatus = nf_inq_dimlen(lun,jxd,nlon)
!    if ( istatus/=nf_noerr ) then
!      write (ilog,'(16x,a)') &
!        '  ERROR Reading JX dimension from '//dir(1:ld)//cfl
!      stop
!    end if
    istatus = nf_inq_dimid(lun,'x',jxd)
    if ( istatus/=nf_noerr ) then
      istatus = nf_inq_dimid(lun,'jx',jxd)
      if ( istatus/=nf_noerr ) then
        write (ilog,'(16x,a)') &
          '  ERROR Reading JX dimension from '//dir(1:ld)//cfl
        stop
      end if
    end if
    istatus = nf_inq_dimlen(lun,jxd,nlon)
    if ( istatus/=nf_noerr ) then
      write (ilog,'(16x,a)') &
        '  ERROR Reading JX dimension from '//dir(1:ld)//cfl
      stop
    end if
    istatus = nf_close(lun)
    if ( istatus/=nf_noerr ) then
      write (ilog,'(16x,a,i6)') '  Error closing file unit ' , lun
      stop
    end if
  end subroutine regcmgetdims

  subroutine cordexnamefile(cdir,cdat,cvar,cfile)
    implicit none
    character(:), allocatable :: cdom, csce, cens, cmod, cver, cfre, cglo
    character(*), intent(in) :: cdir,cdat,cvar
    character(:), allocatable :: cdor
    character(:), allocatable, intent(out) :: cfile
    integer :: ind,indlast,i
    integer, dimension(20) :: indarr
    
    i = 1
    ind = 1
    cdor = trim(cdir)
    indarr = 0
    do while (ind /= 0) 
      indlast = ind
      ind = index( cdor, '/' )
      if (i == 1) then
      indarr(i) = ind
      else
      indarr(i) = indarr(i-1) + ind-1
      end if
      cdor = 'R'//cdor(ind+1:)
      i = i + 1
    end do
    cfre = cdir(indarr(i-3)+1:indarr(i-2)-1)
    cver = cdir(indarr(i-4)+1:indarr(i-3)-1)
    cmod = cdir(indarr(i-5)+1:indarr(i-4)-1)
    cens = cdir(indarr(i-6)+1:indarr(i-5)-1)
    csce = cdir(indarr(i-7)+1:indarr(i-6)-1)
    cglo = cdir(indarr(i-8)+1:indarr(i-7)-1)
    cdom = cdir(indarr(i-10)+1:indarr(i-9)-1)
    cfile = trim(cdir)//trim(cvar)//"/"//trim(cvar)//"_"//cdom//"_"//cglo//"_"//csce//"_"//cens//"_" &
       //cmod//"_"//cver//"_"//cfre//"_"//trim(cdat)
  end subroutine cordexnamefile
  
     
  subroutine acqwascenc1(ora,giorno,mese,anno,flag,pi,la,lo,n)
    implicit none
    integer , parameter :: nlon = 214 , nlat = 198
    integer :: anno , flag , giorno , mese , n , ora
    real , dimension(:) :: la , lo , pi
    intent (in) anno , flag , giorno , mese , ora
    intent (out) la , lo
    intent (inout) n , pi
    character(80) :: cfl , dir , pofl , tofl
    integer :: crec , i , j , lcfl , ld , ilog , lun , plun , prec , tlun ,  &
               trec
    logical :: first
    real , dimension(nlon,nlat) :: lat , lon , rain , temp
    data first/.true./
    data dir/'museo/acqwascenc1/'/
    data tofl , pofl/'nofile' , 'nofile'/
    data trec , prec/ - 1 , -1/
    data tlun , plun , lun/ - 1 , -1 , -1/
    if ( first ) then
      ld = len_trim(dir)
      call mvgetiflags(70,ilog)
      if ( ilog<=0 .or. ilog>=100 ) ilog = 6
      call getlun(lun)
      open (lun,file=dir(1:ld)//'latlongrid.dat',status='old',               &
            form='unformatted')
      read (lun) lat , lon
      close (lun)
      first = .false.
    end if
    if ( flag==1 ) then
      write (cfl,'(a,i6,a)') 'prec/SRF.' , anno*100 + mese , '0100-pre'
      lcfl = len_trim(cfl)
      if ( cfl(1:lcfl)/=pofl(1:lcfl) ) then
        call openregcmfile(plun,dir(1:ld)//cfl(1:lcfl),nlon*nlat*4,ilog)
        pofl = cfl
        prec = -1
      end if
      crec = (giorno-1)*8 + ora/3 + 1
      if ( crec/=prec ) then
        write (ilog,'(16x,a,i3,a)') '  Reading rain record ' , crec ,         &
                                    ' from RegCM C1 scenario'
        read (plun,rec=crec) rain
        prec = crec
      end if
      n = 0
      do i = 1 , nlon
        do j = 1 , nlat
          n = n + 1
          pi(n) = rain(i,j)/24.0
          if ( pi(n)<1.0 ) pi(n) = 0.0
          la(n) = lat(i,j)
          lo(n) = lon(i,j)
        end do
      end do
    else if ( flag==2 ) then
      write (cfl,'(a,i6,a)') 'temp/SRF.' , anno*100 + mese , '0100-tm'
      lcfl = len_trim(cfl)
      if ( cfl(1:lcfl)/=tofl(1:lcfl) ) then
        call openregcmfile(tlun,dir(1:ld)//cfl(1:lcfl),nlon*nlat*4,ilog)
        tofl = cfl
        trec = -1
      end if
      n = 0
      crec = (giorno-1)*8 + ora/3 + 1
      if ( crec/=trec ) then
        write (ilog,'(16x,a,i3,a)') '  Reading temp record ' , crec ,         &
                                    ' from RegCM C1 scenario'
        read (tlun,rec=crec) temp
        trec = crec
      else
        return
      end if
      do i = 1 , nlon
        do j = 1 , nlat
          n = n + 1
          pi(n) = temp(i,j) - 273.15
          la(n) = lat(i,j)
          lo(n) = lon(i,j)
        end do
      end do
    ! call chymerror(8,flag,0.0,'acqwascenc1')
    end if
  end subroutine acqwascenc1
 
  subroutine openregcmfile(lun,cfile,irecl,logfile)
    implicit none
    character(len=*) :: cfile
    integer :: logfile , lun , irecl
    intent (in) logfile , irecl
    if ( lun<0 ) then
      call getlun(lun)
    else
      close (lun)
    end if
    write (logfile,'(16x,a,i2)') '  Opening '//trim(cfile)  &
                               //' on unit ' , lun
    open (lun,file=cfile,status='old',access='direct',action='read',         &
          form='unformatted',convert='big_endian',recl=irecl)
  end subroutine openregcmfile
 
  subroutine openncregcmfile(lun,cfile,logfile)
    implicit none
    include 'netcdf.inc'
    character(len=*) :: cfile
    integer :: logfile , lun
    intent (in) logfile
    integer :: istatus
    istatus = nf_open(cfile,nf_nowrite,lun)
    if ( istatus/=nf_noerr ) then
      write (logfile,'(16x,a)') '  Error opening file '//trim(cfile)
      stop
    end if
    return
  end subroutine openncregcmfile
 
  subroutine acqwascen01(ora,giorno,mese,anno,flag,pi,la,lo,n)
    implicit none
    integer , parameter :: nlon = 18 , nlat = 18
    integer :: anno , flag , giorno , mese , n , ora
    real , dimension(:) :: la , lo , pi
    intent (in) flag
    intent (out) la , lo
    intent (inout) n , pi
    integer :: annop , annot , i , j , ld , ilog , lun , plun , prec , tlun ,&
           trec
    character(80) :: cfl , dir
    real , dimension(nlon,nlat) :: field , lat , lon
    logical :: first
    data first/.true./
    data dir/'museo/acqwascen01/'/
    data annot , annop/ - 1 , -1/
    data trec , prec/ - 1 , -1/
    data tlun , plun , lun/ - 1 , -1 , -1/
    if ( first ) then
      ld = len_trim(dir)
      call mvgetiflags(70,ilog)
      if ( ilog<=0 .or. ilog>=100 ) ilog = 6
      call getlun(lun)
      open (lun,file=dir(1:ld)//'geofile',status='old',form='unformatted')
      read (lun) lat
      read (lun) lon
      close (lun)
      first = .false.
    end if
    if ( flag==1 ) then
      if ( anno/=annop ) then
        if ( plun>0 ) then
          close (plun)
        else
          call getlun(plun)
        end if
        write (cfl,'(a,i4)') dir(1:ld)//'prec/prec_' , anno
        write (ilog,'(18x,a)') 'Opening '//trim(cfl)
        open (plun,file=cfl,status='old',form='unformatted')
        prec = 0
        annop = anno
      end if
      call skipandreadgame(plun,index1h(ora,giorno,mese,anno),prec,0)
      read (plun) ((field(i,j),j=1,nlat),i=1,nlon)
      n = 0
      do i = 1 , nlon
        do j = 1 , nlat
          n = n + 1
          pi(n) = field(i,j)
          if ( pi(n)<1.0 ) pi(n) = 0.0
          la(n) = lat(i,j)
          lo(n) = lon(i,j)
        end do
      end do
    else if ( flag==2 ) then
      if ( anno/=annot ) then
        if ( tlun>0 ) then
          close (tlun)
        else
          call getlun(tlun)
        end if
        write (cfl,'(a,i4)') dir(1:ld)//'temp/t2_' , anno
        write (ilog,'(18x,a)') 'Opening '//trim(cfl)
        open (tlun,file=cfl,status='old',form='unformatted')
        trec = 0
        annot = anno
      end if
      call skipandreadgame(tlun,index1h(ora,giorno,mese,anno),trec,0)
      read (tlun) ((field(i,j),j=1,nlat),i=1,nlon)
      n = 0
      do i = 1 , nlon
        do j = 1 , nlat
          n = n + 1
          pi(n) = field(i,j) - 273.15
          la(n) = lat(i,j)
          lo(n) = lon(i,j)
        end do
      end do
    end if
  end subroutine acqwascen01
 
  subroutine regcmreadnc(nlon,nlat,ora,giorno,mese,anno,iflag,pi,la,lo,n,simname)
    implicit none
    include 'netcdf.inc'
    integer , parameter :: n2dim = 1000000
    integer :: anno , giorno , iflag , mese , n , nlat , nlon , ora
    integer :: anno2 , mese2
    real , dimension(nlon*nlat) :: la , lo , pi
    intent (in) anno , giorno , iflag , mese , ora
    intent (out) la , lo
    intent (inout) n
    character(300) :: cfl , dir , oldcfl , oldcflt , oldcflp
    character(*) :: simname
    integer :: code , i , istatus , itime , j , ld , ilog , lun , pcode ,    &
           irec , scode , tcode  ,  lunp , lunt
    logical :: first , firstp, firstt, firsts
    integer , dimension(4) , save :: icount , istart
    real , dimension(n2dim) , save :: loclat , loclon
    real , save :: pfact
    character(64) , save :: rainmis
    real , dimension(nlon,nlat) :: tmpm
    double precision , dimension(nlon,nlat) :: tmpm_db
    data firstt/.true./
    data firstp/.true./
    data firsts/.true./
    data first/.true./
    integer :: xmese, irecold , irecoldt
    save oldcfl,pcode,lunp,lunt,irecold,tcode
    save oldcflt,irecoldt,oldcflp
    save first , firstp, firstt, firsts, ilog
    xmese=0
    if (mese.gt.1) xmese=julianday(monthlen(mese-1,anno), &
         mese-1,anno)
    irec = xmese+giorno
!    print*,'irec',irec
!    irec=(xmese*24)+1 ; 
    if (irec.eq.0) irec=1
!    if (anno == 2000) then
!      xmese=julianday(monthlen(mese-1,anno), &
!        mese-1,anno)-(31+29)
!      if (mese <= 2) then
!        write(6,'(/,14x,a)')'Persiann Dataset begin at '// &
!          '2000/03/01 00:00'
!        call exit(1)
!      end if
!    end if
    if ( nlon*nlat>n2dim ) then
      write (6,'(/,14x,a)') 'Too big RegCM domain,'//                        &
                   'increase n2dim parameter inside regcmreadnc'
      call exit(1)
    end if
    if ( trim(schym(2)) == 'AFR44-ERAI') then
      if (anno >= 1979 .and. anno <= 1980) then
        write (cfl,'(a,a)') trim(simname)//'1979010112-1980123112.nc'
      else if (anno >= 1981 .and. anno <=1985) then
        write (cfl,'(a,a)') trim(simname)//'1981010112-1985123112.nc'
      else if (anno >= 1986 .and. anno <=1990) then
        write (cfl,'(a,a)') trim(simname)//'1986010112-1990123112.nc'
      else if (anno >= 1991 .and. anno <=1995) then
        write (cfl,'(a,a)') trim(simname)//'1991010112-1995123112.nc'
      else if (anno >= 1996 .and. anno <=2000) then
        write (cfl,'(a,a)') trim(simname)//'1996010112-2000123112.nc'
      else if (anno >= 2001 .and. anno <=2005) then
        write (cfl,'(a,a)') trim(simname)//'2001010112-2005123112.nc'
      else if (anno >= 2006 .and. anno <=2009) then
        write (cfl,'(a,a)') trim(simname)//'2006010112-2009123112.nc'
      end if
    else if ( trim(schym(2)) == 'SAM44-ERAI') then
      if (anno >= 1979 .and. anno <= 1980) then
        write (cfl,'(a,a)') trim(simname)//'1979010112-1980123112.nc'
      else if (anno >= 1981 .and. anno <=1985) then
        write (cfl,'(a,a)') trim(simname)//'1981010112-1985123112.nc'
      else if (anno >= 1986 .and. anno <=1990) then
        write (cfl,'(a,a)') trim(simname)//'1986010112-1990123112.nc'
      else if (anno >= 1991 .and. anno <=1995) then
        write (cfl,'(a,a)') trim(simname)//'1991010112-1995123112.nc'
      else if (anno >= 1996 .and. anno <=2000) then
        write (cfl,'(a,a)') trim(simname)//'1996010112-2000123112.nc'
      else if (anno >= 2001 .and. anno <=2005) then
        write (cfl,'(a,a)') trim(simname)//'2001010112-2005123112.nc'
      else if (anno >= 2006 .and. anno <=2008) then
        write (cfl,'(a,a)') trim(simname)//'2006010112-2008123112.nc'
      end if
    else if ( trim(schym(2)) == 'MED44-ERAIv4') then
      if (anno >= 1979 .and. anno <= 1980) then
        write (cfl,'(a,a)') trim(simname)//'1979010112-1980123112.nc'
      else if (anno >= 1981 .and. anno <=1985) then
        write (cfl,'(a,a)') trim(simname)//'1981010112-1985123112.nc'
      else if (anno >= 1986 .and. anno <=1990) then
        write (cfl,'(a,a)') trim(simname)//'1986010112-1990123112.nc'
      else if (anno >= 1991 .and. anno <=1995) then
        write (cfl,'(a,a)') trim(simname)//'1991010112-1995123112.nc'
      else if (anno >= 1996 .and. anno <=2000) then
        write (cfl,'(a,a)') trim(simname)//'1996010112-2000123112.nc'
      else if (anno >= 2001 .and. anno <=2005) then
        write (cfl,'(a,a)') trim(simname)//'2001010112-2005123112.nc'
      else if (anno >= 2006 .and. anno <=2010) then
        write (cfl,'(a,a)') trim(simname)//'2006010112-2010123112.nc'
      else if (anno >= 2011 .and. anno <=2012) then
        write (cfl,'(a,a)') trim(simname)//'2011010112-2012123112.nc'
      end if
    else if ( trim(schym(2)) == 'MED11-ERAI') then
      if (anno >= 1979 .and. anno <= 1980) then
        write (cfl,'(a,a)') trim(simname)//'1979010112-1980123112.nc'
      else if (anno >= 1981 .and. anno <=1985) then
        write (cfl,'(a,a)') trim(simname)//'1981010112-1985123112.nc'
      else if (anno >= 1986 .and. anno <=1990) then
        write (cfl,'(a,a)') trim(simname)//'1986010112-1990123112.nc'
      else if (anno >= 1991 .and. anno <=1995) then
        write (cfl,'(a,a)') trim(simname)//'1991010112-1995123112.nc'
      else if (anno >= 1996 .and. anno <=2000) then
        write (cfl,'(a,a)') trim(simname)//'1996010112-2000123112.nc'
      else if (anno >= 2001 .and. anno <=2005) then
        write (cfl,'(a,a)') trim(simname)//'2001010112-2005123112.nc'
      else if (anno >= 2006 .and. anno <=2010) then
        write (cfl,'(a,a)') trim(simname)//'2006010112-2010123112.nc'
      end if
    else if ( trim(schym(2)) == 'EUR44-ERAI') then
      if (anno >= 1979 .and. anno <= 1980) then
        write (cfl,'(a,a)') trim(simname)//'1979010112-1980123112.nc'
      else if (anno >= 1981 .and. anno <=1985) then
        write (cfl,'(a,a)') trim(simname)//'1981010112-1985123112.nc'
      else if (anno >= 1986 .and. anno <=1990) then
        write (cfl,'(a,a)') trim(simname)//'1986010112-1990123112.nc'
      else if (anno >= 1991 .and. anno <=1995) then
        write (cfl,'(a,a)') trim(simname)//'1991010112-1995123112.nc'
      else if (anno >= 1996 .and. anno <=2000) then
        write (cfl,'(a,a)') trim(simname)//'1996010112-2000123112.nc'
      else if (anno >= 2001 .and. anno <=2005) then
        write (cfl,'(a,a)') trim(simname)//'2001010112-2005123112.nc'
      else if (anno >= 2006 .and. anno <=2010) then
        write (cfl,'(a,a)') trim(simname)//'2006010112-2010123112.nc'
      else if (anno >= 2011 .and. anno <=2012) then
        write (cfl,'(a,a)') trim(simname)//'2011010112-2012123112.nc'
      end if
    else if ( trim(schym(2)) == 'EUR44-MOHC-his') then
      if (anno >= 1950 .and. anno <=1950) then
        write (cfl,'(a,a)') trim(simname)//'1950010112-1950123012.nc'
      else if (anno >= 1951 .and. anno <=1955) then
        write (cfl,'(a,a)') trim(simname)//'1951010112-1955123012.nc'
      else if (anno >= 1956 .and. anno <=1960) then
        write (cfl,'(a,a)') trim(simname)//'1956010112-1960123012.nc'
      else if (anno >= 1961 .and. anno <=1965) then
        write (cfl,'(a,a)') trim(simname)//'1961010112-1965123012.nc'
      else if (anno >= 1966 .and. anno <=1970) then
        write (cfl,'(a,a)') trim(simname)//'1966010112-1970123012.nc'
      else if (anno >= 1971 .and. anno <=1975) then
        write (cfl,'(a,a)') trim(simname)//'1971010112-1975123012.nc'
      else if (anno >= 1976 .and. anno <=1980) then
        write (cfl,'(a,a)') trim(simname)//'1976010112-1980123012.nc'
      else if (anno >= 1981 .and. anno <=1985) then
        write (cfl,'(a,a)') trim(simname)//'1981010112-1985123012.nc'
      else if (anno >= 1986 .and. anno <=1990) then
        write (cfl,'(a,a)') trim(simname)//'1986010112-1990123012.nc'
      else if (anno >= 1991 .and. anno <=1995) then
        write (cfl,'(a,a)') trim(simname)//'1991010112-1995123012.nc'
      else if (anno >= 1996 .and. anno <=2000) then
        write (cfl,'(a,a)') trim(simname)//'1996010112-2000123012.nc'
      else if (anno >= 2001 .and. anno <=2005) then
        write (cfl,'(a,a)') trim(simname)//'2001010112-2005123012.nc'
      end if
    else if ( trim(schym(2)) == 'CAM44-MPI-his') then
      if (anno >= 1970 .and. anno <=1970) then
        write (cfl,'(a,a)') trim(simname)//'1970010112-1970123112.nc'
      else if (anno >= 1971 .and. anno <=1975) then
        write (cfl,'(a,a)') trim(simname)//'1971010112-1975123112.nc'
      else if (anno >= 1976 .and. anno <=1980) then
        write (cfl,'(a,a)') trim(simname)//'1976010112-1980123112.nc'
      else if (anno >= 1981 .and. anno <=1985) then
        write (cfl,'(a,a)') trim(simname)//'1981010112-1985123112.nc'
      else if (anno >= 1986 .and. anno <=1990) then
        write (cfl,'(a,a)') trim(simname)//'1986010112-1990123112.nc'
      else if (anno >= 1991 .and. anno <=1995) then
        write (cfl,'(a,a)') trim(simname)//'1991010112-1995123112.nc'
      else if (anno >= 1996 .and. anno <=2000) then
        write (cfl,'(a,a)') trim(simname)//'1996010112-2000123112.nc'
      else if (anno >= 2001 .and. anno <=2005) then
        write (cfl,'(a,a)') trim(simname)//'2001010112-2005123112.nc'
      end if
    else if ( trim(schym(2)) == 'CAM44-MPI-rcp85') then
      if (anno >= 2006 .and. anno <=2010) then
        write (cfl,'(a,a)') trim(simname)//'2006010112-2010123112.nc'
      else if (anno >= 2011 .and. anno <=2015) then
        write (cfl,'(a,a)') trim(simname)//'2011010112-2015123112.nc'
      else if (anno >= 2016 .and. anno <=2020) then
        write (cfl,'(a,a)') trim(simname)//'2016010112-2020123112.nc'
      else if (anno >= 2021 .and. anno <=2025) then
        write (cfl,'(a,a)') trim(simname)//'2021010112-2025123112.nc'
      else if (anno >= 2026 .and. anno <=2030) then
        write (cfl,'(a,a)') trim(simname)//'2026010112-2030123112.nc'
      else if (anno >= 2031 .and. anno <=2035) then
        write (cfl,'(a,a)') trim(simname)//'2031010112-2035123112.nc'
      else if (anno >= 2036 .and. anno <=2040) then
        write (cfl,'(a,a)') trim(simname)//'2036010112-2040123112.nc'
      else if (anno >= 2041 .and. anno <=2045) then
        write (cfl,'(a,a)') trim(simname)//'2041010112-2045123112.nc'
      else if (anno >= 2046 .and. anno <=2050) then
        write (cfl,'(a,a)') trim(simname)//'2046010112-2050123112.nc'
      else if (anno >= 2051 .and. anno <=2055) then
        write (cfl,'(a,a)') trim(simname)//'2051010112-2055123112.nc'
      else if (anno >= 2056 .and. anno <=2060) then
        write (cfl,'(a,a)') trim(simname)//'2056010112-2060123112.nc'
      else if (anno >= 2061 .and. anno <=2065) then
        write (cfl,'(a,a)') trim(simname)//'2061010112-2065123112.nc'
      else if (anno >= 2066 .and. anno <=2070) then
        write (cfl,'(a,a)') trim(simname)//'2066010112-2070123112.nc'
      else if (anno >= 2071 .and. anno <=2075) then
        write (cfl,'(a,a)') trim(simname)//'2071010112-2075123112.nc'
      else if (anno >= 2076 .and. anno <=2080) then
        write (cfl,'(a,a)') trim(simname)//'2076010112-2080123112.nc'
      else if (anno >= 2081 .and. anno <=2085) then
        write (cfl,'(a,a)') trim(simname)//'2081010112-2085123112.nc'
      else if (anno >= 2086 .and. anno <=2090) then
        write (cfl,'(a,a)') trim(simname)//'2086010112-2090123112.nc'
      else if (anno >= 2091 .and. anno <=2095) then
        write (cfl,'(a,a)') trim(simname)//'2091010112-2095123112.nc'
      else if (anno >= 2096 .and. anno <=2099) then
        write (cfl,'(a,a)') trim(simname)//'2096010112-2099123012.nc'
      end if
    else if ( trim(schym(2)) == 'AFR44-MPI-his' .or. trim(schym(2)) == 'AFR44-MPIMR-his') then
      if (anno >= 1970 .and. anno <=1970) then
        write (cfl,'(a,a)') trim(simname)//'19700101-19701231.nc'
      else if (anno >= 1971 .and. anno <=1975) then
        write (cfl,'(a,a)') trim(simname)//'19710101-19751231.nc'
      else if (anno >= 1976 .and. anno <=1980) then
        write (cfl,'(a,a)') trim(simname)//'19760101-19801231.nc'
      else if (anno >= 1981 .and. anno <=1985) then
        write (cfl,'(a,a)') trim(simname)//'19810101-19851231.nc'
      else if (anno >= 1986 .and. anno <=1990) then
        write (cfl,'(a,a)') trim(simname)//'19860101-19901231.nc'
      else if (anno >= 1991 .and. anno <=1995) then
        write (cfl,'(a,a)') trim(simname)//'19910101-19951231.nc'
      else if (anno >= 1996 .and. anno <=2000) then
        write (cfl,'(a,a)') trim(simname)//'19960101-20001231.nc'
      else if (anno >= 2001 .and. anno <=2005) then
        write (cfl,'(a,a)') trim(simname)//'20010101-20051130.nc'
      end if
    else if ( trim(schym(2)) == 'AFR44-MPI-rcp85' .or. trim(schym(2)) == 'AFR44-MPIMR-rcp85') then
      if (anno >= 2007 .and. anno <=2010) then
        write (cfl,'(a,a)') trim(simname)//'20070101-20101231.nc'
      else if (anno >= 2011 .and. anno <=2015) then
        write (cfl,'(a,a)') trim(simname)//'20110101-20151231.nc'
      else if (anno >= 2016 .and. anno <=2020) then
        write (cfl,'(a,a)') trim(simname)//'20160101-20201231.nc'
      else if (anno >= 2021 .and. anno <=2025) then
        write (cfl,'(a,a)') trim(simname)//'20210101-20251231.nc'
      else if (anno >= 2026 .and. anno <=2030) then
        write (cfl,'(a,a)') trim(simname)//'20260101-20301231.nc'
      else if (anno >= 2031 .and. anno <=2035) then
        write (cfl,'(a,a)') trim(simname)//'20310101-20351231.nc'
      else if (anno >= 2036 .and. anno <=2040) then
        write (cfl,'(a,a)') trim(simname)//'20360101-20401231.nc'
      else if (anno >= 2041 .and. anno <=2045) then
        write (cfl,'(a,a)') trim(simname)//'20410101-20451231.nc'
      else if (anno >= 2046 .and. anno <=2050) then
        write (cfl,'(a,a)') trim(simname)//'20460101-20501231.nc'
      else if (anno >= 2051 .and. anno <=2055) then
        write (cfl,'(a,a)') trim(simname)//'20510101-20551231.nc'
      else if (anno >= 2056 .and. anno <=2060) then
        write (cfl,'(a,a)') trim(simname)//'20560101-20601231.nc'
      else if (anno >= 2061 .and. anno <=2065) then
        write (cfl,'(a,a)') trim(simname)//'20610101-20651231.nc'
      else if (anno >= 2066 .and. anno <=2070) then
        write (cfl,'(a,a)') trim(simname)//'20660101-20701231.nc'
      else if (anno >= 2071 .and. anno <=2075) then
        write (cfl,'(a,a)') trim(simname)//'20710101-20751231.nc'
      else if (anno >= 2076 .and. anno <=2080) then
        write (cfl,'(a,a)') trim(simname)//'20760101-20801231.nc'
      else if (anno >= 2081 .and. anno <=2085) then
        write (cfl,'(a,a)') trim(simname)//'20810101-20851231.nc'
      else if (anno >= 2086 .and. anno <=2090) then
        write (cfl,'(a,a)') trim(simname)//'20860101-20901231.nc'
      else if (anno >= 2091 .and. anno <=2095) then
        write (cfl,'(a,a)') trim(simname)//'20910101-20951231.nc'
      else if (anno >= 2096 .and. anno <=2099) then
        write (cfl,'(a,a)') trim(simname)//'20960101-20991231.nc'
      end if
    else if ( trim(schym(2)(1:6)) == 'CORDEX') then
!      anno2 = anno ; mese2 = mese + 1
!      if ( mese == 12 ) then
!         anno2=anno+1 ; mese2 = 01
!      end if
!      write (cfl,'(a,i6,a,i6,a)') trim(simname),anno*100 + mese,'010300-',anno2*100 + mese2, &
!             '010000.nc'
      irec = (giorno-1)*8 + (ora/3) + 1
      if (irec.eq.0) irec=1
      cfl = simname
    else if ( trim(schym(2)) == 'standard' ) then
!      irec=(xmese*24)+1 
      write (cfl,'(a,i6,a)') trim(schym(16))//'SRF.' , anno*100 + mese , '0100.nc'
      ! irec = (giorno-1)*24 + ora + 1
      irec = (giorno-1)*8 + (ora/3) + 1
      if (irec.eq.0) irec=1
    end if

    if ( firstt .or. firstp ) then
!      call chymregcmoutput(dir)
!      ld = len_trim(dir)
!      call mvgetiflags(70,ilog)
      if ( ilog<=0 .or. ilog>=100 ) ilog = 6
!      write (cfl,'(a,i6,a)') 'SRF.' , anno*100 + mese , '0100.nc'
!      call getlun(lun)
!      write (ilog,'(18x,a)') 'Opening RegCM file '//trim(cfl)
!      istatus = nf_open(dir(1:ld)//cfl,nf_nowrite,lun)
!      write (ilog,'(18x,a)') 'Opening RegCM file '//trim(cfl)
!      istatus = nf_open(trim(cfl),nf_nowrite,lun)
!      call chkncstatus(istatus,nf_noerr,'open',dir,cfl)
 !     call openncregcmfile(lun,trim(cfl),ilog)
      if ( iflag==1 ) then
        call getlun(lunp)
        oldcflp = cfl
        irecold = irec
        write (ilog,'(18x,a)') 'Opening RegCM file '//trim(cfl)
        istatus = nf_open(trim(cfl),nf_nowrite,lunp)
        call chkncstatus(istatus,nf_noerr,'open',dir,cfl)
!      print*,oldcfl
!      print*,cfl
        istatus = nf_inq_varid(lunp,'lat',code)
        if ( istatus/=nf_noerr ) then
          istatus = nf_inq_varid(lunp,'xlat',code)
          if ( istatus/=nf_noerr ) then
            write (ilog,'(16x,a)') &
              '  ERROR Reading LAT dimension from '//dir(1:ld)//cfl
            stop
          end if
        end if  
!        call chkncstatus(istatus,nf_noerr,'lat',dir,cfl)
        if ( trim(schym(2)(1:6)) == 'CORDEX') then
      !  if ( trim(schym(2)) == 'CORDEX_EUR-11_HadGEM_hist_MARCONI') then
          istatus = nf_get_var_double(lunp,code,tmpm_db)
          call chkncstatus(istatus,nf_noerr,'lat',dir,cfl)
          tmpm = real(tmpm_db)
        else
          istatus = nf_get_var_real(lunp,code,tmpm)
          call chkncstatus(istatus,nf_noerr,'lat',dir,cfl)
        end if
        write (logfile,'(16x,a,F10.2)') &
          ' first LAT dimension ',tmpm(1,1)
        if (rchym(11) == -1) then
          lati = tmpm
        end if
        call regcmmat2vec(tmpm,nlon,nlat,1.0,loclat)
        istatus = nf_inq_varid(lunp,'lon',code)
        if ( istatus/=nf_noerr ) then
          istatus = nf_inq_varid(lunp,'xlon',code)
          if ( istatus/=nf_noerr ) then
            write (ilog,'(16x,a)') &
              '  ERROR Reading LON dimension from '//dir(1:ld)//cfl
            stop
          end if
        end if
!        call chkncstatus(istatus,nf_noerr,'lon',dir,cfl)
        if ( trim(schym(2)(1:6)) == 'CORDEX') then
     !   if ( trim(schym(2)) == 'CORDEX_EUR-11_HadGEM_hist_MARCONI') then
          istatus = nf_get_var_double(lunp,code,tmpm_db)
          call chkncstatus(istatus,nf_noerr,'lon',dir,cfl)
          tmpm = real(tmpm_db)
        else
          istatus = nf_get_var_real(lunp,code,tmpm)
          call chkncstatus(istatus,nf_noerr,'lon',dir,cfl)
        end if
        if (rchym(11) == -1) then
          loni = tmpm
        end if
        write (logfile,'(16x,a,F10.2)') &
          ' first LON dimension ',tmpm(1,1)
        call regcmmat2vec(tmpm,nlon,nlat,1.0,loclon)
        istatus = nf_inq_varid(lunp,'pr',pcode)
        if ( istatus/=nf_noerr ) istatus = nf_inq_varid(lunp,'pr',pcode)
        call chkncstatus(istatus,nf_noerr,'rain',dir,cfl)
        rainmis = ' '
        istatus = nf_get_att_text(lunp,pcode,'units',rainmis)
        call chkncstatus(istatus,nf_noerr,'rain units',dir,cfl)
        pfact = 3600.0
        if ( trim(schym(2)) == 'AFR44-MPI-rcp85' .or. trim(schym(2)) == 'AFR44-MPIMR-rcp85' .or. &
           trim(schym(2)) == 'AFR44-MPI-his' .or. trim(schym(2)) == 'AFR44-MPIMR-his') pfact = 86400
        if ( trim(rainmis)=='kg m-2 day-1' ) pfact = 1.0/24.0
        istatus = nf_inq_varid(lunp,'time',itime)
        call chkncstatus(istatus,nf_noerr,'time',dir,cfl)
        rainmis = ' '
        istatus = nf_get_att_text(lunp,itime,'calendar',rainmis)
        call chkncstatus(istatus,nf_noerr,'calendar','',cfl)
        if ( trim(rainmis)=='noleap' ) then
          call mvsetflags('Calendar',1.0)
          write (6,'(18x,a)') 'No-leap calendar will be used'
        else if ( trim(rainmis)=='360_days' .or. trim(rainmis)=='360_day' ) then
          call mvsetflags('Calendar',2.0)
          write (6,'(18x,a)') '360 days calendar will be used'
        else
          call mvsetflags('Calendar',0.0)
          write (6,'(18x,a)') 'Gregorian calendar will be used'
        end if
        istatus = nf_get_att_real(lunp,NF_GLOBAL,'grid_size_in_meters',regcm_res)
        call chkncstatus(istatus,nf_noerr,'grid_size_in_meters','',cfl)
        if ( trim(schym(2)) == 'standard' .or. trim(schym(2)(1:6)) == 'CORDEX' ) then
          firstp = .false.
        end if
      else if ( iflag==2 ) then
        call getlun(lunt)
        oldcflt = cfl
        irecoldt = irec
        write (logfile,'(18x,a)') 'Opening RegCM file '//trim(cfl)//' for Temperature'
        istatus = nf_open(trim(cfl),nf_nowrite,lunt)
        call chkncstatus(istatus,nf_noerr,'open',dir,cfl)
        istatus = nf_inq_varid(lunt,'lat',code)
        if ( istatus/=nf_noerr ) then
          istatus = nf_inq_varid(lunt,'xlat',code)
          if ( istatus/=nf_noerr ) then
            write (ilog,'(16x,a)') &
              '  ERROR Reading LAT dimension from '//dir(1:ld)//cfl
            stop
          end if
        end if
!        call chkncstatus(istatus,nf_noerr,'lat',dir,cfl)
!        istatus = nf_get_var_real(lunt,code,tmpm)
        if ( trim(schym(2)(1:6)) == 'CORDEX') then
          istatus = nf_get_var_double(lunt,code,tmpm_db)
          call chkncstatus(istatus,nf_noerr,'lat',dir,cfl)
          tmpm = real(tmpm_db)
        else
          istatus = nf_get_var_real(lunt,code,tmpm)
          call chkncstatus(istatus,nf_noerr,'lat',dir,cfl)
        end if
        if (rchym(11) == -1) then
          lati = real(tmpm)
        end if
        write (logfile,'(16x,a,F10.2)') &
          ' first LAT dimension ',tmpm(1,1)
        call regcmmat2vec(tmpm,nlon,nlat,1.0,loclat)
        istatus = nf_inq_varid(lunt,'lon',code)
        if ( istatus/=nf_noerr ) then
          istatus = nf_inq_varid(lunt,'xlon',code)
          if ( istatus/=nf_noerr ) then
            write (ilog,'(16x,a)') &
              '  ERROR Reading LON dimension from '//dir(1:ld)//cfl
            stop
          end if
        end if
!        call chkncstatus(istatus,nf_noerr,'lon',dir,cfl)
!        istatus = nf_get_var_real(lunt,code,tmpm)
        if ( trim(schym(2)(1:6)) == 'CORDEX') then
          istatus = nf_get_var_double(lunt,code,tmpm_db)
          call chkncstatus(istatus,nf_noerr,'lon',dir,cfl)
          tmpm = real(tmpm_db)
        else
          istatus = nf_get_var_real(lunt,code,tmpm)
          call chkncstatus(istatus,nf_noerr,'lon',dir,cfl)
        end if
        write (logfile,'(16x,a,F10.2)') &
          ' first LON dimension ',tmpm(1,1)
        if (rchym(11) == -1) then
          loni = real(tmpm)
        end if
        call regcmmat2vec(tmpm,nlon,nlat,1.0,loclon)
        istatus = nf_inq_varid(lunt,'t2m',tcode)
        if ( istatus/=nf_noerr ) istatus = nf_inq_varid(lunt,'tas',tcode)
        call chkncstatus(istatus,nf_noerr,'temp',dir,cfl)
        istatus = nf_inq_varid(lunt,'time',itime)
        call chkncstatus(istatus,nf_noerr,'time',dir,cfl)
        rainmis = ' '
        istatus = nf_get_att_text(lunt,itime,'calendar',rainmis)
        call chkncstatus(istatus,nf_noerr,'calendar','',cfl)
        if ( trim(rainmis)=='noleap' ) then
          call mvsetflags('Calendar',1.0)
          write (6,'(18x,a)') 'No-leap calendar will be used'
        else if ( trim(rainmis)=='360_days' .or. trim(rainmis)=='360_day' ) then
          call mvsetflags('Calendar',2.0)
          write (6,'(18x,a)') '360 days calendar will be used'
        else
          call mvsetflags('Calendar',0.0)
          write (6,'(18x,a)') 'Gregorian calendar will be used'
        end if
        istatus = nf_get_att_real(lunt,NF_GLOBAL,'grid_size_in_meters',regcm_res)
        call chkncstatus(istatus,nf_noerr,'grid_size_in_meters','',cfl)
      else if ( iflag==3 ) then
        istatus = nf_inq_varid(lun,'scv',scode)
        call chkncstatus(istatus,nf_noerr,'snow',dir,cfl)
      end if
!      istatus = nf_inq_varid(lun,'t2m',tcode)
!      if ( istatus/=nf_noerr ) istatus = nf_inq_varid(lun,'tas',tcode)
!      call chkncstatus(istatus,nf_noerr,'temp',dir,cfl)
!      istatus = nf_inq_varid(lun,'scv',scode)
!      call chkncstatus(istatus,nf_noerr,'snow',dir,cfl)
      istart(1) = 1
      istart(2) = 1
      icount(1) = nlon
      icount(2) = nlat
!      first = .false.
      if ( trim(schym(2)) == 'standard' .or. trim(schym(2)(1:6)) == 'CORDEX' ) then
        firstt = .false.
      else
        firstp = .false.
        firstt = .false.
      end if
    end if
!    write (cfl,'(a,i6,a)') 'SRF.' , anno*100 + mese , '0100.nc'
    if ( iflag == 1 .and. trim(cfl)/=trim(oldcflp) ) then
      istatus = nf_close(lunp)
      call chkncstatus(istatus,nf_noerr,'close',dir,oldcflp)
      call getlun(lunp)
      write (ilog,'(18x,a)') 'Opening RegCM file '//trim(cfl)
!      istatus = nf_open(dir(1:ld)//cfl,nf_nowrite,lun)
      istatus = nf_open(trim(cfl),nf_nowrite,lunp)
!      call openncregcmfile(lun,trim(cfl),ilog)
      call chkncstatus(istatus,nf_noerr,'open',dir,cfl)
      oldcflp = cfl
    else if ( iflag == 2 .and. trim(cfl)/=trim(oldcflt) ) then
            istatus = nf_close(lunt)
      call chkncstatus(istatus,nf_noerr,'close',dir,oldcflt)
      call getlun(lunt)
      write (ilog,'(18x,a)') 'Opening RegCM file '//trim(cfl)
!      istatus = nf_open(dir(1:ld)//cfl,nf_nowrite,lun)
      istatus = nf_open(trim(cfl),nf_nowrite,lunt)
!      call openncregcmfile(lun,trim(cfl),ilog)
      call chkncstatus(istatus,nf_noerr,'open',dir,cfl)
      oldcflt = cfl
    end if
!    irec = (giorno-1)*8 + (ora/3)
!    if ( irec==0 ) irec = 1
    if (  iflag==1 ) then !.and. (irecold/=irec .or. firstp)) then
        istart(3) = irec
        icount(3) = 1
        
        istatus = nf_get_vara_real(lunp,pcode,istart(1:3),icount(1:3),tmpm)
        call chkncstatus(istatus,nf_noerr,'rain filed','',cfl)
        if (rchym(11) == -1) then
          rani = tmpm*pfact
        end if
        write (logfile,'(16x,a,i4,a,F10.2)') &
          ' first RegCM PRECIPITATION value:irec ; tmpm(1,1)*pfact ',irec ," ; ", tmpm(1,1)*pfact
        call regcmmat2vec(tmpm,nlon,nlat,pfact,pi)
!        firstp = .false.
    else if ( iflag==2 ) then !.and. (irecoldt/=irec .or. firstt)) then
        if ( trim(schym(2)) == 'standard' ) then
         istart(3) = 1
         icount(3) = 1
         istart(4) = irec
         icount(4) = 1
         istatus = nf_get_vara_real(lunt,tcode,istart(1:4),icount(1:4),tmpm)
        else
         print*,"irec",irec
         print*,"lunt",lunt
         print*,"tcode",tcode
         istart(3) = irec
         icount(3) = 1
         istatus = nf_get_vara_real(lunt,tcode,istart(1:3),icount(1:3),tmpm)
        end if
!        istart(3) = 1
!        icount(3) = 1
!        istart(4) = irec
!        icount(4) = 1
!        istatus = nf_get_vara_real(lun,tcode,istart(1:4),icount(1:4),tmpm)
        call chkncstatus(istatus,nf_noerr,'temp filed','',cfl)
        call kelvin2centi(tmpm,nlon,nlat,1)
        write (logfile,'(16x,a,F10.2)') &
          ' first TEMPERATURE RegCM ',tmpm(1,1)
        if (rchym(11) == -1) then
          temi = tmpm
        end if
        write (logfile,'(16x,a,i4,a,F10.2)') &
          ' first RegCM TEMPERATURE value:irec ; tmpm(1,1) ',irec ," ; ", tmpm(1,1)
        call regcmmat2vec(tmpm,nlon,nlat,1.0,pi)
!        firstt = .false.
      else if ( iflag==3 ) then
        istart(3) = irec
        icount(3) = 1
        istatus = nf_get_vara_real(lun,scode,istart(1:3),icount(1:3),tmpm)
        call chkncstatus(istatus,nf_noerr,'snow filed',dir,cfl)
        call regcmmat2vec(tmpm,nlon,nlat,1.0,pi)
        firsts = .false.
      else
        write (6,'(16x,a)') '  Bad flag passed to regcmreadnc. Exiting...'
        call exit(1)
    end if
    n = nlon*nlat
    la(1:n) = loclat(1:n)
    lo(1:n) = loclon(1:n)
  end subroutine regcmreadnc
 
  subroutine regcmmat2vec(mat,nlon,nlat,cfact,vec)
    implicit none
    real :: cfact
    integer :: nlat , nlon
    real , dimension(nlon,nlat) :: mat
    real , dimension(nlon*nlat) :: vec
    intent (in) cfact , mat , nlat , nlon
    intent (out) vec
    integer :: i , j , n
    n = 0
    do i = 1 , nlon
      do j = 1 , nlat
        n = n + 1
        vec(n) = mat(i,j)*cfact
      end do
    end do
  end subroutine regcmmat2vec
 
  subroutine chkncstatus(istatus,nf_noerr,field,dir,cfl)
    implicit none
    character(*) :: cfl , dir , field
    integer :: istatus , nf_noerr
    intent (in) istatus , nf_noerr
    if ( istatus==nf_noerr ) return
    if ( len_trim(field)==0 .or. field(1:len_trim(field))=='open' ) then
      write (6,'(/,10x,a)') 'NetCDF Library error while opening file'
      write (6,'(10x,a)') trim(dir)//trim(cfl)
    else if ( field(1:len_trim(field))=='close' ) then
      write (6,'(/,10x,a)') 'NetCDF Library error while closing file'
      write (6,'(10x,a)') trim(dir)//trim(cfl)
    else
      write (6,'(/,10x,a)') 'NetCDF Library error while reading '// &
                      trim(field)//' field from'
      write (6,'(10x,a)') trim(dir)//trim(cfl)
    end if
    write (6,'(10x,a)') 'Exiting with status = 1'
    call exit(1)
  end subroutine chkncstatus
 
  subroutine chymregcmoutput(cfile)
    implicit none
    character(len=*) :: cfile
    intent (out) cfile
    cfile = trim(schym(16))//trim(schym(17))
  end subroutine chymregcmoutput

  subroutine museodata(ora,giorno,mese,anno,pi,la,lo,n)
    implicit none
    real pi(1),la(1),lo(1)
    integer ih,ixmm,jxmm ; parameter(ixmm=70,jxmm=70)
    real , allocatable , dimension(:,:,:) :: p
    real , allocatable , dimension(:,:) :: xlat,xlon
    integer ora,giorno,mese,anno,n
    integer oldgiorno,oldmese,oldanno,ix,jx
    data oldgiorno,oldmese,oldanno /-1,-1,-1/
    save oldgiorno,oldmese,oldanno,ix,jx,p,xlat,xlon
    integer i,j,len_trim
    logical start ; data start /.true./ ; save start

    if (start) then
       allocate(p(ixmm,jxmm,25))
       allocate(xlat(ixmm,jxmm))
       allocate(xlon(ixmm,jxmm))
       start=.false.
    endif
!   write (logun(2),'(16x,a)') '> Using MuSEO archive data'
    if (anno.ne.oldanno.or.mese.ne.oldmese.or.giorno.ne.oldgiorno) then
!      write (logun(2),'(18x,a)') 'Retrieving MM5 data from MuSEO for '//
!   2       now(1:len_trim(now))
       xlat=0 ; xlon=0.0
       call museomm52d('rai',anno,mese,giorno,-25,mchym(12), &
           p,xlat,xlon,ixmm,jxmm,ix,jx)
       oldgiorno=giorno ; oldmese=mese ; oldanno=anno
    endif
    n=0
    if (ix.le.1.or.jx.le.1) then
       write (logfile,'(20x,a)') ' MuSEO data not available for this day.'
       return
    endif
    do i=1,ix ; do j=1,jx
       n=n+1
       pi(n)=(p(i,j,ora+2)-p(i,j,ora+1))*10.0
       if (pi(n).lt.0.0) pi(n)=0.0
       la(n)=xlat(i,j)
       lo(n)=xlon(i,j)
    enddo ; enddo
    return
  end subroutine museodata

end module mod_museo
