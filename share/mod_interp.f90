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

module mod_interp

  use mod_libmv
  use mod_internal
  use mod_vector
  use mod_phys
  use mod_cellaut
  use mod_museo

  integer , private , parameter :: clarke_1866_datum = 1
  integer , private , parameter :: grs_80_datum = 2
  integer , private , parameter :: wgs_84_datum = 3
  real(8) , private , parameter :: lower_eps_limit = 1.0E-14

  integer :: iread
  integer , parameter :: ixm1 = 100
  integer , parameter :: jxm1 = 100
  integer , dimension(ixm1,jxm1) :: land
  integer :: isave , jsave

  contains

  integer function indextm(gph,xlon,xlat,nlon,nlat)
    implicit none
    integer :: nlat , nlon
    real , dimension(nlon,nlat) :: gph , xlat , xlon
    real :: delta , ghgn , ghgs , phi0 , phin , phis , xlo , &
            zphi0 , zphin , zphis
    integer :: ichk1 , idelta , ilonblk , lon
    !          Il range delle longitudini va da 30.0W a 45E e
    !     non corrisponde a quanto riportato sull'articolo
    !     (Verdecchia et Al, 1995). Ma questa e' la scelta che da
    !      meno errori rispetto al file di Tibaldi (47 su 15431).
    indextm = 0
    ilonblk = 0
    do lon = 1 , 21
      xlo = -30.00 + 3.75*(lon-1)
      ichk1 = 0
      do idelta = 1 , 3
        delta = -3.75 + 3.75*(idelta-1)
        phin = 78.75 + delta
        phi0 = 60.00 + delta
        phis = 41.25 + delta
        zphi0 = finterp(gph,xlat,xlon,nlon,nlat,phi0,xlo)
        zphis = finterp(gph,xlat,xlon,nlon,nlat,phis,xlo)
        zphin = finterp(gph,xlat,xlon,nlon,nlat,phin,xlo)
        ghgs = (zphi0-zphis)/(phi0-phis)
        ghgn = (zphin-zphi0)/(phin-phi0)
        if ( ghgs>0.0 .and. ghgn<-10.0 ) ichk1 = ichk1 + 1
      end do
      if ( ichk1<=0 ) then
        ilonblk = 0
      else if ( ilonblk==2 ) then
        indextm = 1
        return
      else if ( ilonblk==0 .or. ilonblk==1 ) then
        ilonblk = ilonblk + 1
      else
        stop ' Flux error inside indextm'
      end if
    end do
  end function indextm

  subroutine chymdownscaling(tm,xlon,xlat,ixm1,jxm1,mt,ca,wk, &
                             nlon,nlat,slon,slat,dij,dji)
    implicit none
    real :: dij , dji , slat , slon
    integer :: ixm1 , jxm1 , nlat , nlon
    integer , dimension(nlon,nlat) :: ca
    real , dimension(nlon,nlat) :: mt , wk
    real , dimension(ixm1,jxm1) :: tm , xlat , xlon
    intent (in) dij , dji , ixm1 , jxm1 , tm
    intent (inout) mt
    intent (out) ca
    real :: dist1 , dist2 , xm
    integer :: lix , ljx
    integer :: i , j , im , jm , np
    np = 0
    xm = 0
    ca = 8
    do im = 1 , ixm1
      do jm = 1 , jxm1
        i = nint((xlon(im,jm)-slon)/dij) + 1
        j = nint((xlat(im,jm)-slat)/dji) + 1
        if ( indexrange(i,1,nlon) .and. indexrange(j,1,nlat) ) then
          mt(i,j) = tm(im,jm)
          ca(i,j) = 0
          xm = xm + tm(im,jm)
          np = np + 1
        end if
      end do
    end do
    do i = 1 , nlon
      do j = 1 , nlat
        if ( ca(i,j)==8 ) mt(i,j) = xm/np
      end do
    end do
    call mvgetiflags(65,i)
    if ( i>0 ) then
      j = i
    else
      dist1 = distance(xlat(1,1),xlon(1,1),xlat(2,2),xlon(2,2))
      dist2 = distance(slat,slon,slat+dji,slon+dij)
      j = (dist1/dist2)*100 + 10
    end if
    do i = 1 , j
      call d2cellcycle(mt,ca,wk,nlon,nlat,0.9)
    end do
  end subroutine chymdownscaling

  subroutine grinterp(fld1,lat1,lon1,ix1,jx1,fld2,lat2,lon2,ix2,jx2)
    implicit none
    integer :: ix1 , ix2 , jx1 , jx2
    real , dimension(ix1,jx1) :: fld1 , lat1 , lon1
    real , dimension(ix2,jx2) :: fld2 , lat2 , lon2
    intent (in) ix1 , jx1 , ix2 , jx2 , fld1 , lat1 , lon1 , lat2 , lon2
    intent (out) fld2
    integer :: i , j
    if ( iflg(23)==0 ) then
      do i = 1 , ix2
        do j = 1 , jx2
          fld2(i,j) = finterp(fld1,lat1,lon1,ix1,jx1,lat2(i,j),lon2(i,j))
        end do
      end do
    else
      do i = 1 , ix2
        do j = 1 , jx2
          fld2(i,j) = ginterp(fld1,lat1,lon1,ix1,jx1,lat2(i,j),lon2(i,j))
        end do
      end do
    end if
  end subroutine grinterp

  real function finterp(field,lat,lon,ix,jx,xlat,xlon)
    implicit none
    integer :: ix , jx
    real :: xlat , xlon
    real , dimension(ix,jx) :: field , lat , lon
    intent (in) field , ix , jx , lat , lon , xlat , xlon
    integer :: i , j
    real :: tw , w1 , w2 , w3 , w4
    do i = 1 , ix - 1
      do j = 1 , jx - 1
        if ( lat(i,j)<=xlat .and. lat(i,j+1)>=xlat .and. &
             lon(i,j)<=xlon .and. lon(i+1,j)>=xlon ) then
          w1 = 1.0/distance(lat(i,j),lon(i,j),xlat,xlon)
          w2 = 1.0/distance(lat(i+1,j),lon(i+1,j),xlat,xlon)
          w3 = 1.0/distance(lat(i,j+1),lon(i,j+1),xlat,xlon)
          w4 = 1.0/distance(lat(i+1,j+1),lon(i+1,j+1),xlat,xlon)
          tw = w1 + w2 + w3 + w4
          finterp = (w1*field(i,j)   + w2*field(i+1,j)+ &
                     w3*field(i,j+1) + w4*field(i+1,j+1))/tw
          return
        end if
      end do
    end do
    write (6,'(a,f7.2,a,f7.2,a)') ' finterp error interpolating ' , xlat , &
                                ',' , xlon , ' on given grid.'
    stop ' abnormal stop inside finterp routine '
  end function finterp

  real function ginterp(field,lat,lon,ix,jx,xlat,xlon)
    implicit none
    integer :: ix , jx
    real :: xlat , xlon
    real , dimension(ix,jx) :: field , lat , lon
    integer :: i , idom , iix , istar , j , jjx , l , lix , ljx , lun , &
               mare , ngood
    real :: tw , w1 , w2 , w3 , w4
    data idom , istar , mare/0 , 0 , 16/
    if ( idom/=iflg(5) .or. istar==0 ) then
      iread = 0
      if ( iflg(5)>=1 .and. iflg(5)<=3 ) then
        iread = -1
        call getlun(lun)
        call openmuseofiles(lun,'mm5terrain.dat',0)
        do l = 1 , iflg(5)
          read (lun,err=20,end=20) iix , jjx , ((land(i,j),i=1,iix),j=1,jjx)
        end do
        iread = 33
 20     close (lun)
        if ( iread<0 ) then
          write (6,'(7x,a,i2)') &
                    'MVlib ginterp: Not avail. land use for domain ', iflg(5)
        else if ( iflg(1)>=10 .and. iflg(7)>=0 ) then
          write (6,'(7x,a,i1,a)') 'MVlib ginterp: landuse info '// &
                    'for domain ' , iflg(5) , ' will be used.'
        end if
      else
        if ( iflg(1)>=10 .and. iflg(7)>=0 ) write (6,'(7x,a)') &
               'MVlib ginterp: Unknown domain, landuse info will not be used'
      end if
      idom = iflg(5)
      istar = 1
    end if
    lix = 1
    do i = 2 , ix
      if ( lat(i,1)>lat(i-1,1) ) lix = i
    end do
    if ( lix<3 .or. lix>ix ) lix = ix
    ljx = 1
    do j = 2 , jx
      if ( lon(1,j)>lon(1,j-1) ) ljx = j
    end do
    if ( ljx<3 .or. ljx>jx ) ljx = jx
    if ( iflg(1)>=200 ) write (6,'(9x,a,i3,a1,i3,a)')      &
        'Mvlib ginterp: Algorithms will work on ' , lix , 'x' , ljx , ' grid.'
    do i = 1 , lix - 1
      do j = 1 , ljx - 1
        if ( lat(i,j)<=xlat .and. lat(i+1,j)>xlat .and. &
             lon(i,j)<=xlon .and. lon(i,j+1)>xlon ) then
          w1 = 1.0/distance(lat(i,j),lon(i,j),xlat,xlon)
          w2 = 1.0/distance(lat(i+1,j),lon(i+1,j),xlat,xlon)
          w3 = 1.0/distance(lat(i,j+1),lon(i,j+1),xlat,xlon)
          w4 = 1.0/distance(lat(i+1,j+1),lon(i+1,j+1),xlat,xlon)
          isave = i
          jsave = j
          ngood = 4
          if ( iread>0 .and. iflg(7)>=0 ) then
            if ( land(i,j)==mare ) then
              w1 = 0.0
              ngood = ngood - 1
            end if
            if ( land(i+1,j)==mare ) then
              w2 = 0.0
              ngood = ngood - 1
            end if
            if ( land(i,j+1)==mare ) then
              w3 = 0.0
              ngood = ngood - 1
            end if
            if ( land(i+1,j+1)==mare ) then
              w4 = 0.0
              ngood = ngood - 1
            end if
          end if
          iflg(50) = ngood
          if ( ngood>0 ) then
            tw = w1 + w2 + w3 + w4
            ginterp = (w1*field(i,j)   + w2*field(i+1,j) + &
                       w3*field(i,j+1) + w4*field(i+1,j+1))/tw
            if ( iflg(1)>=100 .and. ngood<4 ) then
              write (6,'(9x,a,i1,a,2f8.4)')  &
                  'MvLib ginterp: only ' , ngood , &
                  ' grid-pts used for point:' , xlat , xlon
            end if
          else
            if ( iflg(1)>=100 ) write (6,'(9x,a,2f8.4,a)')  &
                  'MvLib ginterp: point' , xlat , xlon , ' is on the see.'
            iflg(50) = 0
            ginterp = closest(field,lat,lon,ix,jx,xlat,xlon)
          end if
          return
        end if
      end do
    end do
    ! differisce dal precedente sul terzo controllo
    do i = 1 , lix - 1
      do j = 1 , ljx - 1
        if ( lat(i,j)<=xlat .and. lat(i+1,j)>xlat .and. &
             lon(i,j)<=xlon .and. lon(i+1,j+1)>xlon ) then
          w1 = 1.0/distance(lat(i,j),lon(i,j),xlat,xlon)
          w2 = 1.0/distance(lat(i+1,j),lon(i+1,j),xlat,xlon)
          w3 = 1.0/distance(lat(i,j+1),lon(i,j+1),xlat,xlon)
          w4 = 1.0/distance(lat(i+1,j+1),lon(i+1,j+1),xlat,xlon)
          ngood = 4
          isave = i
          jsave = j
          if ( iread>0 .and. iflg(7)>=0 ) then
            if ( land(i,j)==mare ) then
              w1 = 0.0
              ngood = ngood - 1
            end if
            if ( land(i+1,j)==mare ) then
              w2 = 0.0
              ngood = ngood - 1
            end if
            if ( land(i,j+1)==mare ) then
              w3 = 0.0
              ngood = ngood - 1
            end if
            if ( land(i+1,j+1)==mare ) then
              w4 = 0.0
              ngood = ngood - 1
            end if
          end if
          iflg(50) = ngood
          if ( ngood>0 ) then
            tw = w1 + w2 + w3 + w4
            ginterp = (w1*field(i,j)   + w2*field(i+1,j) + &
                       w3*field(i,j+1) + w4*field(i+1,j+1))/tw
          else
            if ( iflg(1)>=100 ) write (6,'(9x,a,2f8.4,a)')  &
                     'MvLib ginterp: point' , xlat , xlon , ' is on the see.'
            ginterp = closest(field,lat,lon,ix,jx,xlat,xlon)
          end if
          return
        end if
      end do
    end do
    if ( iflg(1)>=100 ) write (6,'(9x,a,2f8.4,a)') 'MvLib ginterp: point' , &
                    xlat , xlon , ' is outside given domain.'
    ginterp = closest(field,lat,lon,ix,jx,xlat,xlon)
    iflg(50) = -1
  end function ginterp

  real function closest(field,lat,lon,ix,jx,xlat,xlon)
    implicit none
    integer :: ix , jx
    real :: xlat , xlon
    real , dimension(ix,jx) :: field , lat , lon
    intent (in) field , ix , jx , lat , lon , xlat , xlon
    integer :: i , is , j , js , lix , ljx , mare
    real :: xmin , xx
    data mare/16/
    lix = 1
    do i = 2 , ix
      if ( lat(i,1)>lat(i-1,1) ) lix = i
    end do
    if ( lix<3 .or. lix>ix ) lix = ix
    ljx = 1
    do j = 2 , jx
      if ( lon(1,j)>lon(1,j-1) ) ljx = j
    end do
    if ( ljx<3 .or. ljx>jx ) ljx = jx
    if ( iflg(1)>=100 ) write (6,'(9x,a,i2,a1,i2,a)') &
                   'Mvlib closest: Algorithms will work on ' , &
                   lix , 'x' , ljx , ' grid.'
    xmin = 1.E35
    do i = 1 , lix
      do j = 1 , ljx
        if ( iflg(7)<0 .or. iread<=0 .or. land(i,j)/=mare ) then
          xx = distance(lat(i,j),lon(i,j),xlat,xlon)
          if ( xx<xmin ) then
            xmin = xx
            is = i
            js = j
          end if
        end if
      end do
    end do
    closest = field(is,js)
    if ( iflg(1)>=100 ) write (6,'(9x,a,i2,x,i2,a,f12.1,a)')          &
       'Mvlib closest: best point is ' , is , js , ', distance is ' , &
       distance(lat(is,js),lon(is,js),xlat,xlon) , ' meters.'
  end function closest

  subroutine ermejopunto(rlat,rlon,gd,nlat,nlon,slat,slon,ns,lat,lon,cd)
    implicit none
    integer , parameter :: npt = 20
    integer :: lat , lon , nlat , nlon , ns
    integer , dimension(npt) :: cd
    integer , dimension(nlat,nlon) :: gd
    real , dimension(nlat,nlon) :: rlat , rlon
    real , dimension(ns) :: slat , slon
    intent (in) gd , nlon , ns
    intent (out) cd , lat , lon
    real :: dmax , tdis
    integer :: i , ii , j , jj , k
    integer , dimension(1) :: locj
    real , dimension(npt) :: is2 , wk2
    lat = 0
    lon = 0
    dmax = 1.E35
    do i = 1 , nlat
      do j = 1 , nlon
        if ( gd(i,j)==0 ) then
          do k = 1 , ns
            x1dum(k) = distance(rlat(i,j),rlon(i,j),slat(k),slon(k))
            if ( x1dum(k)<1.0 ) then ! i punti coincidono
              lat = i
              lon = j
              cd(1) = k
              do ii = 2 , npt
                cd(ii) = 0
              end do
              return
            end if
            i1dum(k) = k
          end do
          k = ns
          do ii = 1 , nlat
            do jj = 1 , nlon
              if ( gd(ii,jj)>0 ) then
                k = k + 1
                x1dum(k) = distance(rlat(i,j),rlon(i,j), &
                                    rlat(ii,jj),rlon(ii,jj))
                i1dum(k) = kfromij(ii,jj,nlat) + ns
              end if
            end do
          end do
          do ii = 1 , npt
            wk2(ii) = minval(x1dum(1:k))
            locj = minloc(x1dum(1:k))
            jj = locj(1)
            is2(ii) = i1dum(jj)
            x1dum(jj) = x1dum(k)
            i1dum(jj) = i1dum(k)
            k = k - 1
          end do
          tdis = 0
          do ii = 1 , npt
            tdis = tdis + wk2(ii)
          end do
          if ( tdis<dmax ) then
            dmax = tdis
            lat = i
            lon = j
            do ii = 1 , npt
              cd(ii) = is2(ii)
            end do
          end if
        end if
      end do
    end do
  end subroutine ermejopunto

  integer function latlonrange(lat,lon,latmin,latmax,lonmin,lonmax)
    implicit none
    real :: lat , latmax , latmin , lon , lonmax , lonmin
    intent (in) lat , latmax , latmin , lon , lonmax , lonmin
    latlonrange = 0
    if ( lat>=latmin .and. lat<=latmax .and. &
         lon>=lonmin .and. lon<=lonmax ) latlonrange = 1
  end function latlonrange

  logical function realrange(n,n1,n2)
    implicit none
    real :: n , n1 , n2
    intent (in) n , n1 , n2
    realrange = .true.
    if ( n<n1 .or. n>n2 ) realrange = .false.
  end function realrange

  logical function indexrange(n,n1,n2)
    implicit none
    integer :: n , n1 , n2
    intent (in) n , n1 , n2
    indexrange = .true.
    if ( n<n1 .or. n>n2 ) indexrange = .false.
  end function indexrange

  subroutine cressman(val,clat,clon,nc,rdfl,rad,mat,lat,lon,ix,jx)
    implicit none
    integer :: ix , jx , nc , rdfl
    real , dimension(nc) :: clat , clon , val
    real , dimension(ix,jx) :: lat , lon , mat
    real :: rad
    intent (in) ix , jx , nc , rad , rdfl , val
    intent (out) mat
    real :: a , a2 , lamn , lamx , lomn , lomx , r , s , w
    integer :: i , j , n , np
    lamn = minval(lat)
    lomn = minval(lon)
    lamx = maxval(lat)
    lomx = maxval(lon)
    if ( rdfl==1 ) then
      a = rad*1000.
      a2 = a*a
      do j = 1 , jx
        do i = 1 , ix
          w = 0.0
          s = 0.0
          np = 0
          do n = 1 , nc
            if ( latlonrange(clat(n),clon(n),lamn,lamx,lomn,lomx)==1 ) then
              r = distance(clat(n),clon(n),lat(i,j),lon(i,j))
              if ( r<a ) then
                s = s + val(n)*(1-r*r/a2)/(1+r*r/a2)
                w = w + (1-r*r/a2)/(1+r*r/a2)
                np = np + 1
              end if
            end if
          end do
          if ( np>0 ) mat(i,j) = s/w
        end do
      end do
    else
      write (6,'(7x,a,i10)') 'MVLib cressman: unimplemented flag: ' , rdfl
    end if
  end subroutine cressman
 
  subroutine utm2latlon(x,y,xlat,xlon)
    implicit none
    real :: x , xlat , xlon , y
    intent (in) x , y
    intent (out) xlon , xlat
    real(8) :: latitude , longitude , utm_x , utm_y
    if ( utm_grid_zone(1)<=0 ) then
      write (6,'(9x,a)') 'Error: latlon2utm must be called before utm2latlon'
      write (6,'(9x,a)') 'to establish utm grid zone.'
      call exit(0)
    end if
    utm_x = x
    utm_y = y
    call utm2ll(utm_x,utm_y,longitude,latitude,utm_grid_zone,3)
    xlon = longitude
    xlat = latitude
  end subroutine utm2latlon
 
  subroutine latlon2utm(xlat,xlon,x,y)
    implicit none
    real :: x , xlat , xlon , y
    intent (in) xlat , xlon
    intent (out) x , y
    real(8) :: utm_x , utm_y , latitude , longitude
    latitude = xlat
    longitude = xlon
    call ll2utm(longitude,latitude,utm_x,utm_y,utm_grid_zone,3)
    x = utm_x
    y = utm_y
  end subroutine latlon2utm

  subroutine get_grid_zone(longitude,latitude,grid_zone,lambda0)
    implicit none
    real(8) :: lambda0 , latitude , longitude
    integer , dimension(2) :: grid_zone
    intent (in) latitude , longitude
    intent (out) grid_zone , lambda0
    real(8) :: m_pi
    integer :: zone_lat , zone_long
    m_pi = acos(-1.0)
    zone_long = int((longitude+180.0)/6.0) + 1
    zone_lat = nint((latitude+80.0)/8.0)
    grid_zone(1) = zone_long
    grid_zone(2) = zone_lat
    if ( (latitude<-80.0) .or. (latitude>84.0) ) then
      lambda0 = 0.0*m_pi/180.0
      return
    end if
    if ( latitude>72.0 .and. longitude>0.0 .and. longitude<42.0 ) then
      if ( longitude<9.0 ) then
        lambda0 = 4.5*m_pi/180.0
      else if ( longitude<21.0 ) then
        lambda0 = 15.0*m_pi/180.0
      else if ( longitude<33.0 ) then
        lambda0 = 27.0*m_pi/180.0
      else if ( longitude<42.0 ) then
        lambda0 = 37.5*m_pi/180.0
      end if
      return
    end if
    if ( latitude>56.0 .and. latitude<64.0 .and. &
         longitude>0.0 .and. longitude<12.0 ) then
      if ( longitude<3.0 ) then
        lambda0 = 1.5*m_pi/180.0
      else if ( longitude<12.0 ) then
        lambda0 = 7.5*m_pi/180.0
      end if
      return
    end if
    lambda0 = (float(zone_long-1)*6.0+(-180.0)+3.0)*m_pi/180.0
  end subroutine get_grid_zone

  subroutine get_lambda0(grid_zone,lambda0,ierr)
    implicit none
    integer :: ierr
    real(8) :: lambda0
    integer , dimension(2) :: grid_zone
    intent (in) grid_zone
    intent (out) ierr , lambda0
    real(8) :: latitude , longitude , m_pi
    integer :: zone_lat , zone_long
    m_pi = acos(-1.0)
    zone_long = grid_zone(1)
    zone_lat = grid_zone(2)
    if ( (zone_long<1) .or. (zone_long>61) ) then
      write (*,*) 'Invalid grid zone format: ' , zone_long , zone_lat
      ierr = -1
      return
    end if
    longitude = (float(zone_long-1)*6.0) - 180.0
    latitude = (float(zone_lat)*8.0) - 80.0
    if ( (latitude<-80.0) .or. (latitude>84.0) ) then
      lambda0 = 0.0
      ierr = 0
      return
    end if
    if ( latitude>56.0 .and. latitude<64.0 .and. &
         longitude>0.0 .and. longitude<12.0 ) then
      if ( longitude<3.0 ) then
        lambda0 = 1.5*m_pi/180.0
      else if ( longitude<12 ) then
        lambda0 = 7.5*m_pi/180.0
      end if
      ierr = 0
      return
    end if
    if ( latitude>72.0 .and. longitude>0.0 .and. longitude<42.0 ) then
      if ( longitude<9.0 ) then
        lambda0 = 4.5*m_pi/180.0
      else if ( longitude<21.0 ) then
        lambda0 = 15.0*m_pi/180.0
      else if ( longitude<33.0 ) then
        lambda0 = 27.0*m_pi/180.0
      else if ( longitude<42.0 ) then
        lambda0 = 37.5*m_pi/180.0
      end if
      ierr = 0
      return
    end if
    lambda0 = (float(zone_long-1)*6.0+(-180.0)+3.0)*m_pi/180.0
    ierr = 0
  end subroutine get_lambda0
 
  subroutine ll2utm(longitude,latitude,utm_x,utm_y,grid_zone,datum)
    implicit none
    integer :: datum
    real(8) :: latitude , longitude , utm_x , utm_y
    integer , dimension(2) :: grid_zone
    intent (in) longitude , latitude , datum
    intent (inout) :: grid_zone
    intent (out) utm_x , utm_y
    real(8) :: a , aa , aa2 , aa3 , aa4 , aa5 , aa6 , b , cc , e , e2 , e4 , &
              e6 , ep2 , f , k0 , lambda , lambda0 , mm , mm0 , m_pi , nn ,  &
             phi , phi0 , rho , t , tt , x , y
    m_pi = acos(-1.0)
    if ( datum==clarke_1866_datum ) then             ! CLARKE_1866_DATUM:
      a = 6378206.4
      b = 6356583.8
    else if ( datum==grs_80_datum ) then            ! GRS_80_DATUM:
      a = 6378137
      b = 6356752.3
    else if ( datum==wgs_84_datum ) then            ! WGS_84_DATUM:
      a = 6378137.0         !/* semimajor axis of ellipsoid (meters) */
      b = 6356752.31425     !/* semiminor axis of ellipsoid (meters) */
    else
      write (*,*) 'Unknown datum: ' , datum
      return
    end if
    f = 1 - (b/a)
    e2 = 2*f - f*f
    e = sqrt(e2)
    e4 = e2*e2
    e6 = e4*e2
    phi = latitude*m_pi/180.0
    lambda = longitude*m_pi/180.0
    call get_grid_zone(longitude,latitude,grid_zone,lambda0)
    phi0 = 0.0
    if ( latitude>84.0 ) then
      k0 = 0.994
      t = sqrt(((1-sin(phi))/(1+sin(phi))) * &
               (((1+e*sin(phi))/(1-e*sin(phi)))**e))
      rho = 2.0*a*k0*t/sqrt(((1.0+e)**(1.0+e))*((1.0-e)**(1.0-e)))
      x = rho*sin(lambda-lambda0)
      y = -rho*cos(lambda-lambda0)
      x = x + 2000000.0
      y = y + 2000000.0
    else if ( latitude<-80.0 ) then
      phi = -phi
      lambda = -lambda
      lambda0 = -lambda0
      k0 = 0.994
      t = sqrt(((1.0-sin(phi))/(1.0+sin(phi))) * &
               (((1.0+e*sin(phi))/(1.0-e*sin(phi))**e)))
      rho = 2.0*a*k0*t/sqrt(((1+e)**(1+e))*((1-e)**(1-e)))
      !!! Not needed (dhg) m = cos (phi) /  &
      !!!           sqrt (1.0 - e2 * sin (phi) * sin (phi))
      x = rho*sin(lambda-lambda0)
      y = -rho*cos(lambda-lambda0)
      !!! Not needed (dhg) k = rho * a * m
      x = -x
      y = -y
      x = x + 2000000.0
      y = y + 2000000.0
    else
      k0 = 0.9996
      mm = a*((1.0-e2/4.0-3.0*e4/64.0-5.0*e6/256.0) *                   &
             phi-(3.0*e2/8.0+3.0*e4/32.0+45.0*e6/1024.0)*sin(2.0*phi) + &
       (15.0*e4/256.0+45.0*e6/1024.0)*sin(4.0*phi)-(35.0*e6/3072.0) *   &
       sin(6.0*phi))
      mm0 = a*((1.0-e2/4.0-3.0*e4/64.0-5.0*e6/256.0) *                   &
            phi0-(3.0*e2/8.0+3.0*e4/32.0+45.0*e6/1024.0)*sin(2.0*phi0) + &
        (15.0*e4/256.0+45.0*e6/1024.0)*sin(4.0*phi0)-(35.0*e6/3072.0) *  &
        sin(6.0*phi0))
      aa = (lambda-lambda0)*cos(phi)
      aa2 = aa*aa
      aa3 = aa2*aa
      aa4 = aa2*aa2
      aa5 = aa4*aa
      aa6 = aa3*aa3
      ep2 = e2/(1.0-e2)
      nn = a/sqrt(1.0-e2*sin(phi)*sin(phi))
      tt = tan(phi)*tan(phi)
      cc = ep2*cos(phi)*cos(phi)
      x = k0*nn*(aa+(1-tt+cc)*aa3/6+(5-18*tt+tt*tt+72*cc-58*ep2)*aa5/120.0)
      y = k0*(mm-mm0+nn*tan(phi) *                &
         (aa2/2+(5-tt+9*cc+4*cc*cc)*aa4/24.0 +    &
         (61-58*tt+tt*tt+600*cc-330*ep2)*aa6/720))
      x = x + 500000.0
      if ( y<0.0 ) y = y + 10000000.0
    end if
    utm_x = x
    utm_y = y
  end subroutine ll2utm
 
  subroutine utm2ll(utm_x,utm_y,longitude,latitude,grid_zone,datum)
    implicit none
    integer :: datum
    real(8) :: latitude , longitude , utm_x , utm_y
    integer , dimension(2) :: grid_zone
    intent (in) datum , utm_x , utm_y
    intent (out) longitude , latitude
    real(8) :: a , b , cc1 , chi , dd , dd2 , dd3 , dd4 , dd5 , dd6 , e ,   &
               e1 , e12 , e13 , e14 , e2 , e4 , e6 , e8 , ep2 , f , k0 ,    &
               lambda , lambda0 , mm , mm0 , mu , m_pi , m_pi_2 , nn1 ,     &
               phi , phi0 , phi1 , phit , rho , rr1 , t , tt1 , x , y
    integer :: ierr
    m_pi = acos(-1.0)
    m_pi_2 = m_pi*2.0
    if ( datum==clarke_1866_datum ) then             ! CLARKE_1866_DATUM:
      a = 6378206.4
      b = 6356583.8
    else if ( datum==grs_80_datum ) then             ! GRS_80_DATUM:
      a = 6378137
      b = 6356752.3
    else if ( datum==wgs_84_datum ) then             ! WGS_84_DATUM:
      a = 6378137.0      !/* semimajor axis of ellipsoid (meters) */
      b = 6356752.31425  !/* semiminor axis of ellipsoid (meters) */
    else
      write (*,*) 'Unknown datum: ' , datum
      return
    end if
    f = 1.0 - (b/a)
    e2 = (2.0*f) - (f*f)
    e = sqrt(e2)
    e4 = e2*e2
    e6 = e4*e2
    e8 = e4*e4
    call get_lambda0(grid_zone,lambda0,ierr)
    if ( ierr/=0 ) then
      write (*,*) 'Unable to translate UTM to LL'
      return
    end if
    latitude = (float(grid_zone(2))*8.0) - 80.0
    if ( latitude>84.0 ) then           !/* north polar aspect */
      x = utm_x - 2000000.0
      y = utm_y - 2000000.0
      k0 = 0.994
      rho = sqrt(x*x+y*y)
      t = rho*sqrt(((1+e)**(1+e))*((1-e)**(1-e)))/(2*a*k0)
      chi = m_pi_2 - 2*atan(t)
      phit = chi + (e2/2+5*e4/24+e6/12+13*e8/360)*sin(2*chi) + &
             (7*e4/48+29*e6/240+811*e8/11520)*sin(4*chi)     + &
          (7*e6/120+81*e8/1120)*sin(6*chi) + (4279*e8/161280)*sin(8*chi)
      do while ( abs(phi-phit)>lower_eps_limit )
        phi = phit
        phit = m_pi_2 - 2*atan(t*(((1-e*sin(phi))/(1+e*sin(phi)))**(e/2)))
      end do
      lambda = lambda0 + atan2(x,-y)
    else if ( latitude<-80.0 ) then          !/* south polar aspect */
      x = -(utm_x-2000000)
      y = -(utm_y-2000000)
      k0 = 0.994
      rho = sqrt(x*x+y*y)
      t = rho*sqrt(((1+e)**(1+e))*((1-e)**(1-e)))/(2*a*k0)
      chi = m_pi_2 - 2*atan(t)
      phit = chi + (e2/2+5*e4/24+e6/12+13*e8/360)*sin(2*chi) + &
             (7*e4/48+29*e6/240+811*e8/11520)*sin(4*chi)     + &
          (7*e6/120+81*e8/1120)*sin(6*chi) + (4279*e8/161280)*sin(8*chi)
      do while ( abs(phi-phit)>lower_eps_limit )
        phi = phit
        phit = m_pi_2 - 2*atan(t*(((1-e*sin(phi))/(1+e*sin(phi)))**(e/2)))
      end do
      phi = -phi
      lambda = -(-lambda0+atan2(x,-y))
    else
      k0 = 0.9996
      x = utm_x - 500000.0
      y = utm_y
      if ( latitude<0.0 ) y = y - 10000000.0 !/* southern hemisphere */
      phi0 = 0.0
      e1 = (1.0-sqrt(1.0-e2))/(1.0+sqrt(1.0-e2))
      e12 = e1*e1
      e13 = e1*e12
      e14 = e12*e12
      mm0 = a*((1.0-e2/4.0-3.0*e4/64.0-5.0*e6/256.0) *                  &
        phi0-(3.0*e2/8.0+3.0*e4/32.0+45.0*e6/1024.0)*sin(2.0*phi0) +    &
        (15.0*e4/256.0+45.0*e6/1024.0)*sin(4.0*phi0)-(35.0*e6/3072.0) * &
        sin(6.0*phi0))
      mm = mm0 + y/k0
      mu = mm/(a*(1.0-e2/4.0-3.0*e4/64.0-5.0*e6/256.0))
      phi1 = mu + (3.0*e1/2.0-27.0*e13/32.0)*sin(2.0*mu) +              &
         (21.0*e12/16.0-55.0*e14/32.0)*sin(4.0*mu) + (151.0*e13/96.0) * &
        sin(6.0*mu) + (1097.0*e14/512.0)*sin(8.0*mu)
      ep2 = e2/(1.0-e2)
      cc1 = ep2*cos(phi1)*cos(phi1)
      tt1 = tan(phi1)*tan(phi1)
      nn1 = a/sqrt(1.0-e2*sin(phi1)*sin(phi1))
      rr1 = a*(1.0-e2)/((1.0-e2*sin(phi1)*sin(phi1))**1.5)
      dd = x/(nn1*k0)
      dd2 = dd*dd
      dd3 = dd*dd2
      dd4 = dd2*dd2
      dd5 = dd3*dd2
      dd6 = dd4*dd2
      phi = phi1 - (nn1*tan(phi1)/rr1) *                        &
        (dd2/2.0-(5.0+3.0*tt1+10.0*cc1-4.0*cc1*cc1-9.0*ep2) *   &
         dd4/24.0+(61.0+90.0*tt1+298.0*cc1+45.0*tt1*tt1-252.0 * &
         ep2-3.0*cc1*cc1)*dd6/720.0)
      lambda = lambda0 +                                             &
        (dd-(1.0+2.0*tt1+cc1)*dd3/6.0+(5.0-2.0*cc1+28.0*tt1-3.0*cc1* &
        cc1+8.0*ep2+24.0*tt1*tt1)*dd5/120.0)/cos(phi1)
    end if
    latitude = phi*180.0/m_pi
    longitude = lambda*180.0/m_pi
  end subroutine utm2ll
subroutine nearest2d(vari, xxi, yyi, varo, xxo, yyo, nb, nogeo, nxi, nyi, nxo, nyo)
    ! Nearest neighbour interpolation of two curvilinear grids

    implicit none

    ! Parameters
    integer,intent(in) :: nxi, nyi, nxo, nyo, nb
    real, intent(in) :: vari(nyi,nxi)
    real, intent(in) :: xxi(nyi,nxi),yyi(nyi,nxi),xxo(nyo,nxo),yyo(nyo,nxo)
    real,intent(out) :: varo(nyo,nxo)
    integer, intent(in), optional :: nogeo

    ! Local variables
    integer :: ixo,iyo, znb, ixlastline,iylastline,znb2,ixlast,iylast
    integer :: ixmin,ixmax,iymin,iymax, imin,jmin
    logical :: geo
    geo = .not. present(nogeo) .or. nogeo==0
    ! Loop on output points
    print*,"dove int 0"
    if(nb==0)then

        ! Scan all input points everytime
        do ixo = 1, nxo
            do iyo = 1, nyo
                call closest2d(xxi,yyi,xxo(iyo,ixo),yyo(iyo,ixo),nxi,nyi,imin,jmin,.not. geo)
                varo(iyo,ixo) = vari(jmin,imin)
            enddo
!             exit
        enddo

     else
        if(nb<0)then
            znb = 10
        else
            znb = max(2,nb)
        endif
        znb2 = znb/2
        ixlast = 1
        iylast = 1
        ixlastline = 1
        iylastline = 1
    print*,"dove int 1"
        do ixo = 1, nxo
            do iyo = 1, nyo

                ! Try a small block
                ixmin = max(1,ixlast-znb2)
                ixmax = min(nxi,ixlast+znb2)
                iymin = max(1,iylast-znb2)
                iymax = min(nyi,iylast+znb2)
    print*,"dove int 2"
                call closest2d(xxi(iymin:iymax,ixmin:ixmax), &
                    & yyi(iymin:iymax,ixmin:ixmax),xxo(iyo,ixo),yyo(iyo,ixo),nxi,nyi,imin,jmin,.not. geo)
    print*,"dove int 3"
                imin = imin+ixmin-1
                jmin = jmin+iymin-1

                ! Fall on bounds so use full block
                if((imin==ixmin.and.ixmin/=1).or.(imin==ixmax.and.ixmax/=nxi).or.&
                    & (jmin==iymin.and.iymin/=1).or.(jmin==iymax.and.iymax/=nyi))&
                    & call closest2d(xxi,yyi,xxo(iyo,ixo),yyo(iyo,ixo),nxi,nyi,imin,jmin,.not. geo)

                ! Store value
                varo(iyo,ixo) = vari(jmin,imin)
                ! Update min/max positions
                if(ixo==nxo)then
                    ixlast = ixlastline
                    iylast = iylastline
                else
                    ixlast = imin
                    iylast = jmin
                    if(ixo==1)then
                        ixlastline = ixlast
                        iylastline = iylast
                    endif
                endif
            enddo

         enddo
    endif
end subroutine nearest2d
subroutine closest2d(yyi,xxi,yo,xo,nxi,nyi,i,j,nogeo)
    ! Find indices of closest point on 2D axes
    implicit none
    real,intent(in) :: xxi(nxi,nyi),yyi(nxi,nyi),xo,yo
    integer, intent(in):: nyi,nxi
    integer,intent(out) :: i,j
    logical,intent(in) :: nogeo
    real :: dx(nxi,nyi)
    integer:: ij(2)

    if(nogeo)then
        dx = xxi-xo
    else
        dx = abs(xxi-xo)
        where(dx>180d0)dx = 360d0-dx
        dx = dx * cos(yyi*3.14159d0/180d0)
    endif
    ij = minloc(dx**2+(yyi-yo)**2)
    i = ij(2)
    j = ij(1)

end subroutine closest2d



end module mod_interp
