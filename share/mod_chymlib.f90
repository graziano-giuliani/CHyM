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
module mod_chymlib

  use mod_statparams
  use mod_libmv
  use mod_runparams
  use mod_museo
  use mod_phys

  contains

  subroutine plotriverbasin
    implicit none
    logical okplot
!    integer :: iriv,nsec,isec(12000000),jsec(12000000)
    character(len=40) :: rivtitle
    real, dimension(nlon,nlat) :: wk2, work
    real :: cut
    integer :: iriver
    return
    iriver = mchym(18)
    if ( iriver == 0 ) iriver = 2
    okplot=.true.
    cut=rchym(6)/20.0
    call chymmouthsfinder(drai,luse,fmap,nlon,nlat,cut,nsec,isec,jsec)
    ibas=isec(iriver) ; jbas=jsec(iriver) ; call salmone(drai,20.0)
    wk2=dem
    call basinpaint(dem,fmap,luse,work,nlon,nlat,isec(iriver),jsec(iriver),lavonc,-1)
    do iriver=1,nsec
!       call chymmouthsfinder(drai,luse,fmap,nlon,nlat,cut,nsec,isec,jsec)
!       ibas=isec(iriver) ; jbas=jsec(iriver) ; call salmone(drai,20.0)
!       if (drai(isec(iriver),jsec(iriver)) >= MAXVAL(drai)/2) then
         call basinpaint(dem,fmap,luse,work,nlon,nlat,isec(iriver),jsec(iriver),allrivb,iriver)
!       end if
    end do
  end subroutine plotriverbasin

  subroutine chymmouthsfinder(w,luse,fmap,nlon,nlat,cut,nsec,isec,jsec)
    implicit none
    integer mare,lago,fiume ; parameter (mare=15,lago=14)
    integer nlon,nlat,nsec,isec(:),jsec(:)
    real w(nlon,nlat),cut
    integer luse(nlon,nlat),fmap(nlon,nlat)
    integer i,j,idir,ir(9),jr(9) ; save ir,jr
    data ir /-1, 0, 1, 1, 1, 0,-1,-1,0/
    data jr / 1, 1, 1, 0,-1,-1,-1, 0,0/
    nsec=0
    do j=nlat-1,2,-1
      do i=2,nlon-1
        idir=fmap(i,j)
        if (idir.ge.1.and.idir.le.8) then
          if (w(i,j).gt.cut.and. &
                (luse(i+ir(idir),j+jr(idir)).eq.mare.or. &
                 i+ir(idir).eq.1.or.i+ir(idir).eq.nlon.or. &
                 j+jr(idir).eq.1.or.j+jr(idir).eq.nlat)) then
             nsec=nsec+1
!            isec(nsec)=i          ! Vedi sul Logbook i commenti del 9/12/10
!            jsec(nsec)=j
             isec(nsec)=i+ir(idir)
             jsec(nsec)=j+jr(idir)
          endif
        endif
      end do
    end do
    return
  end subroutine chymmouthsfinder

  subroutine salmone(w,cut)            ! risale la corrente come i salmoni
    implicit none
    real, dimension(nlon,nlat) :: w
    real :: cut,xmax,dis,xla,xlo
    integer, dimension(nlon+nlat) :: jlat,ilon
    integer :: i,j,idir,jdir,k
!    integer :: seqi(1000),seqj(1000),seqx(1000),seqy(1000)
    character(len=60) :: river,fiume
    character(len=25) :: test
    logical :: scorri
    integer :: iriver
    iriver = mchym(18)
    nrivp=1
    jlat(nrivp)=jbas
    ilon(nrivp)=ibas
    scorri=.true.
    if (ibas.le.0.or.jbas.le.0) return
    do while (scorri)
       scorri=.false.
       i=ilon(nrivp)
       j=jlat(nrivp)
       xmax=-1.0e36
       jdir=0
       do idir=1,8
          if    (w(i+ir(idir),j+jr(idir)).gt.xmax.and. &
             w(i+ir(idir),j+jr(idir)).lt.w(i,j).and. &
             w(i+ir(idir),j+jr(idir)).gt.cut.and. &
          luse(i+ir(idir),j+jr(idir)).ne.mare.and. &
          luse(i+ir(idir),j+jr(idir)).ne.lago.and. &
          fmap(i+ir(idir),j+jr(idir)).ne.0) then
              xmax=w(i+ir(idir),j+jr(idir))
              jdir=idir
              scorri=.true.
           endif
       enddo
       do idir=1,8
          if    (idir.ne.jdir.and. &
             w(i+ir(idir),j+jr(idir)).lt.w(i,j).and. &
             w(i+ir(idir),j+jr(idir)).gt.cut.and. &
          luse(i+ir(idir),j+jr(idir)).ne.mare.and. &
          luse(i+ir(idir),j+jr(idir)).ne.lago.and. &
          fmap(i+ir(idir),j+jr(idir)).ne.0) then
          nsec=nsec+1
          isec(nsec)=i+ir(idir)
          jsec(nsec)=j+jr(idir)
          endif
       enddo
       if (scorri) then
          nrivp=nrivp+1
          ilon(nrivp)=ilon(nrivp-1)+ir(jdir)
          jlat(nrivp)=jlat(nrivp-1)+jr(jdir)
       endif
    enddo
    do i=1,nrivp
       seqi(i)=ilon(nrivp+1-i)
       seqj(i)=jlat(nrivp+1-i)
       seqx(i)=lon(seqi(i),seqj(i))
       seqy(i)=lat(seqi(i),seqj(i))
    enddo
    call openmuseofiles(10,'river.basin',1) ; fiume='Unknown' ; i=iriver
    do k=1,1000
       read (10,'(a25,25x,2f9.5)',end=100) test,xla,xlo
       dis=distance(lat(isec(i),jsec(i)),lon(isec(i),jsec(i)),xla,xlo)
       if (dis.le.2000.0) fiume=test
    enddo
100 close(10)
    return
  end subroutine salmone

  subroutine basinpaint(dem,fmap,luse,work,nlon,nlat,isa,jsa,wk1,flg)
    implicit none
    integer :: isa , jsa , nlat , nlon , flg
    real , dimension(nlon,nlat) :: dem , wk1 , work
    integer , dimension(nlon,nlat) :: fmap , luse
    intent (in) dem , fmap , isa , jsa , luse , nlat , nlon , flg
    intent (out) wk1
    intent (inout) work
    integer :: i , idir , ii , iind , j , jj
    logical :: notyet
    if ( iflg(1)>=10.0 ) &
      write (6,'(12x,a,i3,a,i3)')                    &
         'Basin Selecion Module called for grid '//'point ' , isa , '-', jsa
    if ( flg==1 .or. flg==0 ) then
      wk1 = -10
      work = 0.0
    endif
    do j = 2 , nlat - 1
      do i = 2 , nlon - 1
        if ( flg==-1 ) then
          if ( drai(i,j) > 80. ) then
!             wk1(i,j)=-5
          end if
        end if
        ii = i
        jj = j
        iind = 0
        notyet = .true.
        do while ( iind<=2*(nlon+nlat) .and. notyet )
          iind = iind + 1
          idir = fmap(ii,jj)
          if ( idir>0 ) then
            ii = ii + ir(idir)
            jj = jj + jr(idir)
            if ( ii==isa .and. jj==jsa ) then
              notyet = .false.
              if ( flg==0 ) then
                wk1(i,j) = dem(i,j)
              else if ( flg==-1 ) then
                wk1(i,j) = -999
              else
                wk1(i,j) = flg
!                if (drai(i,j) > 80.) wk1(i,j)=-5
              end if
            else if ( ii==1 .or. ii==nlon .or. jj==1 .or. jj==nlat ) then
              notyet = .false.                     ! Fuori dominio
            end if
!           else if (luse(ii,jj).eq.mare.or.luse(ii,jj).eq.lago) then
          else if ( luse(ii,jj)==mare ) then              ! In mare
            notyet = .false.
          else                               ! noflow
            wk1(i,j) = -5.0
            if (flg.ne.0) wk1(i,j)=-10
            work(ii,jj) = work(ii,jj) + 1
            notyet = .false.
          end if
        end do
      end do
    end do
    if ( iflg(1)>=10 ) then
      do i = 2 , nlon - 1
        do j = 2 , nlat - 1
          if ( work(i,j)>10. ) write (6,'(/,3x,a,i3,a,i3,a,i4,a,/)') &
               'Severe warning: The no flow point '&
               , i , '-' , j , ' drains ' ,        &
               nint(work(i,j)) , ' cells.'
        end do
      end do
    end if
  end subroutine basinpaint

  subroutine rollingstones2(rntime,fmap,cout,nlon,nlat)
    implicit none
    integer :: nlat , nlon
    integer , dimension(nlon,nlat) :: fmap
    real , dimension(nlon,nlat) :: cout , rntime
    intent (in) fmap , nlat , nlon , rntime
    intent (inout) cout
    integer :: di , dj , i , idir , ist , j , jst , n , ncycle
    logical :: exceed
    data ncycle/10000/
    save ncycle
    exceed = .false.
    cout = rntime
    do i = 1 , nlon
      do j = 1 , nlat
        ist = i
        jst = j
        idir = fmap(ist,jst)
        do n = 1 , ncycle
!          if ( i<=1 .or. i>=nlon .or. j<=1 .or. j>=nlat .or. idir<1 .or. &
!               idir>8 ) exit
          if (i.gt.1.and.i.lt.nlon.and.j.gt.1.and.j.lt.nlat.and.idir.ge. &
               1.and.idir.le.8) then
            call neighborhood(idir,di,dj)
            ist = ist + di
            jst = jst + dj
            cout(ist,jst) = cout(ist,jst) + rntime(i,j)
            idir = fmap(ist,jst)
          else
            exit
          end if
        end do
        if ( n>=ncycle ) exceed = .true.
      end do
    end do
    if ( exceed ) write (6,'(/,10x,a/)')                                &
                          'CHyM severe warning: rollingstones2'//       &
                         ' maximum number of cycles exceeded'
  end subroutine rollingstones2

  subroutine runofftime(rntime,fmap,rout,nlon,nlat)
    implicit none
    integer :: nlat , nlon
    integer , dimension(nlon,nlat) :: fmap
    real , dimension(nlon,nlat) :: rout , rntime
    intent (in) fmap , nlat , nlon , rntime
    intent (inout) rout
    logical :: exceed
    integer :: i , idir , ist , j , jst , n , ncycle
    data ncycle/10000/
    save ncycle
    exceed = .false.
    do i = 1 , nlon
      do j = 1 , nlat
        rout(i,j) = rntime(i,j)/2
        ist = i
        jst = j
        idir = fmap(ist,jst)
        do n = 1 , ncycle
          if ( i>1 .and. i<nlon .and. j>1 .and. j<nlat .and. &
               idir>=1 .and. idir<=8 ) then
            call chymmoveahead(idir,ist,jst)
            rout(i,j) = rout(i,j) + rntime(ist,jst)
            idir = fmap(ist,jst)
          else
            rout(i,j) = rout(i,j)/3600
            exit
          end if
        end do
        if ( n>=ncycle ) exceed = .true.
      end do
    end do
    if ( exceed ) &
      write (6,'(/,10x,a/)') 'CHyM severe warning: runofftime'// &
               ' maximum number of cycles exceeded'
  end subroutine runofftime

  subroutine chymmoveahead(idir,i,j)
    implicit none
    integer :: i , idir , j
    intent (inout) i , j
    integer :: di , dj
    call neighborhood(idir,di,dj)
    i = i + di
    j = j + dj
  end subroutine chymmoveahead

  subroutine neighborhood(id,i,j)
    implicit none
    integer :: i , id , j
    intent (in) id
    intent (out) i , j
    integer , dimension(0:48) :: ir , jr
    data ir(0:8)/0 , -1 , 0 , 1 , 1 , 1 , 0 , -1 , -1/
    data jr(0:8)/0 , 1 , 1 , 1 , 0 , -1 , -1 , -1 , 0/
    data ir(9:24)/ - 2 , -1 , 0 , 1 , 2 , 2 , 2 , 2 , &
                     2 , 1 , 0 , -1 , -2 , -2 , -2 , -2/
    data jr(9:24)/2 , 2 , 2 , 2 , 2 , 1 , 0 , -1 , -2 , &
                  -2 , -2 , -2 , -2 , -1 , 0 , 1/
    data ir(25:48)/ - 3 , -2 , -1 , 0 , 1 , 2 , 3 , 3 , &
                    3 , 3 , 3 , 3 , 3 , 2 , 1 , 0 , -1 , &
                    -2 , -3 , -3 , -3 , -3 , -3 , -3/
    data jr(25:48)/3 , 3 , 3 , 3 , 3 , 3 , 3 , 2 , 1 , 0 , &
                  -1 , -2 , -3 , -3 , -3 , -3 , -3 , -3 ,  &
                  -3 , -2 , -1 , 0 , 1 , 2/
    if ( id<0 .or. id>48 ) then
      write (6,'(10x,a,i4,a)') &
         'Invalid id inside subroutine neighborhood: ', id , '. Exiting...'
      call exit(0)
    end if
    i = ir(id)
    j = jr(id)
  end subroutine neighborhood

  subroutine smoothHydroDEM
    implicit none
    integer :: idir,iii,ii,jj,j,k,numpoints,run,counter
    integer :: iimaxs,jjmaxs,iimax,jjmax
    real    :: maxdra,maxsubdra,demmin,increment
    real, dimension(10000000) :: demval,demval_run,demval_int
    integer, dimension(10000000) :: iriv, jriv
    real, dimension(nlon,nlat) :: drai_t
    drai_t = drai
    do k = 1,numrivdr
    numpoints = 0
    i = 0
    maxdra = 0
    do jj = 1,nlon
       do ii = 1,nlat
          if (drai_t(jj,ii) > maxdra) then
             maxdra = drai_t(jj,ii)
             iimax = ii
             jjmax = jj
          end if
       end do
    end do
    if (maxdra<threshdr) run = 0
    jjmaxs = jjmax
    iimaxs = iimax
    do while(drai_t(jjmax,iimax)>threshdr)
      i=i+1
      drai_t(jjmax,iimax) = 0
      maxsubdra = 0
      jjmaxs = jjmax
      iimaxs = iimax
      do idir = 1,8
        if (fmap(jjmaxs+jr(idir),iimaxs+ir(idir)).ne.0.and.  &
         iimaxs+ir(idir).gt.1.and.iimaxs+ir(idir).lt.nlat.and. &
         jjmaxs+jr(idir).gt.1.and.jjmaxs+jr(idir).lt.nlon) then
         if (drai_t(jjmaxs+jr(idir),iimaxs+ir(idir)) > maxsubdra) then
           maxsubdra = drai_t(jjmaxs+jr(idir),iimaxs+ir(idir))
           iimax = iimaxs+ir(idir)
           jjmax = jjmaxs+jr(idir)
         end if
        end if
      end do
      numpoints = numpoints + 1
      demval(numpoints) = dem(jjmax,iimax)
      iriv(numpoints) = iimax
      jriv(numpoints) = jjmax
    end do
    do iii = 1,numrunave
    do i = numpoints, 1, -1
      if (i < numpoints - uphill) then
      if (i <= 3) then
        dem(jriv(i),iriv(i)) = (dem(jriv(1),iriv(1))+dem(jriv(2),iriv(2))+dem(jriv(3),iriv(3))+dem(jriv(4),iriv(4))    &
          +dem(jriv(5),iriv(5))+dem(jriv(6),iriv(6))+dem(jriv(7),iriv(7))+dem(jriv(8),iriv(8))+dem(jriv(9),iriv(9)))/9
      else if (i >= numpoints-3) then
        dem(jriv(i),iriv(i)) = (dem(jriv(numpoints-1),iriv(numpoints-1))+dem(jriv(numpoints-2),iriv(numpoints-2))      &
          +dem(jriv(numpoints-3),iriv(numpoints-3))+dem(jriv(numpoints-4),iriv(numpoints-4))                           &
          +dem(jriv(numpoints-5),iriv(numpoints-5))+dem(jriv(numpoints-6),iriv(numpoints-6))                           &
          +dem(jriv(numpoints-7),iriv(numpoints-7))+dem(jriv(numpoints-8),iriv(numpoints-8))                           &
          +dem(jriv(numpoints-9),iriv(numpoints-9)))/9
      else
        dem(jriv(i),iriv(i)) = (dem(jriv(i-4),iriv(i-4))+dem(jriv(i-3),iriv(i-3))+dem(jriv(i-2),iriv(i-2))             &
          +dem(jriv(i-1),iriv(i-1))+dem(jriv(i),iriv(i))+dem(jriv(i+1),iriv(i+1))+dem(jriv(i+2),iriv(i+2))             &
          +dem(jriv(i+3),iriv(i+3))+dem(jriv(i+4),iriv(i+4)))/9
      end if
      end if
      demval_run(i) = dem(jriv(i),iriv(i))
    end do
    end do
    counter = 0
    demmin = 0
    do i = numpoints, 2, -1
      if (dem(jriv(i),iriv(i))<dem(jriv(i-1),iriv(i-1)) .and. counter == 0) then
        counter = counter + 1
        demmin = dem(jriv(i),iriv(i))
        cycle
      end if
      if (counter > 0) then
        if (demmin <= dem(jriv(i-1),iriv(i-1)) .and. i > 2) then
          counter = counter + 1
        else if (demmin <= dem(jriv(i-1),iriv(i-1)) .and. i == 2) then
          dem(jriv(i-1),iriv(i-1)) = demmin - 0.5
          increment = (demmin - dem(jriv(i-1),iriv(i-1))) / counter
          do j = 0, counter-1
            dem(jriv(i+j),iriv(i+j))=dem(jriv(i-1),iriv(i-1))+(increment*(j+1))
          end do
          counter = 0
          demmin = 0
        else
          increment = (demmin - dem(jriv(i-1),iriv(i-1))) / counter
          do j = 0, counter-1
            dem(jriv(i+j),iriv(i+j))=dem(jriv(i-1),iriv(i-1))+(increment*(j+1))
          end do
          counter = 0
          demmin = 0
        end if
      end if
    end do
    do i = numpoints, 1, -1
      demval_int(i) = dem(jriv(i),iriv(i))
    end do
    end do
    return
  end subroutine smoothHydroDEM

  subroutine smoothHydroDEM1
    implicit none
    integer :: idir,iii,ii,jj,is,js,k,kk,numpoints,numrivers,numrunave1,run,counter
    integer :: iimaxs,jjmaxs,iimax,jjmax
    real    :: maxdra,maxsubdra,treshdr,demmin,increment,demsum
    real, dimension(10000000) :: demval,demval_run,demval_int
    integer, dimension(10000000) :: iriv, jriv
    real, dimension(nlon,nlat) :: drai_t
    treshdr = 4
    numrunave1 = 1
    do kk = 1, numrunave1
    do jj=1,nlat
      do ii=1,nlon
        numpoints = 1
        demsum = 0
        if (drai(ii,jj)>treshdr) then
          demsum = dem(ii,jj)
          idir = fmap(ii,jj)
          if (idir>8 .or. idir<1) cycle
          is = ii+ir(idir) ; js = jj+jr(idir)
          if (is == 0 .or. is == nlon+1 .or. js == 0 .or. js == nlat+1) cycle
          do k = 1 , 3
          if(luse(is,js).ne.mare.and. &
            luse(is,js).ne.lago.and. &
            fmap(is,js).ne.0) then
              demsum = demsum + dem(is,js)
              numpoints = numpoints+1
              idir = fmap(is,js)
              if (idir>8 .or. idir<1) exit
              is = is+ir(idir) ; js = js+jr(idir)
              if (is == 0 .or. is == nlon+1 .or. js == 0 .or. js == nlat+1) exit
          else
            exit
          end if
          end do
          dem(ii,jj) = demsum/numpoints
        end if
      end do
    end do
    end do
  end subroutine smoothHydroDEM1
  subroutine smoothHydroDEM2
    implicit none
    integer :: idir,iii,ii,jj,is,js,k,kk,numpoints,numrivers,numrunave1,run,counter
    integer :: iimaxs,jjmaxs,iimax,jjmax
    real    :: maxdra,maxsubdra,treshdr,demmin,increment,demsum,weight
    real, dimension(10000000) :: demval,demval_run,demval_int
    integer, dimension(10000000) :: iriv, jriv
    real, dimension(nlon,nlat) :: drai_t
    treshdr = 150
    numrunave1 = 1
    weight = 0.5
    do kk = 1, numrunave
    do jj=1,nlat
      do ii=1,nlon
        numpoints = 1
        demsum = 0
        if (drai(ii,jj)>treshdr .and. accl(ii,jj)>0.002) then
          demsum = dem(ii,jj)
          idir = fmap(ii,jj)
          if (idir>8 .or. idir<1) cycle
          is = ii+ir(idir) ; js = jj+jr(idir)
          if (is == 0 .or. is == nlon+1 .or. js == 0 .or. js == nlat+1) cycle
          do k = 1 , 5
          if(luse(is,js).ne.mare.and. &
            luse(is,js).ne.lago.and. &
            fmap(is,js).ne.0) then
              demsum = demsum + dem(is,js) * weight
              numpoints = numpoints+1
              idir = fmap(is,js)
              if (idir>8 .or. idir<1) exit
              is = is+ir(idir) ; js = js+jr(idir)
              if (is == 0 .or. is == nlon+1 .or. js == 0 .or. js == nlat+1) exit
          else
            exit
          end if
          end do
          dem(ii,jj) = demsum/numpoints
        end if
      end do
    end do
    end do
  end subroutine smoothHydroDEM2

end module mod_chymlib
