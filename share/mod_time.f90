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

module mod_time

  use mod_libmv
  use mod_strman
  use mod_phys
  use mod_internal

  contains

  integer function nslicetoskip(idate1,idate2)
    implicit none
    integer :: idate1 , idate2
    intent(in) :: idate1 , idate2
    integer :: censave , i , idate
    nslicetoskip = -1
    if ( indtimeal(idate1,idate2)==0 ) then
      nslicetoskip = 0
    else if ( indtimeal(idate1,idate2)==1 ) then
      idate = idate1
      call mvgetiflags(57,censave)
      do i = 1 , 1000000
        idate = increasemm5index(idate)
        if ( idate==idate2 ) then
          nslicetoskip = i
          call mvsetflags('Century',float(censave))
          return
        end if
      end do
      call mvsetflags('Century',float(censave))
    end if
  end function nslicetoskip
 
  integer function indtimeal(date1,date2)
    implicit none
    integer :: date1 , date2
    intent(in) :: date1 , date2
    integer :: anno1 , anno2 , giorno1 , giorno2 , &
               mese1 , mese2 , ora1 , ora2
    call gmafrommm5index(date1,ora1,giorno1,mese1,anno1)
    call gmafrommm5index(date2,ora2,giorno2,mese2,anno2)
    if ( anno2>anno1 ) then
      indtimeal = 1
    else if ( anno2<anno1 ) then
      indtimeal = -1
    else if ( mese2>mese1 ) then
      indtimeal = 1
    else if ( mese2<mese1 ) then
      indtimeal = -1
    else if ( giorno2>giorno1 ) then
      indtimeal = 1
    else if ( giorno2<giorno1 ) then
      indtimeal = -1
    else if ( ora2>ora1 ) then
      indtimeal = 1
    else if ( ora2<ora1 ) then
      indtimeal = -1
    else
      indtimeal = 0
    end if
  end function indtimeal
 
  integer function i4digityear(anno)
    implicit none
    integer :: anno
    intent (in) anno
    i4digityear = anno
    if ( i4digityear<1000 ) then
      if ( iflg(57)>0 ) then
        i4digityear = i4digityear + iflg(57)
      else if ( i4digityear<=50 ) then
        i4digityear = i4digityear + 2000
      else
        i4digityear = i4digityear + 1900
      end if
    end if
  end function i4digityear
 
  integer function julianday(id,im,iy)
    implicit none
    integer :: id , im , iy
    intent(in) :: id , im , iy
    julianday = index1d(id,im,iy)
  end function julianday
 
  integer function index1d(id,im,iy)
    implicit none
    integer :: id , im , iy
    intent (in) id , im , iy
    if ( di33/=33 .or. iy/=oldyear ) then
      oldyear = iy
      call definizionedimvtime
    end if
    index1d = dfc(im) + id
  end function index1d
 
  integer function index3h(ora,giorno,mese,anno)
    implicit none
    integer :: anno , giorno , mese , ora
    intent (in) anno , giorno , mese , ora
    index3h = (index1d(giorno,mese,anno)-1)*8 + ora/3 + 1
  end function index3h
 
  integer function index1h(ora,giorno,mese,anno)
    implicit none
    integer :: anno , giorno , mese , ora
    intent (in) anno , giorno , mese , ora
    index1h = (index1d(giorno,mese,anno)-1)*24 + ora + 1
  end function index1h
 
  integer function index15m(minuti,ora,giorno,mese,anno)
    implicit none
    integer :: anno , giorno , mese , minuti , ora
    intent (in) anno , giorno , mese , minuti , ora
    index15m = (index1d(giorno,mese,anno)-1)*4*24 + ora*4 + minuti/15 + 1
  end function index15m
 
  subroutine datafromindex15(dindex,anno,cdata)
    implicit none
    integer :: anno , dindex
    character(len=*) :: cdata
    intent(in) dindex
    intent(inout) :: anno
    intent(out) :: cdata
    integer :: giorno , mese , minuto , ora
    call dayfromindex15(dindex,minuto,ora,giorno,mese,anno)
    call datafrommin(minuto,ora,giorno,mese,anno,cdata)
  end subroutine datafromindex15
 
  subroutine dayfromindex1h(irec,ora,giorno,mese,anno)
    implicit none
    integer :: anno , giorno , mese , ora , irec
    intent (in) irec
    intent (out) anno , giorno , mese , ora
    if ( mod(irec,24)==0 ) then
      call dayfromindex(irec/24,giorno,mese,anno)
      ora = 23
    else
      call dayfromindex(irec/24+1,giorno,mese,anno)
      ora = irec - 24*(irec/24) - 1
    end if
  end subroutine dayfromindex1h
 
  subroutine dayfromindex15(dindex,minuto,ora,giorno,mese,anno)
    implicit none
    integer :: anno , giorno , dindex , mese , minuto , ora
    intent (in) dindex
    intent (out) minuto , ora , giorno , mese , anno
    integer :: nqu
    if ( dindex<=1 ) then
      minuto = 0
      ora = 0
      giorno = 1
      mese = 1
    else if ( mod(dindex-1,96)==0 ) then
      call dayfromindex(1+(dindex-1)/96,giorno,mese,anno)
      minuto = 0
      ora = 0
    else
      call dayfromindex(1+(dindex-1)/96,giorno,mese,anno)
      nqu = dindex - 1 - (index1d(giorno,mese,anno)-1)*96
      ora = nqu/4
      minuto = (nqu-ora*4)*15
    end if
  end subroutine dayfromindex15
 
  integer function indexofyear(id,im,iy)
    implicit none
    integer :: id , im , iy
    intent(in) :: id , im , iy
    indexofyear = index1d(id,im,iy)
  end function indexofyear

  integer function index15ofyear(minuto,ora,gior,mese,anno)
    implicit none
    integer :: anno , gior , mese , minuto , ora
    intent(in) :: anno , gior , mese , minuto , ora
    index15ofyear = index15m(minuto,ora,gior,mese,anno)
  end function index15ofyear

  integer function indexhourofyear(ora,giorno,mese,anno)
    implicit none
    integer :: anno , giorno , mese , ora
    intent(in) anno , giorno , mese , ora
    indexhourofyear = index1h(ora,giorno,mese,anno)
  end function indexhourofyear

  integer function monthlen(im,iy)
    implicit none
    integer :: im , iy
    intent (in) im , iy
    if ( di33/=33 .or. iy/=oldyear ) then
      oldyear = iy
      call definizionedimvtime
    end if
    monthlen = mesi(im)
  end function monthlen
 
  subroutine datafromday(id,im,iy,cdata)
    implicit none
    character(len=*) :: cdata
    integer :: id , im , iy
    call datafromindex(indexofyear(id,im,iy),iy,cdata)
  end subroutine datafromday
 
  integer function decreasemm5index(dindex)
    implicit none
    integer :: dindex
    intent (in) dindex
    integer :: id , ih , im , iy
    character(len=10) :: tmps
    write (tmps,'(i10)') dindex
    read (tmps,'(i4,3i2)') iy , im , id , ih
    ih = ih - 1
    if ( ih==-1 ) then
      ih = 23
      id = id - 1
      if ( id<1 ) then
        oldyear = iy
        call definizionedimvtime
        im = im - 1
        id = mesi(im)
        if ( im<1 ) then
          im = 12
          id = 31
          iy = iy - 1
          oldyear = iy
          call definizionedimvtime
        end if
      end if
    end if
    decreasemm5index = ih + id*100 + im*10000 + iy*1000000
  end function decreasemm5index
 
  subroutine tomorrow(id,im,iy)
    implicit none
    integer :: id , im , iy
    intent (inout) id , im , iy
    oldyear = iy
    call definizionedimvtime
    id = id + 1
    if ( id>mesi(im) ) then
      id = 1
      im = im + 1
      if ( im==13 ) then
        im = 1
        iy = iy + 1
      end if
    end if
  end subroutine tomorrow
 
  integer function increasemm5index(dindex)
    implicit none
    integer :: dindex
    intent (in) dindex
    integer :: id , ih , im , iy
    character(len=8) :: tmps
    if ( dindex==99123123 ) then
      increasemm5index = 00010100
      return
    end if
    write (tmps,'(i8)') dindex
    read (tmps,'(4i2)') iy , im , id , ih
    ih = ih + 1
    if ( ih==24 ) then
      ih = 0
      id = id + 1
      oldyear = iy
      call definizionedimvtime
      if ( id>mesi(im) ) then
        id = 1
        im = im + 1
        if ( im==13 ) then
          im = 1
          iy = iy + 1
        end if
      end if
    end if
    if ( ih<0 .or. ih>23 .or. im<1 .or. im>12 .or. id<1 .or. &
         id>monthlen(im,iy) .or. iy<0 .or. iy>99 ) then
      increasemm5index = -9999
    else
      increasemm5index = ih + id*100 + im*10000 + iy*1000000
    end if
  end function increasemm5index
 
  integer function increasetime(dindex)
    implicit none
    integer :: dindex
    intent (in) dindex
    integer :: ical , id , ih , iii , ilm , im , iy
    character(10) :: tmps
    call mvgetiflags(72,ical)
    if ( dindex==1999123123 .or. (ical==2 .and. dindex==1999123023) ) then
      increasetime = 2000010100
      call mvgetiflags(57,iii)
      if ( iii>0 ) call mvsetflags('Century',float(iii+100))
      return
    end if
    write (tmps,'(i10)',err=100) dindex
    read (tmps,'(i4,3i2)',err=100) iy , im , id , ih
    if ( ical==1 ) then
      ilm = monthlen(im,2001)
    else if ( ical==2 ) then
      ilm = 30
    else
      ilm = monthlen(im,iy)
    end if
    ih = ih + 1
    if ( ih==24 ) then
      ih = 0
      id = id + 1
      if ( id>ilm ) then
        id = 1
        im = im + 1
        if ( im==13 ) then
          im = 1
          iy = iy + 1
        end if
      end if
    end if
    if ( ih<0 .or. ih>23 .or. im<1 .or. im>12 .or. id<1 .or. id>ilm .or.    &
         iy<1900 .or. iy>2099 ) then
      increasetime = -9999
    else
      increasetime = ih + id*100 + im*10000 + iy*1000000
    end if
    return
 100  increasetime = -9999
  end function increasetime
 
  integer function mm5index(ih,id,im,iy)
    implicit none
    integer :: id , ih , im , iy
    intent (in) id , ih , im , iy
    integer :: lyp
    if ( iy>2000 ) then
      lyp = iy - 2000
    else if ( iy==2000 ) then
      lyp = 0
    else if ( iy>1900 ) then
      lyp = iy - 1900
    else
      lyp = iy
    end if
    mm5index = ih + id*100 + im*10000 + lyp*1000000
  end function mm5index
 
  subroutine datafromhour(ih,id,im,iy,cdata)
    implicit none
    character(len=*) :: cdata
    integer :: id , ih , im , iy
    intent(in) :: id , ih , im , iy
    intent(out) :: cdata
    call dataorafromday(ih,id,im,iy,cdata)
  end subroutine datafromhour
 
  subroutine datafrommin(imm,ih,id,im,iy,cdata)
    implicit none
    character(len=*) :: cdata
    integer :: id , ih , im , imm , iy
    intent (in) id , ih , im , imm , iy
    intent (out) cdata
    call dataorafromday(ih,id,im,iy,cdata)
    if ( imm<10 ) then
      write (cdata,'(a,i1)') trim(cdata)//'.0' , imm
    else
      write (cdata,'(a,i2)') trim(cdata)//'.' , imm
    end if
    if ( iflg(10)==4 ) cdata = trim(cdata)//' UGT'
    if ( iflg(10)==5 ) cdata = trim(cdata)//' CEST'
  end subroutine datafrommin
 
  subroutine dataorafromday(ih,id,im,iy,cdata)
    implicit none
    character(len=*) :: cdata
    integer :: id , ih , im , iy
    intent (in) id , ih , im , iy
    intent (out) cdata
    call datafromindex(indexofyear(id,im,iy),iy,cdata)
    if ( ih>=10 ) then
      write (cdata,'(a,i2)') trim(cdata)//' h: ' , ih
    else
      write (cdata,'(a,i1)') trim(cdata)//' h: 0' , ih
    end if
    call no2space(cdata)
    if ( iflg(10)==4 ) cdata = trim(cdata)//' UGT'
    if ( iflg(10)==5 ) cdata = trim(cdata)//' CEST'
  end subroutine dataorafromday
 
  subroutine datafrommm5index(dindex,cdata)
    implicit none
    character(len=*) :: cdata
    integer :: dindex
    intent (in) dindex
    intent (out) cdata
    integer :: id , ih , im , iy
    character(len=10) :: tmps
    write (tmps,'(i0.8)') dindex
    read (tmps,'(4i2)') iy , im , id , ih
    if ( iflg(57)<0 ) then
      iy = iy + 1900
      if ( iy<1950 ) iy = iy + 100
    else
      iy = iy + iflg(57)
    end if
    call datafromindex(indexofyear(id,im,iy),iy,cdata)
    if ( iflg(10)==0 .or. iflg(10)==1 ) then
      write (tmps,'(a5,i3)') ' ore: ' , ih
      cdata = trim(cdata)//tmps
    else if ( iflg(10)==2 ) then
      write (cdata,'(a,1x,i2,a)') trim(cdata)//' h: ' , ih , '.00'
    else if ( iflg(10)==3 ) then
      write (tmps,'(a5,i3)') ' h: ' , ih
      cdata = trim(cdata)//tmps
    else if ( iflg(10)==4 ) then
      write (cdata,'(a,1x,i2,a)') trim(cdata) , ih , ':00 GMT'
    else if ( iflg(10)==5 ) then
      write (cdata,'(a,1x,i2,a)') trim(cdata) , ih , ':00 CEST'
    else if ( iflg(10)==7 ) then
      write (cdata,'(a,i2)') emonths(im)(1:3)//' ' , id
    else if ( iflg(10)==8 ) then
      write (cdata,'(a,i2)') emonths(im)(1:3) , id
      call nospace(cdata)
    else if ( iflg(10)==9 ) then
      call integer2string(id,2,tmps)
      write (cdata,'(a,i5)') tmps(1:2)//' '//emonths(im)(1:3) , iy
    end if
    call no2space(cdata)
    call noinspace(cdata)
  end subroutine datafrommm5index

  subroutine datafromidx(dindex,cdata)
    implicit none
    character(len=*) :: cdata
    integer :: dindex
    intent (in) dindex
    intent (out) cdata
    integer :: id , ih , im , iy
    character(len=10) :: tmps
    character(len=3), dimension(12) :: months
    data months/'JAN','FEB','MAR','APR','MAJ','JUN','JUL','AUG','SEP','OCT','NOV','DEC'/
    write (tmps,'(i0.10)') dindex
    read (tmps,'(i4,3i2)') iy , im , id , ih
    if ( iflg(10)==0 .or. iflg(10)==1 ) then
      write (tmps,'(a5,i3)') ' ore: ' , ih
      cdata = trim(cdata)//tmps
    else if ( iflg(10)==2 ) then
      write (cdata,'(a,1x,i2,1x,i4)')months(im), id, iy
      write (cdata,'(a,1x,i2,a)') trim(cdata)//' h: ' , ih , '.00'
    else if ( iflg(10)==3 ) then
      write (tmps,'(a5,i3)') ' h: ' , ih
      cdata = trim(cdata)//tmps
    else if ( iflg(10)==4 ) then
      write (cdata,'(a,1x,i2,a)') trim(cdata) , ih , ':00 GMT'
    else if ( iflg(10)==5 ) then
      write (cdata,'(a,1x,i2,a)') trim(cdata) , ih , ':00 CEST'
    else if ( iflg(10)==7 ) then
      write (cdata,'(a,i2)') emonths(im)(1:3)//' ' , id
    else if ( iflg(10)==8 ) then
      write (cdata,'(a,i2)') emonths(im)(1:3) , id
      call nospace(cdata)
    else if ( iflg(10)==9 ) then
      call integer2string(id,2,tmps)
      write (cdata,'(a,i5)') tmps(1:2)//' '//emonths(im)(1:3) , iy
    end if
    call no2space(cdata)
    call noinspace(cdata)
  end subroutine datafromidx
 
  integer function indexfrommm5index(dindex)
    implicit none
    integer :: dindex
    intent (in) dindex
    integer :: id , ih , im , iy
    character(8) :: tmps
    write (tmps,'(i8)') dindex
    read (tmps,'(4i2)') iy , im , id , ih
    if ( iflg(57)<0 ) then
      iy = iy + 1900
      if ( iy<1950 ) iy = iy + 100
    else
      iy = iy + iflg(57)
    end if
    indexfrommm5index = indexofyear(id,im,iy)
  end function indexfrommm5index
 
  subroutine gmafrommm5index(dindex,ih,id,im,iy)
    implicit none
    integer :: id , ih , im , dindex , iy
    intent (in) dindex
    intent (out) id , ih , im , iy
    character(8) :: tmps
    write (tmps,'(i8)',err=100) dindex
    read (tmps,'(4i2)',err=100) iy , im , id , ih
    if ( iflg(57)<0 ) then
      iy = iy + 1900
      if ( iy<1950 ) iy = iy + 100
    else
      iy = iy + iflg(57)
    end if
    return
 100  iy = -9999
      im = -9999
      id = -9999
      ih = -9999
  end subroutine gmafrommm5index

  subroutine gmafromindex(dindex,ih,id,im,iy)
    implicit none
    integer :: id , ih , im , dindex , iy
    intent (in) dindex
    intent (out) id , ih , im , iy
    character(10) :: tmps
    write (tmps,'(i10)',err=100) dindex
    read (tmps,'(i4,3i2)',err=100) iy , im , id , ih
!    if ( iflg(57)<0 ) then
!      iy = iy + 100
!      if ( iy<1950 ) iy = iy + 100
!    else
!      iy = iy + iflg(57)
!    end if
    return
 100  iy = -9999
      im = -9999
      id = -9999
      ih = -9999
  end subroutine gmafromindex
 
  subroutine monthofyear(im,iy,cdata)
    implicit none
    character(len=*) :: cdata
    integer :: im , iy
    intent (in) im , iy
    intent (out) cdata
    if ( iflg(10)>=2 .or. iflg(10)<=5 ) then
      write (cdata,'(a)') emonths(im)
    else if ( iflg(10)<2 ) then
      write (cdata,'(a)') months(im)
    else
      write (cdata,'(a)') emonths(im)(1:3)
    end if
    if ( iy>0 ) write (cdata,'(a,1x,i4)') trim(cdata) , iy
  end subroutine monthofyear
 
  subroutine nomedelmese(im,nome)
    implicit none
    integer :: im
    character(len=*) :: nome
    intent (in) im
    intent (out) nome
    if ( iflg(10)>=2 .and. iflg(10)<=5 ) then
      nome = emonths(im)
    else if ( iflg(10)<2 ) then
      nome = months(im)
    else
      nome = emonths(im)(1:3)
    end if
  end subroutine nomedelmese
 
  subroutine datafromindex(dindex,anno,cdata)
    implicit none
    integer :: anno , dindex
    character(len=*) :: cdata
    integer :: id , im , iy
    iy = i4digityear(anno)
    id = -1
    im = -1
    call dayfromindex(dindex,id,im,iy)
    if ( id>0 .and. im>0 ) then
      if ( iflg(10)==1 ) then
        write (cdata,'(i2,1x,a9,i5)') id , months(im) , iy
      else if ( iflg(10)==2 ) then
        write (cdata,'(a,i2,a,i5)') trim(emonths(im))//' ' , id , ', ' , iy
      else if ( iflg(10)>=3 .and. iflg(10)<=5 ) then
        write (cdata,'(a11,1x,a9,1x,i2,i5)') eweekd(iweekday(id,im,iy)) , &
                                             emonths(im) , id , iy
      else if ( iflg(10)==6 ) then
        write (cdata,'(a,i2,a,i5)') emonths(im)(1:3)//' ' , id , ', ' , iy
      else if ( iflg(10)==7 ) then
        write (cdata,'(a4,i2)') emonths(im)(1:3)//' ' , id
      else if ( iflg(10)==8 ) then
        write (cdata,'(a3,i2)') emonths(im)(1:3) , id
      else
        write (cdata,'(a11,1x,i2,1x,a9,i5)') weekd(iweekday(id,im,iy)) , id , &
                                             months(im) , iy
      end if
      call no2space(cdata)
      call noinspace(cdata)
    end if
  end subroutine datafromindex
 
  subroutine dayfromindex(dindex,id,im,iy)
    implicit none
    integer :: id , im , dindex , iy
    intent (in) dindex , iy
    intent (out) id , im
    integer :: i
    if ( di33/=33 .or. iy/=oldyear ) then
      oldyear = iy
      call definizionedimvtime
    end if
    do i = 1 , 12
      if ( dindex>dfc(i) .and. dindex<=dfc(i+1) ) then
        im = i
        id = dindex - dfc(i)
        return
      end if
    end do
    write (6,'(a,i10)') ' Bad index passed to dayfromindex ' , dindex
  end subroutine dayfromindex
 
  integer function iweekday(id,im,iy)
    implicit none
    integer :: id , im , iy
    intent (in) id , im , iy
    integer :: igiorni
    if ( di33/=33 .or. iy/=oldyear ) then
      oldyear = iy
      call definizionedimvtime
    end if
    if ( mod(iy,4)==0 ) then
      igiorni = id + dfc(im) + iy*365 + iy/4 + 3
    else
      igiorni = id + dfc(im) + iy*365 + iy/4 + 4
    end if
    iweekday = mod(igiorni,7) + 1
  end function iweekday
 
  integer function notteogiorno(lat,lon,ore,minuti,giorno,mese)
    implicit none
    integer :: giorno , mese , minuti , ore
    real :: lat , lon
    intent (in) giorno , mese , minuti , ore , lat , lon
    real :: pi , zlat , zlon
    pi = 4.*atan(1.0)
    zlat = (sessa2centi(23,45,0)) * &
      sin((indexofyear(giorno,mese,1999)-indexofyear(21,3,1999))*pi/180.)
    zlon = ((ore+minuti/60.)-12.)*180./12.
    if ( distance(lat,lon,zlat,zlon)<0.5*pi*6371000. ) then
      notteogiorno = 1
    else
      notteogiorno = 0
    end if
  end function notteogiorno
 
  subroutine localtime2ugt(ora,giorno,mese,anno)
    implicit none
    integer :: anno , giorno , mese , ora
    intent (inout) anno , giorno , mese , ora
    integer :: dh , iw
    if ( mod(anno,4)==0 ) then
      mesi(2) = 29
    else
      mesi(2) = 28
    end if
    if ( mese<=2 .or. mese>=11 ) then
      dh = 1
    else if ( mese==3 ) then
      if ( giorno<=24 ) then
        dh = 1
      else
        iw = iweekday(giorno,mese,anno)
        if ( iw==7 ) then
          dh = 2
        else if ( giorno-iw<=24 ) then
          dh = 1
        else
          dh = 2
        end if
      end if
    else if ( mese==10 ) then
      if ( giorno<=24 ) then
        dh = 2
      else
        iw = iweekday(giorno,mese,anno)
        if ( iw==7 ) then
          dh = 1
        else if ( giorno-iw<=24 ) then
          dh = 2
        else
          dh = 1
        end if
      end if
    else
      dh = 2
    end if
    ora = ora - dh
    if ( ora<0 ) then
      ora = ora + 24
      giorno = giorno - 1
      if ( giorno<1 ) then
        mese = mese - 1
        if ( mese==0 ) then
          mese = 12
          anno = anno - 1
        end if
        giorno = mesi(mese)
      end if
    end if
  end subroutine localtime2ugt
 
  subroutine ugt2localtime(ora,giorno,mese,anno)
    implicit none
    integer :: anno , giorno , mese , ora
    intent (inout) anno , giorno , mese , ora
    integer :: dh , iw
    if ( mod(anno,4)==0 ) then
      mesi(2) = 29
    else
      mesi(2) = 28
    end if
    if ( mese<=2 .or. mese>=11 ) then
      dh = 1
    else if ( mese==3 ) then
      if ( giorno<=24 ) then
        dh = 1
      else
        iw = iweekday(giorno,mese,anno)
        if ( iw==7 ) then
          if ( ora<=3 ) then
            dh = 1
          else
            dh = 2
          end if
        else if ( giorno-iw<=24 ) then
          dh = 1
        else
          dh = 2
        end if
      end if
    else if ( mese==10 ) then
      if ( giorno<=24 ) then
        dh = 2
      else
        iw = iweekday(giorno,mese,anno)
        if ( iw==7 ) then
          if ( ora<=3 ) then
            dh = 2
          else
            dh = 1
          end if
        else if ( giorno-iw<=24 ) then
          dh = 2
        else
          dh = 1
        end if
      end if
    else
      dh = 2
    end if
    ora = ora + dh
    if ( ora>23 ) then
      ora = ora - 24
      giorno = giorno + 1
      if ( giorno>mesi(mese) ) then
        giorno = 1
        mese = mese + 1
        if ( mese==13 ) then
          mese = 1
          anno = anno + 1
        end if
      end if
    end if
  end subroutine ugt2localtime
 
  subroutine whattimeisit(minuti,ora,giorno,mese,anno)
    implicit none
    integer :: anno , giorno , mese , minuti , ora
    intent (out) anno , giorno , mese , minuti , ora
    call cheorae(minuti,ora,giorno,mese,anno)
  end subroutine whattimeisit
 
  subroutine cheorae(minuti,ora,giorno,mese,anno)
    implicit none
    integer :: anno , giorno , mese , minuti , ora
    intent (out) anno , giorno , mese , minuti , ora
    integer(4) , dimension(8) :: tarray
    call date_and_time(values=tarray)
    anno = tarray(1)
    mese = tarray(2)
    giorno = tarray(3)
    ora = tarray(5)
    minuti = tarray(6)
  end subroutine cheorae

  real function sinphi(lat,lon,minuti,ore,giorno,mese)
    implicit none
    integer :: giorno , mese , minuti , ore
    real :: lat , lon
    intent (in) lat , lon , minuti , ore , giorno , mese
    real :: ds , dt , phir , pi , xlat , xlon
    data pi , phir /3.1415927,0.409/
    save pi , phir
    ds = phir*cos(2.0*pi*(indexofyear(giorno,mese,1999)-173.)/365.0)
    xlat = deg2rad(lat)
    xlon = deg2rad(lon)
    dt = float((ore+1)*60+minuti)
    sinphi = sin(xlat)*sin(ds) - cos(xlat)*cos(ds) * &
             cos(2.0*pi*dt/(24.*60.)-xlon)
  end function sinphi
 
  real function zenith(lat,lon,minuti,ore,giorno,mese)
    implicit none
    integer :: giorno , mese , minuti , ore
    real :: lat , lon
    intent (in) lat , lon , minuti , ore , giorno , mese
    real :: dst , pi , rad , zlat , zlon
    data rad , pi /6371000.,3.1415927/
    save rad , pi
    zlat = (sessa2centi(23,45,0)) * sin((indexofyear(giorno,mese,1999) - &
      indexofyear(21,3,1999))*pi/180.)
    zlon = ((ore+minuti/60.)-12.)*180./12.
    dst = distance(lat,lon,zlat,zlon)
    if ( dst>0.5*pi*rad ) then
      zenith = pi/2.
    else
      zenith = dst/rad
    end if
  end function zenith

  subroutine convert2hourlyres(y,isen)
    implicit none
    real , dimension(:) :: y
    intent (inout) y
    integer :: isen , mese , giorno , ora , i15m , i1h , n , i
    real :: xm
    do mese = 1 , 12
      do giorno = 1 , monthlen(mese,2008)
        do ora = 0 , 23
          i1h = index1h(ora,giorno,mese,2008)
          i15m = index15m(00,ora,giorno,mese,2008)
          xm = 0.0
          n = 0
          do i = 0 , 3
            if ( nint(y(i15m+i))/=-9999 ) then
              xm = xm + y(i15m+i)
              n = n + 1
            end if
          end do
          if ( n>0 ) then
            y(i1h) = xm
            if (isen/=1) y(i1h) = y(i1h) / n
          else
            y(i1h) = -9999
          end if
        end do
      end do
    end do
  end subroutine convert2hourlyres

  character(len=36) function get_nowdate() result(cdate)
    implicit none
    integer , dimension(8) :: tvals
    character(len=5) :: tzone
    call date_and_time(zone=tzone,values=tvals)
    write(cdate,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,a)') &
      tvals(1), '-', tvals(2), '-', tvals(3), ' ', &
      tvals(5), ':', tvals(6), ':', tvals(7), ' ', tzone
  end function get_nowdate

  subroutine chkcentury(string)
    implicit none 
    integer anno 
    character(len=4) :: yst
    character(len=*) :: string

    do anno=1900,1999
       write(yst,'(i4)') anno
       if (index(string(1:4),yst).ne.0) then
          call mvsetflags('Century',1900.)
          write(6,'(15x,a)') 'Assuming dates refer to XX century.'
       endif
    enddo
    do anno=2000,2099
       write(yst,'(i4)') anno
       if (index(string(1:4),yst).ne.0) then
          call mvsetflags('Century',2000.)
          write(6,'(15x,a)') 'Assuming dates refer to XXI century.'
       endif
    enddo
    return
  end subroutine chkcentury

  subroutine index_to_date(indx,year,month,day,hour)
    implicit none
    integer , intent(in) :: indx
    integer , intent(out) :: year , month , day , hour
    year = indx/1000000
    month = (indx-year*1000000)/10000
    day = (indx-year*1000000-month*10000)/100
    hour = indx-year*1000000-month*10000-day*100
    if ( year > 50 ) then
      year = year + 1900
    else
      year = year + 2000
    end if
  end subroutine index_to_date

end module mod_time
