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
module mod_internal

  use mod_libmv
  use mod_stats
  use mod_strman

  integer , private :: icall
  integer , private , dimension(100) :: iflgs
  real , private , dimension(100) :: rflgs

  contains

  subroutine getlun(lun)
    implicit none
    integer :: lun
    intent (out) lun
    integer :: i
    logical :: test
    do i = 99 , 10 , -1
      inquire(unit=i,opened=test)
      if ( .not.test ) then
        lun = i
        if ( iflg(1)>=10 ) then
          write (6,'(10x,a,i2,a)') 'MVLib: getlun find ' ,  &
                           i , ' as free unit.'
        end if
        return
      end if
    end do
    write (6,'(10x,a)') 'MVLib error: igetlun search failed. Exiting...'
    call exit(0)
  end subroutine getlun

  subroutine mvsetflags(flag,rvalue)
    implicit none
    character(len=*) :: flag
    real :: rvalue
    intent (in) rvalue
    integer :: i , nn
    real :: xx
    character(60) :: iflag
    character(len=32) , dimension(100) :: istrings
    data istrings(1:39)/39*'Noval'/
    data istrings(40)/'plotta color'/
    data istrings(41)/'tipo di carattere'/
    data istrings(42)/'colore freccetelle'/
    data istrings(43)/'colore sfondo plot'/
    data istrings(45)/'boundaries behaviour'/
    data istrings(46:52)/7*'Noval'/
    data istrings(53)/'histogram color'/
    data istrings(54)/'Settata da makepalette'/
    data istrings(55)/'Settata da plottaframe'/
    data istrings(56)/'size delle label'/
    data istrings(57)/'century'/
    data istrings(58)/'video control'/
    data istrings(59)/'displayplot behaviour'/
    data istrings(60)/'som neural init'/
    data istrings(61)/'som neural distance'/
    data istrings(62)/'openmuseofiles behaviour'/
    data istrings(63)/'Usata insieme alla 62'/
    data istrings(64)/'boxplot style'/
    data istrings(65)/'chymdownscaling cycles'/
    data istrings(66)/'html table color'/
    data istrings(67)/'histogram background'/
    data istrings(68)/'histogram statistic'/
    data istrings(69)/'plotta fill'/
    data istrings(70)/'chym log unit'/
    data istrings(71)/'chym standard domain'/
    data istrings(72)/'calendar'/
    data istrings(73:100)/28*' '/

    call cv2lower(flag,iflag)
    call noinspace(iflag)
    call no2space(iflag)
    nn = len_trim(iflag)
    do i = 1 , 100
      if ( trim(iflag)==trim(istrings(i)) ) then
        iflg(i) = nint(rvalue)
        return
      end if
    end do

    if ( iflag(1:nn)=='log level' ) then
      iflg(1) = nint(rvalue)
      write (6,'(7x,a,i4)') 'MVLib: log level setted to ' , iflg(1)
!!      if ( iflg(1)>0 ) call mvlibversion
    else if ( iflag(1:nn)=='x size' ) then
      iflg(2) = nint(rvalue)
      if ( iflg(2)<100 ) iflg(2) = 100
      if ( iflg(2)>2000 ) iflg(2) = 2000
    else if ( iflag(1:nn)=='plot 1d' ) then
      iflg(4) = nint(rvalue)
    else if ( iflag(1:nn)=='domain' .or. iflag(1:nn)=='dominio' ) then
      if ( nint(rvalue)>=0 .and. nint(rvalue)<=3 ) then
        iflg(5) = nint(rvalue)
      else
        write (6,'(7x,a)') 'MVLib: Bad value passwd for '//flag
      end if
    else if ( iflag(1:nn)=='current month' ) then
      iflg(6) = nint(rvalue)
      call definizionedimvtime
    else if ( iflag(1:nn)=='check landuse' ) then
      iflg(7) = nint(rvalue)
    else if ( iflag(1:nn)=='logo' ) then
      iflg(8) = nint(rvalue)
      if ( iflg(8)==2 ) iflg(11) = 1
    else if ( iflag(1:nn)=='html style' ) then
      iflg(9) = nint(rvalue)
    else if ( iflag(1:nn)=='data style' ) then
      iflg(10) = nint(rvalue)
    else if ( iflag(1:nn)=='confini politici' ) then
      write (6,'(7x,a)') 'MVLib: flag confini politici is now obselete, '&
                   //'see flag plot italy.'
    else if ( iflag(1:nn)=='colore titoli' ) then
      iflg(12) = nint(rvalue)
    else if ( iflag(1:nn)=='colore assi' ) then
      iflg(13) = nint(rvalue)
    else if ( iflag(1:nn)=='colore freccette' ) then
      iflg(14) = nint(rvalue)
    else if ( iflag(1:nn)=='reset asse x' ) then
      iflg(15) = nint(rvalue)
    else if ( iflag(1:nn)=='colore confini' ) then
      iflg(16) = nint(rvalue)
    else if ( iflag(1:nn)=='confini ncar' ) then
      iflg(17) = nint(rvalue)
    else if ( iflag(1:nn)=='misura vento' ) then
      iflg(18) = nint(rvalue)
    else if ( iflag(1:nn)=='colore isolinee' ) then
      iflg(19) = nint(rvalue)
    else if ( iflag(1:nn)=='curva di livello' ) then
      iflg(20) = iflg(20) + 1
      if ( iflg(20)>3 .or. iflg(20)<1 ) iflg(20) = 1
      rflg(iflg(20)+4) = rvalue
    else if ( iflag(1:nn)=='numero di curve' ) then
      iflg(20) = nint(rvalue)
    else if ( iflag(1:nn)=='size delle scritte' ) then
      iflg(21) = nint(rvalue)
    else if ( iflag(1:nn)=='coordinate grafiche' ) then
      iflg(22) = nint(rvalue)
    else if ( iflag(1:nn)=='indice latlon' ) then
      iflg(23) = nint(rvalue)
    else if ( iflag(1:nn)=='cellular step' ) then
      iflg(24) = nint(rvalue)
    else if ( iflag(1:nn)=='x left' ) then
      rflg(1) = rvalue
    else if ( iflag(1:nn)=='x right' ) then
      rflg(2) = rvalue
    else if ( iflag(1:nn)=='y bottom' ) then
      rflg(3) = rvalue
    else if ( iflag(1:nn)=='y top' ) then
      rflg(4) = rvalue
    else if ( iflag(1:nn)=='colore sfondo' ) then
      iflg(25) = nint(rvalue)
    else if ( iflag(1:nn)=='numero di colori' ) then
      iflg(26) = nint(rvalue)
    else if ( iflag(1:nn)=='palette di colori' ) then
      iflg(27) = nint(rvalue)
    else if ( iflag(1:nn)=='contorno ai pallocchi' ) then
      iflg(28) = nint(rvalue)
    else if ( iflag(1:nn)=='label bar' ) then
      iflg(29) = nint(rvalue)
    else if ( iflag(1:nn)=='distruggi fotogrammi' ) then
      iflg(30) = nint(rvalue)
    else if ( iflag(1:nn)=='asse x' ) then
      iflg(31) = nint(rvalue)
      call definizionedimvtime
    else if ( iflag(1:nn)=='numero colori plot2d' ) then
      lmeth = nint(rvalue)
    else if ( iflag(1:nn)=='colori plot2d' ) then
      if ( nint(rvalue)==-9999 ) then
        iflg(32) = 0
        return
      end if
      if ( iflg(32)>=100 ) iflg(32) = 0
      iflg(32) = iflg(32) + 1
      rlevels(iflg(32)) = rvalue
    else if ( iflag(1:nn)=='step freccette' ) then
      iflg(33) = nint(rvalue)
    else if ( iflag(1:nn)=='random seed' ) then
      iflg(34) = nint(rvalue)
      if ( iflg(34)==0 ) write (6,'(//,5x,a,//)') &
         'mvsetflags severe error: Random seed setted to zero cause a ' &
          //'zero sequence.'
      xx = setgauss(iflg(34))
      xx = setacaso(iflg(34))
    else if ( iflag(1:nn)=='qnorm number of bin' ) then
      iflg(35) = nint(rvalue)
    else if ( iflag(1:nn)=='asse y' ) then
      iflg(36) = nint(rvalue)
    else if ( iflag(1:nn)=='primo colore' ) then
      iflg(37) = nint(rvalue)
    else if ( iflag(1:nn)=='ogni quanti colori' ) then
      iflg(38) = nint(rvalue)
    else if ( iflag(1:nn)=='scrivi max-min' ) then
      iflg(39) = nint(rvalue)
    else if ( iflag(1:nn)=='colore plotta' ) then
      iflg(40) = nint(rvalue)
    else if ( iflag(1:nn)=='size delle linee' ) then
      rflg(8) = rvalue
    else if ( iflag(1:16)=='cellular weights' ) then
      if ( iflag(nn-1:nn)==' n' ) then
        rflg(9) = rvalue
      else if ( iflag(nn-1:nn)==' e' ) then
        rflg(10) = rvalue
      else if ( iflag(nn-1:nn)==' s' ) then
        rflg(11) = rvalue
      else if ( iflag(nn-1:nn)==' o' ) then
        rflg(12) = rvalue
      else if ( iflag(nn-1:nn)=='ne' ) then
        rflg(13) = rvalue
      else if ( iflag(nn-1:nn)=='se' ) then
        rflg(14) = rvalue
      else if ( iflag(nn-1:nn)=='so' ) then
        rflg(15) = rvalue
      else if ( iflag(nn-1:nn)=='no' ) then
        rflg(16) = rvalue
      else
        write (6,'(7x,a)') 'MVLib: Bad libmv parameter: '//iflag(1:nn)
      end if
    else if ( iflag(1:nn)=='size del carattere' ) then
      rflg(17) = rvalue
    else if ( iflag(1:nn)=='visualizzaplot format' .or. &
              iflag(1:nn)=='displayplot format' ) then
      iflg(46) = nint(rvalue)
    else if ( iflag(1:nn)=='java window size' ) then
      iflg(47) = nint(rvalue)
    else if ( iflag(1:nn)=='plot italy' .or. &
              iflag(1:nn)=='plot confini' ) then
      iflg(48) = nint(rvalue)
    else if ( iflag(1:nn)=='label direction' ) then
      iflg(51) = nint(rvalue)
    else if ( iflag(1:nn)=='stile freccette' ) then
      iflg(52) = nint(rvalue)
    else if ( iflag(1:nn)=='colore istogramma' ) then
      iflg(53) = nint(rvalue)
    else if ( iflag(1:nn)=='visualizzaplot behaviour' ) then
      iflg(59) = nint(rvalue)
    else
      write (6,'(7x,a)') 'MVLib: Bad MVLib parameter: '//iflag(1:nn)
    end if
  end subroutine mvsetflags

  subroutine definizionedimvtime ! internal mvlib use
    implicit none
    integer :: i
    if ( mod(oldyear,4)==0 ) then
      mesi(2) = 29
    else
      mesi(2) = 28
    end if
    dfc(1) = 0
    do i = 2 , 13
      dfc(i) = dfc(i-1) + mesi(i-1)
    end do
    di33 = 33
  end subroutine definizionedimvtime

  subroutine mvgetiflags(comp,ivalue)
    implicit none
    integer :: comp , ivalue
    intent (in) comp
    intent (out) ivalue
    ivalue = iflg(comp)
  end subroutine mvgetiflags

  integer function mvlibmagicnum(i)
    implicit none
    include '../doc/settings.inc'
    integer , intent(in) :: i
    mvlibmagicnum=-9999
    if ( i==1 ) mvlibmagicnum=N01
    if ( i==2 ) mvlibmagicnum=1997
    if ( i==3 ) mvlibmagicnum=AANNO
    if ( i==4 ) mvlibmagicnum=1990
    if ( i==5 ) mvlibmagicnum=2009
    if ( i==6 ) mvlibmagicnum=AANNO
    if ( i==7 ) mvlibmagicnum=2012
    return
  end function mvlibmagicnum

  subroutine locateij(xla,xlo,slat,slon,dlat,dlon,nlat,nlon,i,j)
    implicit none
    real :: dlat , dlon , slat , slon , xla , xlo
    integer :: i , j , nlat , nlon
    intent (in) dlat , dlon , nlat , nlon , slat , slon , xla , xlo
    intent (inout) i , j
    i = nint((xlo-slon)/dlon) + 1
    j = nint((xla-slat)/dlat) + 1
    if ( i<=0 .or. j<=0 .or. i>nlon .or. j>nlat ) i = -1
  end subroutine locateij

  subroutine museopath
    implicit none
    include '../doc/settings.inc'
    pathm = MPATH
    npathm = len_trim(pathm)
  end subroutine museopath

  subroutine mvliberror(caller,msg,icode,rcode)
    implicit none
    character(len=*) :: caller , msg
    integer :: icode
    real :: rcode
    intent (in) icode , rcode , caller , msg
    character(len=78) :: str
    write (6,99001)
    write (6,'(11x,a)') 'Severe error from '//trim(caller) &
                     //' subroutine.'
    if ( len_trim(msg)>0 ) write (6,'(11x,a)') msg
    if ( icode/=-9999 ) then
      write (str,'(a,i12)') 'Value of integer parameter was ' , icode
      call no2space(str)
      write (6,'(11x,a)') trim(str)
    end if
    if ( nint(rcode)/=-9999 ) then
      if ( abs(rcode)<1000.0 ) then
        write (str,'(a,f10.6)') 'Value of real parameter is ' , rcode
      else if ( abs(rcode)<1000000.0 ) then
        write (str,'(a,f10.3)') 'Value of real parameter is ' , rcode
      else
        write (str,'(a,e14.6)') 'Value of real parameter is ' , rcode
      end if
      call no2space(str)
      write (6,'(11x,a)') trim(str)
    end if
    write (6,'(11x,a)') 'Exiting from fortran code with status=1'
    write (6,99001)
    call exit(1)
99001 format (11x,50('*'))
  end subroutine mvliberror

  subroutine mvsaveflags
    implicit none
    icall = 100
    iflgs = iflg
    rflgs = rflg
  end subroutine mvsaveflags

  subroutine mvrestflags
    implicit none
    if ( icall/=100 ) then
      write (6,'(7x,a)') 'MVLib Error: mvrestflags called before mvsaveflags'
    else
      iflg = iflgs
      rflg = rflgs
    end if
  end subroutine mvrestflags

  subroutine resetmvlibint
    implicit none
    isub = -9999
    rsub = -9999.0
  end subroutine resetmvlibint

  subroutine setmvlibintr(i,r)
    implicit none
    integer :: i
    real :: r
    intent (in) i , r
    if ( i>=1 .and. i<=nintpar ) rsub(i) = r
  end subroutine setmvlibintr

  subroutine setmvlibinti(i,j)
    implicit none
    integer :: i , j
    intent (in) i , j
    if ( i>=1 .and. i<=nintpar ) isub(i) = j
  end subroutine setmvlibinti

  real function getmvlibint(i)
    implicit none
    integer i
    if (i.ge.1.and.i.le.nintpar) then
       getmvlibint=rsub(i)
    else
       getmvlibint=-9999.0
    endif
    return
  end function getmvlibint

  integer function igetmvlibint(i)
    implicit none
    integer i
    if (i.ge.1.and.i.le.nintpar) then
       igetmvlibint=isub(i)
    else
       igetmvlibint=-9999.0
    endif
    return
  end function igetmvlibint

end module
