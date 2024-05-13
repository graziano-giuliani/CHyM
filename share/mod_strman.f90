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

module mod_strman

  contains

  subroutine noinspace(a)
    implicit none
    character(len=*) , intent(inout) :: a
    a = adjustl(a)
  end subroutine noinspace

  subroutine nospace(a)
    implicit none
    character(len=*) :: a
    intent (inout) a
    integer :: i , n
    nospace_loop: &
    do
      n = len_trim(a)
      do i = 1 , n
        if ( ichar(a(i:i))==32 ) then
          a = a(1:i-1)//a(i+1:)
          cycle nospace_loop
        end if
      end do
      exit
    end do nospace_loop
  end subroutine nospace

  subroutine space2undescore(str)
    implicit none
    character(len=*) :: str
    intent (inout) str
    integer :: i , ln
    character(132) :: lstr
    ln = len_trim(str)
    do i = 1 , ln
      if ( str(i:i)==' ' ) then
        lstr = str(1:i-1)//'_'//str(i+1:)
        str = lstr(1:ln)
      end if
    end do
  end subroutine space2undescore
 
  subroutine tab2space(b)
    implicit none
    character(len=*) :: b
    intent (inout) b
    integer :: i
    do i = 1 , len_trim(b)
      if ( ichar(b(i:i))==9 ) b = b(1:i-1)//' '//b(i+1:)
    end do
  end subroutine tab2space
 
  subroutine integer2string(ii,ilen,cfile)
    implicit none
    character(len=*) :: cfile
    integer :: ii , ilen
    intent (in) ii , ilen
    intent (inout) :: cfile
    character(80) :: lfile
    integer :: lt
    write (lfile,'(i29)') ii
    call noinspace(cfile)
    lfile = '000000000000000000000'//lfile(1:len_trim(lfile))
    call nospace(lfile)
    lt = len_trim(lfile)
    cfile = lfile(lt-ilen+1:lt)
  end subroutine integer2string
 
  subroutine latlon2str2(lat,lon,str)
    implicit none
    real :: lat , lon
    character(len=*) :: str
    intent (in) lat , lon
    intent (out) str
    if ( lat>0 ) then
      write (str,'(f7.4,a)') lat , 'N-'
    else if ( lat<0 ) then
      write (str,'(f7.4,a)') abs(lat) , 'S-'
    else
      write (str,'(f7.4,a)') 'Eq-'
    end if
    if ( lon>0 ) then
      write (str,'(a,f7.4,a)') str(1:len_trim(str)) , lon , 'E'
    else
      write (str,'(a,f7.4,a)') str(1:len_trim(str)) , abs(lon) , 'W'
    end if
    call nospace(str)
  end subroutine latlon2str2
 
  subroutine latlon2str(lat,lon,str)
    implicit none
    real :: lat , lon
    character(len=*) :: str
    intent (in) lat , lon
    intent (out) str
    if ( lat>0 ) then
      write (str,'(f7.4,a)') lat , '\S\o\N\N-'
    else
      write (str,'(f7.4,a)') abs(lat) , '\S\o\N\S-'
    end if
    if ( lon>0 ) then
      write (str,'(a,f7.4,a)') str(1:len_trim(str)) , lon , '\S\o\N\E'
    else
      write (str,'(a,f7.4,a)') str(1:len_trim(str)) , abs(lon) , '\S\o\N\W'
    end if
    call nospace(str)
  end subroutine latlon2str
 
  subroutine lon2str(lon,str)
    implicit none
    real :: lon
    character(len=*) :: str
    intent (in) lon
    intent (out) str
    if ( lon>0 .and. lon<180.01 ) then
      write (str,'(f5.1,a)') lon , 'E'
    else if ( abs(lon)<180.01 ) then
      write (str,'(f5.1,a)') abs(lon) , 'W'
    else
      str = '*****'
    end if
    call nospace(str)
  end subroutine lon2str
 
  subroutine lat2str(lat,str)
    implicit none
    real :: lat
    character(len=*) :: str
    intent (in) lat
    intent (out) str
    if ( lat>0 .and. lat<90.01 ) then
      write (str,'(f4.1,a)') lat , 'N'
    else if ( abs(lat)<90.01 ) then
      write (str,'(f4.1,a)') abs(lat) , 'S'
    else
      str = '*****'
    end if
    call nospace(str)
  end subroutine lat2str
 
  subroutine strsub(a,b,c)
    implicit none
    character(len=*) :: a , b , c
    intent (in) :: a , b
    intent (inout) c
    if ( index(c,a)==0 ) return
    if ( len_trim(b)==0 .and. ichar(c(1:1))==48 ) then
      c = c(1:index(c,a)-1)//b(1:1)//c(index(c,a)+len_trim(a):)
    else
      c = c(1:index(c,a)-1)//b(1:len_trim(b))//c(index(c,a)+len_trim(a):)
    end if
  end subroutine strsub
 
  subroutine no2space(a)
    implicit none
    character(len=*) :: a
    intent (inout) a
    integer :: i
    no2space_loop: &
    do
      do i = 1 , len_trim(a) - 1
        if ( ichar(a(i:i))==32 .and. ichar(a(i+1:i+1))==32 ) then
          a = a(1:i)//a(i+2:)
          cycle no2space_loop
        end if
      end do
      exit
    end do no2space_loop
  end subroutine no2space
 
  subroutine cv2lower(inpstr,cout)
    implicit none
    character(len=*) :: inpstr , cout
    intent (in) :: inpstr
    intent (out) cout
    character(1) :: c
    integer :: i
    character(1024) :: outstr
    do i = 1 , len_trim(inpstr)
      if ( ichar(inpstr(i:i))>=65 .and. ichar(inpstr(i:i))<=90 ) then
        c = char(ichar(inpstr(i:i))+32)
      else
        c = inpstr(i:i)
      end if
      if ( i==1 ) then
        outstr = c
      else
        outstr = outstr(1:i-1)//c
      end if
    end do
    if ( len_trim(inpstr)>0 ) then
      cout = outstr
    else
      cout = inpstr
    end if
  end subroutine cv2lower
 
  subroutine cv2upper(inpstr,cout)
    implicit none
    character(len=*) :: inpstr , cout
    intent (in) :: inpstr
    intent (out) cout
    character(1) :: c
    integer :: i
    character(1024) :: outstr
    do i = 1 , len_trim(inpstr)
      if ( ichar(inpstr(i:i))>=97 .and. ichar(inpstr(i:i))<=122 ) then
        c = char(ichar(inpstr(i:i))-32)
      else
        c = inpstr(i:i)
      end if
      if ( i==1 ) then
        outstr = c
      else
        outstr = outstr(1:i-1)//c
      end if
    end do
    if ( len_trim(inpstr)>0 ) then
      cout = outstr
    else
      cout = inpstr
    end if
  end subroutine cv2upper
 
  subroutine makefilename(cstart,n,cend,cfile)
    implicit none
    character(len=*) :: cend , cfile , cstart
    integer :: n
    intent (in) cstart , cend , n
    intent (out) cfile
    if ( n>999 ) then
      write (cfile,'(a,i4,a)') cstart(1:len_trim(cstart)) , n , cend
    else if ( n>99 ) then
      write (cfile,'(a,i3,a)') cstart(1:len_trim(cstart)) , n , cend
    else if ( n>9 ) then
      write (cfile,'(a,i2,a)') cstart(1:len_trim(cstart))//'0' , n , cend
    else
      write (cfile,'(a,i1,a)') cstart(1:len_trim(cstart))//'00' , n , cend
    end if
  end subroutine makefilename
 
  subroutine getcommandline(prom,chin,chlen)
    implicit none
    character(len=*) :: chin , prom
    integer :: chlen
    intent (in) prom
    intent (out) chlen
    intent (inout) chin
    character(12) :: cform
    integer :: i , linp , lpro
    lpro = len(prom)
    linp = len(chin)
    write (cform,'(a1,i5)') 'a' , lpro
    call nospace(cform)
    write (cform,'(a2,i5,a2)') '(a' , lpro , '$)'
    call nospace(cform)
    write (6,cform) prom
    write (cform,'(a2,i5,a1)') '(a' , linp , ')'
    call nospace(cform)
    read (5,cform) chin
    chlen = 0
    do i = 1 , linp
      if ( chin(i:i)/=' ' ) chlen = i       ! si puo' sostuire con len_trim?
    end do
99001 format (a)                   ! Cosi' vine brutto chymlab
99002 format (a)                   ! Cosi' vine brutto chymlab
  end subroutine getcommandline
 
  subroutine getargsfromcommand(comm,args,n)
    implicit none
    character(len=*) :: comm
    integer :: n
    character(len=*) , dimension(:) :: args
    intent (in) comm
    intent (out) args
    intent (inout) n
    integer :: i , is1 , llco
    character(80) :: lco
    n = 0
    lco = comm
    call no2space(lco)
    call noinspace(lco)
    llco = len_trim(lco)
    if ( llco==0 ) return
    is1 = 1
    do i = 2 , llco
      if ( ichar(lco(i:i))==32 ) then
        n = n + 1
        args(n) = lco(is1:i-1)
        is1 = i + 1
      end if
    end do
    if ( is1<=llco ) then
      n = n + 1
      args(n) = lco(is1:llco)
    end if
  end subroutine getargsfromcommand
 
  subroutine readiarg(i,j)
    implicit none
    integer :: i , j
    intent (in) i
    intent (out) j
    character(60) :: str
    call getarg(i,str)
    j = -9999
    if ( len_trim(str)==0 ) then
      return
    else
      read (str,'(i60)',err=99999) j
    end if
99999 continue
  end subroutine readiarg
 
  subroutine readrarg(i,r)
    implicit none
    integer :: i
    real :: r
    intent (in) i
    intent (out) r
    character(60) :: str
    call getarg(i,str)
    r = -9999.0
    if ( len_trim(str)==0 ) then
      return
    else
      read (str,*,err=99999) r
    end if
99999 continue
  end subroutine readrarg
 
  subroutine readsarg(i,str)
    implicit none
    integer :: i
    character(len=*) :: str
    intent (in) i
    intent (out) str
    call getarg(i,str)
  end subroutine readsarg
 
  subroutine videocontrol(row,col,str)
    use mod_libmv
    implicit none
    integer :: col , row
    character(len=*) :: str
    intent (in) col , row
    intent (out) str
    character(1) :: esc
    if ( iflg(58)==0 ) then
      str = ' '
      return
    end if
    esc = char(27)
    if ( row==1 .and. col==0 ) then              ! reverse
      str = esc//'[7m'
    else if ( row==2 .and. col==0 ) then         ! clear screen
      str = esc//'[2J'
    else if ( col>0 ) then                       ! position
      write (str,'(a,i2,a,i2,a)') esc//'[' , row , ';' , col , 'H'
      call nospace(str)
    else                                         ! normal
      str = esc//'[0m'
    end if
  end subroutine videocontrol
 
  subroutine takefields(si,c,field,n)
    implicit none
    character(len=*) :: c , si
    integer :: n
    character(len=*) , dimension(:) :: field
    intent (inout) n
    integer :: iend
    character(6000) :: s
    if ( len_trim(si)>6000 ) then
      write (6,'(20x,a)') 'MVLib sever error: takefields cannot manipulate'
      write (6,'(20x,a)') 'strings longer than 6000 characters. Exiting.'
      call exit(0)
    end if
    s = si
    call tab2space(s)
    call noinspace(s)
    call no2space(s)
    n = 0
    if ( s(len_trim(s)-len_trim(c)+1:len_trim(s))/=c(1:len_trim(c)) ) then
      s = s(1:len_trim(s))//c
    end if
    do while ( len_trim(s)>0 )
      iend = index(s,c)
      n = n + 1
      field(n) = s(1:iend-1)
      s = s(iend+1:)
    end do
    if ( len_trim(field(n))==0 ) n = n - 1
  end subroutine takefields
 
  subroutine extractval(str,v)
    implicit none
    character(len=*) :: str
    real :: v
    intent (in) str
    intent (out) v
    integer :: i , lend , ll
    ll = 0
    lend = 0
    v = -9999.0
    do i = 1 , len_trim(str)
      if ( (ichar(str(i:i))>=48 .and. ichar(str(i:i))<=57) .or. &
            ichar(str(i:i))==46 .or. ichar(str(i:i))==45 ) then
        if ( lend==0 ) ll = i
      else if ( ll>0 ) then
        lend = 1
      end if
    end do
    if ( ll==0 ) return
    read (str(1:ll),*,err=100) v
    return
 100  v = -9999.0
  end subroutine extractval
 
  subroutine extractint(str,v)
    implicit none
    character(len=*) :: str
    integer :: v
    intent (in) str
    intent (out) v
    integer :: i , lend , ll
    ll = 0
    lend = 0
    v = -9999
    do i = 1 , len_trim(str)
      if ( (ichar(str(i:i))>=48 .and. ichar(str(i:i))<=57) .or.              &
            ichar(str(i:i))==45 ) then
        if ( lend==0 ) ll = i
      else if ( ll>0 ) then
        lend = 1
      end if
    end do
    if ( ll==0 ) return
      read (str(1:ll),*,err=100) v
      return
 100  v = -9999
  end subroutine extractint
 
  subroutine htmlspecialchar(str)
    implicit none
    character(len=*) :: str
    intent (inout) str
    integer :: i , n
    character(1024) :: lstr
    lstr = str
    htmlspecialchar_loop: &
    do
      n = len_trim(lstr)
      do i = 1 , n - 1
        if ( ichar(lstr(i:i))==195 .and. ichar(lstr(i+1:i+1))==164 ) then
          lstr = lstr(1:i-1)//'&auml;'//lstr(i+2:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==195 .and. ichar(lstr(i+1:i+1))==182 ) then
          lstr = lstr(1:i-1)//'&ouml;'//lstr(i+2:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==195 .and. ichar(lstr(i+1:i+1))==156 ) then
          lstr = lstr(1:i-1)//'&Uuml;'//lstr(i+2:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==195 .and. ichar(lstr(i+1:i+1))==188 ) then
          lstr = lstr(1:i-1)//'&uuml;'//lstr(i+2:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==223 ) then
          lstr = lstr(1:i-1)//'&szlig;'//lstr(i+1:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==195 .and. ichar(lstr(i+1:i+1))==159 ) then
          lstr = lstr(1:i-1)//'&szlig;'//lstr(i+2:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==194 .and. ichar(lstr(i+1:i+1))==178 ) then
          lstr = lstr(1:i-1)//'&sup2;'//lstr(i+2:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==228 ) then
          lstr = lstr(1:i-1)//'&auml;'//lstr(i+1:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==246 ) then
          lstr = lstr(1:i-1)//'&ouml;'//lstr(i+1:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==252 ) then
          lstr = lstr(1:i-1)//'&uuml;'//lstr(i+1:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==191 ) then
          lstr = lstr(1:i-1)//'&euro;'//lstr(i+1:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==194 .and. ichar(lstr(i+1:i+1))==191 ) then
          lstr = lstr(1:i-1)//'"'//lstr(i+2:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==195 .and. ichar(lstr(i+1:i+1))==160 ) then
          lstr = lstr(1:i-1)//'&agrave;'//lstr(i+2:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==195 .and. ichar(lstr(i+1:i+1))==162 ) then
          lstr = lstr(1:i-1)//'&acirc;'//lstr(i+2:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==195 .and. ichar(lstr(i+1:i+1))==168 ) then
          lstr = lstr(1:i-1)//'&egrave;'//lstr(i+2:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==195 .and. ichar(lstr(i+1:i+1))==169 ) then
          lstr = lstr(1:i-1)//'&eacute;'//lstr(i+2:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==195 .and. ichar(lstr(i+1:i+1))==170 ) then
          lstr = lstr(1:i-1)//'&ecirc;'//lstr(i+2:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==195 .and. ichar(lstr(i+1:i+1))==185 ) then
          lstr = lstr(1:i-1)//'&ugrave;'//lstr(i+2:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==195 .and. ichar(lstr(i+1:i+1))==178 ) then
          lstr = lstr(1:i-1)//'&ograve;'//lstr(i+2:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==197 .and. ichar(lstr(i+1:i+1))==147 ) then
          lstr = lstr(1:i-1)//'&oelig;'//lstr(i+2:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==195 .and. ichar(lstr(i+1:i+1))==180 ) then
          lstr = lstr(1:i-1)//'&ocirc;'//lstr(i+2:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==194 .and. ichar(lstr(i+1:i+1))==171 ) then
          lstr = lstr(1:i-1)//'"'//lstr(i+2:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==194 .and. ichar(lstr(i+1:i+1))==187 ) then
          lstr = lstr(1:i-1)//'"'//lstr(i+2:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==226 .and. ichar(lstr(i+1:i+1))==128 .and. &
             ichar(lstr(i+2:i+2))==157 ) then
          lstr = lstr(1:i-1)//'"'//lstr(i+3:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==226 .and. ichar(lstr(i+1:i+1))==128 .and. &
             ichar(lstr(i+2:i+2))==153 ) then
          lstr = lstr(1:i-1)//''''//lstr(i+3:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==226 .and. ichar(lstr(i+1:i+1))==130 .and. &
             ichar(lstr(i+2:i+2))==172 ) then
          lstr = lstr(1:i-1)//'&euro;'//lstr(i+3:)
          cycle htmlspecialchar_loop
        end if
        if ( ichar(lstr(i:i))==226 .and. ichar(lstr(i+1:i+1))==128 .and. &
            (ichar(lstr(i+2:i+2))==158 .or. ichar(lstr(i+2:i+2))==156) ) then
          lstr = lstr(1:i-1)//'"'//lstr(i+3:)
          cycle htmlspecialchar_loop
        end if
      end do
      exit
    end do htmlspecialchar_loop
    str = lstr
  end subroutine htmlspecialchar
 
  subroutine nameadjust(str)
    implicit none
    character(len=*) :: str
    intent (inout) str
    integer :: i , ln
    character(1024) :: lstr
    call cv2lower(str,lstr)
    call noinspace(lstr)
    call no2space(lstr)
    ln = len_trim(lstr)
    if ( ln<1 ) return
    do i = 1 , ln
      if ( lstr(i:i)=='_' ) then
        str = lstr(1:i-1)//' '//lstr(i+1:)
        lstr = str
      end if
    end do
    if ( ichar(lstr(1:1))>=97 .and. ichar(lstr(1:1))<=122 ) then
      str = char(ichar(lstr(1:1))-32)//lstr(2:)
      lstr = str
    end if
    do i = 2 , ln - 1
      if ( lstr(i:i)==' ' .and. ichar(lstr(i+1:i+1))>=97 .and.  &
                                ichar(lstr(i+1:i+1))<=122 ) then
        str = lstr(1:i)//char(ichar(lstr(i+1:i+1))-32)//lstr(i+2:)
        lstr = str
      end if
      if ( ichar(lstr(i:i))==39 .and. ichar(lstr(i+1:i+1))>=97 .and. &
           ichar(lstr(i+1:i+1))<=122 ) then
        str = lstr(1:i)//char(ichar(lstr(i+1:i+1))-32)//lstr(i+2:)
        lstr = str
      end if
    end do
    str = lstr
  end subroutine nameadjust

end module mod_strman
