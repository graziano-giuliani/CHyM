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
module mod_mssg

  use mod_strman
  use mod_html
  
  contains
  subroutine chymerror(ic,i,r,s)
    implicit none
    integer :: i , ic
    real :: r
    character(len=*) :: s
    intent (in) ic , r
    character(len=24) :: chyms
    character(len=12) :: e
    character(len=132) :: stringona
    data chyms/'CHyM model severe error.'/
    data e/', exiting...'/
    if ( ic==1 ) then
      write (6,'(6x,a)') chyms
      write (6,'(6x,a)') 'Maximum secondary inside '//trim(s)//e
      call exit(1)
    else if ( ic==2 ) then
      write (6,'(6x,a)') chyms
      write (6,'(6x,a)') 'Unkown precipitation sources: '//trim(s)//e
      call exit(1)
    else if ( ic==3 ) then
      write (6,'(/,6x,a)') chyms
      write (6,'(6x,a)') 'Maximum vector lenght reached while processing '// &
                         trim(s)//' data.'
      write (6,'(6x,a,i9,a)') 'Inside buildrainfield parameter mxnicorw '//  &
                              'must be at least' , i , e
      call exit(1)
    else if ( ic==4 ) then
      write (6,'(/,6x,a)') chyms
      write (6,'(6x,a)') 'Maximum vector lenght reached while processing '// &
                         trim(s)//' data.'
      write (6,'(6x,a,i9,a)') 'parameter n1dv must be at least' , i , e
      call exit(1)
    else if ( ic==5 ) then
      write (6,'(/,6x,a)') chyms
      write (6,'(6x,a)') 'Error reading previous '//trim(s)//' fields.'
      call writetime(trim(s)//' time series start at',i)
      call writetime('          This run should start at',nint(r))
      write (6,'(6x,a)') e(3:)
      call exit(1)
    else if ( ic==6 ) then
      write (6,'(/,6x,a)') chyms
      write (6,'(6x,a)') 'Insufficient number of time slices in '//trim(s) &
                         //' file.'
      write (stringona,'(i5,a,i5,a)') i , ' needed, ' , nint(r) , ' available.'
      call noinspace(stringona)
      call no2space(stringona)
      write (6,'(6x,a)') trim(stringona)
      write (6,'(6x,a)') e(3:)
      call exit(1)
    else if ( ic==7 ) then
      write (6,'(/,6x,a)') chyms
      write (6,'(6x,a)') 'Cannot open MODIS data file: '//trim(s)
      write (6,'(6x,a)') 'Check the file or set MODIS = 0 in CHyM script'
      write (6,'(6x,a)') e(3:)
      call exit(1)
    else if ( ic==8 ) then
      write (6,'(/,6x,a)') chyms
      write (6,'(6x,a)') 'Flux error inside '//trim(s)//' routine.'
      write (6,'(6x,a)') e(3:)
      call exit(1)
    else if ( ic==9 ) then
      write (6,'(/,6x,a)') chyms
      write (6,'(6x,a)') trim(s)//' routine cannot open input file.'
      write (6,'(6x,a)') e(3:)
      call exit(1)
    else if ( ic==20 ) then
      write (6,'(6x,a)') chyms
      write (6,'(6x,a)') trim(s)//' rainfall source module not '// &
                         'yet implemented'//e
      call exit(1)
    else if ( ic==21 ) then
      write (6,'(6x,a)') chyms
      write (6,'(6x,a,i3,a)') 'Temperature source code for acqwa '//   &
                              trim(s)//' scenario not yet implemented.'
      write (6,'(6x,a)') e(3:)
      call exit(1)
    else if ( ic==22 ) then
      write (6,'(6x,a)') chyms
      write (6,'(6x,a,i3,a)') 'Temperature source code ' , i ,  &
                              ' not yet implemented'//e
      call exit(1)
    else if ( ic==23 ) then
      write (6,'(6x,a)') chyms
      write (6,'(6x,a,i3,a)') 'Integration type flag ' , i ,  &
                              ' not existent, please insert integrflag=0 for old' &
                              // ' and =1 for new in the namelist file'//e
      call exit(1)
    else
      write (stringona,'(6x,a,i10)') 'chymerror, unknown errot code: ' , ic
      call no2space(stringona)
      write (6,'(a)') stringona
    end if
  end subroutine chymerror

  subroutine writetime(str,idate)
    implicit none
    integer :: idate
    character(len=*) :: str
    character(len=8) :: sstr1
    call integer2string(idate,8,sstr1)
    write (6,'(12x,a)') trim(str)//' '//sstr1
  end subroutine writetime

  subroutine writestatus(string)
    implicit none
    character(len=*) :: string
    intent (in) string
    integer :: file , lun
    data file/0/
    if ( file<0 ) return
    call getlun(lun)
    file = -1
    open (lun,file = 'html/status.html',status = 'unknown',err = 100)
    file = 1
    call htmlheader('CHyM Run Status',' ',lun)
    write (lun,99002) trim(schym(4))
    write (lun,99003) trim(schym(6))
    write (lun,99001) trim(schym(5)), string
    close (lun)
   100  return
99001 format ('<BR>CHyM running on : ',a,a)
99002 format ('<BR>The simulation started ',a)
99003 format ('<BR>Title of experiment is ',a)
  end subroutine writestatus

end module mod_mssg
