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

module mod_phys

  contains

  subroutine kelvin2centi(t,ix,jx,kx)
    implicit none
    integer :: ix , jx , kx
    real , dimension(ix,jx,kx) :: t
    intent (in) ix , jx , kx
    intent (inout) t
    integer :: i , j , k
    do i = 1 , ix
      do j = 1 , jx
        do k = 1 , kx
          t(i,j,k) = t(i,j,k) - 273.15
        end do
      end do
    end do
  end subroutine kelvin2centi

  real function centi2fahrenheit(xf)
    implicit none
    real :: xf
    intent (in) xf
    centi2fahrenheit = (5./9.)*(xf-32.)
  end function centi2fahrenheit

  real function fahrenheit2centi(xc)
    implicit none
    real :: xc
    intent (in) xc
    fahrenheit2centi = (9./5.)*xc + 32.
  end function fahrenheit2centi

  real function winddirection(ux,vx)
    implicit none
    real :: ux , vx
    intent (in) ux , vx
    if ( abs(ux)<1.0E-30 .and. vx>=0.0 ) then
      winddirection = 0.0
    else if ( abs(ux)<1.0E-30 .and. vx<=0.0 ) then
      winddirection = 180.0
    else if ( ux>=0.0 .and. vx>=0.0 ) then
      winddirection = rad2deg(atan(vx/ux))
    else if ( ux<=0.0 .and. vx>=0.0 ) then
      winddirection = rad2deg(atan(vx/ux)) + 180
    else if ( ux<=0.0 .and. vx<0.0 ) then
      winddirection = rad2deg(atan(vx/ux)) + 180
    else if ( ux>=0.0 .and. vx<0.0 ) then
      winddirection = rad2deg(atan(vx/ux)) + 360
    else
      write (6,'(10x,a)') 'MVLib: flux error inside winddirection function.'
      write (6,'(10x,a,2f20.10)') 'Input values were: ' , ux , vx
      winddirection = 0.0
      return
    end if
  end function winddirection

  real function sessa2centi(gradi,primi,secondi)
    implicit none
    integer :: gradi , primi , secondi
    intent (in) gradi , primi , secondi
    sessa2centi = gradi + (primi*60+secondi)/3600.
  end function sessa2centi

  real function distance(lat1,lon1,lat2,lon2)
    implicit none
    real rad,dpi ; parameter(rad=6371000.0,dpi=6.2831855)
    real lat1,lon1,lat2,lon2,lt1,lt2,ln1,ln2,x,y
    lt1=lat1*dpi/360. ; lt2=lat2*dpi/360.
    ln1=lon1*dpi/360. ; ln2=lon2*dpi/360.
    if (abs(lat1-lat2).lt.0.2.and.abs(lon1-lon2).lt.0.2) then
        x=(rad*cos(lt1)*(ln1-ln2))*(rad*cos(lt2)*(ln1-ln2))
        y=(rad*(lt1-lt2))**2
        distance=sqrt(x+y)
    else
       x=sin(lt1)*sin(lt2)+cos(lt1)*cos(lt2)*cos((ln1)-(ln2))
       if (x.gt.1) x=1.0
       distance=acos(x)*rad
    endif
    if (distance.lt.0.1) distance=0.1
    return
  end function distance

  real function olddistance(lat1,lon1,lat2,lon2)
    implicit none
    real :: lat1 , lat2 , lon1 , lon2
    intent (in) lat1 , lat2 , lon1 , lon2
    logical , save :: convert , ifirst
    real , save :: dpi , ln1 , ln2 , lt1 , lt2 , rad , teta , x , x1 , x2
    data convert/.true./                 ! if true input data deg else rad
    data ifirst/.true./
    data rad/6371000./
    if ( ifirst ) then
      ifirst = .false.
      dpi = 8.0*atan(1.0)
    end if
    if ( convert ) then
      lt1 = lat1*dpi/360.0
      lt2 = lat2*dpi/360.0
      ln1 = lon1*dpi/360.0
      ln2 = lon2*dpi/360.0
    else
      lt1 = lat1
      lt2 = lat2
      ln1 = lon1
      ln2 = lon2
    end if
    x = sin(lt1)*sin(lt2) + cos(lt1)*cos(lt2)*cos((ln1)-(ln2))
    if ( x>=1 ) then
      olddistance = 0.0
    else
      if ( x<-1.0 ) x = -1.0
      teta = acos(x)
      olddistance = teta*rad
    end if
    if ( olddistance<20000. ) then
      x1 = rad*sin(lt1-lt2)
      x2 = rad*cos(ln2)*sin(ln1-ln2)
      olddistance = sqrt(x1*x1+x2*x2)
    end if
    if ( olddistance<0.1 ) olddistance = 0.1
  end function olddistance

  real function deg2rad(alpha)
    implicit none
    real :: alpha
    intent (in) alpha
    deg2rad = 4.*atan(1.)*alpha/180.
  end function deg2rad

  real function rad2deg(alpha)
    implicit none
    real :: alpha
    intent (in) alpha
    rad2deg = 45.*alpha/atan(1.)
  end function rad2deg

end module mod_phys
