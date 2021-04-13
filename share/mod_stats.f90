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
module mod_stats

  use mod_vector
  use mod_libmv

  integer , private :: nbin , ndata
  integer , private , dimension(101) :: nbxch
  real , private :: xmax , xmin
  integer , private :: iy , jsum , jy , loop

  contains

  real function average(x,n)
    implicit none
    integer :: n
    real , dimension(n) :: x
    intent (in) n , x
    integer :: i
    average = 0.0
    do i = 1 , n
      average = average + x(i)/n
    end do
  end function average
 
  real function sigma2(x,n)
    implicit none
    integer :: n
    real , dimension(n) :: x
    intent (in) n , x
    integer :: i
    real :: xm
    sigma2 = 0
    xm = average(x,n)
    do i = 1 , n
      sigma2 = sigma2 + (x(i)-xm)*(x(i)-xm)/n
    end do
  end function sigma2
 
  real function wcorr(x,y,n)
    implicit none
    integer :: n
    real , dimension(n) :: x , y
    intent (in) n , x , y
    integer :: i
    real :: s1 , s2 , xm , ym
    wcorr = 0.0
    if ( n<=0 ) return
    xm = average(x,n)
    ym = average(y,n)
    do i = 1 , n
      wcorr = wcorr + (x(i)-xm)*(y(i)-ym)/n
    end do
    s1 = sqrt(sigma2(x,n))
    if ( s1<1E-30 ) then
      ! write (6,'(7x,a)')'MVlib wcorr: ''ndondolo sigma del vettore 1 e'' 0!'
      wcorr = 0.0
      return
    end if
    s2 = sqrt(sigma2(y,n))
    if ( s2<1E-30 ) then
      ! write (6,'(7x,a)')'MVlib wcorr: ''ndondolo sigma del vettore 2 e'' 0!'
      wcorr = 0.0
      return
    end if
    wcorr = wcorr/(s1*s2)
    ! wcorr=wcorr/(sqrt(sigma2(x,n))*sqrt(sigma2(y,n)))
  end function wcorr
 
  real function vcorr(x,y,n)
    implicit none
    integer :: n
    real , dimension(n) :: x , y
    intent (in) n , x , y
    vcorr = abs(wcorr(x,y,n))
  end function vcorr
 
  real function qanormal(id,x)
    implicit none
    integer :: id
    real :: x
    intent (in) x , id
    integer :: i
    real :: q1 , q2 , xdiv
    if ( x<=xmin ) then
      qanormal = xmin
    else if ( x>=xmax ) then
      qanormal = xmax
    end if
    xdiv = (xmax-xmin)/nbin
    do i = 1 , 100
      q1 = float(nbxch(i))/ndata
      q2 = float(nbxch(i+1))/ndata
      if ( x>=q1 .and. x<=q2 ) then
        qanormal = xmin + xdiv*(i-1) + xdiv*(x-q1)/(q2-q1)
        return
      end if
    end do
    stop ' Flux error inside qanormal'
  end function qanormal
 
  real function qnormal(id,x)
    implicit none
    integer :: id
    real :: x
    intent (in) x , id
    integer :: ig
    real :: xdiv , xleft , xright
    if ( x<=xmin ) then
      qnormal = 0.0
      return
    else if ( x>=xmax ) then
      qnormal = 1.0
      return
    end if
    xdiv = (xmax-xmin)/nbin
    ig = int((x-xmin)/xdiv) + 1
    if ( ig==nbin+1 ) ig = nbin
    if ( ig==0 ) ig = 1
    if ( ig<1 .or. ig>nbin ) stop ' Flux error inside qnormal'
    xleft = xmin + xdiv*(ig-1)
    xright = xmin + xdiv*ig
    qnormal = (nbxch(ig+1)-nbxch(ig))*(x-xleft)/(xright-xleft) + nbxch(ig)
    qnormal = qnormal/ndata
  end function qnormal
 
  subroutine qnormdef(id,x,n)
    implicit none
    integer :: id , n
    real , dimension(n) :: x
    intent (in) id , n
    integer :: i , ig , nnbin
    logical :: searchagain
    real :: xdiv
    real , dimension(101) :: y
    if ( iflg(35)>100 ) iflg(35) = 100
    if ( iflg(35)>0 ) then
      nbin = iflg(35)
    else
      nbin = n/20
    end if
    if ( nbin>100 ) nbin = 100
    xmax = maxval(x(1:n))
    xmin = minval(x(1:n))
    searchagain = .true.
    do while ( nbin>0 .and. searchagain )
      y(1:nbin) = 0.0
      xdiv = (xmax-xmin)/nbin
      do i = 1 , n
        ig = int((x(i)-xmin)/xdiv) + 1
        if ( ig==nbin+1 ) ig = nbin
        if ( ig==0 ) ig = 1
        if ( ig<1 .or. ig>nbin ) stop ' Flux error inside qnormdef'
        y(ig) = y(ig) + 1
      end do
      searchagain = .false.
      if ( iflg(35)<=0 ) then
        nnbin = nbin
        do i = 1 , nbin
          if ( nint(y(i))==0 ) then
            searchagain = .true.
            nnbin = nnbin - 3
          end if
        end do
        nbin = nnbin
      end if
    end do
    if ( nbin<10 ) write (6,'(7x,a,i1,a)') 'MVLib qnormdef: This is '//     &
                                           'a poor distribution - only ' ,  &
                                           nbin , ' defined.'
    nbxch(1) = 0
    do i = 2 , nbin + 1
      nbxch(i) = y(i-1) + nbxch(i-1)
    end do
    ndata = n
    if ( iflg(1)>=100 ) then
      write (6,'(7x,a,i2)') 'MVlib QNormDef: defined normalization ' , id
      write (6,'(8(i3,f7.1))') (i,y(i),i=1,nbin)
    end if
  end subroutine qnormdef
 
  subroutine rainscore(mm5,arssa,ns,alfa,resu)
    implicit none
    real :: alfa
    integer :: ns
    real , dimension(ns) :: arssa , mm5
    real , dimension(3) :: resu
    intent (in) alfa , arssa , mm5 , ns
    intent (out) resu
    real :: bias , far , thsco
    real , intrinsic :: float
    integer :: ia , ib , ic , id , inostat1 , inostat2 , inostat3 , n
    ia = 0
    ib = 0
    ic = 0
    id = 0
    do n = 1 , ns
      if ( mm5(n)>=alfa .and. arssa(n)>=alfa ) ia = ia + 1
      if ( mm5(n)>=alfa .and. arssa(n)<alfa ) ib = ib + 1
      if ( mm5(n)<alfa .and. arssa(n)>=alfa ) ic = ic + 1
      if ( mm5(n)<alfa .and. arssa(n)<alfa ) id = id + 1
    end do
    inostat1 = ia + ib + ic
    inostat2 = ia + ib
    inostat3 = ia + ic
    if ( inostat1==0 ) then
      thsco = 99.99
    else
      thsco = (float(ia))/(float(inostat1))
    end if
    if ( inostat2==0 ) then
      far = 99.99
    else
      far = (float(ib))/(float(inostat2))
    end if
    if ( inostat3==0 ) then
      bias = 99.99
    else
      bias = float(ia+ib)/float(inostat3)
    end if
 
    ! calcola l'indice di rousseau
    ! rouss=float(4*ia*id-(ib+ic))/float((2*ia+ib+ic)*(2*id+ib+ic))
 
    resu(1) = thsco
    resu(2) = far
    resu(3) = bias
  end subroutine rainscore
 
  real function gauss(aver,sigma)
    implicit none
    real , parameter :: cons = 1.1920928955078E-07
    integer , parameter :: mask31 = 2147483647
    real :: aver , sigma
    intent (in) aver , sigma
    logical , save :: first
    data first/.true./
    if ( first ) then
      if ( iflg(34)==0 ) iflg(34) = 875949887
      iy = iflg(34)
      first = .false.
    end if
    jsum = 0
    do loop = 1 , 12
      iy = iy*69069
      iy = iand(iy,mask31)
      jy = iy/256
      jsum = jsum + jy
    end do
    jsum = (jsum+134)/256*256
    gauss = cons*jsum - 6.0
    gauss = aver + sigma*gauss
  end function gauss

  real function setgauss(idummy)
    implicit none
    integer :: idummy
    intent (in) idummy
    iy = idummy
    setgauss = 0.0
  end function setgauss
 
  real function acaso(x)
    implicit none
    real :: x
    intent (in) x
    real , parameter :: cons = 4.6566128730774E-10
    integer , parameter :: mask31 = 2147483647
    logical , save :: first
    data first/.true./
    if ( first ) then
      if ( iflg(34)==0 ) iflg(34) = 875949887
      iy = iflg(34)
      first = .false.
    end if
    iy = iy*69069
    iy = iand(iy,mask31)
    jy = iy/256*256
    acaso = cons*jy
  end function acaso

  real function setacaso(idummy)
    implicit none
    integer :: idummy
    intent(in) :: idummy
    iy = idummy
    setacaso = 0.0
  end function setacaso
 
  subroutine mvft(dati,n,a,b,ns)
    implicit none
    integer :: n , ns
    real , dimension(0:ns-1) :: a , b
    real , dimension(n) :: dati
    intent (in) dati , n , ns
    intent (inout) a , b
    real :: dx , x
    integer :: i , k
    a(0) = 0
    b(0) = 0
    do i = 1 , n
      a(0) = a(0) + dati(i)
    end do
    a(0) = a(0)/n
    dx = 4*atan(1.0)/n
    do k = 1 , ns - 1
      a(k) = 0
      b(k) = 0
      do i = 1 , n
        x = dx*i
        a(k) = a(k) + dati(i)*sin(k*x)
        b(k) = b(k) + dati(i)*cos(k*x)
      end do
      a(k) = a(k)/n
      b(k) = b(k)/n
    end do
  end subroutine mvft
 
  subroutine mvftinv(dati,n,a,b,ns)
    implicit none
    integer :: n , ns
    real , dimension(0:ns-1) :: a , b
    real , dimension(n) :: dati
    intent (in) a , b , n , ns
    intent (inout) dati
    real :: dx , x
    integer :: i , k
    dx = 4*atan(1.0)/n
    do i = 1 , n
      x = dx*i
      dati(i) = a(0)/2
      do k = 1 , ns - 1
        dati(i) = dati(i) + a(k)*sin(k*x) + b(k)*cos(k*x)
      end do
    end do
  end subroutine mvftinv
 
  subroutine mvftfilt(inp,ou,n,spett1,spett2)
    implicit none
    integer :: n , spett1 , spett2
    real , dimension(n) :: inp , ou
    intent (in) spett1
    real , dimension(0:1000) :: a , b
    integer :: i
    call mvft(inp,n,a,b,spett2)
    do i = 1 , spett1 - 1
      a(i) = 0.0
      b(i) = 0.0
    end do
    call mvftinv(ou,n,a,b,spett2)
  end subroutine mvftfilt
 
  subroutine mvfthffilter(y1,y2,ncomp,a,b,scent,swide)
    implicit none
    integer :: ncomp
    real :: scent , swide
    real , dimension(0:ncomp-1) :: a , b
    real , dimension(ncomp) :: y1 , y2
    intent (inout) a , b
    integer :: i
    call mvft(y1,ncomp,a,b,ncomp)
    do i = 0 , ncomp - 1
      a(i) = a(i)*(1-sigmoide(float(i),scent,swide))
      b(i) = b(i)*(1-sigmoide(float(i),scent,swide))
    end do
    call mvftinv(y2,ncomp,a,b,ncomp)
  end subroutine mvfthffilter
 
  subroutine mvftlffilter(y1,y2,ncomp,a,b,scent,swide)
    implicit none
    integer :: ncomp
    real :: scent , swide
    real , dimension(0:ncomp-1) :: a , b
    real , dimension(ncomp) :: y1 , y2
    intent (inout) a , b
    integer :: i
    call mvft(y1,ncomp,a,b,ncomp)
    do i = 0 , ncomp - 1
      a(i) = a(i)*sigmoide(float(i),scent,swide)
      b(i) = b(i)*sigmoide(float(i),scent,swide)
    end do
    call mvftinv(y2,ncomp,a,b,ncomp)
  end subroutine mvftlffilter
 
  real function sigmoide(x,a,b) ! a = centro della sigmoide
    implicit none
    real :: a , b , x
    intent (in) a , b , x
    sigmoide = 1.0/(1.0+exp(-(x-a)/(b/2.19723))) ! b = distanza da a in cui la
  end function sigmoide
 
  subroutine naverage(x,y,n,nav)
    implicit none
    integer :: n , nav
    real , dimension(n) :: x , y
    intent (in) n , nav , x
    intent (inout) y
    integer :: i , i1 , i2 , j , j1 , j2
    if ( mod(nav,2)==0 ) then
      i1 = nav/2 + 1
      i2 = n - nav/2 + 1
      j1 = -nav/2
      j2 = nav/2 - 1
    else
      i1 = nav/2 + 1
      i2 = n - nav/2
      j1 = -nav/2
      j2 = nav/2
    end if
    do i = 1 , i1 - 1
      y(i) = 0.0
      do j = 1 , i1 + i - 1
        y(i) = y(i) + x(j)
      end do
      y(i) = y(i)/(i1+i-1)
    end do
    do i = i1 , i2
      y(i) = 0.0
      do j = j1 , j2
        y(i) = y(i) + x(i+j)
      end do
      y(i) = y(i)/nav
    end do
    do i = i2 + 1 , n
      y(i) = 0.0
      do j = i - nav/2 , n
        y(i) = y(i) + x(j)
      end do
      y(i) = y(i)/(1+n-i+nav/2)
    end do
  end subroutine naverage
 
  real function bias(x1,x2,n)
    implicit none
    integer :: n
    real , dimension(n) :: x1 , x2
    intent (in) n , x1 , x2
    integer :: i
    bias = 0.0
    do i = 1 , n
      bias = bias + (x1(i)-x2(i))
    end do
    if ( n>0 ) bias = -bias/n
  end function bias
 
  real function rmst(x1,x2,n)
    implicit none
    integer :: n
    real , dimension(n) :: x1 , x2
    intent (in) n , x1 , x2
    integer :: i
    rmst = 0.0
    do i = 1 , n
      rmst = rmst + (x1(i)-x2(i))*(x1(i)-x2(i))
    end do
    if ( n>0 ) rmst = sqrt(rmst/n)
  end function rmst
 
  real function rmstn(x1,x2,n)
    implicit none
    integer :: n
    real , dimension(n) :: x1 , x2
    intent (in) n , x1 , x2
    rmstn = rmst(x1,x2,n)/sqrt(sigma2(x1,n))
  end function rmstn
 
  real function cioppskill(x,y,n)
    implicit none
    integer :: n
    real , dimension(n) :: x , y
    intent (in) n , x , y
    integer :: i
    real :: xden , xnum , ymean
    xnum = 0.0
    xden = 0.0
    ymean = average(y,n)
    do i = 1 , n
      xnum = xnum + ((x(i)-y(i))**2)
      xden = xden + ((ymean-y(i))**2)
    end do
    cioppskill = 1.0 - (xnum/xden)
  end function cioppskill
 
  subroutine autocorrelation(x,n,y,maxlags)
    implicit none
    integer :: maxlags , n
    real , dimension(n) :: x
    real , dimension(maxlags) :: y
    intent (in) x , maxlags , n
    intent (out) y
    integer :: i , lag
    do lag = 1 , maxlags
      do i = 1 , n - lag
        x1dum(i) = x(i+lag)
      end do
      y(lag) = wcorr(x,x1dum,n-lag)
    end do
  end subroutine autocorrelation
 
  subroutine o1derivative(x,y,n)
    implicit none
    integer :: n
    real , dimension(n) :: x , y
    intent (in) n , x
    intent (out) y
    integer :: i
    if ( n>=1 ) y(1) = x(1)
    do i = 2 , n
      y(i) = x(i) - x(i-1)
    end do
  end subroutine o1derivative
 
  subroutine medianfilter(x,y,n,m)
    implicit none
    integer :: m , n
    real , dimension(n) :: x , y
    intent (in) n , x
    intent (out) y
    integer :: i , j , n1 , n2
    real , dimension(100) :: w
    n2 = m/2
    if ( mod(m,2)/=0 ) then
      n1 = m/2
    else
      n1 = m/2 - 1
    end if
    do i = n1 + 1 , n - n2
      do j = -n1 , n2
        w(j+1+n1) = x(i+j)
      end do
      call vec_ordinad(w,1,m)
      y(i) = w(n2)
    end do
    do i = 1 , n1
      y(i) = x(i)
    end do
    do i = n - n2 + 1 , n
      y(i) = x(i)
    end do
  end subroutine medianfilter

  subroutine linearreg(x,y,n,a,b)
    implicit none
    integer , parameter :: nmaxstep = 10000
    real :: a , b
    integer :: n
    real , dimension(n) :: x , y
    intent (in) :: n , x , y
    intent (inout) a , b
    real :: d , err , sterr
    integer :: ii , istep
    sterr = erronlinearreg(x,y,n,a,b)
    err = sterr
    do ii = 1 , 20
      d = 10000.0/(ii*ii*nmaxstep)
      d = 1.0/(ii*ii)
      istep = 0
      do while ( istep<nmaxstep )
        istep = istep + 1
        err = erronlinearreg(x,y,n,a,b)
        if ( err>erronlinearreg(x,y,n,a+d,b) ) then
          a = a + d
        else
          a = a - d
        end if
        err = erronlinearreg(x,y,n,a,b)
        if ( err>erronlinearreg(x,y,n,a,b+d) ) then
          b = b + d
        else
          b = b - d
        end if
      end do
    end do
  end subroutine linearreg
 
  real function erronlinearreg(x,y,n,a,b)
    implicit none
    real :: a , b
    integer :: n
    real , dimension(n) :: x , y
    intent (in) a , b , n , x , y
    integer :: i
    erronlinearreg = 0.0
    do i = 1 , n
      erronlinearreg = erronlinearreg + (y(i)-(a+b*x(i)))**2
    end do
  end function erronlinearreg

end module mod_stats
