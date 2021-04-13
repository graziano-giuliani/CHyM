subroutine basinpaint(thres,rlon,rlat,dem,fmap,luse,nlat,nlon,wk1)
  implicit none
  integer nlat,nlon,latriv,lonriv,pdrain,thres
  integer, dimension(nlat,nlon) :: fmap, luse
  real, dimension(nlat,nlon) :: dem, work
  real, intent(out),dimension(nlat,nlon) ::  wk1
  integer, dimension(9) :: ir, jr
  save ir,jr
  integer i,j,indx,ii,jj,idir,iflg,rlon,rlat
  data jr /-1, 0, 1, 1, 1, 0,-1,-1,0/
  data ir / 1, 1, 1, 0,-1,-1,-1, 0,0/
  integer, parameter :: mare=15, lago=14
  logical :: notyet
  latriv = rlat+1
  lonriv = rlon+1
  iflg = 11
  pdrain = 0
  if (iflg.ge.10)                                                            & 
      write (6,'(12x,a,i3,a,i3)') 'Basin Selecion Module called for grid '//      &
      'point ',latriv,'-',lonriv
  wk1=-10 ; work=0.0
  do i=2,nlat-1
     do j=2,nlon-1
        ii=i
        jj=j
        indx=0
        notyet=.true.
        do while (indx.le.2*(nlon+nlat).and.notyet)
           indx=indx+1
           idir=fmap(ii,jj)
           if (idir.gt.0) then
              ii=ii+ir(idir)
              jj=jj+jr(idir)
              if (ii.eq.latriv.and.jj.eq.lonriv.and.dem(i,j)>thres) then
                 notyet=.false.
                 if (thres>0) then
                   wk1(i,j)=thres
                 else
                   wk1(i,j)=dem(i,j)
                 end if
                 pdrain = pdrain +1
              else if (ii.eq.1.or.ii.eq.nlat.or.jj.eq.1.or.jj.eq.nlon)then
                 notyet=.false.                           ! Fuori dominio
              endif
           else if (luse(ii,jj).eq.mare) then             ! In mare
              notyet=.false.
           else                                           ! noflow
              wk1(i,j)=-5.0
              work(ii,jj)=work(ii,jj)+1
              notyet=.false.
           endif
        enddo
     enddo
  enddo
  write(6,'(/,3x,a,i7,/)') 'Number of drained points by the river = ',pdrain
  if (iflg.ge.10) then
     do i=2,nlat-1 ; do j=2,nlon-1
        if (work(i,j).gt.10.) write(6,'(/,3x,a,i3,a,i3,a,i4,a,/)')                &
           'Severe warning: The no flow point ',                                  &
           i,'-',j,' drains ',nint(work(i,j)),' cells.'
     enddo ; enddo
  endif
  return
end
