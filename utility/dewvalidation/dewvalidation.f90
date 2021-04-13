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
    program dewvalidation
        use mod_museo
        use mod_time
        use mod_ncio
        use mod_mpimess
        use mod_varandtypes
        use mod_dewfun
	implicit none
	integer j,n1,annoi,latlonrange,ndata,nrec,lunl,n,isource,irec,last,ii
!	parameter (nrec=3338)
	character com*64,sens*32,title*128,file*60,cfile*60
        real , dimension(700000) :: late , lone , pe
	real ylat,ylon,x1(366*24*4),y1(366*24*4)
	integer plot,hist
        integer, parameter :: mxnicorw=700000
        real, dimension(mxnicorw) :: la,lo,pi
	logical log
        save lunl
!	anno 2001 - Numri al lotto per i record 0041 0126 0459 0548 0760 0966 
!                   1145 1565 1726 1751 2106 2009
	annoi=2001
	plot=1 ; log = .true.
	nrec=mvlibmagicnum(1)
	call mvsetflags('X size',900.0)
	call tsvalidation('dewrain',annoi,0.05,2.5)
        call exit(0)
        print*,"unoooo"
        do irec=1,366*24
           last = 0
	   y1=-9999
	  call museoanagrafica('dewrain',plot,com,sens,ylat,ylon)
          print*,"sensore",sens,com,ylat,ylon,plot
          plot=2
	  call museoanagrafica('dewrain',plot,com,sens,ylat,ylon)
          print*,"sensore",sens,com,ylat,ylon,plot
           call exit(0)
	  ! call museotimeseries('dewrain',plot,annoi,x1,y1,n1)
	  ! call convert2hourlyres(y1,1) ; ndata=8760
           write (cfile,'(a,i4,a)') 'pioggiaoraria' , annoi , '.dat'
           call openmuseofiles(lunl,cfile,0)
           if ( irec>=1 .and. irec<=366*24 ) then
             if ( irec<=last ) then
               rewind (lunl)
               last = 0
             end if
             do ii = last + 1 , irec - 1
               read (lunl,err=100,end=100) n
             end do
             n = 0
             read (lunl,err=100,end=100) n , (late(i),lone(i),isource,pe(i),i=1,n)
             print*,"n",n
             last = irec
           end if
           if (irec==50) then 
           print*,"precipitazione",pe(1:n)
           print*,"lat",late(1:n)
           print*,"lon",lone(1:n)
           print*,"isource",isource
           call exit(0)
           end if
        end do
	   do i=1,ndata 
              x1(i)=i/24.
           enddo
	   call extractgooddata(x1,y1,ndata,-9999)
	   write(title,'(a,i4,a,i4,a)') ' Record ',plot,' - year ',annoi, &
		     ': '//com(1:len_trim(com))//' - '//sens(1:len_trim(sens))
!	   call plottaborder(0.10,0.90,0.52,0.90)
!	   call plottaframe(0.0,365.0,0.0,100.0) 
!	   if (ndata.gt.2) call plottavw(x1,y1,ndata,' ')
	call tsvalidation('dewrain',annoi,0.05,2.5)
!F2	if (plot.ge.1.and.plot.le.nrec) then
!F2	   write(file,'(a,i4,a)') 'dewrain',annoi,'.dat'
!F2	   call getlun(lun)
!F2	   open(lun,file=file,status='old',form='unformatted')
!F2!	   open(lun,file='dewrain2001.dat-ori',status='old',form='unformatted')
!F2	   do j=1,plot
!F2	      y1=-9999 ; read(lun) n1,(y1(i),i=1,n1)
!F2	   enddo
!F2	   close(lun)
!F2	   call convert2hourlyres(y1,1) ; ndata=8760
!F2	   do i=1,ndata 
!F2             x1(i)=i/24.
!F2           enddo
!F	   call extractgooddata(x1,y1,ndata,-9999)
!F	   call plottaborder(0.10,0.90,0.12,0.50)
!F	   call plottaframe(0.0,365.0,0.0,100.0) 
!F	   if (ndata.gt.2) call plottavw(x1,y1,ndata,' ')
!F	   call displayplot(' ',title,' ')
!F2	endif
  100  close (lun)
	end program dewvalidation

