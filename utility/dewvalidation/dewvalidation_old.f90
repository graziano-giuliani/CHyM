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
	integer j,n1,annoi,latlonrange,ndata,nrec
!	parameter (nrec=3338)
	character com*64,sens*32,title*128,file*60
	real ylat,ylon,x1(366*24*4),y1(366*24*4)
	integer plot,hist
	logical log
!	anno 2001 - Numri al lotto per i record 0041 0126 0459 0548 0760 0966 
!                   1145 1565 1726 1751 2106 2009
	annoi=2001
	plot=0 ; log = .true.
	nrec=mvlibmagicnum(1)
	call mvsetflags('X size',900.0)
	if (plot.ge.1.and.plot.le.nrec) then
	   y1=-9999
	   call museoanagrafica('dewrain',plot,com,sens,ylat,ylon)
	   call museotimeseries('dewrain',plot,annoi,x1,y1,n1)
	   call convert2hourlyres(y1,1) ; ndata=8760
	   do i=1,ndata 
              x1(i)=i/24.
           enddo
	   call extractgooddata(x1,y1,ndata,-9999)
	   write(title,'(a,i4,a,i4,a)') ' Record ',plot,' - year ',annoi, &
		     ': '//com(1:len_trim(com))//' - '//sens(1:len_trim(sens))
	   call plottaborder(0.10,0.90,0.52,0.90)
	   call plottaframe(0.0,365.0,0.0,100.0) 
	   if (ndata.gt.2) call plottavw(x1,y1,ndata,' ')
	endif
	call tsvalidation('dewrain',annoi,0.05,2.5)
	if (plot.ge.1.and.plot.le.nrec) then
	   write(file,'(a,i4,a)') 'dewrain',annoi,'.dat'
	   call getlun(lun)
	   open(lun,file=file,status='old',form='unformatted')
!	   open(lun,file='dewrain2001.dat-ori',status='old',form='unformatted')
	   do j=1,plot
	      y1=-9999 ; read(lun) n1,(y1(i),i=1,n1)
	   enddo
	   close(lun)
	   call convert2hourlyres(y1,1) ; ndata=8760
	   do i=1,ndata 
             x1(i)=i/24.
           enddo
!F	   call extractgooddata(x1,y1,ndata,-9999)
!F	   call plottaborder(0.10,0.90,0.12,0.50)
!F	   call plottaframe(0.0,365.0,0.0,100.0) 
!F	   if (ndata.gt.2) call plottavw(x1,y1,ndata,' ')
!F	   call displayplot(' ',title,' ')
	endif
	end program dewvalidation

