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
program dranet
  use netcdf
  use mod_statparams
  use mod_param
  use mod_museo
  use mod_mpimess
  use mod_varandtypes
  implicit none

  integer :: nlon1,nlat1,nchymp1,ii,jj,icode,istatus
  real :: dij1,maxdra,maxsubdra,treshdr
  real,parameter :: thdiv=250
  integer :: jjmax,iimax,idir,ios
  integer :: jjsmax,iismax
  character(len=256) :: arg
  logical :: run = .true.
  call mpi_init(mpierr)
  call mpi_comm_dup(mpi_comm_world,mycomm,mpierr)
  call mpi_comm_size(mycomm, nproc, mpierr)
  call mpi_comm_rank(mycomm, myid, mpierr)

  maxdra  = 0
  iimax   = 0
  jjmax   = 0
  ii      = 1

  call setparam
  call acquirescriptpar

  call getlun(lun)
  istatus = nf90_open(trim(schym(12)),nf90_nowrite,lun)
  call chkncstatus(istatus,nf90_noerr,'open',trim(schym(12)),'')
  
  istatus = nf90_inq_dimid(lun,'lat',icode)
  call chkncstatus(istatus,nf90_noerr,'nlat code',trim(schym(12)),'')
  istatus = nf90_inquire_dimension(lun,icode,len=nlat1)
  call chkncstatus(istatus,nf90_noerr,'nlat dim',trim(schym(12)),'')
  istatus = nf90_inq_dimid(lun,'lon',icode)
  call chkncstatus(istatus,nf90_noerr,'nlon code',trim(schym(12)),'')
  istatus = nf90_inquire_dimension(lun,icode,len=nlon1)
  call chkncstatus(istatus,nf90_noerr,'nlon dim',trim(schym(12)),'')
  istatus = nf90_get_att(lun,nf90_global,'Start_longitude',slon)
  call chkncstatus(istatus,nf90_noerr,'Start_longitude',    &
    trim(schym(12)),'')
  istatus = nf90_get_att(lun,nf90_global,'Start_latitude',slat)
  call chkncstatus(istatus,nf90_noerr,'Start_latitude',trim(schym(12)),'')
  istatus = nf90_get_att(lun,nf90_global, &
        'approximate_lat-lon_resolution_in_meters',dij1)
  call chkncstatus(istatus,nf90_noerr,   &
    'approximate_lat-lon_resolution_in_meters',trim(schym(12)),'')
!  istatus = nf90_inq_varid(lun,'dem',icode)
!  call chkncstatus(istatus,nf90_noerr,'dem',trim(schym(12)),'')
!  istatus = nf90_get_var(lun,icode,dem)
!  call chkncstatus(istatus,nf90_noerr,'dem',trim(schym(12)),'')
  istatus = nf90_inq_varid(lun,'fdm',icode)
  call chkncstatus(istatus,nf90_noerr,'fdm',trim(schym(12)),'')
  istatus = nf90_get_var(lun,icode,fmap)
  call chkncstatus(istatus,nf90_noerr,'fdm',trim(schym(12)),'')
!  istatus = nf90_inq_varid(lun,'acc',icode)
!  call chkncstatus(istatus,nf90_noerr,'acc',trim(schym(12)),'')
!  istatus = nf90_get_var(lun,icode,accl)
!  call chkncstatus(istatus,nf90_noerr,'acc',trim(schym(12)),'')
  istatus = nf90_inq_varid(lun,'lus',icode)
  call chkncstatus(istatus,nf90_noerr,'lus',trim(schym(12)),'')
  istatus = nf90_get_var(lun,icode,luse)
  call chkncstatus(istatus,nf90_noerr,'lus',trim(schym(12)),'')
!  istatus = nf90_inq_varid(lun,'aer',icode)
!  call chkncstatus(istatus,nf90_noerr,'aer',trim(schym(12)),'')
!  istatus = nf90_get_var(lun,icode,area)
!  call chkncstatus(istatus,nf90_noerr,'aer',trim(schym(12)),'')
  istatus = nf90_inq_varid(lun,'dra',icode)
  call chkncstatus(istatus,nf90_noerr,'dra',trim(schym(12)),'')
  istatus = nf90_get_var(lun,icode,drai)
  call chkncstatus(istatus,nf90_noerr,'dra',trim(schym(12)),'')
!  istatus = nf90_inq_varid(lun,'run',icode)
!  call chkncstatus(istatus,nf90_noerr,'run',trim(schym(12)),'')
!  istatus = nf90_get_var(lun,icode,runt)
!  call chkncstatus(istatus,nf90_noerr,'run',trim(schym(12)),'')
  istatus = nf90_close(lun)
  call chkncstatus(istatus,nf90_noerr,'close',trim(schym(12)),'')
  call getlun(lun)
  OPEN(UNIT=lun, IOSTAT=ios, FILE=trim(adjustl(schym(11)))//'.drainage_net.csv', &
        STATUS='replace', ACTION='write')!,FORM='formatted',POSITION="append")
  write(lun,'(a)')"line_id,point_id,x,y,drainage_net"
  treshdr = MAXVAL(drai)/thdiv
  do while(run)
   call catchriver(nlat1,nlon1,ii,treshdr,run)
   ii = ii + 1
  end do
  close(lun)
    
end program dranet

subroutine catchriver(nlat1,nlon1,river,treshdr,run)
  use mod_statparams
  use mod_param
  use mod_museo
  implicit none
  integer :: nlon1,nlat1,ii,jj,river
  real :: maxdra,maxsubdra,treshdr
  integer :: jjmax,iimax,idir
  integer :: jjsmax,iismax
  character(100000) :: a,b,c,d,e
  logical :: run
  i = 0
  maxdra = 0
  do ii = 1,nlat1
     do jj = 1,nlon1
        if (drai(jj,ii) > maxdra) then
           maxdra = drai(jj,ii)
           iimax = ii
           jjmax = jj
        end if
     end do
  end do
  if (maxdra<treshdr) run = .false.
  jjsmax = jjmax
  iismax = iimax
  do while(drai(jjmax,iimax)>treshdr)
    i=i+1
    write(a,'(i5)')river
    write(b,'(i5)')i
    write(c,'(f8.3)')lon(jjmax,iimax)
    write(d,'(f8.3)')lat(jjmax,iimax)
    write(e,'(f11.3)')drai(jjmax,iimax)
    write(lun,'(A,A,A,A,A,A,A,A,A)')trim(adjustl(a)),",",trim(adjustl(b)),",",trim(adjustl(c)),",",&
        trim(adjustl(d)),",",trim(adjustl(e))
!    WRITE(lun,'(I3,A,I4,A,f8.3,a,f8.3)') river,",",i,",",lon(jjmax,iimax),",",lat(jjmax,iimax)
    drai(jjmax,iimax) = 0
    maxsubdra = 0
    jjsmax = jjmax
    iismax = iimax
    do idir = 1,8
       if (luse(jjsmax+jr(idir),iismax+ir(idir)).ne.mare.and. &
          fmap(jjsmax+jr(idir),iismax+ir(idir)).ne.0.and.  &
          iismax+ir(idir).gt.1.and.iismax+ir(idir).lt.nlat1.and. &
          jjsmax+jr(idir).gt.1.and.jjsmax+jr(idir).lt.nlon1) then
         if (drai(jjsmax+jr(idir),iismax+ir(idir)) > maxsubdra) then
            maxsubdra = drai(jjsmax+jr(idir),iismax+ir(idir))
            iimax = iismax+ir(idir)
            jjmax = jjsmax+jr(idir)
         end if
       end if
    end do
  end do
  return
end subroutine catchriver
