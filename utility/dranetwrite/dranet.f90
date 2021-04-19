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
!  real,parameter :: thdiv=500
  real,parameter :: thdiv=50
  integer :: jjmax,iimax,idir,ios
  integer :: jjsmax,iismax
  character(len=256) :: arg
  logical :: run = .true.
  integer, allocatable, dimension(:,:) :: dummy
  integer, dimension(2) :: dimensions
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
!  istatus = nf90_get_att(lun,nf90_global,'Start_longitude',slon)
!  call chkncstatus(istatus,nf90_noerr,'Start_longitude',    &
!    trim(schym(12)),'')
!  istatus = nf90_get_att(lun,nf90_global,'Start_latitude',slat)
!  call chkncstatus(istatus,nf90_noerr,'Start_latitude',trim(schym(12)),'')
!  istatus = nf90_get_att(lun,nf90_global, &
!        'approximate_lat-lon_resolution_in_meters',dij1)
!  call chkncstatus(istatus,nf90_noerr,   &
!    'approximate_lat-lon_resolution_in_meters',trim(schym(12)),'')
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
  istatus = nf90_inq_varid(lun,'Q100',icode)
  if (istatus .ne. 0) then
     istatus = nf90_inq_varid(lun,'QRP',icode)
    if (istatus .ne. 0) then
       istatus = nf90_inq_varid(lun,'qmax',icode)
      if (istatus .ne. 0) then
         istatus = nf90_inq_varid(lun,'dis',icode)
      end if
    end if
  end if
  call chkncstatus(istatus,nf90_noerr,'Q100',trim(schym(12)),'')
  istatus = nf90_get_var(lun,icode,area)
  call chkncstatus(istatus,nf90_noerr,'Q100',trim(schym(12)),'')
!  istatus = nf90_inq_varid(lun,'run',icode)
!  call chkncstatus(istatus,nf90_noerr,'run',trim(schym(12)),'')
!  istatus = nf90_get_var(lun,icode,runt)
!  call chkncstatus(istatus,nf90_noerr,'run',trim(schym(12)),'')
  istatus = nf90_close(lun)
  call chkncstatus(istatus,nf90_noerr,'close',trim(schym(12)),'')
  call getlun(lun)
  !  dij1 = 0.0083333
  !  slat = 27.00417 
  !  slon = -25.38744 
  OPEN(UNIT=lun, IOSTAT=ios, FILE=trim(adjustl(schym(11)))//'.drainage_net_ths50.csv', &
        STATUS='replace', ACTION='write')!,FORM='formatted',POSITION="append")
  write(lun,'(a)')"line_id,point_id,x,y,drainage_net,q100ch"
!  treshdr = MAXVAL(drai)/thdiv
  dimensions = shape(drai)
  print*,'DIMENSIONS',dimensions
  allocate(dummy(dimensions(1),dimensions(2)))
  dummy = drai
  treshdr = thdiv
  do while(run)
   call catchriver(nlat1,nlon1,ii,treshdr,run,dummy,dimensions(1),dimensions(2))
   ii = ii + 1
  end do
  close(lun)
    
end program dranet

subroutine catchriver(nlat1,nlon1,river,treshdr,run,dummy,n1,n2)
  use mod_statparams
  use mod_param
  use mod_museo
  implicit none
  integer :: nlon1,nlat1,ii,jj,river,n1,n2
  real :: maxdra,maxsubdra,treshdr
  integer :: jjmax,iimax,idir
  integer :: jjsmax,iismax
  integer, dimension(n1,n2) :: dummy
  character(100) :: a,b,c,d,e,f
  logical :: run
  if (mod(river,1000) == 0) print*,river
  i = 0
  maxdra = 0
  do jj = 1,nlat1
     do ii = 1,nlon1
        if (dummy(ii,jj) > maxdra ) then
           maxdra = dummy(ii,jj)
           iimax = ii
           jjmax = jj
        end if
     end do
  end do
  if (maxdra<=treshdr) run = .false.
  jjsmax = jjmax
  iismax = iimax
  ! add the linkin point between the two river chunks
  if(fmap(iimax,jjmax) .gt.0) then
!    print*,fmap(jjmax,iimax)
!    print*,iimax+ir(fmap(jjmax,iimax)),jjmax+jr(fmap(jjmax,iimax))
!    print*,luse(iimax+ir(fmap(jjmax,iimax)),jjmax+jr(fmap(jjmax,iimax)))
!    print*,area(iimax+ir(fmap(jjmax,iimax)),jjmax+jr(fmap(jjmax,iimax)))
    if (luse(iimax+ir(fmap(iimax,jjmax)),jjmax+jr(fmap(iimax,jjmax))).ne.mare.and. &
        iimax+ir(fmap(iimax,jjmax)).gt.1.and.iimax+ir(fmap(iimax,jjmax)).lt.nlon1.and. &
        jjmax+jr(fmap(iimax,jjmax)).gt.1.and.jjmax+jr(fmap(iimax,jjmax)).lt.nlat1.and. &
        fmap(iimax+ir(fmap(iimax,jjmax)),jjmax+jr(fmap(iimax,jjmax))).ne.0.and. &
        drai(iimax,jjmax).lt.drai(iimax+ir(fmap(iimax,jjmax)),jjmax+jr(fmap(iimax,jjmax)))) then
      write(a,'(i5)')river
!!      write(b,'(i5)')i
      i=i+1
      write(b,'(i5)')i
      write(c,'(f8.3)')lon(iimax+ir(fmap(iimax,jjmax)),jjmax+jr(fmap(iimax,jjmax)))
      write(d,'(f8.3)')lat(iimax+ir(fmap(iimax,jjmax)),jjmax+jr(fmap(iimax,jjmax)))
      write(e,'(f11.3)')drai(iimax,jjmax)
      write(f,*)area(iimax+ir(fmap(iimax,jjmax)),jjmax+jr(fmap(iimax,jjmax)))
      if (area(iimax+ir(fmap(iimax,jjmax)),jjmax+jr(fmap(iimax,jjmax))).lt.1e20) then
          write(lun,'(A,A,A,A,A,A,A,A,A,A,A)')trim(adjustl(a)),",",trim(adjustl(b)),",",trim(adjustl(c)),",",&
            trim(adjustl(d)),",",trim(adjustl(e)),",",trim(adjustl(f))
      end if
    end if
  end if
  do while(dummy(iimax,jjmax)>treshdr)
    iismax = iimax
    jjsmax = jjmax
    if (area(iismax,jjsmax).lt.1e20) then
      i=i+1
      write(a,'(i5)')river
      write(b,'(i5)')i
      write(c,'(f8.3)')lon(iimax,jjmax)
      write(d,'(f8.3)')lat(iimax,jjmax)
      write(e,'(f11.3)')drai(iimax,jjmax)
!      write(f,'(f15.5)')area(jjmax,iimax)
      write(f,*)area(iimax,jjmax)
      write(lun,'(A,A,A,A,A,A,A,A,A,A,A)')trim(adjustl(a)),",",trim(adjustl(b)),",",trim(adjustl(c)),",",&
        trim(adjustl(d)),",",trim(adjustl(e)),",",trim(adjustl(f))
    end if
!    WRITE(lun,'(I3,A,I4,A,f8.3,a,f8.3)') river,",",i,",",lon(jjmax,iimax),",",lat(jjmax,iimax)
    dummy(iimax,jjmax) = 0
    maxsubdra = 0
    do idir = 1,8
       !if (luse(iismax+ir(idir),jjsmax+jr(idir)).ne.mare.and. &
       if(fmap(iismax+ir(idir),jjsmax+jr(idir)).ne.0.and.  &
          jjsmax+jr(idir).gt.1.and.jjsmax+jr(idir).lt.nlat1.and. &
          iismax+ir(idir).gt.1.and.iismax+ir(idir).lt.nlon1) then
         if (dummy(iismax+ir(idir),jjsmax+jr(idir)) > maxsubdra) then
            maxsubdra = dummy(iismax+ir(idir),jjsmax+jr(idir))
            iimax = iismax+ir(idir)
            jjmax = jjsmax+jr(idir)
         end if
       end if
    end do
  end do
  return
end subroutine catchriver
