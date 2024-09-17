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
module mod_crtdyn

use mod_ncio
use mod_runparams
use mod_libmv
use mod_time
use mod_strman
use mod_param
use mod_internal
use mod_crtstatic
use mod_museo
use mod_statparams
use mod_interp
use mod_varandtypes
use mod_mpimess
use mod_mm5

implicit none
contains

  subroutine createnc1(hourstep)
    implicit none
    integer ora,giorno,mese,anno,yy,mm
    integer :: oldmese             ! = 01
    integer , intent(in) :: hourstep
    character(len=10) :: cfl
    character(len=135) :: mess
    logical :: first
    data first /.true./
    save first
    save mm,yy,oldmese

    call gmafromindex(time,ora,giorno,mese,anno)
!    if (hourstep == mchym(14)) then
    if ( first ) then
       read(schym(3)(5:6),'(i2)') mm
       read(schym(3)(1:4),'(i4)') yy
!       yy = ($CENTURY - 1) * 100 + yy
       oldmese = mese
       first = .false.
    end if
    if (mese /= oldmese .and. hourstep > 1) then
      mm = mm + 1
      if (mm >= 13) then
         mm =  1
           yy = yy + 1
      end if
      if (mm < 10) then
         write(cfl,'(i4,a,i1,a)') yy,'0',mm,'0100'
      else
         write(cfl,'(i4,i2,a)') yy,mm,'0100'
      end if
      call closefile
      write (mess,'(15x,a,a)') 'Creating new output file:  '// &
         trim(schym(11)),'_'//trim(cfl)//'.nc'
      write (6,'(a)') mess(1:len_trim(mess))
!      call createfile(trim(schym(11))//'_'//trim(cfl)//'.nc', &
!          mchym,rchym,schym,time)
      call createfile(trim(schym(11))//'_'//trim(cfl)//'.nc',time)
    end if
    oldmese = mese
    return
  end subroutine createnc1

  subroutine createnc2(hourstep)
    implicit none
    integer ora,giorno,mese,anno,yy,mm
    integer :: oldmese             ! = 01
    integer , intent(in) :: hourstep
    character(len=10) :: cfl
    character(len=135) :: mess
    logical :: first
    data first /.true./
    save first
    save mm,yy,oldmese

    call gmafromindex(time,ora,giorno,mese,anno)
!    if (hourstep == mchym(26)) then
    if ( first ) then
       read(schym(3)(5:6),'(i2)') mm
       read(schym(3)(1:4),'(i4)') yy
!       yy = ($CENTURY - 1) * 100 + yy
       oldmese = mm
       first = .false.
       call closerstfile
    end if
    if (mese /= oldmese .and. hourstep > 1) then
      mm = mm + 1
      if (mm >= 13) then
         mm =  1
           yy = yy + 1
      end if
      if (mm < 10) then
         write(cfl,'(i4,a,i1,a)') yy,'0',mm,'0100'
      else
         write(cfl,'(i4,i2,a)') yy,mm,'0100'
      end if
!      call closerstfile
      write (mess,'(15x,a)') 'Creating new restart file:  '// &
          trim(schym(11))//'_'//trim(cfl)//'_rst.nc'
      write (6,'(a)') mess(1:len_trim(mess))
      call createfile_rst(trim(schym(11))//'_'//trim(cfl)//'_rst.nc',time)
      call addrst_timestep
      call writerstfile
      call closerstfile
    end if
    oldmese = mese
    return
  end subroutine createnc2

  subroutine writeoutfile(flag)
    implicit none
    integer flag,len_trim,n,i,j
    character(len=80) :: scst,sstr1,string
    integer(long) , allocatable , dimension(:,:) :: avgport_int
    if (flag.eq.0) then
       call integer2string(time,10,sstr1)
       write (string,'(15x,a,i8,a)') 'CHyM integration step number ', &
           hourstep,': '//sstr1(1:10)
       write (logfile,'(a)') string(1:len_trim(string))
       write (6,'(a)') string(1:len_trim(string))
       return
    endif
    if (flag.eq.1) then
       call field2save(mchym,schym,'rai',j)
       if (j.eq.1) then
          avgrain = avgrain / mchym(14)
          if (mchym(20).gt.0) write(lun1) avgrain
          call write_dynvar('rai',avgrain)
       endif
       call field2save(mchym,schym,'rsr',j)
       if (j.eq.1) then
          avgrsrm = avgrsrm / mchym(14)
          if (mchym(20).gt.0) write(lun1) avgrsrm
          call write_dynvar('rsr',avgrsrm)
       endif
    else if (flag.eq.2) then
    else if (flag.eq.3) then
       write (logfile,'(45x,a)') 'Saving data at current step'
       call field2save(mchym,schym,'ara',j)
!       avgarai = avgarai / mchym(14)
       if (integrflag == 0) then
!         call runoff(avgarai/nsave)
!         call MPI_GATHER(port_sub(1,3),nlon*(nlat/nproc),MPI_REAL, &
!           port,nlon*(nlat/nproc),MPI_REAL,0,comm,ierr)
       else if(integrflag == 1) then
         call runoff1(avgarai/(nsave*3600))
       else
         call chymerror(23,integrflag,0.0,'integrflag')
       end if
       if (j.eq.1) then
          if (mchym(20).gt.0) write(lun1) avgarai
          call write_dynvar('ara',avgarai)
       endif
       call field2save(mchym,schym,'por',j)
       if (j.eq.1) then
          avgport = avgport / mchym(14)
          if (mchym(20).gt.0) write(lun1) avgport
          if (mchym(29) > 0) then
            if (.not. allocated(avgport_int)) allocate(avgport_int(nlon,nlat))
            avgport_int = nint((avgport - add_offset_por)/scale_factor_por)
            call write_dynvar('por',avgport_int)
          else
            call write_dynvar('por',avgport)
          end if
!          call write_dynvar('por',port)
       endif
       call field2save(mchym,schym,'wet',j)
       if (j.eq.1) then
          avgbwet = avgbwet / mchym(14)
          if (mchym(20).gt.0) write(lun1) avgbwet
          call write_dynvar('wet',avgbwet)
       endif
       call field2save(mchym,schym,'gwt',j)
       if (j.eq.1) then
          avggh2o = avggh2o / mchym(14)
          if (mchym(20).gt.0) write(lun1) avggh2o
          call write_dynvar('gwt',avggh2o)
       endif
       call field2save(mchym,schym,'gwt',j)
       if (j.eq.1) then
          avgh2o = avgh2o / mchym(14)
          if (mchym(20).gt.0) write(lun1) avgh2o
          call write_dynvar('h2o',avgh2o)
       endif
       call field2save(mchym,schym,'evp',j)
       if (j.eq.1) then
          avgevap = avgevap / mchym(14)
          if (mchym(20).gt.0) write(lun1) avgevap
          call write_dynvar('evp',avgevap)
       endif
       call field2save(mchym,schym,'sno',j)
       if (j.eq.1) then
          avgsnow = avgsnow / mchym(14)
          if (mchym(20).gt.0) write(lun1) avgsnow
          call write_dynvar('sno',avgsnow)
       endif
       call field2save(mchym,schym,'tem',j)
       if (j.eq.1) then
          avgtemp = avgtemp / mchym(14)
          if (mchym(20).gt.0) write(lun1) avgtemp
          call write_dynvar('tem',avgtemp)
       endif
       call field2save(mchym,schym,'dgw',j)
       if (j.eq.1) then
          avgdeepw = avgdeepw / mchym(14)
          if (mchym(20).gt.0) write(lun1) avgdeepw
          call write_dynvar('dgw',avgdeepw)
       endif
       call field2save(mchym,schym,'ddw',j)
       if (j.eq.1) then
          avgddeepw = avgddeepw / mchym(14)
          if (mchym(20).gt.0) write(lun1) avgddeepw
          call write_dynvar('ddw',avgddeepw)
       endif
!       call field2save(mchym,schym,'alf',j)
!       if (j.eq.1) then
!          if (mchym(20).gt.0) write(lun1) alfa
!          call write_dynvar('alf',alfa)
!       endif
      ! write (scst,'(i10)') time
      ! if (scst(5:10).eq.'123100') then
      !    call writetime('   Creating restart file at',time)
      !    call getlun(lun)
      !    open(lun,file='tmp/restart.chym',status='unknown', &
      !      form='unformatted')
      !    write(lun) nlon,nlat,rchym(1),rchym(2),rchym(9), &
      !      rchym(8),time
      !    write(lun) port
      !    write(lun) gh2o
      !    write(lun) snow
      !    write(lun) deepw
!     !    write(lun) ddeepw
      !    close(lun)
      ! endif
    endif
  end subroutine writeoutfile

  subroutine writerstfile
    implicit none
    character(len=80) :: sstr

    call integer2string(time,10,sstr1)
    write (string,'(15x,a,i5,a)') 'Restart file written at: '//sstr1(1:10)
    write (logfile,'(a)') string(1:len_trim(string))
    write (6,'(a)') string(1:len_trim(string))

    call writerst_dynvar_real('rgwt',gh2o)
    call writerst_dynvar_real('rh2o',h2o)
    call writerst_dynvar_real('rsno',snow)
    call writerst_dynvar_real('rpor',port)
    call writerst_dynvar_real('rdgw',deepw)
  end subroutine writerstfile

  subroutine runoff(rainloc)
    implicit none
    logical :: first
    data first /.true./
    save first
    integer ::  i,j,idir,imax,imin,jmax,jmin,tmin,step
    integer,save ::  i1,i2,j1,j2
    real dm,rainload,deltat
    save deltat,step
    real, dimension(nlon,nlat) :: cum , rainloc
    real :: cumulo
    cumulo = 0
    cum = 0
    if (myid==0) then
      write (logfile,'(16x,a)') '> Calculating Runoff and Flow Discharge'
    endif
    if (first) then
       call mpi_barrier(mycomm,mpierr)
       if (myid == 0) then
         call runoffspeed
       end if
       call mpi_bcast(alfa(1,1),nlat*nlon,MPI_REAL, 0,mycomm,mpierr)
       call mpi_bcast(dx(1,1),nlat*nlon,MPI_REAL, 0,mycomm,mpierr)
       call mpi_barrier(mycomm,mpierr)
       deltat=3600/chym_steps
!!     deltat=3600/mchym(22)
       first=.false.
!       print*,"maxval(dx)",maxval(dx)
       call mpi_barrier(mycomm,mpierr)
       i1 = ide1gb+1
       i2 = ide2gb-1
       j1 = jde1gb+1
       j2 = jde2gb-1
       if (ma%has_bdytop) i2 = ide2gb
       if (ma%has_bdybottom) i1 = ide1gb
       if (ma%has_bdyleft) j1 = jde1gb
       if (ma%has_bdyright) j2 = jde2gb
    endif
    do tmin=1,chym_steps          !Fmchym(22)
       wkm1sub=0.0
!       call exchange_bdy(h2osub,2,jde1ga,jde2ga,ide1ga,ide2ga,jde1gb,jde2gb,ide1gb,ide2gb)
!       call exchange_bdy(bwetsub,2,jde1ga,jde2ga,ide1ga,ide2ga,jde1gb,jde2gb,ide1gb,ide2gb)
!       call exchange_bdy(portsub,2,jde1ga,jde2ga,ide1ga,ide2ga,jde1gb,jde2gb,ide1gb,ide2gb)
       call exchange_bdy(h2osub(jde1gb:jde2gb,ide1gb:ide2gb),2,jde1,jde2,ide1,ide2,jde1gb,jde2gb,ide1gb,ide2gb)
       call exchange_bdy(bwetsub(jde1gb:jde2gb,ide1gb:ide2gb),2,jde1,jde2,ide1,ide2,jde1gb,jde2gb,ide1gb,ide2gb)
       call exchange_bdy(portsub(jde1gb:jde2gb,ide1gb:ide2gb),2,jde1,jde2,ide1,ide2,jde1gb,jde2gb,ide1gb,ide2gb)
       do i=i1,i2
          do j=j1,j2
!             if (i>=3 .and. i<=nlat-3) then
!             if (j>=3 .and. j<=nlon-3) then
             idir=fmap(j,i)
             if (luse(j,i).ne.mare.and.idir.ge.1.and.idir.le.8) then
                dm=portsub(j,i)*deltat
!                write(6,'(12x,2i4,2f9.4)') i,j,dm,port(i,j)
                if (dm.gt.h2osub(j,i)) dm=h2osub(j,i)
                wkm1sub(j,i)=wkm1sub(j,i)-dm
!                write(6,'(12x,2i4,2f9.4)') i,j,wrk1(i,j),port(i,j)
                wkm1sub(j+ir(idir),i+jr(idir))=wkm1sub(j+ir(idir),i+jr(idir))+dm
!             endif
!             endif
             endif
          enddo
       enddo
       call exchange_bdy(wkm1sub(jde1gb:jde2gb,ide1gb:ide2gb),2,jde1,jde2,ide1,ide2,jde1gb,jde2gb,ide1gb,ide2gb)
       do i=i1,i2
          do j=j1,j2
!             if (i>=3 .and. i<=nlat-3) then
!             if (j>=3 .and. j<=nlon-3) then
             idir=fmap(j,i)
             if (luse(j,i).ne.mare.and.idir.ge.1.and.idir.le.8) then
!!FF                rainload=area(i,j)*1.0e+03*(rainloc(i,j)/chym_steps)!/mchym(22))
!!F                rainload=area(i,j)*1.0e+03*(rainloc(i,j)/(3600*nsave))*deltat
                rainload=area(j,i)*1.0e+03*             &
                   rainloc(j,i)/chym_steps!F(3600*nsave))*deltat
                cumulo=cumulo + rainload
!               rainload=area(i,j)*1.0e+03*rundt(i,j)*deltat
!Fabio          !rainload=area(i,j)*1.0e+03*rain(i,j)/(mchym(22)*rchym(5)*100)
                h2osub(j,i)=h2osub(j,i)+wkm1sub(j,i)+rainload
                bwetsub(j,i)=h2osub(j,i)/dx(j,i)
                portsub(j,i)=alfa(j,i)*bwetsub(j,i)
!             endif
!             endif
             endif
          enddo
       enddo
    enddo
    return
  end subroutine runoff

  subroutine runoff1(rainloc)
    implicit none
    logical :: first
    data first /.true./
    save first
    integer ::  i,j,idir,imin,step,jmin
    real dm,rainload,deltat
    save deltat,step
    real, dimension(nlon,nlat) :: cum , rainloc
    cum = 0
    write (logfile,'(16x,a)') '> Calculating Runoff and Flow Discharge'
    if (first) then
       call runoffspeed
       deltat=(nsave*3600)/chym_steps
!       deltat=3600/mchym(22)
       first=.false.
    endif
    do imin=1,chym_steps          !Fmchym(22)
       wrk1=0.0
       do i=2,nlon-1
          do j=2,nlat-1
             idir=fmap(i,j)
             if (luse(i,j).ne.mare.and.idir.ge.1.and.idir.le.8) then
                dm=port(i,j)*deltat
!                write(6,'(12x,2i4,2f9.4)') i,j,dm,port(i,j)
                if (dm.gt.h2o(i,j)) dm=h2o(i,j)
                wrk1(i,j)=wrk1(i,j)-dm
!                write(6,'(12x,2i4,2f9.4)') i,j,wrk1(i,j),port(i,j)
                wrk1(i+ir(idir),j+jr(idir))=wrk1(i+ir(idir),j+jr(idir))+dm
             endif
          enddo
       enddo
       do i=2,nlon-1
          do j=2,nlat-1
             idir=fmap(i,j)
             if (luse(i,j).ne.mare.and.idir.ge.1.and.idir.le.8) then
                rainload=area(i,j)*1.0e+03*rainloc(i,j)*deltat
!                rainload=area(i,j)*1.0e+03*rainloc(i,j)/chym_steps
                cum(i,j)=rainload
!Fabio          !rainload=area(i,j)*1.0e+03*rain(i,j)/(mchym(22)*rchym(5)*100)
                h2o(i,j)=h2o(i,j)+wrk1(i,j)+rainload
                bwet(i,j)=h2o(i,j)/(dx(i,j)*(rchym(5)*100))
!                bwet(i,j)=h2o(i,j)/dx(i,j)
                port(i,j)=alfa(i,j)*bwet(i,j)
             endif
          enddo
       enddo
    enddo
    print*,'maxvaldx',maxval(dx)
    print*,'minvaldx',minval(dx)
    print*,'maxvalrainloc',maxval(rainloc)
    print*,"maxvalBWET:        ",maxval(bwet)
    print*,"RAINLOADDDDDDDDDDDDDDDDDDDDDDDDDDDD:        ",maxval(cum)
!    print*,"H2OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO:        ",sum(h2o)
    print*,"H2OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO:        ",maxval(h2o)
    print*,"POOOOOOOOOOOOOOORRRRRRRRRRRRRTTTTTT:        ",maxval(port)
  end subroutine runoff1

  subroutine temperature
    implicit none
    integer :: oldmonth,oldday
    data oldmonth,oldday /-1,-1/
    save oldmonth,oldday
    logical :: first
    data first /.true./
    save first
    real :: dl
    character(len=20) :: str
    character(len=40) :: now
    integer :: ora,giorno,mese,anno,i,j,ii,nn,nskip,jul
    integer :: len_trim

    call gmafromindex(time,ora,giorno,mese,anno)
    call datafromidx(time,now)
    if (first) then
       call getlun(lun)
       if (mchym(23).eq.2) then
          write(6,'(12x,a)') 'Temperature field will be read '// &
           'from '//schym(14)(1:len_trim(schym(14)))
          open(lun,file=schym(14),status='old',form='unformatted')
          read(lun) ii,nn
          if (ii.ne.time) write(6,'(12x,a)') 'Unaligned sequence of '// &
              'temperature in '//schym(14)(1:len_trim(schym(14)))
          if (indtimeal(ii,time).eq.1) then
             write(6,'(15x,a,i10,a,i10)') 'Skipping from ',ii,' to ',time
             nskip=nslicetoskip(ii,time)
             if (mchym(21).gt.(nn-nskip)) &
                  call chymerror(6,mchym(21),float(nn-nskip),'temperature')
             do i=1,nskip ; read(lun) ; enddo
          else if (indtimeal(ii,time).eq.-1) then
             call chymerror(5,ii,float(time),'temperature')
          endif
       else if (mchym(23).eq.1) then
          open(lun,file=schym(14),status='unknown',form='unformatted')
          write(lun) mchym(4),mchym(16)
       endif
       first=.false.
    endif
    if (mchym(23).eq.2) then
       read(lun) temp
    else if (mchym(15).eq.1) then
       if (mese.eq.oldmonth) then
          if (mchym(23).eq.1) write(lun) temp
          return
       endif
       call tempfromera40
       if (mchym(23).eq.1) write(lun) temp
       oldmonth=mese
    else if (mchym(15).eq.2) then
       if (giorno.eq.oldday) then
          if (mchym(23).eq.1) write(lun) temp
          return
       endif
       call tempfrommuseo(mchym(15))
       if (mchym(23).eq.1) write(lun) temp
       oldday=giorno
    else if (mchym(15).eq.3) then
       call tempfrommuseo(mchym(15))
    else if (mchym(15).eq.4) then
       call temperaturefield(4)
    else if (mchym(15).eq.5) then
       call temperaturefield(1)
    else if (mchym(15).eq.6) then
       call temperaturefield(6)
    else if (mchym(15).eq.7) then
       call temperaturefield(7)
    else if (mchym(15).eq.8) then
       if (giorno.eq.oldday) then
          if (mchym(23).eq.1) write(lun) temp
          return
       endif
       call temperaturefield(8)
       oldday=giorno
    else if (mchym(15).eq.21) then
       call temperaturefield(21)
    else if (mchym(15).eq.22) then
       call temperaturefield(22)
    else if (mchym(15).eq.23) then
       call temperaturefield(23)
    else if (mchym(15).eq.24) then
       call temperaturefield(24)
    else if (mchym(15).eq.25) then
       call temperaturefield(25)
    else if (mchym(15).eq.26) then
       call temperaturefield(26)
    else if (mchym(15).eq.27) then
       call temperaturefield(27)
    else if (mchym(15).eq.28) then
       call temperaturefield(28)
    else if (mchym(15).eq.29) then
       call temperaturefield(29)
    endif
    if (mchym(23).eq.1) write(lun) temp
    jul=indexofyear(giorno,mese,anno)
    evap=0.0
    do j=1,nlat ; dl=daylenght(lat(1,j),jul) ; do i=1,nlon
       if (temp(i,j).gt.2.0.and.luse(i,j).ne.mare.and.luse(i,j).ne.lago.and. &
                               luse(i,j).ne.fiume) then
          call evapot(temp(i,j),lat(i,j),lon(i,j),dem(i,j),jul,dl,evap(i,j))
          evap(i,j)=evap(i,j)*evc(luse(i,j),mese)*(1./1.)
          if (evap(i,j).lt.0.0) evap(i,j)=0.0
       endif
    enddo ; enddo
    evap=evap/24
  end subroutine temperature

  subroutine temperaturefield(iflag)
    implicit none
    character(len=20) :: lsource
    integer, dimension(nlon,nlat) :: lca
    integer :: ora,giorno,mese,anno
    integer :: i,j,k,i1,i2,j1,j2,ngoo,len_trim,np,ncacyc,n,iflag,ii,jj
    real :: avg,rad
    real, dimension(nlon,nlat) :: wk
    integer, parameter :: n1dv=1100000
    real, dimension(n1dv) :: xwk1, xwk2, xwk3
    character(len=80) :: scst
    character(len=300) :: simname
    logical :: useca = .true.
    character(:), allocatable :: cname
    character(:), allocatable :: cvar
    character(30) :: cdat
    integer :: anno2 , mese2
    call gmafromindex(time,ora,giorno,mese,anno)
    if (iflag.eq.1) then
       lsource='IntDB'
       rad=rchym(10)*1000
       ncacyc=10
       n=0
       call hourlytemp(ora,giorno,mese,anno,xwk1,xwk2,xwk3,n)
       call worldtemp (ora,giorno,mese,anno,xwk1,xwk2,xwk3,n)
    else if (iflag.eq.4) then
      lsource='RegCM' ; rad=50000.0 ; ncacyc=10
      if ( trim(schym(2)) == 'AFR44-ERAI') then
        simname = './museo/RegCM/AFR-44/ICTP/ECMWF-ERAINT/evaluation/r1i1p1/ICTP-RegCM4-3'// &
                  '/v4/day/tas/tas_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_ICTP-RegCM4-3_v4_day_'
      else if ( trim(schym(2)) == 'CORDEX_EUR-11_HadGEM_hist_MARCONI') then
        simname =trim(schym(16))//'/ICTP/HadGEM/historical/r1i1p1/' // &
                  'ICTP-RegCM4.6/v1/3hr/tas/tas_EURO-CORDEX_HadGEM_historical_r1i1p1_'//  &
                  'ICTP-RegCM4.6_v1_3hr_'
      else if ( trim(schym(2)) == 'CORDEX_EUR-44_HadGEM_rcp85_MARCONI') then
        simname =trim(schym(16))//'/ICTP/HadGEM/rcp85/r1i1p1/ICTP-RegCM4.6/v1/3hr/tas/' // &
                 'tas_EURO-CORDEX_HadGEM_rcp85_r1i1p1_ICTP-RegCM4.6_v1_3hr_'
      else if ( trim(schym(2)) == 'standard') then
!        rad=12000.0
        simname = trim(schym(16))
      else if ( trim(schym(2)) == 'CORDEX' ) then
        anno2 = anno ; mese2 = mese + 1
        if ( mese == 12 ) then
           anno2=anno+1 ; mese2 = 01
        end if
        write (cdat,'(i6,a,i6,a)') anno*100 + mese ,'010300-',anno2*100 + mese2,'010000.nc'
        cvar = "tas"
        call cordexnamefile(trim(schym(16)),cdat,cvar,cname)
        simname = cname
      end if
      if ( nreg1==0 ) then
        call regcmgetdims(ora,giorno,mese,anno,nreg1,nreg2,simname)
        if ( rchym(11) == -1) then
          if ( .not.allocated(lati) ) allocate(lati(nreg1,nreg2))
          if ( .not.allocated(loni) ) allocate(loni(nreg1,nreg2))
          if ( .not.allocated(rani) ) allocate(rani(nreg1,nreg2))
          if ( .not.allocated(temi) ) allocate(temi(nreg1,nreg2))
          if ( .not.allocated(mapi) ) allocate(mapi(nlon,nlat))
          if ( .not.allocated(mapj) ) allocate(mapj(nlon,nlat))
        end if
      end if
      call regcmreadnc(nreg1,nreg2,ora,giorno,mese,anno,2,xwk1,xwk2,xwk3,n,simname)
      if ( trim(schym(2)) == 'standard' .or. trim(schym(2)(1:6)) == "CORDEX") then
        if (rchym(11) > 0) then
          rad = rchym(11) * 1000.0
        else if (rchym(11) == 0.) then
          rad = ceiling(((regcm_res/2.)*sqrt(2.))/1000)*1000
        else if (rchym(11) ==-1.) then
          !no radius NN interpolation
        end if
      end if
      if (rchym(11) == -1) then
       if ( nnflag ) then
       do ii = 1, nlat
         do jj = 1, nlon
           call closest2d(lati,loni,lat(jj,ii),lon(jj,ii),nreg1,nreg2,mapi(jj,ii),mapj(jj,ii),.true.)
         enddo
       enddo
       nnflag = .false.
       end if
       do i=1,nlat
         do j=1,nlon
            temp(j,i) = temi(mapj(j,i),mapi(j,i))
         end do
       end do
       return
      end if
!      call regcmreadnc(ora,giorno,mese,anno,2,xwk1,xwk2,xwk3,n)
    else if (iflag.eq.6) then
!      lsource='EIN75' ; rad=75000.0 ; ncacyc=10
       lsource='EIN75' ; rad=90000.0 ; ncacyc=30 ; useca=.true.
       call ein75readnc(480,241,ora,giorno,mese,anno,1,xwk1,xwk2,xwk3,n)
    else if (iflag.eq.7) then
!      lsource='EIN75' ; rad=75000.0 ; ncacyc=10
       lsource='ERA5' ; rad=30000.0 ; ncacyc=30 ; useca=.true.
       call ein75readnc(1440,721,ora,giorno,mese,anno,1,xwk1,xwk2,xwk3,n)
    else if (iflag.eq.8) then
!      lsource='EIN75' ; rad=75000.0 ; ncacyc=10
       lsource='ERA5' ; rad=30000.0 ; ncacyc=30 ; useca=.true.
       call ein75readnc(129,120,ora,giorno,mese,anno,1,xwk1,xwk2,xwk3,n)
!n    else if (iflag.eq.21) then
!n       lsource='REMO25x25' ; rad=30000.0 ; ncacyc=10
!n       call acqwascen01(ora,giorno,mese,anno,2,xwk1,xwk2,xwk3,n)
!n    else if (iflag.eq.22) then
!n       lsource='RegCM25x25' ; rad=30000.0 ; ncacyc=30
!n       call acqwascen02(ora,giorno,mese,anno,6,xwk1,xwk2,xwk3,n)
!n    else if (iflag.eq.23) then
!n       lsource='REMO10x10' ; rad=15000.0 ; ncacyc=10
!n       call acqwascen03(ora,giorno,mese,anno,2,xwk1,xwk2,xwk3,n)
!n    else if (iflag.eq.24) then
!n       lsource='RegCM3x3'
!n       rad=3000.0
!n       ncacyc=30
!n       call acqwascen04(ora,giorno,mese,anno,2,xwk1,xwk2,xwk3,n)
!n    else if (iflag.eq.25) then
!n       lsource='PPREMO25x25'
!n       rad=30000.0
!n       ncacyc=30
!n       call acqwascen5a(ora,giorno,mese,anno,2,xwk1,xwk2,xwk3,n)
!n    else if (iflag.eq.26) then
!n       lsource='PP REMO10'
!n       rad=2000.0
!n       ncacyc=10
!n       call acqwascen5b(ora,giorno,mese,anno,2,xwk1,xwk2,xwk3,n)
!n    else if (iflag.eq.27) then
!n       lsource='PPRegCM25x25'
!n       rad=30000.0
!n       ncacyc=30
!n       call acqwascen6a(ora,giorno,mese,anno,2,xwk1,xwk2,xwk3,n)
!n    else if (iflag.eq.28) then
!n       lsource='PPRegCM3x3'
!n       rad=2000.0
!n       ncacyc=10
!n       call acqwascen6b(ora,giorno,mese,anno,2,xwk1,xwk2,xwk3,n)
!n    else if (iflag.eq.29) then
!n       lsource='PPRegCM50x50'
!n       rad=50000.0
!n       ncacyc=10
!n       call acqwascenc1(ora,giorno,mese,anno,2,xwk1,xwk2,xwk3,n)
    else
       call chymerror(22,iflag,0.0,lsource)
    endif
    if (n.gt.n1dv) call chymerror(4,n,0.0,lsource)
    write (logfile,'(16x,a)')'> Building temperature field using '//lsource
    lca=-1 ; np=nint((rad)/rchym(6)) ; ngoo=0 ; avg=0
    do k=1,n
       i=nint((xwk3(k)-rchym(1))/rchym(8))+1
       j=nint((xwk2(k)-rchym(2))/rchym(9))+1
       if (inside(i,1,nlon).and.inside(j,1,nlat)) then
          lca(i,j)=0 ; wk(i,j)=xwk1(k) ;avg=avg+xwk1(k) ; ngoo=ngoo+1
          i1=i-np ; i2=i+np ; j1=j-np ; j2=j+np
          if (i1.lt.1) i1=1 ; if (i2.gt.nlon) i2=nlon
          if (j1.lt.1) j1=1 ; if (j2.gt.nlat) j2=nlat
          do i=i1,i2 ; do j=j1,j2
             if (distance(lat(i,j),lon(i,j),xwk2(k),xwk3(k)).lt.rad) then
                if (lca(i,j).eq.-1) then
                   lca(i,j)=1
                   wk(i,j)=xwk1(k)
                else if (lca(i,j).gt.0) then
                   lca(i,j)=lca(i,j)+1
                   wk(i,j)=wk(i,j)+xwk1(k)
                endif
             endif
          enddo ; enddo
       endif
    enddo
    if (ngoo.gt.5) then
       write (scst,'(i8,a,i8,a)') ngoo,' of ',n,' data inside domain'
       call no2space(scst) ; call noinspace(scst)
       write (logfile,'(18x,a)') scst(1:len_trim(scst))
       avg=avg/ngoo
    else
       write (logfile,'(18x,a)') 'Too low number of observations, '// &
          'temperature field unchanged.'
       return
    endif
    do i=1,nlon ; do j=1,nlat
       if (lca(i,j).gt.0) then
          temp(i,j)=wk(i,j)/lca(i,j)
          lca(i,j)=8
!         wrk1(i,j)=10.
       else if (lca(i,j).eq.0) then
          temp(i,j)=wk(i,j)
!         wrk1(i,j)=10.
       else
          temp(i,j)=avg
          lca(i,j)=-1
!         wrk1(i,j)=20.
       endif
    enddo ; enddo
    if ((mchym(15).eq.5.and.schym(10)(1:10).eq.'acqwarhone').or. &
       (mchym(15).eq.27.and.schym(10)(1:7).eq.'acqwapo')) then
        call cressman(xwk1,xwk2,xwk3,n,1,30.0,temp,lat,lon,nlon,nlat)
    endif
!    if (useca) then
      do i=1,ncacyc
        call d2cellcycle(temp,lca,wk,nlon,nlat,0.9)
      enddo
!    endif
  end subroutine temperaturefield

  subroutine rainfall
    implicit none
    logical :: first
    integer :: i , ih , ii , j , lun , nn , nskip
    integer :: ora,giorno,mese,anno
    integer :: oldday
    data first/.true./
    data oldday/-1/
    save first,oldday
    call gmafromindex(time,ora,giorno,mese,anno)
    if ( first ) then
      call getlun(lun)
      if ( mchym(17)==2 ) then
        write (6,'(12x,a)') 'Precipitation field will be read '//'from '// &
                            trim(schym(13))
        open (lun,file = schym(13),status = 'old',form = 'unformatted')
        read (lun) ii , nn
        if ( ii/=idate ) write (6,'(12x,a)') 'Unaligned sequence of '//  &
                                'precipitation in '//trim(schym(13))
        if ( indtimeal(ii,idate)==1 ) then
          write (6,'(15x,a,i10,a,i10)') 'Skipping from ' , ii , ' to ' , idate
          nskip = nslicetoskip(ii,idate)
          if ( mchym(21)>(nn-nskip) ) then
            call chymerror(6,mchym(21),float(nn-nskip),'rainfall')
          end if
          do i = 1 , nskip
            read (lun)
          end do
        else if ( indtimeal(ii,idate)==-1 ) then
          call chymerror(5,ii,float(idate),'precipitation')
        end if
      else if ( mchym(17)==1 ) then
        open (lun,file = schym(13),status = 'unknown',form = 'unformatted')
        write (lun) mchym(4) , mchym(16)
      end if
      first = .false.
    end if
    if ( mchym(17)==2 ) then
      read (lun) rain
      deepw(1:nlon,1:nlat) = deepw(1:nlon,1:nlat) + rain(1:nlon,1:nlat)
      return
    end if
    call writestatus('Building rain fields for '//now)
! CA Legend
! 8  point should be calculated in the current module
! 0  point is already calculated
! -1  point cannot be calculated in the current module
! -2  point must not be considered (for ex. sea point)
    ih = idhour
    caplot = 0
    rain = 0.0
    rsrm = 0.0
    ca = 8
    do i = 1 , nlon
      do j = 1 , nlat
        if ( luse(i,j)==mare ) ca(i,j) = -2
      end do
    end do
    actsource = ' '
    if (srcflag(15)) call buildrainfield(15,ca)     ! raingauges from text file
    if (srcflag(01)) call buildrainfield(1,ca)      ! raingauges
    if (srcflag(30)) call buildrainfield(30,ca)     ! filtered raingauges
    if (srcflag(31)) call buildrainfield(31,ca)     ! filtered raingauges nc
    if (srcflag(02)) call buildrainfield(2,ca)      ! radar
    if (srcflag(03)) call buildrainfield(3,ca)      ! ERA-Interim
    if (srcflag(05)) call buildrainfield(5,ca)      ! WRF Domain 2
    if (srcflag(04)) call buildrainfield(4,ca)      ! WRF Domain 1
    if (srcflag(06)) call buildrainfield(6,ca)      ! tmp file
    if (srcflag(13)) call buildrainfield(13,ca)     ! museo
    if (srcflag(12)) call buildrainfield(12,ca)     ! mm5
    if (srcflag(07)) call buildrainfield( 7,ca)     ! RegCM
    if (srcflag(08)) call buildrainfield(8,ca)      ! PERSIANN
    if (srcflag(09)) call buildrainfield(9,ca)      ! TRIMM
    if (srcflag(10)) call buildrainfield(10,ca)     ! EIN75
    if (srcflag(11)) call buildrainfield(11,ca)     ! EIN15
    if (srcflag(14)) call buildrainfield(14,ca)     ! ERA5
    if (srcflag(16)) then
       if (giorno.eq.oldday) then
          return
       endif
       call buildrainfield(16,ca)     ! ERA5 (limited area)
       oldday=giorno
     end if
!    if (srcflag(21)) call buildrainfield(21,ca)     ! ACQWA-ReMo 25x25
!    if (srcflag(22)) call buildrainfield(22,ca)     ! ACQWA-RegCM 25x25
!    if (srcflag(23)) call buildrainfield(23,ca)     ! ACQWA-ReMo 10x10
!    if (srcflag(24)) call buildrainfield(24,ca)     ! ACQWA-RegCM 3x3
!    if (srcflag(25)) call buildrainfield(25,ca)     ! ACQWA-PP REMO 25x25
!    if (srcflag(26)) call buildrainfield(26,ca)     ! ACQWA-PP REMO 10x10
!    if (srcflag(27)) call buildrainfield(27,ca)     ! ACQWA-PP RegCM 25x25
!    if (srcflag(28)) call buildrainfield(28,ca)     ! ACQWA-PP RegCM 10x10
!    if (srcflag(29)) call buildrainfield(29,ca)     ! ACQWA-C1 RegCM 50x50
    do i = 1 , nlon
      rain(i,1) = 0.0
      rain(i,nlat) = 0.0
    end do
    do j = 1 , nlat
      rain(1,j) = 0.0
      rain(nlon,j) = 0.0
    end do
    deepw(1:nlon,1:nlat) = deepw(1:nlon,1:nlat) + rain(1:nlon,1:nlat)
    if ( mchym(17)==1 ) write (lun) rain
  end subroutine rainfall

  subroutine buildrainfield(idata,ca)     ! former datisparsi
!  Poiche' si passa ca a d2cellcycle viene aggiornata tutta la mappa,
!  anche nei punti distanti piu' di radius. conviene sempre?
    implicit none
    integer :: idata
    integer, dimension(nlon,nlat) :: ca, canc
    integer :: ora,giorno,mese,anno,i,j,ic,index,n,ncacyc,ngreat,srccode,nin
    integer nls,len_trim,ngood
    character(len=40) :: data
    character(len=24) :: lsource
    character(len=132) :: scst
    character(len=300) :: simname
    character(len=80) :: now,sstr1,sstr2
    real :: radius
    logical :: useca,piove
!Fabio-persiann         integer, parameter :: mxnicorw=200000
!Fabio-Era5    integer, parameter :: mxnicorw=700000
    integer, parameter :: mxnicorw=1100000
    real, dimension(mxnicorw) :: la,lo,pi,pi2val
!    real, allocatable, dimension(:) :: la,lo,pi,pi2val
    logical :: first
    data first /.true./
    save first
    real :: dij
    save pi,radius
    real :: distaic,minlat,minlon,regcm_res_deg
    integer :: ii,jj
    real :: tallint,tallintf,tallinti
    character(:), allocatable :: cname
    character(:), allocatable :: cvar
    character(30) :: cdat
    integer :: anno2 , mese2


    call gmafromindex(time,ora,giorno,mese,anno)
    index=indexofyear(giorno,mese,anno)
    call dataorafromday (ora,giorno,mese,anno,data)
    now=data
    if ( idata==1 .or. idata==30 .or. idata==31 .or. idata==15) then
      if ( idata==1 ) then
        lsource = 'RainGauges'
      else if ( idata==30 ) then
        lsource = 'Fil-RainGauges'
      else if ( idata==31 ) then
        lsource = 'Fil-RainGauges-nc'
      else if ( idata==15 ) then
        lsource = 'Text-RainGauges'
      end if
      nls = len_trim(lsource)
      radius = rchym(10)
      useca = .true.
      ncacyc = 100
      srccode = 1
      if ( idata==1 ) then
        write (logfile,'(16x,a)') '> Using Rain gauges by '//lsource
      else if ( idata==30 ) then
        write (logfile,'(16x,a)') '> Using filtered Rain gauges by '//lsource
      else if ( idata==31 ) then
        write (logfile,'(16x,a)') '> Using NetCDF filtered Rain gauges by '//lsource
      else if ( idata==15 ) then
        write (logfile,'(16x,a)') '> Using Rain gauges by '//lsource
      end if
      n = 0
      if ( idata==1 .or. idata==30) then
        if (year < 2001 .and. first) then
          write (6,'(15x,a)') '>WARNING< The '//lsource //' dataset is avalaible only from 2001'
          first = .false.
        else if (year >=2001) then
          call hourlyrain(ora,giorno,mese,anno,pi,la,lo,n)
        end if
      else if ( idata==31 ) then
        if (year < 2001 .and. first) then
          write (6,'(15x,a)') '>WARNING< The '//lsource //' dataset is avalaible only from 2001'
          first = .false.
        else if (year >=2001) then
          call hourlyrainnc(ora,giorno,mese,anno,pi,la,lo,n)
        end if
      else if ( idata==15 ) then
        call read_stations(ora,giorno,mese,anno,pi,la,lo,n)
      end if
!      call worldrain (ora,giorno,mese,anno,pi,la,lo,n)
      pi2val(1:n)=pi(1:n)
      call raingaugesval(pi2val,pi,la,lo,n)
!     call marcoverdecchia(ora,giorno,mese,anno,pi,la,lo,n)  !Was in chymf90 branch
    else if ( idata==2 ) then
      lsource = 'Radar'
      nls = len_trim(lsource)
      useca = .true.
      ncacyc = 5
      radius = 2.0
      srccode = 2
      write (logfile,'(16x,a)') '> Using rain rate estimates by '//lsource
      call montemidiavec(index1h(ora,giorno,mese,anno),anno,pi,la,lo,n)
    else if ( idata==3 ) then
      lsource = 'ERA-Interim'
      nls = len_trim(lsource)
      radius = 65.0
      useca = .true.
      ncacyc = 30
      srccode = 3
      write (logfile,'(16x,a)') '> Using re-analysis by '//lsource(1:nls) &
                                //' project.'
      call erainterim(ora,giorno,mese,anno,pi,la,lo,n)
    else if ( idata==7 ) then
      useca = .false.
      ncacyc = 0
      srccode = 7
      lsource = 'RegCM'
      nls = len_trim(lsource)
      write (logfile,'(16x,a)') '> Using Rain Rate prediction by '//lsource
      if ( trim(schym(2)) == 'AFR44-ERAI') then
         radius = 50.0
         simname = './museo/RegCM/AFR-44/ICTP/ECMWF-ERAINT/evaluation/r1i1p1/ICTP-RegCM4-3'// &
                   '/v4/day/pr/pr_AFR-44_ECMWF-ERAINT_evaluation_r1i1p1_ICTP-RegCM4-3_v4_day_'
      else if ( trim(schym(2)) == 'SAM44-ERAI') then
         radius = 50.0
         simname = './museo/RegCM/SAM-44/ICTP/ECMWF-ERAINT/evaluation/r1i1p1/ICTP-RegCM4-3'// &
                   '/v4/day/pr/pr_SAM-44_ECMWF-ERAINT_evaluation_r1i1p1_ICTP-RegCM4-3_v4_day_'
      else if ( trim(schym(2)) == 'MED44-ERAIv4') then
         radius = 50.0
         simname = './museo/RegCM/MED-44/ICTP/ECMWF-ERAINT/evaluation/r1i1p1/ICTP-RegCM4-3'// &
                   '/v4/day/pr/pr_MED-44_ECMWF-ERAINT_evaluation_r1i1p1_ICTP-RegCM4-3_v4_day_'
      else if ( trim(schym(2)) == 'MED11-ERAI') then
         radius = 50.0
         simname = './museo/RegCM/MED-11/ICTP/ECMWF-ERAINT/evaluation/r1i1p1/ICTP-RegCM4-3'// &
                   '/v4/day/pr/pr_MED-11_ECMWF-ERAINT_evaluation_r1i1p1_ICTP-RegCM4-3_v4_day_'
      else if ( trim(schym(2)) == 'EUR44-ERAI') then
         radius = 50.0
         simname = './museo/RegCM/EUR-44/ICTP/ECMWF-ERAINT/evaluation/r1i1p1/ICTP-RegCM4-3'// &
                   '/v1/day/pr/pr_EUR-44_ECMWF-ERAINT_evaluation_r1i1p1_ICTP-RegCM4-3_v1_day_'
      else if ( trim(schym(2)) == 'EUR44-MOHC-his') then
         radius = 50.0
         simname = './museo/RegCM/EUR-44/ICTP/MOHC-HadGEM2-ES/historical/r1i1p1/ICTP-RegCM4-3'// &
                   '/v1/day/pr/pr_EUR-44_MOHC-HadGEM2-ES_historical_r1i1p1_ICTP-RegCM4-3_v1_day_'
      else if ( trim(schym(2)) == 'CAM44-MPI-his') then
         radius = 50.0
         simname = './museo/RegCM/CAM-44/ICTP/MPI-M-MPI-ESM-MR/historical/r1i1p1/ICTP-RegCM4-3'// &
                   '/v4/day/pr/pr_CAM-44_MPI-M-MPI-ESM-MR_historical_r1i1p1_ICTP-RegCM4-3_v4_day_'
      else if ( trim(schym(2)) == 'CAM44-MPI-rcp85') then
         radius = 50.0
         simname = './museo/RegCM/CAM-44/ICTP/MPI-M-MPI-ESM-MR/rcp85/r1i1p1/ICTP-RegCM4-3'// &
                   '/v4/day/pr/pr_CAM-44_MPI-M-MPI-ESM-MR_rcp85_r1i1p1_ICTP-RegCM4-3_v4_day_'
      else if ( trim(schym(2)) == 'AFR44-MPI-his') then
         radius = 50.0
         simname = './museo/RegCM/AFR-44/ICTP/MPI-ESM-MR/historical/r1i1p1/ICTP-RegCM4-3'// &
                   '/4.3-rc15/day/pr/pr_AFR-44_MPI-ESM-MR_historical_r1i1p1_ICTP-RegCM4-3_4.3-rc15_day_'
      else if ( trim(schym(2)) == 'AFR44-MPI-rcp85') then
         radius = 50.0
         simname = './museo/RegCM/AFR-44/ICTP/MPI-ESM-MR/rcp85/r1i1p1/ICTP-RegCM4-3'// &
                   '/4.3-rc15/day/pr/pr_AFR-44_MPI-ESM-MR_rcp85_r1i1p1_ICTP-RegCM4-3_4.3-rc15_day_'
      else if ( trim(schym(2)) == 'AFR44-MPIMR-his') then
         radius = 50.0
         simname = './museo/RegCM/AFR-44/ICTP/MPI-M-MPI-ESM-MR/historical/r1i1p1/ICTP-RegCM4-3'// &
                   '/v4/day/pr/pr_AFR-44_MPI-M-MPI-ESM-MR_historical_r1i1p1_ICTP-RegCM4-3_v4_day_'
      else if ( trim(schym(2)) == 'AFR44-MPIMR-rcp85') then
         radius = 50.0
         simname = './museo/RegCM/AFR-44/ICTP/MPI-M-MPI-ESM-MR/rcp85/r1i1p1/ICTP-RegCM4-3'// &
                   '/v4/day/pr/pr_AFR-44_MPI-M-MPI-ESM-MR_rcp85_r1i1p1_ICTP-RegCM4-3_v4_day_'
      else if ( trim(schym(2)) == 'CORDEX_EUR-11_HadGEM_hist_MARCONI') then
        simname =trim(schym(16))//'/ICTP/HadGEM/historical/r1i1p1/' // &
                  'ICTP-RegCM4.6/v1/3hr/pr/pr_EURO-CORDEX_HadGEM_historical_r1i1p1_'//  &
                  'ICTP-RegCM4.6_v1_3hr_'
      else if ( trim(schym(2)) == 'CORDEX_EUR-44_HadGEM_rcp85_MARCONI') then
        simname =trim(schym(16))//'/ICTP/HadGEM/rcp85/r1i1p1/ICTP-RegCM4.6/v1/3hr/pr/' // &
                 'pr_EURO-CORDEX_HadGEM_rcp85_r1i1p1_ICTP-RegCM4.6_v1_3hr_'
      else if ( trim(schym(2)) == 'standard') then
         simname = trim(schym(16))
      else if ( trim(schym(2)) == 'CORDEX' ) then
        anno2 = anno ; mese2 = mese + 1
        if ( mese == 12 ) then
           anno2=anno+1 ; mese2 = 01
        end if
        write (cdat,'(i6,a,i6,a)') anno*100 + mese ,'010300-',anno2*100 +mese2,'010000.nc'
        cvar = "pr"
        call cordexnamefile(trim(schym(16)),cdat,cvar,cname)
        simname = cname
      else
         write (6,'(15x,a)')'Bad flag for RegCM data. Please check chym_ifile1. Exiting...'
         call exit(0)
      end if
      if ( nreg1==0 ) then
        call regcmgetdims(ora,giorno,mese,anno,nreg1,nreg2,simname)
        if (rchym(11) == -1) then
          if ( .not.allocated(lati) ) allocate(lati(nreg1,nreg2))
          if ( .not.allocated(loni) ) allocate(loni(nreg1,nreg2))
          if ( .not.allocated(rani) ) allocate(rani(nreg1,nreg2))
          if ( .not.allocated(mapi) ) allocate(mapi(nlon,nlat))
          if ( .not.allocated(mapj) ) allocate(mapj(nlon,nlat))
        end if
      end if

      call regcmreadnc(nreg1,nreg2,ora,giorno,mese,anno,1,pi,la,lo,n,simname)
      if ( trim(schym(2)) == 'standard' .or. trim(schym(2)(1:6)) == 'CORDEX') then
        if (rchym(11) > 0.) then
          radius = rchym(11)
        else if (rchym(11) == 0.) then
          radius = ceiling((regcm_res/2.)*sqrt(2.)/1000.)
        else if (rchym(11) ==-1.) then
          !no radius NN interpolation
        else
          write (6,'(15x,a)')'Bad flag for RegCM int Radius. Please check chym_regcm_rad. Exiting...'
          call exit(1)
        end if
      end if

    else if (idata.eq.8) then
      useca=.false. ; ncacyc=30 ; srccode=30 ; radius=27.5
      lsource='persiann' ; nls=len_trim(lsource)
      write (logfile,'(16x,a)') '> Using Rain estimates by '//lsource
      call persiannreadncll(1440,480,ora,giorno,mese,anno,1,pi,la,lo,n)
    else if (idata.eq.9) then
!      useca=.true. ; ncacyc=30 ; srccode=50 ; radius=27.5
      useca=.false. ; ncacyc=30 ; srccode=50 ; radius=27.5    !Fabio
      lsource='trimm' ; nls=len_trim(lsource)
      write (logfile,'(16x,a)') '> Using Rain estimates by '//lsource
      call trimmreadncll(1440,400,ora,giorno,mese,anno,1,pi,la,lo,n)
    else if (idata.eq.10) then
!      useca=.true. ; ncacyc=10 ; srccode=70 ; radius=70.0
      useca=.false. ; ncacyc=10 ; srccode=70 ; radius=80.0
      lsource='ein75' ; nls=len_trim(lsource)
      write (logfile,'(16x,a)') '> Using re-analysis by '//lsource
!      call erareadncll(480,241,ora,giorno,mese,anno,1,n,10)
!      call erareadncll(480,241,ora,giorno,mese,anno,1,pi,la,lo,n,10)
      call erareadncll3(480,241,ora,giorno,mese,anno,1,pi,la,lo,n,10)
    else if (idata.eq.11) then
      useca=.true. ; ncacyc=10 ; srccode=80 ; radius=165.0
      lsource='ein15' ; nls=len_trim(lsource)
      write (logfile,'(16x,a)') '> Using re-analysis by '//lsource
      call erareadncll(240,121,ora,giorno,mese,anno,1,pi,la,lo,n,8)
!FF      call erareadncll(30,20,ora,giorno,mese,anno,1,pi,la,lo,n,8)
    else if (idata.eq.12) then
      useca=.true. ; ncacyc=30 ; srccode=10 ; radius=3**(4-mchym(12))
      lsource='MM5' ; nls=len_trim(lsource)
      call mm5data(ora,giorno,mese,anno,pi,la,lo,n)
    else if (idata.eq.13) then
      useca=.true. ; ncacyc=30 ; srccode=11 ; radius=3**(4-mchym(12))
      lsource='MuSEO' ; nls=len_trim(lsource)
      call museodata(ora,giorno,mese,anno,pi,la,lo,n)
    else if (idata.eq.14) then
      useca=.false. ; ncacyc=10 ; srccode=90 ; radius=30.0
      lsource='era5' ; nls=len_trim(lsource)
      write (logfile,'(16x,a)') '> Using re-analysis by '//lsource
      call erareadncll3(1440,721,ora,giorno,mese,anno,1,pi,la,lo,n,14)
    else if (idata.eq.16) then
      useca=.false. ; ncacyc=10 ; srccode=90 ; radius=30.0
      lsource='era5' ; nls=len_trim(lsource)
      write (logfile,'(16x,a)') '> Using re-analysis by '//lsource
      call erareadncll3(129,120,ora,giorno,mese,anno,1,pi,la,lo,n,16)
!    else if (idata.eq.15) then
!      useca=.false. ; ncacyc=30 ; srccode=99 ; radius=27.5
!      lsource='persiann' ; nls=len_trim(lsource)
!      write (logfile,'(16x,a)') '> Using Rain estimates by '//lsource
!      call readnetcdf(1440,480,ora,giorno,mese,anno,1,pi,la,lo,n)
    else
       return
    endif
    nls=len_trim(lsource)
    write (logfile,'(16x,a)') '> Building rainfall field from '//lsource
!   write (*,'(16x,a)') '> Building rainfall field from '//lsource !Fabio-nc
    if (n.gt.mxnicorw) call chymerror(3,n,0.0,lsource)
    nin=0 ; ngood=0 ; piove=.false. ; ngreat=0 ; dij=rchym(5)

    do ic=1,n
      i=nint((lo(ic)-rchym(1))/dij)+1
      j=nint((la(ic)-rchym(2))/dij)+1
      if (inside(i,1,nlon).and.inside(j,1,nlat)) then
        nin=nin+1
        if (ca(i,j).eq.8.and.pi(ic).lt.100.0) then
          if (pi(ic) <= -1.0) then
             ca(i,j) = 8
          else
          ngood=ngood+1
          rain(i,j)=pi(ic)
!          if (idata.eq.1) then
!           ca(i,j)=0
!          else
!            ca  (i,j)=8 !F Per risolvere il problema nell\'interpolazione (punti singolari)
!            ca  (i,j)=1 !Originale
            ca  (i,j)=0 !Originale
!          endif
          rsrm(i,j)=srccode
          if (rain(i,j).gt.1.0) ngreat=ngreat+1
          end if
        endif
      else
!Fabio        pi(ic) = -2
      endif

    end do
    if (ngreat.gt.1) piove=.true.
    write (scst,'(i6,a,i5,a,i5,a)') n,' data available, ',nin, &
         ' inside domain, ',ngood,' actually used.'
    call noinspace(scst)
    call no2space(scst)
    write (logfile,'(18x,a)') scst(1:len_trim(scst))
    if (ngood.lt.1) then
      write (logfile,'(21x,a)') 'too low number of data, '// &
          lsource(1:nls)//' not used at current step.'
!      if (allocated(pi)) deallocate(pi)
!      if (allocated(la)) deallocate(la)
!      if (allocated(lo)) deallocate(lo)
      return
    else
      call addsource(lsource,rainsrc)
    endif
!    CALL CPU_TIME(tallinti)
    if ( idata == 7 ) then
      if (rchym(11) == -1) then
        write (logfile,'(18x,a)') "Nearest Neighbour Interpolation method for RegCM precipitation applyed"
        ! if ( first ) then
        !   minlat=minval(lati)
        !   minlon=minval(loni)
        !   regcm_res_deg = regcm_res/111111
        !   do j=1,nlat
        !     do i=1,nlon
        !       do ic = 1,n
        !       end do
        !
        !       mapi(i,j) = ii
        !       mapj(i,j) = jj
        !       !distaic = 100000000
        !       !do ic = 1,n
        !       !  if (inside(i,1,nlon).and.inside(j,1,nlat)) then
        !       !  if (distance(la(ic),lo(ic),lat(i,j),lon(i,j)) < distaic) then
        !       !    distaic = distance(la(ic),lo(ic),lat(i,j),lon(i,j))
        !       !    mapic(i,j) = ic
        !       !  end if
        !       !  end if
        !       !end do
        !     end do
        !   end do
        !   first = .false.
        ! end if
        ! do j=1,nlat
        !   do i=1,nlon
        !      rain(i,j) = rani(mapi(i,j),mapj(i,j))
        !      rsrm(i,j)=srccode
        !   end do
        ! end do
!        call nearest2d(rani, lati, loni, rain, lat, lon, 0, 1, nreg2, nreg1, nlat, nlon)
        if ( nnflag ) then
        !!call closest2d(lati,loni,lat,lon,nreg1,nreg2,mapi(jj,ii),mapj(jj,ii),.true.)
!        rani = rani(:,nreg2:1:-1)
        do ii = 1, nlat
          do jj = 1, nlon
            call closest2d(lati,loni,lat(jj,ii),lon(jj,ii),nreg1,nreg2,mapi(jj,ii),mapj(jj,ii),.true.)

!            mapj(jj,ii) = jmin
!            mapi(jj,ii) = imin
          enddo
        enddo
        nnflag = .false.
        end if
        do i=1,nlat
          do j=1,nlon
!             rain(i,j) = rani(mapi(nlon-i+1,nlat-j+1),mapj(nlon-i+1,nlat-j+1))
             rain(j,i) = rani(mapj(j,i),mapi(j,i))

             rsrm(j,i)=srccode
          end do
        end do
!         CALL CPU_TIME(tallintf)
!         tallint = tallintf - tallinti
        return
      end if
    end if
!   write (logfile,'(18x,a)') 'Cressman Algorithm to build rain field.'
!    CALL CPU_TIME(tcressi)
    call ncressman(pi,la,lo,n,rchym,radius,rain,ca,wrk2,lat,lon, &
       canc,nlon,nlat)
!    CALL CPU_TIME(tcressf)
!    tcresst = tcresst + tcressf - tcressi
!    call cpu_time(tautoi)
    if (useca.and.piove) then
!     write (logfile,'(18x,a,i3,a)') 'Calling CA Algorithm to build '//
!   2      'rain field ',ncacyc,' cycles.'
      do i=1,ncacyc
         call d2cellcycle(rain,ca,wrk2,nlon,nlat,0.9)
      end do
    else if (useca) then
      write (logfile,'(18x,a,i3,a)') 'No relevant rain on the domain, '// &
           'CA Algorithm not used'
    endif
!    call cpu_time(tautof)
!    tautot = tautot + tautof - tautoi
    do j=1,nlat
      do i=1,nlon
        if (canc(i,j).eq.1.and.ca(i,j).eq.8) then
          ca(i,j)=0
          rsrm(i,j)=srccode
        end if
      end do
    end do
!    CALL CPU_TIME(tallintf)
!    tallint = tallintf - tallinti
!    print*,tallint
  end subroutine buildrainfield

  subroutine snowcover
    implicit none
    integer :: ora,giorno,mese,anno
    logical :: first
    integer :: i , j , lastday
    data lastday/ - 1/
    data first/.true./
    save lastday,first
    call gmafromindex(time,ora,giorno,mese,anno)
!    if ( first ) then
!      modis = 0.0
!      first = .false.
!    end if
    if ( mchym(24)/=0 ) write (logfile,'(16x,a)') '> Estimating snow cover'
    if ( abs(mchym(24))==21 ) then
      if ( giorno/=lastday ) then
        mchym(24) = 21
        call acqwamodis(giorno,mese,anno,modis)
        lastday = giorno
      end if
      do i = 1 , nlon
        do j = 1 , nlat
          if ( nint(modis(i,j))==1 ) rain(i,j) = 0.0
        end do
      end do
    else if ( abs(mchym(24))==1 ) then
      if ( giorno/=lastday ) then
        call modisitaly(giorno,mese,anno,modis,nlon,nlat)
        lastday = giorno
      end if
      do i = 1 , nlon
        do j = 1 , nlat
          if ( nint(modis(i,j))==1 ) rain(i,j) = 0.0
        end do
      end do
    end if
  end subroutine snowcover

  subroutine snowacc
    implicit none
    integer :: i , j
    do i = 1 , nlon
      do j = 1 , nlat
        if ( temp(i,j)<0.0 ) then
          snow(i,j) = snow(i,j) + rain(i,j)
          rain(i,j) = 0.0
        end if
      end do
    end do
  end subroutine snowacc

  subroutine evapotranspiration
    implicit none
    integer :: i , j , l2
    write (logfile,'(16x,a)') '> Calculating Evaporatranspiration term'
    do i = 2 , nlon - 1
      do j = 2 , nlat - 1
!        l2 = luse(i,j)
!        if ( l2/=mare .and. l2/=lago .and. l2>=1 .and. l2<=lntypes ) then
        if (luse(i,j).ne.mare.and.luse(i,j).ne.lago) then
          gh2o(i,j) = gh2o(i,j) - cvmm2m3(evap(i,j),i,j)
          if ( gh2o(i,j)<0.0 ) gh2o(i,j) = 0.0
        end if
      end do
    end do
  end subroutine evapotranspiration

  subroutine melting
    implicit none
    real :: alb , s , sc , sigh , sigl , sigm , tr
    integer :: i , j
    integer ora,giorno,mese,anno,ltime
    data alb/0.80/   ! Albedo
    data sc/1368.0/  ! Solar Constant
    data sigh , sigm , sigl/0.2 , 0.2 , 0.2/  ! Cloud cover fraction for
                                              ! high, middle and low atmosph.
    save alb,sc,sigh,sigm,sigl
    write (logfile,'(16x,a)') '> Calculating melting term'
    if (time.eq.2000010100) then
      ltime=1999123123
    else
      ltime=decreasemm5index(time)
    end if
    call gmafromindex(ltime,ora,giorno,mese,anno)
    wrk1 = 0.0
    do i = 2 , nlon - 1
      do j = 2 , nlat - 1
        if ( luse(i,j)/=mare .and. snowcv(i,j)==1 ) then
          if ( temp(i,j)>1.0 ) then
            wrk1(i,j) = cpar(4)*temp(i,j)
            s = sinphi(lat(i,j),lon(i,j),30,ora,giorno,mese)
            if ( s>0.0 ) then
              tr = (0.6+0.2*s)*(1.0-0.4*sigh)*(1.0-0.7*sigm)*(1.0-0.4*sigl)
              wrk1(i,j) = wrk1(i,j) + cpar(5)*(1.-alb)*sc*tr*s
            end if
          end if
        end if
      end do
    end do
    do i = 1 , nlon
      do j = 1 , nlat
        if ( wrk1(i,j)>0.0 ) then
          snow(i,j) = snow(i,j) - wrk1(i,j)
          if ( snow(i,j)<0.0 ) snow(i,j) = 0.0
        end if
      end do
    end do
    rain = rain + wrk1
  end subroutine melting

  integer function snowcv(i,j)
    implicit none
    integer :: i , j
    intent (in) i , j
    real , dimension(12) :: snowl
    data snowl/1800. , 2000. , 2000. , 2100. , 2200. , 2500. , 2500. , &
               2200. , 2000. , 2000. , 2000. , 1900./
    save snowl
    snowcv = 0
    if ( mchym(24)==1 .or. mchym(24)==21 ) then
      if ( nint(modis(i,j))==1 ) snowcv = 1
    else if ( snow(i,j)>0.1 .or. luse(i,j)==12 ) then
      snowcv = 1
    end if
  end function snowcv

  subroutine groundwater
    implicit none
    real :: actinf , relh , xinfl
    integer :: i , j , l2
    write (logfile,'(16x,a)') '> Calculating percolation/filtration (2) term'
    do i = 1 , nlon
      do j = 1 , nlat
!        l2 = luse(i,j)
!        if ( l2/=mare .and. l2/=lago .and. l2>=1 .and. l2<=lntypes ) then
        if (luse(i,j).ne.mare.and.luse(i,j).ne.lago) then
          xinfl = cvmm2m3(infi(luse(i,j)),i,j)
          if ( xinfl>100 ) then
            gh2o(i,j) = gh2o(i,j) - cvmm2m3(perc(luse(i,j)),i,j)
!           deepw(i,j)=deepw(i,j)+perc(luse(i,j))
            if ( gh2o(i,j)<0.0 ) gh2o(i,j) = 0.0
            relh = gh2o(i,j)/xinfl
            if ( relh<0.999 ) then
              actinf = cvm32mm(h2o(i,j),i,j)*sigmoide(1-relh,0.5,0.165)
              actinf = cvmm2m3(actinf,i,j)
              gh2o(i,j) = gh2o(i,j) + actinf
              h2o(i,j) = h2o(i,j) - actinf
            end if
          end if
        end if
      end do
    end do
  end subroutine groundwater

  subroutine returnflow
    implicit none
    integer :: i , j , ora , giorno , mese , anno, ltime
    if (time.eq.2000010100) then
       ltime=1999123123
    else
       ltime=decreasemm5index(time)
    endif
    call gmafromindex(ltime,ora,giorno,mese,anno)
    deepw(1:nlon,1:nlat) = deepw(1:nlon,1:nlat) * &
                          ((cpar(7)*24.0-1.0)/(cpar(7)*24.0))
    if ( ora==0 ) then
      wk(1:nlon,1:nlat) = 1000.0*deepw(1:nlon,1:nlat)*area(1:nlon,1:nlat)
      call rollingstones2(wk,fmap,ddeepw,nlon,nlat)
      do i = 2 , nlon - 1
        do j = 2, nlat - 1
          ddeepw(i,j) = .001*(ddeepw(i,j)/max(drai(i,j),0.1))
        end do
      end do
    end if
    do i = 2 , nlon - 1
      do j = 2 , nlat - 1
        if ( drai(i,j)>=cpar(9) ) then
          rain(i,j) = rain(i,j) + cpar(1)*ddeepw(i,j)*drai(i,j)
        end if
      end do
    end do
  end subroutine returnflow

  subroutine tempfromera40
    implicit none
    real, dimension(360,181) :: tmpt
    real :: w1,w2,w3,w4,wt
    character(len=20) :: str
    integer :: i,j,mese,lun,ii,jj
    call gmafromindex(time,i,j,mese,lun)
    call monthofyear(mese,-1,str)
    write (logfile,'(16x,a)') 'Acquiring ERA40 monthly average '// &
          'temperature for '//str
    lun=-1
    call openmuseofiles(lun,'worldtemp.dat',0)
    do i=1,mese
      read(lun) tmpt
    enddo
    close(lun)
    call kelvin2centi(tmpt,360,181,1)
    do i=1,nlon
      do j=1,nlat
        ii=int(lon(i,j))+181
        jj=int(lat(i,j))+91
        w1=1.0/(distance(lat(i,j),lon(i,j),-91.0+jj+0,-181.0+ii+0))**2
        w2=1.0/(distance(lat(i,j),lon(i,j),-91.0+jj+0,-181.0+ii+1))**2
        w3=1.0/(distance(lat(i,j),lon(i,j),-91.0+jj+1,-181.0+ii+0))**2
        w4=1.0/(distance(lat(i,j),lon(i,j),-91.0+jj+1,-181.0+ii+1))**2
        wt=w1+w2+w3+w4
        temp(i,j)=(tmpt(ii,jj)*w1+tmpt(ii+1,jj)*w2+ &
                   tmpt(ii,jj+1)*w3+tmpt(ii+1,jj+1)*w4)/wt
      enddo
    enddo
  end subroutine tempfromera40

  subroutine tempfrommuseo(iflg)
    implicit none
    integer, parameter :: ixm = 70
    integer, parameter :: jxm = 70
    real, dimension(ixm,jxm,25) :: p
    real, dimension(ixm,jxm,24) :: t
    real, dimension(ixm,jxm) :: xlat, xlon
    integer :: r2read
    integer, save :: oldgiorno,oldmese,oldanno,ix,jx,lun,rec,tslice
    integer :: ora,giorno,mese,anno,lca(nlon,nlat),iflg,i,j,ih,ii,jj,im,jm
    integer :: lkx,lns,n,np,len_trim
    character(len=30) :: ifile
    character(len=30), save :: ofile
    character(len=60) :: now
    data ofile,lun /'notyetread',-1/
    call gmafromindex(time,ora,giorno,mese,anno)
    call datafromidx(time,now)
    write(ifile,'(a,i4,i1,a)') 'mm5t2m',anno,mchym(12),'.dat'
    tslice=ora+1
    if (ifile.ne.ofile) then
       if (lun.gt.0) close(lun)
       call openmuseofiles(lun,ifile,0) ; ofile=ifile ; rec=0
    endif
    r2read=indexofyear(giorno,mese,anno)
    if (r2read.ne.rec+1.and.r2read.ne.rec) then
       do i=rec+1,r2read-1 ; read(lun) ; enddo
       rec=r2read-1
    endif
    if (r2read.ne.rec) then
       write (logfile,'(16x,a)') '> Acquiring museo temperatures for ' &
           //now(1:len_trim(now))
       read(lun) ix,jx,lkx,lns, &
        ((xlat(ii,jj),ii=1,ix),jj=1,jx),((xlon(ii,jj),ii=1,ix),jj=1,jx), &
        (((t(ii,jj,n),ii=1,ix),jj=1,jx),n=1,24)
       if (ix.le.1.or.jx.le.1) then
          write (logfile,'(20x,a)') 'MuSEO data not available for this day'
          write (logfile,'(20x,a)') 'Temperature field unchanged.'
       else if (iflg.eq.2) then
          do im=1,ix ; do jm=1,jx
             do ih=2,24
                t(im,jm,1)=t(im,jm,ih)+t(im,jm,1)
             enddo
             t(im,jm,1)=t(im,jm,1)/24.
          enddo ; enddo
          tslice=1
       endif
       rec=r2read
    endif
    lca=8 ; temp=5
    do im=1,ix ; do jm=1,jx
       i=nint((xlon(im,jm)-rchym(1))/rchym(8))+1
       j=nint((xlat(im,jm)-rchym(2))/rchym(9))+1
       if (inside(i,1,nlon).and.inside(j,1,nlat)) then
          temp(i,j)=t(im,jm,tslice)
          lca(i,j)=0
       endif
    enddo ; enddo
    np=0.8*float(3**(4-mchym(12)))/(0.001*rchym(6))
    do ii=1,nlon ; do jj=1,nlat
       if (lca(ii,jj).eq.0) then
          do i=ii-np,ii+np  ; do j=jj-np,jj+np
             if(inside(i,1,nlon).and.inside(j,1,nlat)) temp(i,j)=temp(ii,jj)
          enddo ; enddo
       endif
    enddo ; enddo
    do i=1,nlon ; do j=1,nlat
       if(luse(i,j).eq.mare) lca(i,j)=0
    enddo ; enddo
    do i=1,30 ; call d2cellcycle(temp,lca,wrk1,nlon,nlat,0.9) ; enddo
    return
  end subroutine tempfrommuseo

  subroutine raingaugesval(pi2val,pi,la,lo,n)
    implicit none
    integer :: n
    real , dimension(:) :: la , lo , pi , pi2val
    intent (in) pi2val
    intent (inout) la , lo , n , pi
    integer :: i , ln
    ln = n
    n = 0
    do i = 1 , ln
      if ( pi(i)<100.0 ) then
        n = n + 1
        pi(n) = pi2val(i)
        la(n) = la(i)
        lo(n) = lo(i)
      end if
    end do
  end subroutine raingaugesval

  subroutine erainterim(ora,giorno,mese,anno,pi,la,lo,n)
    implicit none
    integer , parameter :: nlon = 480 , nlat = 241
    integer :: anno , giorno , mese , n , ora
    real , dimension(:) :: la , lo , pi
    intent (out) la , lo , pi
    intent (inout) n
    integer :: i , j , lun
    real , dimension(nlon,nlat) :: rain
    data lun/ - 1/
    call eradata(lun,ora,giorno,mese,anno,rain)
    n = 0
    do i = 1 , nlon
      do j = 1 , nlat
        n = n + 1
        pi(n) = rain(i,j)
        la(n) = eralat(i,j)
        lo(n) = eralon(i,j)
      end do
    end do
  end subroutine erainterim

  subroutine evapot(temp,lat,lon,quota,j,n,eto)
    implicit none
    integer , parameter :: ord = 1 , nord = ord + 1 , &
                    nr = 12 , llat = 121 , nlon = 360
    real :: eto , lat , lon , n , quota , temp
    integer :: j
    intent (in) lat , lon , n
    intent (out) eto
    intent (inout) quota , temp
    real , dimension(nlon,llat,nr) :: ba
    real , dimension(2) :: coeff , e , wq
    real :: den , ilat , jlon , num , wta , x
    real , dimension(4) :: eh , h
    logical :: first
    integer , dimension(2) :: ih , iq , klat , klon
    real , dimension(3) :: xabc
    integer :: il , jl , k , l , lun
    data first/.true./   ! ordine eq., n=numero di coeff.
    save first
    if ( first ) then
      lun = -1
      call openmuseofiles(lun,'BAtot.dat',0)
      read (lun) ba
      close (lun)
      first = .false.
    end if
    if ( lat>60.0 .or. lat<-60.0 ) then
      eto = 0.0
      return
    end if
    if ( temp==0.0 ) temp = 0.0001
    if ( quota<=0.0 ) quota = 0.01
    call id_quota(quota,iq,ih)
    ilat = 61 + lat        !indice latitudine per BA
    jlon = 181 + lon       !indice longitudine per BA
    klat(1) = int(ilat)
    klat(2) = int(ilat) + 1
    klon(1) = int(jlon)
    klon(2) = int(jlon) + 1
    do k = 1 , 2
      xabc(1) = abc(iq(k),1)
      xabc(2) = abc(iq(k),2)
      xabc(3) = abc(iq(k),3)
      wta = polyval(temp,3,xabc)       ! Wta,eq: Wta=A*T^2+B*T+C
      x = temp*n*wta ! polyval per ETo= b*(N*T*Wta)+a
      num = 0
      den = 0
      l = 0
      do jl = 1 , 2
        do il = 1 , 2
          l = l + 1
          coeff(1) = ba(klon(jl),klat(il),iq(k))
          coeff(2) = ba(klon(jl),klat(il),iq(k)+6)
          eh(l) = polyval(x,nord,coeff)
          h(l) = 1./distance(ilat,jlon,float(klat(il)),float(klon(jl)))
          num = num + eh(l)*h(l)
          den = den + h(l)
        end do
      end do
      e(k) = num/den
      if ( abs(quota-ih(k))>0.00001 ) then
        wq(k) = 1/abs(quota-ih(k))
      else
        wq(k) = 1
      end if
    end do           !k
    eto = (e(1)*wq(1)+e(2)*wq(2))/(wq(1)+wq(2))     ! Evapotrasp. di Riferimento
  end subroutine evapot

  subroutine id_quota(quota,iq,ih)
    implicit none
    real :: quota
    integer , dimension(2) :: ih , iq
    intent (in) quota
    intent (out) ih , iq
    integer :: i , k
    integer , dimension(6) :: valq
    data (valq(i),i = 1,6)/0 , 500 , 1000 , 2000 , 3000 , 4000/

    k = 0
    do while ( k<5 )
      k = k + 1
      if ( quota>=valq(k) .and. quota<valq(k+1) ) then
        iq(1) = k
        iq(2) = k + 1
        ih(1) = valq(k)
        ih(2) = valq(k+1)
        k = 6
      else if ( quota>=4000 ) then
        iq(1) = 6
        iq(2) = 6
        ih(1) = 4000
        ih(2) = 4000
      end if
    end do
  end subroutine id_quota

  real function polyval(x,nord,c)
    implicit none
    integer :: nord
    real :: x
    real , dimension(nord) :: c
    intent (in) c , nord , x
    integer :: i
    real :: y
    y = 0.
    do i = 1 , nord
      y = y + c(i)*(x**(nord-i))
    end do
    polyval = y
  end function polyval

  real function daylenght(l,j)           !L=latitudine,J=giorno giuliano
    implicit none
    integer :: j
    real :: l
    intent (in) j , l
    real :: p , pi , ws
    pi = 4.*atan(1.)
    p = asin(.39795*cos(.2163108+2*atan(0.9671396*tan(.00860*(j-186)))))
    ws = (sin(0.8333*pi/180)+sin(l*pi/180)*sin(p))/(cos(l*pi/180)*cos(p))
    daylenght = 24 - (24/pi)*acos(ws)
  end function daylenght

  subroutine addsource(lsource,actsource)
    implicit none
    character(len=*) :: lsource,actsource
    integer len_trim
    if (len_trim(actsource).eq.0) then
      actsource=lsource
    else
      actsource=actsource(1:len_trim(actsource))//','//lsource
    endif
    return
  end subroutine addsource

  subroutine ncressman(val,clat,clon,nc,rchym,rad,mat,iwk,w, &
                       lat,lon,canc,nlon,nlat)
    implicit none
    integer :: nc , nlat , nlon
    integer , dimension(nlon,nlat) :: canc , iwk
    real , dimension(nc) :: clat , clon , val
    real , dimension(nlon,nlat) :: lat , lon , mat , w
    real :: rad
    real , dimension(:) :: rchym
    intent (in) iwk , nc , nlat , nlon , rad , rchym , val
    intent (out) canc
    intent (inout) mat , w
    real :: a , a2 , avr , dij , r , r2 , slat , slon , weight
    integer :: i , i1 , i2 , j , j1 , j2 , n , nij

    avr = rchym(6)
    dij = rchym(5)
    slon = rchym(1)
    slat = rchym(2)
    w = 0.0
    do j = 1 , nlat
      do i = 1 , nlon
        if ( iwk(i,j)>0 ) mat(i,j) = 0.0
      end do
    end do
    nij = nint((rad*1000.)/avr)
    a = rad*1000.
    a2 = a*a
    do n = 1 , nc
      if (val(n)>-1.0 ) then
      i = nint((clon(n)-slon)/dij) + 1
      j = nint((clat(n)-slat)/dij) + 1
      i1 = i - nij
      i2 = i + nij
      if ( i1<1 ) i1 = 1
      if ( i2>nlon ) i2 = nlon
      if ( i2<1 ) i2 = 0
      j1 = j - nij
      if ( j1<1 ) j1 = 1
      j2 = j + nij
      if ( j2>nlat ) j2 = nlat
      if ( j2<1 ) j2 = 0
      do j = j1 , j2
        do i = i1 , i2
          if ( iwk(i,j)>0 ) then
            r = distance(clat(n),clon(n),lat(i,j),lon(i,j))
            if ( r<a  .and. val(n)>-1.0 ) then
              r2 = r*r
              weight = (1-r2/a2)/(1+r2/a2)
              mat(i,j) = mat(i,j) + val(n)*weight
              w(i,j) = w(i,j) + weight
            end if
          end if
        end do
      end do
      end if
    end do
    do j = 1 , nlat
      do i = 1 , nlon
        if ( iwk(i,j)>0 .and. abs(w(i,j))>1E-10 ) then
          mat(i,j) = mat(i,j)/w(i,j)
          canc(i,j) = 1
        else
          canc(i,j) = 0
        end if
      end do
    end do
  end subroutine

  subroutine acqwamodis(giorno,mese,anno,mat)
    implicit none
    integer :: anno , giorno , mese
    real , dimension(mchym(2),mchym(3)) :: mat
    intent (out) mat
    real , dimension(580,878) :: dati
    real  :: xx,yy,xla,xlo
    real , save :: dlat , dlon , dx , dy , slat , slon , x1 ,y1
    logical :: first , nolonger
    integer :: i , il , j , jl , logfile , lun , nx , ny , irec
    integer , save :: oldrec , oldyear
    data nolonger/.false./
    data first/.true./
    data oldyear/ - 1/
    save nolonger,first,logfile,lun , nx , ny
    if ( first ) then
      call getlun(lun)
      call mvgetiflags(70,logfile)
      if ( logfile<=0 .or. logfile>=100 ) logfile = 6
      if ( schym(10)(1:7)=='acqwapo' ) then
        nx = 580
        ny = 848
        slon = 6.1805
        slat = 43.4679
        dlat = 0.0042
        dlon = 0.0058
        call openmuseofiles(lun,'0000/scPo00-09.dat',0)
        call latlon2utm(45.0,8.0,xx,yy)
        x1 = 271930.435
        y1 = 4816643.679
        dy = 463.3127
        dx = 463.3127
      else if ( schym(10)(1:10)=='acqwarhone' ) then
        nx = 410
        ny = 366
        slon = 6.5800
        slat = 45.2016
        dlat = 0.0050
        dlon = 0.0059
        call openmuseofiles(lun,'0000/scRhone00-09.dat',0)
        x1 = 309940.197
        y1 = 5038625.126
        dy = 463.3127
        dx = 463.3127
        call latlon2utm(45.0,8.0,xx,yy)
      else
        write (6,'(12x,a)') 'Unknown basin: '//trim(schym(10)) &
                           //', ACQWA modis module no longer used.'
        mchym(24) = 0
        nolonger = .true.
      end if
      mat = 0.0
      first = .false.
    end if
    if ( nolonger ) return
    irec = 0
    if ( anno==2000 ) then
      irec = (julianday(giorno,mese,anno)-46)/8
    else if ( anno>=2001 .and. anno<=2009 ) then
      irec = julianday(giorno,mese,anno)/8 + (anno-2000)*46 - 6
    end if
    if ( irec<1 .or. irec>453 ) then
      mchym(24) = -21
      irec = mm5index(0,giorno,mese,anno)
      if ( oldyear/=anno ) write (6,'(15x,a,i8)')  &
                     'ACQWA Modis data not available for ' , irec
      oldyear = anno
      mat = 0
      return
    else if ( irec==oldrec ) then
      return
    else if ( irec==oldrec+1 ) then
      write (logfile,'(18x,a,i5,a)') &
        'Read record ' , irec , ' of ACQWA modis db'
      read (lun) dati(1:nx,1:ny)
      oldrec = irec
    else
      write (logfile,'(18x,a,i5,a)') 'Skip ' , irec - 1 , &
                                  ' records of ACQWA modis db'
      rewind (lun)
      do i = 1 , irec - 1
        read (lun)
      end do
      write (logfile,'(18x,a,i5,a)') &
        'Read record ' , irec , ' of ACQWA modis db'
      read (lun) dati(1:nx,1:ny)
      oldrec = irec
    end if
    mat = 0.0
    do jl = 1 , ny
      do il = 1 , nx
        xx = (il-1)*dx + x1
        yy = (jl-1)*dy + y1
        call utm2latlon(xx,yy,xla,xlo)
        call locateij(xla,xlo,rchym(2),rchym(1),rchym(9),rchym(8), &
                      mchym(3),mchym(2),i,j)
        if ( i>=1 .and. i<=mchym(2) .and. j>=1 .and. j<=mchym(3) ) then
          if ( nint(dati(il,jl))==200 .or. nint(dati(il,jl))==100 ) mat(i,j) = 1
        end if
        i = i + 1
        if ( i>=1 .and. i<=mchym(2) .and. j>=1 .and. j<=mchym(3) ) then
          if ( nint(dati(il,jl))==200 .or. nint(dati(il,jl))==100 ) mat(i,j) = 1
        end if
      end do
    end do
    return
  end subroutine acqwamodis

  subroutine modisitaly(giorno,mese,anno,modis,nlon,nlat)
    implicit none
    integer , parameter :: n2 = 2529 , n1 = 3458
    real , dimension(n1,n2) :: xmod
    integer :: anno , giorno , mese , nlat , nlon
    real , dimension(nlon,nlat) :: modis
    intent (out) modis
    real :: dla , dlat , dlo , dlon , lat , lon , r1 , r2 , &
            r3 , r4 , slat , slon
    integer :: i , j , k , l , logfile , m , n3 , n4 , nt , rec
    logical :: start
    data start/.true./
    save start,slat,slon,dlat,dlon,logfile
    if (anno<mvlibmagicnum(6).or.anno>mvlibmagicnum(7)) then
      return
    endif
    if ( start ) then
      call chymgetinfo(1,logfile,slon,' ')
      call chymgetinfo(2,k,slon,' ')
      call chymgetinfo(3,k,slat,' ')
      call chymgetinfo(4,k,dlon,' ')
      call chymgetinfo(5,k,dlat,' ')
      start = .false.
    end if
    rec = julianday(giorno,mese,anno)
    call snowfrommodis(rec,anno,xmod,n1,n2,n3,n4,r1,r2,r3,r4)
    dla = (r4-r3)/(n4-1)
    dlo = (r2-r1)/(n3-1)
    nt = 0
    do k = 1 , n3
      do l = 1 , n4
        lon = r1 + (k-1)*dlo
        lat = r3 + (l-1)*dla
        call locateij(lat,lon,slat,slon,dlat,dlon,nlat,nlon,i,j)
        if ( i>0 ) then
          m = nint(xmod(k,l))
          if ( m==200 ) then
            modis(i,j) = 1
            nt = nt + 1
          else if ( m==25 .or. m==37 ) then
            modis(i,j) = 0
            nt = nt + 1
          end if
        end if
      end do
    end do
    write (logfile,'(18x,i6,a)') nt , ' grid points defined by ModIS module.'
  end subroutine modisitaly

  subroutine chymgetinfo(flag,intv,rea,str)
    implicit none
    integer :: flag , intv
    real :: rea
    character(len=*) :: str
    intent (in) flag
    intent (out) intv , rea
    if ( flag==1 ) then
      intv = logfile
    else if ( flag==2 ) then
      rea = rchym(1)
    else if ( flag==3 ) then
      rea = rchym(2)
    else if ( flag==4 ) then
      rea = rchym(8)
    else if ( flag==5 ) then
      rea = rchym(9)
    end if
  end subroutine chymgetinfo

  real function cvmm2m3(x,i,j)
    implicit none
    integer :: i , j
    real :: x
    intent (in) i , j , x
    cvmm2m3 = area(i,j)*1.0E+03*x
  end function cvmm2m3

  real function cvm32mm(x,i,j)
    implicit none
    integer :: i , j
    real :: x
    intent (in) i , j , x
    cvm32mm = x*1.0E-03/area(i,j)
  end function cvm32mm

        subroutine mm5data(hour,day,month,year,pi,la,lo,n)
        implicit none
        integer hour,day,month,year ; real pi(1),la(1),lo(1)
        integer ora,giorno,mese,anno,n,len_trim,im,jm,islice,iphase,istart
        integer ix,jx,kx,nn
        integer u,uu,ih
        integer ixm,jxm ; parameter (ixm=100,jxm=100)
        real rmm5(ixm,jxm),smm5(ixm,jxm),mmlat(ixm,jxm),mmlon(ixm,jxm)
        data islice,iphase,istart /0,0,0/
        save islice,iphase,istart,rmm5,smm5,mmlat,mmlon
        integer ildate,idatemm5 ; save ildate,idatemm5
        uu=logfile ; n=0 ; nn=len_trim(schym(2)) ; u=lu(4)
!------------------------------------------------------------------------------
!  Following lines are needed to use more than 1 mm5 files. - To check
!
!       if (istart.ne.0.and.schym(2)(1:nn).ne.oldfile(1:len_trim(oldfile)) then
!          write (6,'(16x,a)')'> mm5data: using now a new MM5 output file.'
!          close(u) ; oldfile=schym(2)
!          open(u,status='old',form='unformatted',file=schym(2),iostat=istat)
!          if (istat.ne.0) then
!             write (6,'(15x,a)')'Cannot open '//schym(2)(1:nn)//' - NOT USED'
!             iphase=-1
!             return
!          endif
!          islice=0 ; iphase=0 ; istart=0
!       endif
!------------------------------------------------------------------------------
        if (iphase.lt.0) return
!------------------------------------------------------------------------------
!    Following lines are needed because of difference between local and UGT time
        ora=hour ; giorno=day ; mese=month ; anno=year
        call localtime2ugt(ora,giorno,mese,anno)
        if (anno.ge.2000) then
           ildate=mm5index(ora,giorno,mese,anno-2000)
        else
           ildate=mm5index(ora,giorno,mese,anno-1900)
        endif
!------------------------------------------------------------------------------
        if (istart.eq.0) then
           istart=1
           call mm5readfld (u,'terrain',islice,mmlat,ixm,jxm,1,ix,jx,kx)
           if (ix.le.1.or.jx.le.1.or.islice.lt.0) then
              call chymerror(12,istart,0.0,schym(2))
              iphase=-1
              return
           endif
           idatemm5=mm5mif(1,mm5mif(1,1))
           call mm5readfld (u,'latitcrs',islice,mmlat,ixm,jxm,1,ix,jx,kx)
           if (ix.le.1.or.jx.le.1.or.islice.lt.0) then
              call chymerror(12,istart,1.0,schym(2))
              iphase=-1
              return
           endif
           call mm5readfld (u,'longicrs',islice,mmlon,ixm,jxm,1,ix,jx,kx)
           if (ix.le.1.or.jx.le.1.or.islice.lt.0) then
              call chymerror(12,istart,2.0,schym(2))
              iphase=-1
              return
           endif
           islice=0             ! Next mm5readfld call will rewind input file
        endif
        if (iphase.eq.0) then
           if (idatemm5.gt.ildate) then
!             write (uu,'(18x,a)') 'First slice is in the future, not yet used.'
              islice=0
              return
           else if (idatemm5.eq.ildate) then
              call mm5rain(u,mm5mif(1,1),islice,rmm5,ixm,jxm,1,ix,jx,kx)
              call mm5rain(u,mm5mif(1,1),islice,smm5,ixm,jxm,1,ix,jx,kx)
              if (ix.le.1.or.jx.le.1.or.islice.le.0) then
                 call chymerror(12,istart,3.0,schym(2))
                 iphase=-1
                 return
              endif
              if (increasemm5index(ildate).ne.mm5mif(1,mm5mif(1,1))) then
                call chymerror(14,increasemm5index(ildate), &
                        float(mm5mif(1,mm5mif(1,1))),schym(2))
                iphase=-1
                return
              endif
              iphase=1
           else if (idatemm5.lt.ildate) then
              do while (islice.ge.0.and.mm5mif(1,mm5mif(1,1)).ne.ildate)
                 call mm5rain(u,mm5mif(1,1),islice,rmm5,ixm,jxm,1,ix,jx,kx)
              end do
              if (ix.le.1.or.jx.le.1.or.islice.le.0) then
                 call chymerror(12,istart,4.0,schym(2))
                 iphase=-1
                 return
              endif
              call mm5rain(u,mm5mif(1,1),islice,smm5,ixm,jxm,1,ix,jx,kx)
              if (ix.le.1.or.jx.le.1.or.islice.le.0) then
                 call chymerror(12,istart,5.0,schym(2))
                 iphase=-1
                 return
              endif
              if (increasemm5index(ildate).ne.mm5mif(1,mm5mif(1,1))) then
                call chymerror(14,increasemm5index(ildate), &
                        float(mm5mif(1,mm5mif(1,1))),schym(2))
                iphase=-1
                return
              endif
              iphase=1
           endif
        else if (iphase.eq.1) then
           call mm5rain(u,mm5mif(1,1),islice,smm5,ixm,jxm,1,ix,jx,kx)
           if (ix.le.1.or.jx.le.1.or.islice.le.0) then
              write(uu,'(18x,a)') 'mm5data module, end of file, no longer used.'
              iphase=-1
              return
           endif
           if (increasemm5index(ildate).ne.mm5mif(1,mm5mif(1,1))) then
              call chymerror(14,increasemm5index(ildate), &
                        float(mm5mif(1,mm5mif(1,1))),schym(2))
           endif
        else if (iphase.lt.0) then
           return
        endif
        do im=1,ix-1 ; do jm=1,jx-1
           n=n+1
           pi(n)=(smm5(im,jm)-rmm5(im,jm))*10.0
           la(n)=mmlat(im,jm)
           lo(n)=mmlon(im,jm)
        enddo ; enddo
        rmm5=smm5
        return
        end subroutine mm5data

  subroutine mm5rain(u,mif,islice,rmm5,ixm,jxm,kxm,ix,jx,kx)
    implicit none
    integer , parameter :: lkxm = 150
    integer :: islice , ix , ixm , jx , jxm , kx , kxm , mif , u
    real , dimension(ixm,jxm,kxm) :: rmm5
    intent (in) kxm , mif
    intent (inout) rmm5
    integer :: i , j
    real , dimension(lkxm,lkxm) :: smm5
    if ( mif==7 ) then
      call mm5readfld(u,'rain tot',islice,rmm5,ixm,jxm,1,ix,jx,kx)
    else if ( mif==6 ) then
      call mm5readfld(u,'rain con',islice,rmm5,ixm,jxm,1,ix,jx,kx)
      call mm5readfld(u,'rain non',islice,smm5,lkxm,lkxm,1,ix,jx,kx)
      do i = 1 , ixm
        do j = 1 , jxm
          rmm5(i,j,1) = rmm5(i,j,1) + smm5(i,j)
        end do
      end do
    else if ( mif==11 .or. mif==8 ) then
      call mm5readfld(u,'rain con',islice,rmm5,ixm,jxm,1,ix,jx,kx)
      call mm5readfld(u,'rain non',islice,smm5,lkxm,lkxm,1,ix,jx,kx)
      do i = 1 , ixm
        do j = 1 , jxm
          rmm5(i,j,1) = rmm5(i,j,1) + smm5(i,j)
        end do
      end do
    else
      write (6,'(2x,a,i2)') 'leggilapioggia flux error, mif was ' , mif
      stop ' '
    end if
  end subroutine mm5rain
  subroutine read_stations(hour,day,month,year,pi,lat,lon,nstat)
      implicit none
      integer :: i,j
      integer , intent(in) :: day,month,year,hour
      integer :: n,nhour,nstat
      integer  :: adate,shour,sdate
      real :: cs
      character(len=50) :: a
      character(len=5) :: d
      character(len=256) :: lfile
!      real , allocatable , dimension(:) :: pi,pi1,lat,lon
      integer, parameter :: mxnicorw=1100000
      real, dimension(mxnicorw) :: lat,lon,pi
      adate = year*10000+month*100+day
      lfile = "doc/stations.txt"  !Name of the file from which you want to read precipitation data from station
      n = 1                    !Number of of records
      open(62,file=trim(lfile),status='old',err=10,action='read')
      read(62,*) a,nstat
      read(62,*) a
      pi=9999; lat=-999; lon=-999
      do j=1,nstat
        read(62,*) a, n, lon(j), lat(j)
        do i=1,n
          read(62,*) sdate, shour,cs
          if (adate == sdate .and. hour == shour) then
            pi(j) = cs
          end if
        enddo
      enddo
      close(62)
      return
10    print*,"error in opening file"
  end subroutine read_stations

!old  subroutine updateBOUND(var)
!old    real,dimension(:,:),intent(inout) :: var
!old
!old    if (myid.eq.0) then
!old!     call MPI_SEND(var(1,(nlat/nproc)+2),nlon,MPI_REAL,rank+1,     &
!old!       rank,comm,mpistatus,ierr)
!old     call MPI_SEND(var(1,(nlat/nproc)+2),nlon,MPI_REAL,rank+1,     &
!old       rank,mycomm,mpierr)
!old     call MPI_RECV(var(1,(nlat/nproc)+3),nlon,MPI_REAL,rank+1,     &
!old       rank+nproc+1,mycomm,mpistatus1,mpierr)
!old    else if(rank.eq.nproc-1) then
!old     call MPI_RECV(var(1,2),nlon,MPI_REAL,rank-1,rank-1,mycomm,     &
!old       mpistatus,mpierr)
!old!     call MPI_SEND(var(1,3),nlon,MPI_REAL,rank-1,rank+nproc,comm, &
!old!       mpistatus1,ierr)
!old     call MPI_SEND(var(1,3),nlon,MPI_REAL,rank-1,rank+nproc,mycomm, &
!old       mpierr)
!old    else
!old     call MPI_SENDRECV(var(1,(nlat/nproc)+2),nlon,MPI_REAL,rank+1, &
!old       rank,var(1,2),nlon,MPI_REAL,rank-1,rank-1,mycomm,mpistatus,   &
!old       mpierr)
!old     call MPI_SENDRECV(var(1,3),nlon,MPI_REAL,rank-1,rank+nproc, &
!old       var(1,(nlat/nproc)+3),nlon,MPI_REAL,rank+1,rank+nproc+1,    &
!old       mycomm,mpistatus1,mpierr)
!old    end if
!old  end subroutine

end module mod_crtdyn
