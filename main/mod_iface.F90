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
    module mod_iface

    use mod_internal
    use mod_runparams
    use mod_param
    use mod_mssg
    use mod_crtstatic
    use mod_ncio
    use mod_crtdyn
    use mod_libmv
    use mod_statparams
    use mod_mpimess
    use mod_varandtypes

    implicit none
    private

    public chym_init
    public chym_run
    public chym_close

    contains

    subroutine chym_init
      implicit none
      include '../doc/settings.inc'
      integer :: i1,j1,iretval,iv
      integer :: i,displacem
      integer :: century
!
!-----------------------------------------------------------------------
!     Print model informations
!-----------------------------------------------------------------------
!
      if (myid == 0) then
        write (6,'(7x,a,1x,a,1x,a)') 'Start of CHYM version', VERSION, &
         'fortran code'
        write (6,99001)  SVN_REV, __DATE__ , __TIME__
99001   format(2x,' SVN Revision: ',a,' compiled at: data : ',a,  &
             '  time: ',a,/)
      end if
      call setparam
      call acquirescriptpar
      jx = nlon
      iy = nlat
      call set_nproc
      call setup_model_indexes
      if (allocated(portsub)) deallocate(portsub)
      allocate(portsub(jde1gb:jde2gb,ide1gb:ide2gb))
      if (allocated(wkm1sub)) deallocate(wkm1sub)
      allocate(wkm1sub(jde1gb:jde2gb,ide1gb:ide2gb))
      if (allocated(h2osub)) deallocate(h2osub)
      allocate(h2osub(jde1gb:jde2gb,ide1gb:ide2gb))
      if (allocated(bwetsub)) deallocate(bwetsub)
      allocate(bwetsub(jde1gb:jde2gb,ide1gb:ide2gb))
      portsub = 0.
      wkm1sub = 0.
      h2osub = 0.
      bwetsub = 0.
      allocate(cartesian_np(nproc), cartesian_dis(nproc))
      allocate(cartesian_npoint(nproc), cartesian_displ(nproc))
      if (myid == 0) then
        allocate(ide1p(nproc), ide2p(nproc), jde1p(nproc), jde2p(nproc))
        allocate(iypp(nproc) , jxpp(nproc))
        allocate(ide1gbp(nproc), ide2gbp(nproc), jde1gbp(nproc), jde2gbp(nproc))
        allocate(iypgbp(nproc) , jxpgbp(nproc))
      end if
      call mpi_reduce(nngb, nngbp, 1, mpi_integer, mpi_sum, 0, mycomm, &
        mpierr)
      if (myid  == 0) then
        call writestatus('Reading static fields.')
        allocate(port1d(nngbp),bwet1d(nngbp),h2o1d(nngbp))
        port1d=0.
        bwet1d=0.
        h2o1d=0.
      end if
      call mpi_gather(iyp,1,mpi_integer,                          &
                       iypp,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(jxp,1,mpi_integer,                          &
                       jxpp,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(ide1,1,mpi_integer,                          &
                       ide1p,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(ide2,1,mpi_integer,                          &
                       ide2p,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(jde1,1,mpi_integer,                          &
                       jde1p,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(jde2,1,mpi_integer,                          &
                       jde2p,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(ide1gb,1,mpi_integer,                          &
                       ide1gbp,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(ide2gb,1,mpi_integer,                          &
                       ide2gbp,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(jde1gb,1,mpi_integer,                          &
                       jde1gbp,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(jde2gb,1,mpi_integer,                          &
                       jde2gbp,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(jde2gb,1,mpi_integer,                          &
                       jde2gbp,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(iypgb,1,mpi_integer,                          &
                       iypgbp,1,mpi_integer,0,mycomm,mpierr)
      call mpi_gather(jxpgb,1,mpi_integer,                          &
                       jxpgbp,1,mpi_integer,0,mycomm,mpierr)
      allocate(h2osub1d(iypgb*jxpgb))
      allocate(portsub1d(iypgb*jxpgb))
      allocate(bwetsub1d(iypgb*jxpgb))
      allocate(portsub1d_n(iyp*jxp))
      allocate(h2osub1d_n(iyp*jxp))
      allocate(port1d_n(jx*iy))
      allocate(h2o1d_n(jx*iy))
      h2osub1d=0.
      portsub1d=0.
      bwetsub1d=0.
      portsub1d_n=0.
      h2osub1d_n=0.
      port1d_n=0.
      h2o1d_n=0.
      call read_restart_stat_NC
      call mpi_bcast(mchym(9),1,MPI_integer, 0,mycomm,mpierr)
!     print*,"number 1, myid", myid
!     call mpi_barrier(mycomm,mpierr)
!     call exit(0)
      if (myid == 0) then
!       call mpi_barrier(comm,ierr)
         cartesian_dis(1) = 0
         cartesian_displ(1) = 0
         do i = 1,nproc
           cartesian_np(i) = iypp(i)*jxpp(i)
           cartesian_npoint(i) = iypgbp(i)*jxpgbp(i)
           if (i == nproc) exit
           cartesian_dis(i+1) = cartesian_dis(i) + cartesian_np(i)
           cartesian_displ(i+1) = cartesian_displ(i) + cartesian_npoint(i)
         end do
      end if
      if (mchym(9).eq.1) then
        if (myid == 0) then
          call read_restart_dyn_NC
!        call mpi_barrier(comm,ierr)
          displacem = 0
          cartesian_dis(1) = 0
          cartesian_displ(1) = 0
          do i = 1,nproc
            call mypack_real_grid(h2o,h2o1d,ide1gbp(i),ide2gbp(i), &
                             jde1gbp(i),jde2gbp(i),displacem)
            call mypack_real_grid(port,port1d,ide1gbp(i),ide2gbp(i), &
                             jde1gbp(i),jde2gbp(i),displacem)
            call mypack_real_grid(bwet,bwet1d,ide1gbp(i),ide2gbp(i), &
                             jde1gbp(i),jde2gbp(i),displacem)
            displacem = displacem + iypgbp(i)*jxpgbp(i)
            cartesian_np(i) = iypp(i)*jxpp(i)
            cartesian_npoint(i) = iypgbp(i)*jxpgbp(i)
            if (i == nproc) exit
            cartesian_dis(i+1) = cartesian_dis(i) + cartesian_np(i)
            cartesian_displ(i+1) = cartesian_displ(i) + cartesian_npoint(i)
          end do
        end if
        call mpi_scatterv(h2o1d, cartesian_npoint, cartesian_displ, &
           mpi_real,h2osub1d,iypgb*jxpgb, mpi_real, 0, mycomm, mpierr)
        call mpi_scatterv(port1d, cartesian_npoint, cartesian_displ, &
           mpi_real,portsub1d,iypgb*jxpgb, mpi_real, 0, mycomm, mpierr)
        call mpi_scatterv(bwet1d, cartesian_npoint, cartesian_displ, &
           mpi_real,bwetsub1d,iypgb*jxpgb, mpi_real, 0, mycomm, mpierr)
        call myunpack_real_grid(h2osub1d,h2osub,ide1gb,ide2gb,jde1gb,jde2gb)
        call myunpack_real_grid(portsub1d,portsub,ide1gb,ide2gb,jde1gb,jde2gb)
        call myunpack_real_grid(bwetsub1d,bwetsub,ide1gb,ide2gb,jde1gb,jde2gb)
      endif
      if (myid == 0) then
        call mvgetiflags(57,century)
        time=mchym(4)
        call gmafromindex(time,hour,day,month,year)
!        call gmafrommm5index(time,hour,day,month,year)
        write(filename,'(a,a)') trim(schym(11))//'_',trim(schym(3))// &
            '.nc'
        write(filenamerst,'(a,a)') trim(schym(11))//'_',trim(schym(3))// &
            '_rst.nc'
!        call createfile(trim(filename),mchym,rchym,schym,time)
        call createfile(trim(filename),time)
        call createfile_rst(trim(filenamerst),time)
      end if
    end subroutine chym_init

    subroutine chym_run()
      implicit none
      real :: tiorun, trunini, trunfin
      real :: tiorai, traiini, traifin
      real :: ti1, tf1, tot1
      real :: ti2, tf2, tot2
      real :: ti3, tf3, tot3
      real :: ti4, tf4, tot4
      real, dimension(nlon,nlat) :: runrai
      integer :: i,j,k,iv,displacem
      tot1 = 0.0
      tot2 = 0.0
      tot3 = 0.0
      tot4 = 0.0
      tiorai = 0.0
      tiorun = 0.0
      if (myid == 0) then
        write (6,'(/12x,a)') 'Start of integration. Hereafter output is'// &
                 ' sent to chymlog file'
        call getlun(logfile)
        call mvsetflags('CHyM log unit',float(logfile))
        call calibration
        open(logfile,file='chymlog',status='unknown')
!         call mpi_barrier(comm,ierr)
        cartesian_dis(1) = 0
        cartesian_displ(1) = 0
        do i = 1,nproc
          cartesian_np(i) = iypp(i)*jxpp(i)
          cartesian_npoint(i) = iypgbp(i)*jxpgbp(i)
          if (i == nproc) exit
          cartesian_dis(i+1) = cartesian_dis(i) + cartesian_np(i)
          cartesian_displ(i+1) = cartesian_displ(i) + cartesian_npoint(i)
        end do

      end if
      call mpi_bcast(cartesian_np,nproc,mpi_integer,0,mycomm,mpierr)
      call mpi_bcast(cartesian_dis,nproc,mpi_integer,0,mycomm,mpierr)
      call mpi_bcast(mchym(16),1,MPI_integer, 0,mycomm,mpierr)
      call mpi_bcast(mchym(14),1,MPI_integer, 0,mycomm,mpierr)
      call mpi_bcast(chym_steps,1,MPI_integer, 0,mycomm,mpierr)
      call mpi_bcast(dx(1,1),nlat*nlon,MPI_REAL, 0,mycomm,mpierr)
      call mpi_bcast(ienddate,1,MPI_integer, 0,mycomm,mpierr)
      call mpi_bcast(time,1,MPI_integer, 0,mycomm,mpierr)
      call mpi_barrier(mycomm,mpierr)
!      do hourstep=1,mchym(16)
      hourstep = 0
      do while (time <= ienddate)
        hourstep = hourstep + 1
        if (myid == 0 ) then
        !  if (hourstep /= 1) then
        !    avgtemp = avgtemp + temp
        !    avgevap = avgevap + evap
        !    avgrsrm = avgrsrm + rsrm
        !    avgbwet = avgbwet + bwet
        !    avgport = avgport + port
        !    avggh2o = avggh2o + gh2o
        !    avgh2o = avgh2o + h2o
        !    avgsnow = avgsnow + snow
        !    avgdeepw = avgdeepw + deepw
        !    avgddeepw = avgddeepw + ddeepw
         ! end if
          call createnc1(hourstep)
        end if
        if ( mod(hourstep,mchym(14)).eq.0 ) then
          if (myid == 0) then
            CALL CPU_TIME(ti1)
         !   call createnc1(hourstep)
            call add_timestep
            call writeoutfile(0)
            call temperature
            CALL CPU_TIME(tf1)
            tot1 = tot1 + tf1 - ti1
            CALL CPU_TIME(traiini)
            call rainfall
            CALL CPU_TIME(traifin)
            tiorai = tiorai + traifin - traiini
            !'rain' has to be updated here because rain is modified by
            !one of the following calls too, before being written as 'ara'
            CALL CPU_TIME(ti2)
            avgrain = avgrain + rain
            call snowcover
            time=increasetime(time)
            call dataorafromday (hour,day,month,year,now) !I'm not sure about that (???)
            avgrsrm = avgrsrm + rsrm
            call writeoutfile(1)
            call snowacc
            call evapotranspiration
            call melting
            CALL CPU_TIME(tf2)
            tot2 = tot2 + tf2 - ti2
            CALL CPU_TIME(ti3)
            call groundwater
            call returnflow
            avgarai = avgarai + rain
            call writeoutfile(2)
            CALL CPU_TIME(tf3)
            tot3 = tot3 + tf3 - ti3
!          call runoff(rain)
          end if
          call mpi_barrier(mycomm,mpierr)
!          call mpi_bcast(evap(1,1),nlat*nlon,MPI_REAL, 0,comm,mpierr)
          call mpi_bcast(time,1,MPI_integer, 0,mycomm,mpierr)
          call mpi_bcast(avgarai(1,1),nlat*nlon,MPI_REAL, 0,mycomm,mpierr)
          call mpi_bcast(rain(1,1),nlat*nlon,MPI_REAL, 0,mycomm,mpierr)
!          call mpi_bcast(deepw(1,1),nlat*nlon,MPI_REAL, 0,comm,mpierr)
!          call mpi_bcast(h2o(1,1),nlat*nlon,MPI_REAL, 0,comm,mpierr)
!          call mpi_bcast(rsrm(1,1),nlat*nlon,MPI_REAL, 0,comm,mpierr)
          call mpi_barrier(mycomm,mpierr)
          if (myid==0)then
            CALL CPU_TIME(trunini)
          end if
!          if (avgarai == 0) then
!            runrai = 0
!          else
!            runrai = avgarai/nsave
!          endif
          if (myid == 0) then
            displacem = 0
            do i = 1,nproc
              call mypack_real_grid(h2o,h2o1d,ide1gbp(i),ide2gbp(i), &
                           jde1gbp(i),jde2gbp(i),displacem)
              displacem = displacem + iypgbp(i)*jxpgbp(i)
              if (i == nproc) exit
            end do
          end if
          call mpi_scatterv(h2o1d, cartesian_npoint, cartesian_displ, &
              mpi_real,h2osub1d,iypgb*jxpgb, mpi_real, 0, cartesian_communicator, mpierr)
!              mpi_real,h2osub1d,iypgb*jxpgb, mpi_real, 0, mycomm, mpierr)
          call myunpack_real_grid(h2osub1d,h2osub,ide1gb,ide2gb,jde1gb,jde2gb)

          call runoff(rain)
          if (myid==0)then
            CALL CPU_TIME(trunfin)
            tiorun = tiorun + trunfin - trunini
          end if
          iv = 1
          do i=ide1,ide2
            do j=jde1,jde2
!            hsub1d(iv)=hsub(j,i)
              portsub1d_n(iv)=portsub(j,i)
              h2osub1d_n(iv)=h2osub(j,i)
              iv = iv+1
            end do
          end do
          iv = cartesian_np(myid+1)
          call mpi_gatherv(portsub1d_n,iv,MPI_REAL,port1d_n,&
             cartesian_np, cartesian_dis, &
             MPI_REAL,0,cartesian_communicator, mpierr)
          call mpi_gatherv(h2osub1d_n,iv,MPI_REAL,h2o1d_n,&
             cartesian_np, cartesian_dis, &
             MPI_REAL,0,cartesian_communicator, mpierr)
          if (myid == 0) then
            iv = 1
            do k=1,nproc
              do i=ide1p(k),ide2p(k)
                do j=jde1p(k),jde2p(k)
                  port(j,i) = port1d_n(iv)
                  h2o(j,i) = h2o1d_n(iv)
                  iv = iv + 1
                end do
              end do
            end do
          end if
          call mpi_barrier(mycomm,mpierr)
          if (myid == 0) then
            CALL CPU_TIME(ti4)
            avgtemp = avgtemp + temp
            avgevap = avgevap + evap
            avgrsrm = avgrsrm + rsrm
            avgbwet = avgbwet + bwet
            avgport = avgport + port
            avggh2o = avggh2o + gh2o
            avgh2o = avgh2o + h2o
            avgsnow = avgsnow + snow
            avgdeepw = avgdeepw + deepw
            avgddeepw = avgddeepw + ddeepw
            call writeoutfile(3)
            CALL CPU_TIME(tf4)
            tot4 = tot4 + tf4 - ti4
            avgtemp = 0.0
            avgevap = 0.0
            avgrain = 0.0
            avgarai = 0.0
            avgrsrm = 0.0
            avgbwet = 0.0
            avgport = 0.0
            avggh2o = 0.0
            avgh2o = 0.0
            avgsnow = 0.0
            avgdeepw = 0.0
            avgddeepw = 0.0
          end if
        else
          if (myid == 0) then
            call temperature
            call rainfall
            !'rain' has to be updated here because rain is modified by
            !one of the following calls too, before being written as 'ara'
            avgrain = avgrain + rain
            call snowcover
            time=increasetime(time)
            call dataorafromday (hour,day,month,year,now)
            call snowacc
            call evapotranspiration
            call melting
            call groundwater
            call returnflow
            avgarai = avgarai + rain
          end if
          call mpi_bcast(time,1,MPI_integer, 0,mycomm,mpierr)
          call mpi_bcast(rain(1,1),nlat*nlon,MPI_REAL, 0,mycomm,mpierr)
          call mpi_barrier(mycomm,mpierr)
          if (myid == 0) then
            displacem = 0
            do i = 1,nproc
              call mypack_real_grid(h2o,h2o1d,ide1gbp(i),ide2gbp(i), &
                           jde1gbp(i),jde2gbp(i),displacem)
              displacem = displacem + iypgbp(i)*jxpgbp(i)
              if (i == nproc) exit
            end do
          end if
          if (allocated(h2osub1d)) deallocate(h2osub1d);allocate(h2osub1d(iypgb*jxpgb))
          call mpi_scatterv(h2o1d, cartesian_npoint, cartesian_displ, &
              mpi_real,h2osub1d,iypgb*jxpgb, mpi_real, 0, mycomm, mpierr)
          call myunpack_real_grid(h2osub1d,h2osub,ide1gb,ide2gb,jde1gb,jde2gb)
!bug          call runoff(rain)
          iv = 1
          do i=ide1,ide2
            do j=jde1,jde2
!            hsub1d(iv)=hsub(j,i)
              portsub1d_n(iv)=portsub(j,i)
              h2osub1d_n(iv)=h2osub(j,i)
              iv = iv+1
            end do
          end do
          iv = cartesian_np(myid+1)
          call mpi_gatherv(portsub1d_n,iv,MPI_REAL,port1d_n,&
             cartesian_np, cartesian_dis, &
             MPI_REAL,0,cartesian_communicator, mpierr)
          call mpi_gatherv(h2osub1d_n,iv,MPI_REAL,h2o1d_n,&
             cartesian_np, cartesian_dis, &
             MPI_REAL,0,cartesian_communicator, mpierr)
          if (myid == 0) then
            iv = 1
            do k=1,nproc
              do i=ide1p(k),ide2p(k)
                do j=jde1p(k),jde2p(k)
                  port(j,i) = port1d_n(iv)
                  h2o(j,i) = h2o1d_n(iv)
                  iv = iv + 1
                end do
              end do
            end do
          end if
          call mpi_barrier(mycomm,mpierr)
          if (myid == 0) then
            avgtemp = avgtemp + temp
            avgevap = avgevap + evap
            avgrsrm = avgrsrm + rsrm
            avgbwet = avgbwet + bwet
            avgport = avgport + port
            avggh2o = avggh2o + gh2o
            avgh2o = avgh2o + h2o
            avgsnow = avgsnow + snow
            avgdeepw = avgdeepw + deepw
            avgddeepw = avgddeepw + ddeepw
          end if
!          call runoff(rain)
        end if
        if ( myid == 0 ) then
          call createnc2(hourstep)
!          if ( mod(hourstep,mchym(26)).eq.0 ) then
!            call createnc2(hourstep)
!            call addrst_timestep
!            call writerstfile
!          end if
        end if
        call mpi_barrier(mycomm,mpierr)
      end do
      if (myid == 0) then
      end if
      call mpi_barrier(mycomm,mpierr)
    end subroutine


    subroutine chym_close()
      if (myid == 0) then
        write (6,'(/12x,a)') 'Closing all files'
        call closefile
        write (6,'(/12x,a)') 'Simulation completed!'
      end if
      call mpi_barrier(mycomm,mpierr)
      call mpi_finalize(mpierr)
    end subroutine


    end module
