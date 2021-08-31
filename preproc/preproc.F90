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
program preproc      

use mod_internal
use mod_runparams
use mod_param
use mod_mssg
use mod_crtstatic
use mod_ncio
use mod_mpimess
use mod_varandtypes
implicit none
include '../doc/settings.inc'
call mpi_init(mpierr)
call mpi_comm_dup(mpi_comm_world,mycomm,mpierr)
call mpi_comm_size(mycomm, nproc, mpierr)
call mpi_comm_rank(mycomm, myid, mpierr)

!
!-----------------------------------------------------------------------
!     Print model informations
!-----------------------------------------------------------------------
!
    write (6,'(7x,a,1x,a,1x,a)') 'Start of CHYM preprocessing version', VERSION, 'fortran code'
    write (6,99001)  SVN_REV, __DATE__ , __TIME__
99001   format(2x,' SVN Revision: ',a,' compiled at: data : ',a,  &
        '  time: ',a,/)
!
!-----------------------------------------------------------------------
!     Read configuration parameters 
!-----------------------------------------------------------------------
!
    call setparam
    call acquirescriptpar
    call calibration
!
!-----------------------------------------------------------------------
!    Create static fields
!-----------------------------------------------------------------------
!
    call writestatus('Creating static fields.')
    if (mchym(13)==10) then
!      call buildhyddirmap
      call hydrodem_1k
    else
      call builddem                      ! Fill DEM matrix
    end if
    call buildlandusemap               ! Fill landuse matrix
    if (.not.mchym(13)==10) then 
      call buildflowdirmap(.True.)               ! Fill fmap and smooth dem matrix
    end if
    call areamatrix                    ! Fill drai and area matrix
    call buildacclivitymap             ! Fill accl matrix
    if (.not.mchym(13)==10) then 
      call reconnectdem                  ! Further Check on dem & fmap
    end if
    call riveronlanduse                ! Fill luse matrix with river code
    call runoffspeed
!!    call buildicemap                   ! Fill the ice matrix
    call write_stat_NC

end program preproc
