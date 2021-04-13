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

module mod_mm5

  use mod_libmv
  use mod_strman
  use mod_museo

  contains

  subroutine mm5colorlevels(lm,levels)
    implicit none
    integer :: lm
    real , dimension(:) :: levels
    intent (in) levels , lm
    integer :: i , n
    lmeth = lm
    if ( lmeth==0 ) then
      n = 0
    else if ( lmeth==-1 ) then
      n = 3
    else if ( lmeth==-2 ) then
      n = 1
    else if ( lmeth>0 ) then
      n = 0
    else if ( lmeth<-2 ) then
      n = abs(lmeth) - 2
    end if
    do i = 1 , n
      rlevels(i) = levels(i)
    end do
  end subroutine mm5colorlevels
 
  real function mm5mrf(i,j)
    implicit none
    integer :: i , j
    intent (in) i , j
    mm5mrf = mrf(i,j)
  end function mm5mrf
 
  integer function mm5mif(i,j)
    implicit none
    integer :: i , j
    intent (in) i , j
    mm5mif = mif(i,j)
  end function mm5mif
 
  subroutine bhi2mif
    implicit none
    integer :: id , ih , im , imm , iy
    mif(1,1) = bhi(1,1)            ! Data type
    mif(2,1) = bhi(5,1)            ! I Coarse Domain Grid Dim
    mif(3,1) = bhi(6,1)            ! J Coarse Domain Grid Dim
    mif(4,1) = bhi(7,1)            ! Map Projection
    mif(8,1) = bhi(11,1)           ! Grid Ofset I Direction
    mif(9,1) = bhi(12,1)           ! Grid Ofset J Direction
    mif(101,1) = bhi(13,1)         ! Domain ID
    mif(104,1) = bhi(16,1)         ! Domain Grid Dimension I
    mif(105,1) = bhi(17,1)         ! Domain Grid Dimension J
    mif(108,1) = bhi(20,1)         ! Grid Size Ratio
    read (current_date,'(i4,1x,i2,1x,i2,1x,i2,1x,i2)') &
      iy , im , id , ih , imm
    if ( iy>=2000 ) then
      iy = iy - 2000
    else
      iy = iy - 1900
    end if
    if ( imm==59 ) ih = ih + 1
    mif(1,mif(1,1)) = ih + id*100 + im*10000 + iy*1000000
    mif(1,7) = ih + id*100 + im*10000 + iy*1000000
 
    mrf(2,1) = bhr(2,1)            ! Coarse Domain Center Latitude
    mrf(3,1) = bhr(3,1)            ! Coarse Domain Center Longitude
    mrf(5,1) = bhr(5,1)            ! True Latitude 1
    mrf(6,1) = bhr(6,1)            ! True Longitude 1
    mrf(101,1) = bhr(9,1)/1000.    ! Domain Grid Distance (km->m)
    mrf(102,1) = bhr(10,1)         ! I Location in Coarse Dom. (1,1)
    mrf(103,1) = bhr(11,1)         ! J Location in Coarse Dom. (1,1)
    mrf(104,1) = bhr(12,1)         ! I Location in Coarse Dom. (IX,JX)
    mrf(105,1) = bhr(13,1)         ! J Location in Coarse Dom. (IX,JX)
 
    ! read (current_date,'(i)') mif(1,mif(1,1))
  end subroutine bhi2mif

  subroutine mm5basefld(fld,dom,dum,ixm,jxm,ix,jx,errore)
    implicit none
    integer :: dom , errore , ix , ixm , jx , jxm
    character(len=*) :: fld
    real , dimension(ixm,jxm) :: dum
    intent (in) dom , ixm , jxm
    intent (out) dum
    intent (inout) errore , ix , jx
    integer :: i , if , irec , j , lun
    character(8) :: lfld
 
    if ( iflg(1)>=100 ) &
      write (6,'(7x,a,i1)') 'MVLib mm5basefld called for field '//fld// &
                            ' on domain ' , dom
    errore = 0
    call cv2lower(fld,lfld)
    if ( lfld(1:8)=='land use' ) then
      if = 1
    else if ( lfld(1:7)=='terrain' ) then
      if = 2
    else if ( lfld(1:8)=='latitdot' ) then
      if = 3
    else if ( lfld(1:8)=='longidot' ) then
      if = 4
    else if ( lfld(1:8)=='latitcrs' ) then
      if = 5
    else if ( lfld(1:8)=='longicrs' ) then
      if = 6
    else
      errore = 1
    end if
 
    if ( dom==0 .and. iflg(5)<=3 .and. iflg(5)>=1 ) then
      irec = (if-1)*3 + iflg(5)
    else if ( dom>=1 .and. dom<=3 ) then
      irec = (if-1)*3 + dom
    else
      errore = 2
    end if
 
    if ( errore==0 ) then
      call getlun(lun)
      call openmuseofiles(lun,'mm5.dat',0)
      do i = 1 , irec - 1
        read (lun,err=200)
      end do
      read (lun,err=200) ix , jx , ((dum(i,j),i=1,ix),j=1,jx)
    end if
 100  close (lun)
    if ( errore/=0 ) &
      write (6,'(7x,a,i1)') 'MVLib mm5basefld: error number ', errore
    return
 200  errore = 4
    go to 100
  end subroutine mm5basefld

  subroutine mm5readfld(u,fld,n,dum,ixm,jxm,kxm,ix,jx,kx)
    implicit none
    character(*) :: fld
    integer :: ix , ixm , jx , jxm , kx , kxm , n , u
    real , dimension(ixm,jxm,kxm) :: dum
    intent (inout) dum , ix , jx , kx , n
    logical :: debug , notfound
    character(8) , save :: fldes
    integer , save :: i , idate , ii , ind , jj , kk , lim1 , lim2 ,   &
                 n2d , n3d
    character(8) , dimension(maxfl) , save :: ich
    integer , dimension(maxfl) , save :: icr
    integer , dimension(100) :: lcount , lix , ljx , lkx , n2 , n3
    character(8) , dimension(100,maxfl) , save :: lich
    data lcount/100*1/
    data debug/.false./
    if ( u<1 .or. u>100 ) then
      write (6,'(a,i10)') ' MM5Readfld: Bad unit number passed ' , u
      stop ' MM5Readfld: program stopped.'
    end if
    if ( n<0 ) then
      write (6,'(7x,a,i3)') 'MVLib MM5Readfld: "strange" value for '//       &
                            'slices counter ' , n
      write (6,'(25x,a)') 'Control returned to calling without read.'
      return
    else if ( n==0 ) then
      lcount(u) = 1
      mm5flag = 99
      rewind (u)
      read (u,err=50) mm5flag
 50   rewind (u)
      if ( mm5flag==0 ) then
        lcount(u) = -3
        if ( iflg(1)>0 ) write (6,'(7x,a)') 'MM5readfld V3 file detected.'
        call mm5v3readfld(u,fld,n,dum,ixm,jxm,kxm,ix,jx,kx)
        return
      end if
    else if ( lcount(u)==-3 ) then
      call mm5v3readfld(u,fld,n,dum,ixm,jxm,kxm,ix,jx,kx)
      return
    end if
 
    notfound = .false.
100 if ( lcount(u)==1 ) then
      if ( n<0 ) return
      nd(1) = 0
      nd(3) = 0
      call mm5readcnt(u,n,ind,n3d,n2d,ix,jx,kx,idate,ich,icr)
      if ( n<0 ) then
        return
      else if ( n==1 ) then
        lix(u) = ix
        ljx(u) = jx
        lkx(u) = kx
        n3(u) = n3d
        n2(u) = n2d
        do i = 1 , n2(u) + n3(u)
          lich(u,i) = ich(i)
        end do
      else
        do i = 1 , nd(1)
          if ( ic1(i)==201 .and. ir1(i)==ind ) then
            if ( debug ) &
              write (6,'(a,i2,a,i2,a,i2)') ' changing n3(' , u , &
                         ') from ' , n3(u) , ' to ' , id1(i)
            n3(u) = id1(i)
          end if
          if ( ic1(i)==202 .and. ir1(i)==ind ) then
            if ( debug ) &
              write (6,'(a,i2,a,i2,a,i2)') ' changing n2(' , u ,  &
                         ') from ' , n2(u) , ' to ' , id1(i)
            n2(u) = id1(i)
          end if
        end do
        lim1 = 204 + 1
        lim2 = 204 + n3(u) + n2(u)
        do i = 1 , nd(3)
          if ( ic3(i)>=lim1 .and. ic3(i)<=lim2 .and. ir3(i)==ind ) &
            lich(u,ic3(i)-lim1+1) = ch1(i)
        end do
      end if
      lcount(u) = 1
    end if
    fldes = fld
    ix = lix(u)
    jx = ljx(u)
    kx = lkx(u)
    call cv2upper(fldes,fldes)
    do i = lcount(u) , n3(u) + n2(u)
      if ( trim(lich(u,i))//' '==trim(fldes)//' ' ) then
        if ( i>n3(u) ) kx = 1
        if ( ix>ixm .or. jx>jxm .or. kx>kxm ) then
          write (6,'(2x,a)') 'MM5readfld error: matrix too small '// &
                   ' to load '//trim(fldes)//' mm5 field'
          write (6,'(12x,a,3i4)') ' Grid dimensions are ' , ix , jx , kx
          write (6,'(12x,a,3i4)') ' Matrix dimensions are ' , ixm , jxm , kxm
          write (6,'(12x,a)') ' program stopped'
          stop ' '
        end if
        read (u,err=200,end=200) (((dum(ii,jj,kk),ii=1,ix),jj=1,jx),kk=1,kx)
        if ( iflg(1)>0 ) &
          write (6,'(7x,a,i2,a,i2,a,f10.2)')        &
            'MM5readfld Read '//fldes(1:8)          &
            //' for slice # ' , n , ' unit ' , u ,  &
            '. Typ. value: ' , dum(ix/2,jx/2,1)
        if ( i==n3(u)+n2(u) ) then
          lcount(u) = 1
        else
          lcount(u) = i + 1
        end if
        return
      else
        read (u,err=200)
        if ( debug ) &
          write (6,'(2x,a,i2,a,i2,a,i2)') 'MM5readfld Skip '//    &
               lich(u,i)(1:8)//' for slice # ' , n ,              &
               ' from unit ' , u , ' i=' , i
      end if
    end do
    if ( notfound .or. lcount(u)==1 ) then
      write (6,'(8x,a,i2)') 'MM5readfld error: not found field '// &
                   trim(fldes)//' for slice # ' , n
      write (6,'(10x,a)') ' The following fields are available:'
      do i = 1 , n3(u) + n2(u)
        write (6,'(12x,a)') lich(u,i)
      end do
      stop ' Stop inside MM5readfld'
    else
      lcount(u) = 1
      notfound = .true.
      go to 100
    end if
200 write (6,'(7x,a,i2)') &
          'MM5readfld error: Abnormal read-error when '&
          //'reading '//fldes(1:8)//' for slice ' , n
    n = -40
  end subroutine mm5readfld

  subroutine mm5readcnt(iunit,n,ind,num3d,num2d,ix,jx,kx,id,ich,icr)
    implicit none
    integer :: id , ind , iunit , ix , jx , kx , n , num2d , num3d
    character(8) , dimension(:) :: ich
    integer , dimension(:) :: icr
    intent (out) ich , icr , id , ix , jx , kx
    intent (inout) ind , num2d , num3d
    integer , save :: i
    call mm5readhd(iunit,n,mif,mrf,mifc,mrfc)
    ind = mif(1,1)
    ix = mif(104,1)
    jx = mif(105,1)
    if ( ind<=3 .or. ind==7 ) then
      kx = mif(101,ind)
    else
      kx = nint(mrf(101,ind))
    end if
    num3d = mif(201,ind)
    num2d = mif(202,ind)
    id = mif(1,6)
    do i = 1 , num3d + num2d
      ich(i) = mifc(204+i,ind)(1:8)
      icr(i) = mif(204+i,ind)
    end do
  end subroutine mm5readcnt

  subroutine mm5readhd(iunit,n,mmif,mmrf,mmifc,mmrfc)
    implicit none
    integer :: iunit , n
    integer , dimension(1000,20) :: mmif
    character(80) , dimension(1000,20) :: mmifc , mmrfc
    real , dimension(1000,20) :: mmrf
    intent (in) iunit
    intent (out) mmifc , mmrf , mmrfc
    intent (inout) mmif , n
    logical :: debug
    data debug/.false./
    if ( n==0 ) then
      rewind (iunit)
      read (iunit,end=100,err=100) mmif , mmrf , mmifc , mmrfc
      if ( debug ) write (6,'(a)') '  MM5readhd Read first initial header '
      lindex(iunit) = mmif(1,1)
      if ( mmif(1,1)>20 ) mmif(1,1) = mmif(1,1) - 20
      n = 1
      ! iflg(5)=mif(101,1)          ! Conviene commentarla ?
    else if ( n>0 ) then
      if ( lindex(iunit)<=20 ) then
        read (iunit,end=200,err=300) mmif , mmrf , mmifc , mmrfc
      else
        write (6,'(7x,a,i2)') 'MVLib MM5ReadHd: Flux error 01'
        call exit(0)
!
!  Questo e' il caso dei files versione 2 compressi
!  Non dovrebbe piu' accadere.
!
!       read (iunit,end=998) nd,(id1(i),ic1(i),ir1(i),i=1,nd(1)), &
!                               (rd1(i),ic2(i),ir2(i),i=1,nd(2)), &
!                               (ch1(i),ic3(i),ir3(i),i=1,nd(3)), &
!                               (ch2(i),ic4(i),ir4(i),i=1,nd(4))
!       do i=1,nd(1)
!         mif(ic1(i),ir1(i))=id1(i)
!       end do
!       do i=1,nd(2)
!         mrf(ic2(i),ir2(i))=rd1(i)
!       end do
!       do i=1,nd(3)
!         mifc(ic3(i),ir3(i))=ch1(i)
!       end do
!       do i=1,nd(4)
!         mrfc(ic4(i),ir4(i))=ch2(i)
!       end do
!       if (debug) then
!         write (6,'(a,i2,a,i2)') ' **** Unita'': ',iunit,' slice ',n+1
!         write (6,'(5x,i2,a)') nd(1),' diff. found for mif array'
!         write (6,'(5x,i2,a)') nd(2),' diff. found for mrf array '
!         write (6,'(5x,i2,a)') nd(3),' diff. found for mifc array '
!         write (6,'(5x,i2,a)') nd(4),' diff. found for mrfc array '
!       endif
!       mif(1,1)=lindex(iunit)-20
      end if
      n = n + 1
      ! iflg(5)=mmif(101,1)          ! Conviene commentarla ?
    else
      write (6,'(7x,a,i3)') 'MVLib MM5ReadHd: "strange" value for '// &
                      'slices counter ' , n
      write (6,'(24x,a)') 'Control returned to calling without read.'
    end if
    return
100 n = -30
    write (6,'(7x,a,i2)') 'MVLib MM5ReadHd: abnormal error while'//         &
                          ' reading first header from unit ' , iunit
    write (6,'(24x,a)') 'Control returned to calling without read.'
    return
200 n = -10
    return
300 write (6,'(7x,a,i3)') 'MVLib MM5ReadHd: abnormal error while'//   &
                          ' reading header for slice: ' , n
    n = -20
  end subroutine mm5readhd
 
  subroutine mm5v3readfld(u,fldes,n,dum,ixm,jxm,kxm,ix,jx,kx)
    implicit none
    character(len=*) :: fldes
    integer :: ix , ixm , jx , jxm , kx , kxm , n , u
    real , dimension(ixm,jxm,kxm) :: dum
    intent (in) fldes , ixm , jxm , kxm , u
    intent (inout) dum , ix , jx , kx , n
    character(len=10) , dimension(100) , save :: campi
    character(len=8) , save :: fld
    integer , save :: i , icont , ixs , j , jxs , k , kxs , mindex , numfield
    logical , save :: nextslice

    fld = fldes
    call cv2upper(fld,fld)
    icont = 0
    if ( n<0 ) then
      write (6,'(7x,a,i3)') 'ReadMM5v3: "strange" value for '// &
                            'slices counter ' , n
      write (6,'(25x,a)') 'Control returned to calling without read.'
      return
    end if
    if ( u<1 .or. u>100 ) then
      write (6,'(7x,a,i10)') 'ReadMM5v3: Bad unit number passed ' , u
      return
    end if
    nextslice = .false.
    if ( n==0 ) rewind (u)
 100  read (u,end=200) mm5flag
    if ( mm5flag==0 ) then
      read (u) bhi , bhr , bhic , bhrc
      mindex = bhi(1,1)
      if ( iflg(1)>=10 ) write (6,'(7x,a,i4,3(x,i2))') 'ReadMM5v3: '//     &
                   trim(bhic(1,1))//' - Beginnig Time: ' , bhi(5,mindex) , &
                   bhi(6,mindex) , bhi(7,mindex) , bhi(8,mindex)
      n = 1
      numfield = 0
      go to 100
    else if ( mm5flag==1 ) then
      read (u) mm5ndim , start_index , end_index , xtime , staggering ,      &
               ordering , current_date , mm5_name , mm5_unit , description
      ixs = end_index(1)
      jxs = end_index(2)
      kxs = end_index(3)
      numfield = numfield + 1
      write (campi(numfield),'(a)') trim(mm5_name)
      if ( trim(mm5_name)==trim(fld) ) then
        if ( ordering=='YX' .or. ordering=='CA' ) then
          read (u) ((dum(i,j,1),i=1,ixs),j=1,jxs)
          if ( ixs>ixm .or. jxs>jxm ) then
            write (6,'(12x,a,3i4)') 'ixm,jxm,kxm=' , ixs , jxs , kx
            stop ' MM5v3Readfld stopped - Please increase ixm,jxm,kxm '
          end if
        else if ( ordering=='P' .or. ordering=='S' ) then
          read (u) (dum(i,1,1),i=1,ixs)
        else
          if ( ixs>ixm .or. jxs>jxm .or. kxs>kxm ) then
            write (6,'(12x,a,3i4)') 'ixm,jxm,kxm=' , ixs , jxs , kxs
            stop ' MM5v3Readfld stopped - Please increase ixm,jxm,kxm '
          end if
          read (u) (((dum(i,j,k),i=1,ixs),j=1,jxs),k=1,kxs)
        end if
        icont = icont + 1
        ix = end_index(1)
        jx = end_index(2)
        kx = end_index(3)
        if ( iflg(1)>=10 ) &
          write (6,'(7x,a,i2,a,i2,a,f10.2)') 'ReadMM5v3:  Reading '// &
                   trim(mm5_name)//' for slice # ' , n ,&
                   ' unit ' , u , '. Typ. value: ' , dum(ix/2,jx/2,1)
        call bhi2mif
        return
      else
        if ( iflg(1)>=100 ) &
          write (6,'(7x,a)') 'ReadMM5v3: Skipping '//mm5_name
        read (u)
        go to 100
      end if
    else if ( mm5flag==2 ) then
      if ( nextslice ) then
        write (6,'(7x,a)') 'ReadMM5v3: field '//fld//' not exist!'
        write (6,'(10x,a)') 'List of good fields follow.'
        write (6,'(12x,a)') (campi(i),i=1,numfield)
        stop '                   MVLib ReadMM5v3 - Abnormal end'
      end if
      numfield = 0
      nextslice = .true.
      if ( iflg(1)>=100 ) write (6,'(7x,a,i2)') &
                  'ReadMM5v3:        End of time step: ' , n
      n = n + 1
      go to 100
    end if
 200  n = -1
  end subroutine mm5v3readfld
 
end module mod_mm5
