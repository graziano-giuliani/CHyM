module mod_dewfun

        use mod_museo
        use mod_time
        contains
        subroutine tsvalidation(db,annoi,pmin,pmax)
        implicit none
        integer nrec,rec,annoi,ndata,n1,counta,ii,countabis
        character com*64,sens*32
        real ylat,ylon
        integer mesei,i1,i2,ndi
        real x1(366*24*4),y1(366*24*4),y2(366*24*4),y3(366*24*4),totp,pmin,pmax
        real, allocatable, dimension(:,:) :: y4
        real, allocatable, dimension(:) :: ylat1,ylon1
        character, allocatable, dimension(:) :: sens1*32
        character db*(*),file*32
        counta = 0
        nrec=mvlibmagicnum(1)                           ! Must depend from db
        write(file,'(a,i4,a)') db(1:len_trim(db)),annoi,'-fil.dat'
        call getlun(lun)
        open(lun,file=file,status='unknown',form='unformatted')
        allocate(y4(366*24,nrec))
        allocate(ylat1(nrec))
        allocate(ylon1(nrec))
        allocate(sens1(nrec))
        do rec=1,nrec
           y1=-9999
           call museotimeseries(db,rec,annoi,x1,y1,n1)
           y2=y1
           call convert2hourlyres(y2,1)         ! Must depend from db
           ndata=8760
           call extractgooddata(x1,y2,ndata,-9999)
           totp=0.0
           do i=1,ndata
              totp=totp+y2(i)
           enddo
           if (totp/ndata.gt.pmax) then
              do mesei=1,12
                 i1=index15m(00,00,01,mesei,annoi)
                 i2=index15m(45,23,monthlen(mesei,annoi),mesei,annoi)
                 totp=0.0 ; ndi=0
                 do i=i1,i2     
                    if (nint(y1(i)).ne.-9999) then
                       totp=totp+y1(i)
                       ndi=ndi+1
                    endif
                 enddo
                 if (ndi.gt.0) then
                    if (totp/ndi.gt.pmax) y1(i1:i2)=-9999.0
                 endif
              enddo
              n1=0
              do i=1,366*24*4
                 if (nint(y1(i)).ne.-9999) n1=i
              enddo
              write(6,'(10x,a,i4,a,i5)')'Rec. ',rec,' filtered '// &
                'by tsvalidation, remaining components: ',n1
           endif
           if (n1.gt.0) call spikefilter(y1,n1,rec)
           if (n1.gt.0) then
           counta = counta + 1
           call museoanagrafica('dewrain',rec,com,sens,ylat,ylon)
           y3=y1
           n1 = 366*24
           call convert2hourlyres(y3,1)         ! Must depend from db
           y4(:,counta) = y3(1:n1)
           ylat1(counta) = ylat
           ylon1(counta) = ylon
           sens1(counta) = sens
!           write(lun) n1,(y3(i),i=1,n1)
           end if
        enddo
        do i=1,366*24
           countabis = 0
           do ii=1,counta
              if (y4(i,ii)>=-99) then
                countabis = countabis + 1 
              end if
           end do
           write(lun) countabis
           do ii=1,counta
             if (y4(i,ii)>=-99) then 
               write(lun) ylat1(ii),ylon1(ii),sens,y4(i,ii)
             end if
!             write(lun) counta,(ylat1(ii),ylon1(ii),sens,y4(i,ii),ii=1,counta)
           end do
        end do
        print*,counta
        print*,"sensore",sens1(counta),com,ylat1(counta),ylon1(counta),rec
        close(lun)
        write(6,'(15x,a)') 'A new version of '//file(1:len_trim(file))//' created'
        return
        end subroutine tsvalidation


        subroutine spikefilter(y,n,rec)
        implicit none
        integer n,j,il,rec,maxd,ncor,iri
        real y(n)
        data maxd /30/ ; save maxd
        logical log ; data log /.false./ ; save log
        ncor=0
        do i=2,n-1
           if (y(i).gt.100) then
              ncor=ncor+1
              y(i)=-9999
           endif
        enddo
        do i=5,n-5
           if (nint(y(i)).ne.-9999) then
              il=-9999 ; iri=-9999
              do j=i-4,i-1
                 if (nint(y(j)).ne.-9999) il=j
              enddo
              do j=i+4,i+1,-1
                 if (nint(y(j)).ne.-9999) iri=j
              enddo
              if (il.ne.-9999.and.iri.ne.-9999) then
                 if (y(i)-y(il).gt.maxd.and.y(i)-y(iri).gt.maxd) then
                    ncor=ncor+1
                    y(i)=-9999
                 endif
              endif
           endif
        enddo
        if (ncor.gt.100) then
           write(6,'(10x,a,i4,a)') 'Rec. ',rec,' discarded by spikefilter'
           n=0
        else if (ncor.gt.15) then
           write(6,'(10x,a,i4,i5,a)') &
                'Rec. ',rec,ncor,' values discarded - spikefilter'
        endif
        return
        end subroutine spikefilter


        subroutine extractgooddata(x,y,n,noval)
        implicit none
        integer n,ln,i,noval ; real x(n),y(n)
        real , allocatable , dimension(:) :: lx,ly
        allocate(lx(n)) ; allocate(ly(n))
        ln=0
        do i=1,n
           if (nint(y(i)).ne.-9999) then
              ln=ln+1 ; lx(ln)=x(i) ; ly(ln)=y(i)
           endif
        enddo
        x(1:ln)=lx(1:ln) ; y(1:ln)=ly(1:ln) ; n=ln
        deallocate(lx) ; deallocate(ly)
        return
        end subroutine extractgooddata
end module mod_dewfun
