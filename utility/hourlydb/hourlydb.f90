        program hourlydb
        use mod_hourdb
        use mod_time 
        use mod_museo
        character museod*42,file*80,data*60,strum*16
        parameter (maxdata=5000)
        real llat(maxdata),llon(maxdata),lp(maxdata)
        integer source(maxdata)
        common /timefl/ lgiorno,lmese,lanno,lidate,logf
        integer lgiorno,lmese,lanno,lidate ; logical logf,dew
        !iflag per scegliere il datatset da convertire 2 per temporaria e 1 per pioggiaoraria
        logf=.true. ;  lanno=2001 ; lanno=i4digityear(lanno) ; iflag=1
        if (iflag.gt.10) then
           iflag=iflag-10
           dew=.false.
        else
           dew=.true.
        endif
        if (iflag.ne.1.and.iflag.ne.2) stop 'Invalid flag passed as 1st arg.'
        call mvsetflags('Data style',7.0) ; call museodir(museod)
        if (iflag.eq.2) then
           write (file,'(a,i4,a,i4,a)') museod(1:len_trim(museod)), &
                       lanno,'/temporaria',lanno,'.dat'
           strum=' temperatures'
        else
           write (file,'(a,i4,a,i4,a)') museod(1:len_trim(museod)), &
                       lanno,'/pioggiaoraria',lanno,'.dat'
           strum=' rain gauges'
        endif
        open (11,file='temporary.t',status='unknown',form='unformatted')
        write(6,'(10x,a,i4,a,i4)') 'Updating record from ',1,' to ', &
               index1h(23,31,12,lanno)
        if (logf) write(6,'(x,a,10x,a,7x,a)') &
               'Rec','Date','Arssa  iira Um/La acqwa Dewet'
        do irec=1,index1h(23,31,12,lanno)
           if (mod(irec,24).eq.0) then
              iday=irec/24
              ihou=24
           else
              iday=irec/24+1
              ihou=mod(irec,24)
           endif
           call dayfromindex(iday,lgiorno,lmese,lanno)
           lidate=mm5index(ihou,lgiorno,lmese,lanno)
           call datafromindex(iday,lanno,data)
           write(data,'(a,i2,a)') data(1:len_trim(data))//' h: ',ihou,'.00'
           n=0 ; n1=0 ; n2=0 ; n3=0 ; n4=0 ; n5=0
!                                                       arssa ---> source=1
!                                                        iira ---> source=2
!                                                       lazio ---> source=3
!                                                       acqwa ---> source=4
!                                                     dewetra ---> source=5
           if (lanno.ge.2000.and.dew) then
              call dewetra(iflag,lanno,irec,llat,llon,source,lp,n)
              n5=n
           endif
!          if (lanno.ge.1997.and.lanno.le.2007.and.n5.eq.0) then
           if (lanno.ge.1997.and.lanno.le.2007) then
              call arssandb (iflag,lanno,irec,llat,llon,source,lp,n)
              n1=n-n5
           endif
!          if (lanno.ge.2003.and.lanno.le.2011.and.n5.eq.0) then
           if (lanno.ge.2003.and.lanno.le.2011) then
              call iira (iflag,lanno,irec,llat,llon,source,lp,n)
              n2=n-(n1+n5)
           endif
!          if (lanno.ge.2007.and.lanno.le.2011.and.n5.eq.0) then
           if (lanno.ge.2007.and.lanno.le.2011) then
              call lazio(iflag,lanno,irec,llat,llon,source,lp,n)
              n3=n-(n2+n1+n5)
           endif
           if (lanno.ge.1999.and.lanno.le.2010.and.n5.eq.0) then
              call acqwa(iflag,lanno,irec,llat,llon,source,lp,n)
              n4=n-(n3+n2+n1+n5)
           endif
           if (logf) write(6,'(i5,a20,4(i5,a),2(i5,a))') irec,' ('// &
                data(1:len_trim(data))//')  ',n1,'+',n2,'+',n3,'+',n4,'+',n5, &
                '=',n,strum(1:len_trim(strum))
!          do i=1,n ; write(6,'(i5,3x,4f8.4)') i,p(i),lat(i),lon(i) ; enddo
           if (iflag.eq.1) call validazione (llat,llon,source,p,n)
           write(11) n,(llat(i),llon(i),source(i),lp(i),i=1,n)
        enddo
        if (logf) write(6,'(x,a,10x,a,7x,a)') &
               'Rec','Date','Arssa  iira Um/La acqwa Dewet'
!------------------------------------------------------------------------------
!                               Probabilmente queste righe sono un residuo della
!                               vecchia versione in cui veniva creato il file 
!                               inizialmente vuoto. Oppure il GL si sta 
!                               rincoglionendo
        do irec=1,366*24
!          write(11) n,(lat(i),lon(i),source(i),p(i),i=1,n)
        enddo
!------------------------------------------------------------------------------
        close(11)
        write(6,'(5x,a)')  'Updating '//file(1:len_trim(file))//'...'
!        call system('/bin/cp ${IO}.t '//file)
!        call system('/bin/rm ${IO}.t')
        write(6,'(5x,a)')  'Done.'
        end program hourlydb
